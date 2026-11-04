"use strict";
// Looking after runs: the machine's state, resuming, the queue, downloads, config
// files and examples, a first-run guide and notifications.

function fmtSize(n) { return n === null || n === undefined ? "" : fmtBytes(n); }

function renderMachine(m) {
  const box = document.getElementById("machine");
  if (!box || !m) return;
  const parts = [`${m.in_use} of ${m.max_cpus} CPUs in use by these runs`];
  if (m.queued) parts.push(`${m.queued} waiting`);
  if (m.load !== null) parts.push(`machine load ${m.load.toFixed(0)}`);
  parts.push(`${fmtSize(m.disk_free)} free`);
  box.textContent = parts.join(" · ");
  box.classList.toggle("warn-text", m.load !== null && m.load > m.cpus);
  box.title = m.load !== null && m.load > m.cpus ? "This machine is busier than it has CPUs, so runs will be slower than the estimate." : "";
}

function runActions(name, d, actions, body) {
  if (d.state === "queued") {
    body.append(el("p", { class: "setup" }, "Waiting for free CPUs. It starts on its own, also if you close this page."));
    actions.append(el("button", { type: "button", class: "danger", onclick: () => cancel(name) }, "Remove from the queue"));
  }
  if (d.resumable) {
    const left = d.jobs_total - d.jobs_finished;
    actions.prepend(el("button", { type: "button", class: "primary", onclick: () => resumeRun(name),
      title: "Keeps the jobs that finished and runs the others again, with the same seeds" },
      d.state === "partial" ? "Finish the run" : `Resume (${left} of ${d.jobs_total} job${d.jobs_total > 1 ? "s" : ""} left)`));
  }
  if (d.mergeable) {
    actions.append(el("button", { type: "button", onclick: () => mergeFinished(name),
      title: "Make output files from the jobs that finished, without running the others" }, "Use the finished jobs"));
  }
  if (d.partial) {
    const text = Object.entries(d.partial).map(([s, p]) => `${STAGE_NAMES[s] || s}: ${fmtInt(p.events)} of ${fmtInt(p.total)} events`).join("; ");
    body.append(el("p", { class: "note note-warn" }, `These outputs hold part of the run: ${text}. Σw / attempts is still right, because each job's attempts are counted.`));
  }
  if (d.size) body.append(el("p", { class: "setup" }, `${fmtSize(d.size)} on disk`));
  if (d.state === "stopped" && !d.orphans && d.resumable) {
    body.append(el("p", { class: "setup" }, "The run stopped without finishing, for example because the machine restarted. Nothing it finished is lost."));
  }
}

async function resumeRun(name) {
  try { await api(`api/runs/${encodeURIComponent(name)}/resume`, {}); toast("Resuming " + name); }
  catch (e) { alert(e.message); }
  refreshDetails(); showDetails(name);
}

async function mergeFinished(name) {
  if (!confirm("Make output files from the jobs that finished? You can still resume the run afterwards.")) return;
  document.getElementById("d-body").prepend(el("p", { class: "muted" }, "Merging…"));
  try { await api(`api/runs/${encodeURIComponent(name)}/merge`, {}); }
  catch (e) { alert(e.message); }
  delete resultsAsked[name + "/generate"]; delete resultsAsked[name + "/fsi"];
  await refreshDetails();
}

function outputButtons(name, stage, o, d) {
  const file = o.path.split("/").pop();
  const out = [el("a", { class: "btn small" + (o.exists ? "" : " hidden"), href: tokenUrl(`api/runs/${encodeURIComponent(name)}/${encodeURIComponent(file)}`) + "&download=1",
    download: file, title: "Download the event file" }, "Download")];
  if (stage === "generate" && o.exists && !["running", "merging", "starting", "queued"].includes(d.state)) {
    out.push(el("button", { type: "button", class: "small", title: "Fill the form to run the GENIE cascade on these events", onclick: () => addFsi(o.path) }, "Add FSI…"));
  }
  return out;
}

function addFsi(path) {
  document.querySelector("input[name=stages][value='fsi']").checked = true;
  updateVisibility();
  document.getElementById(fieldId("fsi", "input")).value = path;
  document.getElementById(fieldId("run", "name")).value = defaultName();
  showFormTab("setup");
  changed();
  document.getElementById("form-card").scrollIntoView({ block: "start" });
  toast("Choose the cascade settings and press Start run");
}

function logLink(name, l) {
  const base = `api/runs/${encodeURIComponent(name)}/` + l.file.split("/").map(encodeURIComponent).join("/");
  return el("a", { href: tokenUrl(base), target: "_blank", rel: "noopener", class: "loglink" }, `full log (${fmtSize(l.size)})`);
}

// -- config files and examples ---------------------------------------------
const EXAMPLE_LABELS = { "3N": "3N on 12C, with FSI", "2N": "2N on 12C, with FSI", MF: "Mean field on 12C", fsi: "FSI on an existing file", dis: "Tagged DIS on 12C (genDIS)" };

async function loadConfigText(text, source, keepName) {
  let cfg;
  try { cfg = await api("api/parse-config", { text }); }
  catch (e) { showMsg(`${source}: ${plain(e)}`, "err"); return false; }
  if (cfg.run) delete cfg.run.output_dir;
  applyConfig(cfg);
  if (!keepName || !(cfg.run && cfg.run.name)) document.getElementById(fieldId("run", "name")).value = defaultName();
  showFormTab("setup");
  changed();
  toast("Form filled from " + source);
  return true;
}

function buildConfigTools(examples) {
  const bar = document.querySelector("#form-card .actions");
  const file = el("input", { type: "file", accept: ".toml,text/plain", class: "hidden", id: "config-file",
    onchange: async () => { if (file.files[0]) await loadConfigText(await file.files[0].text(), file.files[0].name, true); file.value = ""; } });
  const sel = el("select", { id: "example-select", "aria-label": "Fill the form from an example",
    onchange: async () => { const k = sel.value; sel.value = ""; if (k) await loadConfigText(examples[k], `the ${k} example`, false); } },
    el("option", { value: "" }, "Examples…"), ...Object.keys(examples).map(k => el("option", { value: k }, EXAMPLE_LABELS[k] || k)));
  bar.append(el("button", { type: "button", onclick: () => file.click(), title: "Fill the form from a TOML config you saved" }, "Load config…"), sel, file);
}

// -- first run --------------------------------------------------------------
function startingGuide() {
  return el("div", { class: "guide" },
    el("h3", {}, "Getting started"),
    el("ol", {},
      el("li", {}, "Press ", el("button", { type: "button", class: "small", onclick: fillTestRun }, "Fill in a small test run"), " to get a run that takes under a minute, then ", el("b", {}, "Start run"), "."),
      el("li", {}, "Click the run here to follow it. When it finishes you get its numbers, plots and output files."),
      el("li", {}, "Change settings in the form, on ", el("b", {}, "Run setup"), " and ", el("b", {}, "Physics inputs"), ". Each ", el("b", {}, "?"), " explains its field, and ", el("b", {}, "Check & show config"), " estimates the time and disk space.")),
    el("p", { class: "setup" }, "Runs carry on when you close this page or log out. Come back any time to the same address. The ", el("a", { href: tokenUrl("docs") , target: "_blank", rel: "noopener" }, "docs"), " explain what the generators compute."));
}

function fillTestRun() {
  resetForm();
  document.getElementById(fieldId("generate", "model")).value = "2N";
  updateVisibility();
  document.getElementById(fieldId("generate", "events")).value = 5000;
  const jobs = Math.min(2, Number(SCHEMA.run.find(o => o.key === "jobs").default) + 1);
  document.getElementById(fieldId("run", "jobs")).value = jobs;
  document.getElementById(fieldId("run", "name")).value = "first-test-" + defaultName().split("-").slice(-1)[0];
  changed();
  document.getElementById("btn-start").focus();
  toast("A small test run is ready: press Start run");
}

// -- problems of this machine -------------------------------------------------
function renderDoctor(results) {
  const box = document.getElementById("doctor");
  const bad = (results || []).filter(c => !c.ok);
  if (!box) return;
  box.classList.toggle("hidden", !bad.length);
  box.replaceChildren(...bad.map(c => el("div", { class: c.required ? "doctor-bad" : "" },
    el("b", {}, c.name), c.message ? `: ${c.message}` : "", c.fix ? " → " : "", c.fix ? el("code", {}, c.fix) : "")));
}

// -- notifications ------------------------------------------------------------
let lastStates = {};
function notifyChanges() {
  const on = (() => { try { return localStorage.getItem("gcf-gui-notify") === "1"; } catch (e) { return false; } })();
  for (const r of runsCache) {
    const before = lastStates[r.name];
    if (on && before && before !== r.state && ["done", "failed", "partial"].includes(r.state) && window.Notification && Notification.permission === "granted") {
      new Notification(`GCF run ${r.state}: ${r.name}`, { body: r.setup || "" });
    }
    lastStates[r.name] = r.state;
  }
}

function buildNotifyToggle() {
  if (!window.Notification) return;
  const box = el("input", { type: "checkbox", id: "notify-toggle" });
  try { box.checked = localStorage.getItem("gcf-gui-notify") === "1"; } catch (e) {}
  box.addEventListener("change", async () => {
    if (box.checked && Notification.permission !== "granted") box.checked = (await Notification.requestPermission()) === "granted";
    try { localStorage.setItem("gcf-gui-notify", box.checked ? "1" : "0"); } catch (e) {}
  });
  document.getElementById("machine-line").append(el("label", { class: "notify", for: "notify-toggle", title: "A desktop notification when a run ends, while this page is open" }, box, "Notify me"));
}
