"use strict";
// The Scan tab of the form, the list of scans, and a scan's table and plot.

let selectedScan = null, scanTimer = null, scansCache = [], scanPreviewTimer = null;

function scanParameters() {
  const st = stages(), out = [];
  if (st.includes("generate")) {
    const m = model();
    for (const o of SCHEMA.generate) {
      if (!o.models.includes(m) || ["model", "target", "phase_space", "events"].includes(o.key)) continue;
      if (o.kind === "float" || o.kind === "int" || o.kind === "str") out.push({ path: `generate.${o.key}`, label: o.label + (o.unit ? ` (${o.unit})` : ""), o, group: "Generator" });
    }
  }
  if (st.includes("fsi")) {
    for (const o of SCHEMA.fsi) {
      if (["input", "keep", "tree", "off", "scale", "target", "overwrite"].includes(o.key)) continue;
      if (o.kind === "float" || o.kind === "str") out.push({ path: `fsi.${o.key}`, label: o.label, o, group: "FSI" });
    }
    const fm = fsiModel();
    for (const c of CHANNELS.filter(c => (fm === "hA" ? c.hA : c.hN) && c.scalable))
      out.push({ path: `fsi.scale:${c.name}`, label: `${c.name} scale factor`, o: { kind: "float" }, group: "FSI channels" });
  }
  out.push({ path: "run.seed", label: "Base seed (repeat with other events)", o: { kind: "int" }, group: "Statistics" });
  return out;
}

function scanAxisRow(path, values) {
  const sel = el("select", { class: "scan-path", "aria-label": "Parameter to scan" });
  const fill = () => {
    const cur = sel.value || path, params = scanParameters();
    const groups = {};
    for (const p of params) (groups[p.group] = groups[p.group] || []).push(p);
    sel.replaceChildren(...Object.entries(groups).map(([g, ps]) => el("optgroup", { label: g }, ...ps.map(p => el("option", { value: p.path }, p.label)))));
    if (params.some(p => p.path === cur)) sel.value = cur;
  };
  fill();
  const input = el("input", { type: "text", class: "scan-values", value: values || "", "aria-label": "Values",
    placeholder: "0.10, 0.15, 0.20   or   0.1:0.2:5" });
  const suggest = () => {
    const p = scanParameters().find(x => x.path === sel.value);
    if (p && p.o.choices && !input.value) input.value = p.o.choices.join(", ");
    input.placeholder = p && p.o.choices ? p.o.choices.join(", ") : "0.10, 0.15, 0.20   or   0.1:0.2:5";
    scheduleScanPreview();
  };
  sel.addEventListener("change", () => { input.value = ""; suggest(); });
  input.addEventListener("input", scheduleScanPreview);
  const row = el("div", { class: "scan-axis" }, sel, input, el("button", { type: "button", class: "small", "aria-label": "Remove this parameter",
    onclick: () => { row.remove(); scheduleScanPreview(); } }, "✕"));
  row.refill = fill;
  suggest();
  return row;
}

function readScanAxes() {
  return [...document.querySelectorAll(".scan-axis")].map(r => ({ path: r.querySelector(".scan-path").value, values: r.querySelector(".scan-values").value }));
}

function buildScanPanel() {
  const panel = document.getElementById("panel-scan");
  panel.append(
    el("p", { class: "setup" }, "Run the setup from the other tabs for several values of a parameter. Each point is an ordinary run that waits its turn in the queue, and all points share a seed."),
    el("fieldset", {}, el("legend", {}, "Parameters"), el("div", { id: "scan-axes" }),
      el("button", { type: "button", class: "small", onclick: () => { document.getElementById("scan-axes").append(scanAxisRow("", "")); } }, "Add a parameter")),
    el("p", { id: "scan-preview", class: "setup", "aria-live": "polite" }),
    el("div", { class: "actions" }, el("button", { type: "button", class: "primary", id: "btn-scan", onclick: startScan }, "Start scan")));
  document.getElementById("scan-axes").append(scanAxisRow("", ""));
}

function refreshScanPanel() {
  document.querySelectorAll(".scan-axis").forEach(r => r.refill && r.refill());
  scheduleScanPreview();
}

function scheduleScanPreview() {
  clearTimeout(scanPreviewTimer);
  scanPreviewTimer = setTimeout(previewScan, 300);
}

async function previewScan() {
  const box = document.getElementById("scan-preview");
  if (!box) return;
  const axes = readScanAxes();
  if (!axes.some(a => a.values.trim())) { box.textContent = "Give the values to run."; box.classList.remove("msg", "err"); return; }
  try {
    const r = await api("api/scans/preview", { config: collect(), axes });
    const ev = r.events ? ` · ${fmtInt(r.events)} events in all` : "";
    box.textContent = `${r.points} run${r.points === 1 ? "" : "s"}${ev}, such as ${r.first}`;
    box.className = "setup";
  } catch (e) { box.textContent = friendly(plain(e)); box.className = "msg err"; }
}

async function startScan() {
  const btn = document.getElementById("btn-scan"), box = document.getElementById("scan-preview");
  btn.disabled = true; box.textContent = "Checking every point and starting…"; box.className = "setup";
  try {
    const r = await api("api/scans", { config: collect(), axes: readScanAxes() });
    box.textContent = `Submitted scan ${r.name} with ${r.points} runs. It carries on if you close this page.`;
    document.getElementById(fieldId("run", "name")).value = defaultName();
    await refreshRuns();
    showScan(r.name);
  } catch (e) { box.textContent = friendly(plain(e)); box.className = "msg err"; }
  finally { btn.disabled = false; }
}

// -- the list of scans ----------------------------------------------------------
function renderScans(scans) {
  scansCache = scans || [];
  const box = document.getElementById("scans");
  if (!box) return;
  box.classList.toggle("hidden", !scansCache.length);
  box.replaceChildren(el("h3", {}, "Scans"), ...scansCache.map(s => {
    const done = s.states.done || 0;
    const active = (s.states.running || 0) + (s.states.queued || 0) + (s.states.starting || 0) + (s.states.merging || 0);
    const failed = (s.states.failed || 0) + (s.states.stopped || 0) + (s.states.cancelled || 0);
    return el("button", { type: "button", class: "scan-row" + (s.name === selectedScan ? " selected" : ""), onclick: () => showScan(s.name) },
      el("span", {}, s.name, el("span", { class: "setup" }, ` · ${s.axes.join(" × ")}`)),
      el("span", { class: "setup" }, `${done} of ${s.points} done` + (active ? `, ${active} active` : "") + (failed ? `, ${failed} not finished` : "")));
  }));
}

// -- a scan's results --------------------------------------------------------------
function showScan(name) {
  selectedScan = name;
  selected = null;
  if (detailsTimer) { clearInterval(detailsTimer); detailsTimer = null; }
  document.getElementById("d-title").textContent = "Scan " + name;
  document.getElementById("d-body").replaceChildren(el("p", { class: "muted" }, "Reading…"));
  refreshScan();
  if (scanTimer) clearInterval(scanTimer);
  scanTimer = setInterval(refreshScan, 3000);
  renderRuns(); renderScans(scansCache);
  if (window.innerWidth <= 980) document.getElementById("details").scrollIntoView({ block: "start" });
}

function leaveScan() {
  selectedScan = null;
  if (scanTimer) { clearInterval(scanTimer); scanTimer = null; }
}

async function refreshScan() {
  const name = selectedScan;
  if (!name) return;
  let s;
  try { s = await api(`api/scans/${encodeURIComponent(name)}`); }
  catch (e) { document.getElementById("d-body").replaceChildren(el("p", { class: "msg err" }, plain(e))); return; }
  if (name !== selectedScan) return;
  const active = s.points.some(p => ["running", "merging", "starting", "queued"].includes(p.state));
  if (!active && scanTimer) { clearInterval(scanTimer); scanTimer = null; }
  const body = el("div");
  const counts = {};
  s.points.forEach(p => { counts[p.state] = (counts[p.state] || 0) + 1; });
  body.append(el("div", { class: "summary" }, Object.entries(counts).map(([k, v]) => `${v} ${k}`).join(" · ")));
  const chart = scanChart(s);
  if (chart) body.append(chart);
  for (const n of s.notes || []) body.append(el("p", { class: "note" }, n));
  const first = s.points.find(p => p.yield !== null);
  body.append(el("table", { class: "compact scan-table" },
    el("thead", {}, el("tr", {}, ...s.axes.map(a => el("th", { scope: "col" }, a.label)), ...["State", "N_eff", "Σw / attempts", first ? "Ratio to first" : ""].map(h => el("th", { scope: "col" }, h)))),
    el("tbody", {}, ...s.points.map(p => el("tr", { class: "run", tabindex: "0", role: "button", title: p.name, "aria-label": `Open run ${p.name}`,
      onclick: () => { leaveScan(); showDetails(p.name); }, onkeydown: e => { if (e.key === "Enter") { leaveScan(); showDetails(p.name); } } },
      ...p.values.map(v => el("td", {}, typeof v === "number" ? fmtNum(v) : String(v))), el("td", {}, badge(p.state)),
      el("td", {}, p.n_eff === null ? "–" : fmtInt(Math.round(p.n_eff))),
      el("td", {}, p.yield === null ? "–" : `${fmtNum(p.yield)} ± ${fmtNum(p.yield_err, 2)}`),
      el("td", {}, p.yield === null || !first ? "–" : fmtNum(p.yield / first.yield)))))));
  const actions = el("div", { class: "actions" },
    el("button", { type: "button", onclick: () => downloadScanTable(s) }, "Download table (CSV)"),
    el("button", { type: "button", disabled: !active ? null : true, class: "danger", onclick: () => deleteScan(s.name) }, "Delete all its runs"));
  body.append(actions);
  document.getElementById("d-body").replaceChildren(...body.childNodes);
}

// Σw/attempts against the parameter, for a scan of one numeric parameter.
function scanChart(s) {
  if (s.axes.length !== 1) return null;
  const pts = s.points.filter(p => p.yield !== null && typeof p.values[0] === "number");
  if (pts.length < 2 || pts.length !== s.points.filter(p => p.yield !== null).length) return null;
  const W = 520, H = 230, L = 62, R = 14, T = 12, B = 40;
  const xs = pts.map(p => p.values[0]), ys = pts.flatMap(p => [p.yield - p.yield_err, p.yield + p.yield_err]);
  const pad = (a, b) => { const d = (b - a) || Math.abs(a) || 1; return [a - 0.08 * d, b + 0.08 * d]; };
  const [x0, x1] = pad(Math.min(...xs), Math.max(...xs)), [y0, y1] = pad(Math.min(...ys), Math.max(...ys));
  const X = v => L + (v - x0) / (x1 - x0) * (W - L - R), Y = v => H - B - (v - y0) / (y1 - y0) * (H - T - B);
  const NS = "http://www.w3.org/2000/svg";
  const svg = document.createElementNS(NS, "svg");
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  svg.setAttribute("class", "scan-chart");
  svg.setAttribute("role", "img");
  svg.setAttribute("aria-label", `Σw per attempt against ${s.axes[0].label}; the table below has the same values`);
  const add = (tag, attrs, text) => { const e = document.createElementNS(NS, tag); for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v); if (text !== undefined) e.textContent = text; svg.append(e); return e; };
  const ticks = (a, b, n) => { const step = Math.pow(10, Math.floor(Math.log10((b - a) / n))); const err = (b - a) / n / step; const k = err >= 7 ? 10 : err >= 3 ? 5 : err >= 1.5 ? 2 : 1;
    const s_ = k * step, out = []; for (let v = Math.ceil(a / s_) * s_; v <= b + 1e-12; v += s_) out.push(Number(v.toPrecision(6))); return out; };
  for (const v of ticks(y0, y1, 4)) {
    add("line", { x1: L, x2: W - R, y1: Y(v), y2: Y(v), class: "scan-grid" });
    add("text", { x: L - 6, y: Y(v) + 4, "text-anchor": "end", class: "scan-tick" }, fmtNum(v, 3));
  }
  for (const v of ticks(x0, x1, 5)) add("text", { x: X(v), y: H - B + 16, "text-anchor": "middle", class: "scan-tick" }, fmtNum(v, 3));
  add("line", { x1: L, x2: W - R, y1: H - B, y2: H - B, class: "scan-axis-line" });
  add("text", { x: (L + W - R) / 2, y: H - 6, "text-anchor": "middle", class: "scan-label" }, s.axes[0].label);
  add("text", { x: 12, y: (T + H - B) / 2, "text-anchor": "middle", class: "scan-label", transform: `rotate(-90 12 ${(T + H - B) / 2})` }, "Σw / attempts");
  const sorted = [...pts].sort((a, b) => a.values[0] - b.values[0]);
  add("polyline", { points: sorted.map(p => `${X(p.values[0])},${Y(p.yield)}`).join(" "), class: "scan-line" });
  for (const p of sorted) {
    const x = X(p.values[0]);
    add("line", { x1: x, x2: x, y1: Y(p.yield - p.yield_err), y2: Y(p.yield + p.yield_err), class: "scan-err" });
    add("circle", { cx: x, cy: Y(p.yield), r: 4.5, class: "scan-dot" });
    const hit = add("circle", { cx: x, cy: Y(p.yield), r: 13, class: "scan-hit" });
    const t = document.createElementNS(NS, "title");
    t.textContent = `${s.axes[0].label} = ${fmtNum(p.values[0])}: ${fmtNum(p.yield)} ± ${fmtNum(p.yield_err, 2)}, N_eff ${fmtInt(Math.round(p.n_eff || 0))}`;
    hit.append(t);
  }
  return el("div", { class: "scan-chart-box" }, svg);
}

function downloadScanTable(s) {
  const q = v => /[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : String(v);
  const head = [...s.axes.map(a => a.label), "run", "state", "events", "attempts", "n_eff", "yield", "yield_err"];
  const rows = s.points.map(p => [...p.values, p.name, p.state, p.events ?? "", p.attempts ?? "", p.n_eff ?? "", p.yield ?? "", p.yield_err ?? ""]);
  const text = [head, ...rows].map(r => r.map(q).join(",")).join("\n") + "\n";
  const a = el("a", { href: URL.createObjectURL(new Blob([text], { type: "text/csv" })), download: `${s.name}.csv` });
  document.body.append(a); a.click(); a.remove();
}

async function deleteScan(name) {
  if (!confirm(`Delete scan ${name} and all its runs and event files? This cannot be undone.`)) return;
  try { await api(`api/scans/${encodeURIComponent(name)}/delete`, {}); } catch (e) { alert(e.message); return; }
  leaveScan();
  document.getElementById("d-title").textContent = "Run details";
  document.getElementById("d-body").replaceChildren(el("p", { class: "placeholder" }, "Select a run to see its progress, outputs, plots and logs."));
  refreshRuns();
}
