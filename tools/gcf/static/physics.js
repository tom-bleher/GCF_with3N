"use strict";
// The "Physics inputs" tab of the form: the nuclear inputs of the chosen setup
// with their table values, read from the generator itself, and what can override them.

let physicsParams = null, physicsKey = null, physicsSeq = 0;
const PHYSICS_ROWS = { cpp0: "contacts", cnn0: "contacts", cpn0: "contacts", cpn1: "contacts" };
// Used by genQE -M. Measured on 12C, so it is only offered there (#153).
const ESTAR_12C = { estar: 0.01732, sigma_e: 0.009571 };

function formTabs() { return [...document.querySelectorAll("[data-tab-panel]")]; }

function showFormTab(name) {
  for (const p of formTabs()) p.classList.toggle("hidden", p.dataset.tabPanel !== name);
  document.querySelectorAll("#form-tabs [role=tab]").forEach(b => {
    const on = b.dataset.tab === name;
    b.setAttribute("aria-selected", on ? "true" : "false");
    b.tabIndex = on ? 0 : -1;
  });
  try { localStorage.setItem("gcf-gui-tab", name); } catch (e) {}
  if (name === "physics") refreshPhysics();
  if (name === "scan" && typeof refreshScanPanel === "function") refreshScanPanel();
}

function buildFormTabs() {
  const strip = el("div", { id: "form-tabs", role: "tablist", "aria-label": "Parts of the form" },
    ...[["setup", "Run setup"], ["physics", "Physics inputs"], ["scan", "Scan"]].map(([name, label]) =>
      el("button", { type: "button", role: "tab", "data-tab": name, id: `tab-${name}`, "aria-controls": `panel-${name}`,
        onclick: () => showFormTab(name),
        onkeydown: e => {
          if (e.key !== "ArrowLeft" && e.key !== "ArrowRight") return;
          const order = ["setup", "physics", "scan"].filter(t => !document.getElementById(`tab-${t}`).disabled);
          const next = order[(order.indexOf(name) + (e.key === "ArrowRight" ? 1 : order.length - 1)) % order.length];
          showFormTab(next); document.getElementById(`tab-${next}`).focus();
        } }, label, el("span", { class: "n", id: `tab-count-${name}` }))));
  document.getElementById("form-tabs-slot").replaceWith(strip);
  let start = "setup";
  try { start = localStorage.getItem("gcf-gui-tab") || "setup"; } catch (e) {}
  showFormTab(["physics", "scan"].includes(start) ? start : "setup");
}

function physicsAvailable() { return stages().includes("generate"); }

function updatePhysicsTab() {
  const tab = document.getElementById("tab-physics");
  if (!tab) return;
  tab.disabled = !physicsAvailable();
  tab.title = tab.disabled ? "These settings are for the Generate stage" : "";
  if (tab.disabled && tab.getAttribute("aria-selected") === "true") showFormTab("setup");
  const scanTab = document.getElementById("tab-scan");
  scanTab.disabled = stages().includes("program") || stages().includes("eic");
  scanTab.title = scanTab.disabled ? "Scans are for the Generate and FSI stages" : "";
  if (scanTab.disabled && scanTab.getAttribute("aria-selected") === "true") showFormTab("setup");
  let n = 0;
  document.querySelectorAll('[data-tab-panel="physics"] .field:not(.hidden)').forEach(f => {
    const o = opt(f.dataset.section, f.dataset.key);
    try { const v = valueOf(f.dataset.section, o); if (v !== null && JSON.stringify(v) !== JSON.stringify(o.default)) n++; } catch (e) { n++; }
  });
  document.getElementById("tab-count-physics").textContent = n ? ` (${n} changed)` : "";
  document.getElementById("physics-note").classList.toggle("hidden", m3NNote());
}
function m3NNote() { return model() !== "3N"; }

function setTableValue(key, value, err) {
  const f = document.querySelector(`.field[data-section="generate"][data-key="${key}"]`);
  if (!f) return;
  const input = f.querySelector("input[type=number]");
  if (input) input.placeholder = value === null || value === undefined ? "default" : String(value);
  let line = f.querySelector(".tableval");
  if (!line) { line = el("span", { class: "tableval" }); f.append(line); }
  line.textContent = value === null || value === undefined ? "" :
    `Table: ${value}` + (err ? ` ± ${err}` : "");
}

async function refreshPhysics() {
  if (!physicsAvailable() || !document.getElementById("panel-physics") ||
      document.getElementById("panel-physics").classList.contains("hidden")) return;
  const m = model();
  const A = document.getElementById(fieldId("generate", "target", "A")).value, Z = document.getElementById(fieldId("generate", "target", "Z")).value;
  const nnField = document.getElementById(fieldId("generate", "nn"));
  const nn = m === "2N" && nnField ? nnField.value : "";
  const key = [m, A, Z, nn].join("/");
  if (key === physicsKey) return;
  physicsKey = key;
  const seq = ++physicsSeq, head = document.getElementById("physics-head");
  head.textContent = "Reading the tables…";
  let p;
  try { p = await api(`api/nuclear?model=${encodeURIComponent(m)}&A=${encodeURIComponent(A)}&Z=${encodeURIComponent(Z)}&nn=${encodeURIComponent(nn)}`); }
  catch (e) { if (seq === physicsSeq) { physicsParams = null; head.textContent = "No table values: " + e.message; } return; }
  if (seq !== physicsSeq) return;
  physicsParams = p;
  const name = ((TARGETS[m] || []).find(t => t.A === p.A && t.Z === p.Z) || {}).name || `A = ${p.A}, Z = ${p.Z}`;
  const parts = [`${name} in the ${m} model`];
  if (p.nn) parts.push(`${p.nn} interaction`);
  if (p.contact_source) parts.push(`contacts: ${p.contact_source}`);
  head.textContent = parts.join(" · ") + ".";
  setTableValue("sigma_cm", p.sigma_cm ? p.sigma_cm[0] : null, p.sigma_cm ? p.sigma_cm[1] : null);
  setTableValue("estar", p.estar, null);
  setTableValue("sigma_e", null);
  const contacts = p.contacts || {};
  for (const k of Object.keys(PHYSICS_ROWS)) setTableValue(k, contacts[k] ? contacts[k][0] : null, contacts[k] ? contacts[k][1] : null);
  setTableValue("kF", p.kF === undefined ? null : p.kF);
  document.getElementById("physics-preset").classList.toggle("hidden", !(m === "2N" && p.A === 12 && p.Z === 6));
  const masses = Object.entries(p.mRes || {}).map(([k, v]) => `${k} ${Number(v).toFixed(4)}`).join(", ");
  document.getElementById("physics-masses").textContent = `Masses (GeV): target ${Number(p.mA).toFixed(4)}; residual ${masses}` +
    (p.estar_max ? `. E* range used when varying the inputs: 0 to ${p.estar_max} GeV.` : ".");
}

function clearPhysics() {
  document.querySelectorAll('[data-tab-panel="physics"] .field').forEach(f => {
    const o = opt(f.dataset.section, f.dataset.key);
    const input = document.getElementById(fieldId(f.dataset.section, f.dataset.key));
    if (!input) return;
    if (o.kind === "bool") input.checked = !!o.default; else input.value = o.default === null ? "" : o.default;
  });
  changed();
}

function fillEstar12C() {
  for (const [k, v] of Object.entries(ESTAR_12C)) document.getElementById(fieldId("generate", k)).value = v;
  changed();
}
