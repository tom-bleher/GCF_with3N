"use strict";
// The numbers of a run's outputs, their warnings and provenance, and the
// selection applied to its plots.

const cutPlots = {};        // run + stage -> result of the last selection
const resultsAsked = {};    // run + stage -> summary requested

function fmtNum(x, digits = 4) {
  if (x === null || x === undefined) return "–";
  return Math.abs(x) >= 1e5 || (x !== 0 && Math.abs(x) < 1e-3) ? x.toExponential(digits - 1) : Number(x.toPrecision(digits)).toString();
}

function resultsRow(label, value, help) {
  return [el("dt", { title: help || null }, label), el("dd", {}, value)];
}

function resultsPanel(name, stage, o, d) {
  const box = el("div", { class: "results" });
  const key = name + "/" + stage;
  const r = (d.results || {})[stage];
  if (!r) {
    if (o.exists && d.state === "done" && !resultsAsked[key]) {
      resultsAsked[key] = true;
      api(`api/runs/${encodeURIComponent(name)}/results`, { stage }).then(() => refreshDetails())
        .catch(e => { box.replaceChildren(el("p", { class: "setup" }, "No summary: " + e.message)); });
    }
    if (o.exists) box.append(el("p", { class: "setup" }, "Reading the summary…"));
    return box;
  }
  const yieldText = r.yield === null ? "–" : `${fmtNum(r.yield)} ± ${fmtNum(r.yield_err, 2)}`;
  box.append(el("dl", { class: "numbers" },
    ...resultsRow("Events", fmtInt(r.events)),
    ...resultsRow("Attempts", r.attempts ? `${fmtInt(r.attempts)} (${fmtNum(100 * r.efficiency, 3)}% accepted)` : "–", "Generation attempts, including rejected ones"),
    ...resultsRow("Σw / attempts", yieldText, "The quantity to compare between samples of the same model. Statistical error only."),
    ...resultsRow("N_eff", r.n_eff === null ? "–" : `${fmtInt(Math.round(r.n_eff))}${r.n_eff > 0 ? ` (±${fmtNum(100 / Math.sqrt(r.n_eff), 2)}% on the yield)` : ""}`, "(Σw)²/Σw²: how many equal-weight events this sample is worth"),
    ...(r.max_w_fraction === null ? [] : resultsRow("Heaviest event", `${fmtNum(100 * r.max_w_fraction, 2)}% of Σw`)),
    ...(r.fsi_ran_fraction === null ? [] : resultsRow("Cascade ran", `${fmtNum(100 * r.fsi_ran_fraction, 4)}% of events`))));
  for (const [level, text] of r.notes || []) box.append(el("p", { class: `note note-${level}` }, text));
  const prov = Object.entries(r.provenance || {});
  if (prov.length) {
    box.append(el("details", {}, el("summary", {}, `Provenance: ${stage === "generate" ? "PWIA" : "with FSI"}`),
      ...prov.map(([s, p]) => el("dl", { class: "numbers" },
        ...resultsRow("Stage", s), ...resultsRow("Code version", p.versions.join(", ")), ...resultsRow("Date", p.date),
        ...resultsRow("Jobs", String(p.seeds.length)), ...resultsRow("Seeds", p.seeds.slice(0, 6).join(", ") + (p.seeds.length > 6 ? ", …" : "")),
        ...resultsRow("Command", el("code", {}, p.command))))));
  }
  return box;
}

const CUT_FIELDS = [["q2", "Q²", "GeV²"], ["xb", "x_B", ""], ["theta_e", "θ_e", "deg"], ["pmiss", "p_miss", "GeV/c"], ["protons", "Protons", "after FSI if run"]];

function selectionForm(name, stage, o) {
  const key = name + "/" + stage;
  const last = cutPlots[key];
  const inputs = {};
  const rows = CUT_FIELDS.map(([k, label, unit]) => {
    const vals = last ? (last.cuts[k] || [null, null]) : [null, null];
    inputs[k] = [0, 1].map(i => el("input", { type: "number", step: "any", placeholder: i ? "max" : "min",
      value: vals[i] === null ? null : vals[i], "aria-label": `${label} ${i ? "maximum" : "minimum"}` }));
    return el("div", { class: "cut" }, el("label", {}, label, unit ? el("span", { class: "unit" }, ` (${unit})`) : null), el("div", { class: "pair" }, ...inputs[k]));
  });
  const msg = el("p", { class: "setup" }, last ? `${last.description}: ${fmtInt(last.events)} events, N_eff = ${fmtInt(Math.round(last.n_eff))}, ${fmtNum(100 * last.fraction, 3)}% of Σw` + (last.yield_ === null ? "" : `, Σw / attempts = ${fmtNum(last.yield_)} ± ${fmtNum(last.yield_err, 2)}`) : "");
  const apply = el("button", { type: "button", class: "small", onclick: async () => {
    const cuts = {};
    for (const [k] of CUT_FIELDS) {
      const vals = inputs[k].map(i => i.value.trim() === "" ? null : Number(i.value));
      if (vals.some(x => x !== null)) cuts[k] = vals;
    }
    if (!Object.keys(cuts).length) { delete cutPlots[key]; refreshDetails(); return; }
    apply.disabled = true; msg.textContent = "Drawing…";
    try {
      const r = await api(`api/runs/${encodeURIComponent(name)}/quicklook`, { stage, cuts });
      cutPlots[key] = Object.assign(r.stats, { plot: r.plot, t: Date.now() });
      await refreshDetails();
    } catch (e) { msg.textContent = e.message; apply.disabled = false; }
  } }, "Apply selection");
  const clear = el("button", { type: "button", class: "small", onclick: () => { delete cutPlots[key]; refreshDetails(); } }, "Clear");
  return el("details", { open: !!last }, el("summary", {}, `Selection: ${stage === "generate" ? "PWIA" : "with FSI"}`), el("div", { class: "cuts" }, ...rows),
    el("div", { class: "actions", style: "margin-top:6px" }, apply, last ? clear : null), msg);
}

function selectedPlot(name, stage) {
  const last = cutPlots[name + "/" + stage];
  if (!last) return null;
  const src = tokenUrl(`api/runs/${encodeURIComponent(name)}/${encodeURIComponent(last.plot)}`);
  return el("a", { href: src, target: "_blank", rel: "noopener", title: "Open full size" },
    el("img", { class: "plot", alt: `Control plots with the selection ${last.description}`, src: src + "&t=" + last.t }));
}
