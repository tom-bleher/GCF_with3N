"use strict";
// Compare the outputs of two runs: tools/validation/compare_samples.py, on the page.

function compareDialog() {
  let dlg = document.getElementById("compare-dialog");
  if (!dlg) {
    dlg = el("dialog", { id: "compare-dialog", class: "wide-dialog", "aria-labelledby": "compare-title" });
    document.body.append(dlg);
  }
  return dlg;
}

const STAGE_NAMES = { generate: "PWIA", fsi: "with FSI" };

function comparableRuns() {
  return runsCache.filter(r => Object.values(r.outputs || {}).some(o => o.exists));
}

function stageOptions(run) {
  return Object.entries(run.outputs || {}).filter(([, o]) => o.exists).map(([s]) => el("option", { value: s }, STAGE_NAMES[s] || s));
}

function openCompare(name) {
  const dlg = compareDialog();
  const runs = comparableRuns();
  const pick = (id, label, chosen) => {
    const sel = el("select", { id: `${id}-run`, "aria-label": `${label} run` }, ...runs.map(r => el("option", { value: r.name, selected: r.name === chosen }, r.name)));
    const stage = el("select", { id: `${id}-stage`, "aria-label": `${label} output` });
    const fill = () => { const cur = stage.value; stage.replaceChildren(...stageOptions(runs.find(r => r.name === sel.value) || {}));
      if ([...stage.options].some(o => o.value === cur)) stage.value = cur; else if ([...stage.options].some(o => o.value === "fsi")) stage.value = "fsi"; };
    sel.addEventListener("change", fill); fill();
    return { sel, stage, row: el("div", { class: "field" }, el("label", { for: `${id}-run` }, label), el("div", { class: "pair" }, sel, stage)) };
  };
  const other = (runs.find(r => r.name !== name) || {}).name || name;
  const a = pick("cmp-a", "Sample A", name), b = pick("cmp-b", "Sample B", other);
  const events = el("input", { type: "number", id: "cmp-events", min: "1000", step: "1000", value: 300000 });
  const out = el("div", { id: "compare-out" });
  const go = el("button", { type: "button", class: "primary", onclick: async () => {
    go.disabled = true; out.replaceChildren(el("p", { class: "muted" }, "Comparing; the weighted test takes up to a few minutes on large samples…"));
    try {
      const r = await api("api/compare", { a: { run: a.sel.value, stage: a.stage.value }, b: { run: b.sel.value, stage: b.stage.value }, max_events: Number(events.value) });
      out.replaceChildren(...compareResult(r));
    } catch (e) { out.replaceChildren(el("p", { class: "msg err" }, plain(e))); }
    finally { go.disabled = false; }
  } }, "Compare");
  dlg.replaceChildren(el("h2", { id: "compare-title" }, "Compare two samples"),
    el("div", { class: "grid" }, a.row, b.row,
      el("div", { class: "field" }, el("label", { for: "cmp-events" }, "Events per sample for the shapes"), events)),
    el("div", { class: "actions" }, go, el("button", { type: "button", onclick: () => dlg.close() }, "Close")), out);
  dlg.showModal();
}

function compareResult(r) {
  const parts = [];
  const flagged = r.observables.filter(o => o.p_weighted !== null && o.p_weighted <= r.alpha);
  const unweighted = r.observables.filter(o => o.p <= r.alpha && !flagged.includes(o));
  const lowest = Math.min(r.a.n_eff, r.b.n_eff);
  let verdict;
  if (flagged.length) verdict = el("p", { class: "note note-warn" }, `The shapes differ in: ${flagged.map(o => o.name).join("; ")}.`);
  else verdict = el("p", { class: "note" }, `No difference in shape found. This test can only see differences larger than the statistics of the samples (smallest N_eff = ${fmtInt(Math.round(lowest))}).`);
  parts.push(verdict);
  if (unweighted.length) parts.push(el("p", { class: "note" }, `The unweighted event distributions also differ in: ${unweighted.map(o => o.name).join("; ")}.`));
  for (const text of r.notes || []) parts.push(el("p", { class: "note" }, text));
  const row = (label, s) => el("tr", {}, el("td", {}, label), el("td", {}, s.name), el("td", {}, fmtInt(s.events)), el("td", {}, `${fmtNum(s.yield_)} ± ${fmtNum(s.yield_err, 2)}`), el("td", {}, fmtInt(Math.round(s.n_eff))));
  parts.push(el("table", { class: "compact" },
    el("thead", {}, el("tr", {}, ...["", "Sample", "Events", "Σw / attempts", "N_eff"].map(h => el("th", { scope: "col" }, h)))),
    el("tbody", {}, row("A", r.a), row("B", r.b))),
    el("p", { class: "setup" }, `Pull of Σw / attempts, A − B: ${r.pull >= 0 ? "+" : ""}${r.pull.toFixed(2)} σ.`));
  const pfmt = p => p === null ? "–" : (p <= 0.0011 ? "≤ 0.001" : p.toFixed(3));
  parts.push(el("table", { class: "compact" },
    el("thead", {}, el("tr", {}, ...["Observable", "Weighted p", "Unweighted p", "Mean A", "Mean B"].map(h => el("th", { scope: "col" }, h)))),
    el("tbody", {}, ...r.observables.map(o => el("tr", { class: o.p_weighted !== null && o.p_weighted <= r.alpha ? "flagged" : null },
      el("td", {}, o.name), el("td", {}, pfmt(o.p_weighted)), el("td", {}, pfmt(o.p)), el("td", {}, fmtNum(o.mean_a)), el("td", {}, fmtNum(o.mean_b)))))));
  if (r.plot) {
    const src = tokenUrl(`api/compare/${encodeURIComponent(r.plot)}`);
    parts.push(el("a", { href: src, target: "_blank", rel: "noopener", title: "Open full size" },
      el("img", { class: "plot", alt: "Overlays of the two samples, weighted and at unit area", src: src + "&t=" + Date.now() })));
  }
  return parts;
}
