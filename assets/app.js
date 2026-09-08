/* Quantum_AgroX — interactive engineering layer (per architecture doc) */
(function () {
  "use strict";

  /* ---------------- tiny utils ---------------- */
  var $ = function (sel) { return document.querySelector(sel); };
  var $$ = function (sel) { return Array.prototype.slice.call(document.querySelectorAll(sel)); };
  function el(tag, cls, html) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (html != null) e.innerHTML = html;
    return e;
  }
  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }
  function num(x) { var v = parseFloat(x); return isFinite(v) ? v : 0; }
  function round(x, d) { var p = Math.pow(10, d == null ? 2 : d); return Math.round(num(x) * p) / p; }
  function clamp01(x) { return Math.max(0, Math.min(1, x)); }
  function mean(arr) { return arr.reduce(function (a, b) { return a + num(b); }, 0) / arr.length; }
  function pctLink(q, label) {
    var s = "https://scholar.google.com/scholar?q=" + encodeURIComponent(q);
    var p = "https://pubmed.ncbi.nlm.nih.gov/?term=" + encodeURIComponent(q);
    return '<a class="lnk" target="_blank" rel="noopener" href="' + s + '">' + (label || "Scholar") + "</a> &nbsp;·&nbsp; " +
           '<a class="lnk" target="_blank" rel="noopener" href="' + p + '">PubMed</a>';
  }
  function svgIcon(path) {
    return '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="' + path + '"/></svg>';
  }
  function nowStamp() { return new Date().toISOString().replace("T", " ").slice(0, 19); }

  var D = null, M = null;         /* dashboard data + model layer */
  var APP = {
    RDKit: null, rdkitReady: false, tab: "home",
    pipe: 0, trail: [], done: {}, targetKey: null, lastLab: null, ledger: [], pendingPipe: null
  };
  window.APP = APP;

  function logTrail(msg) {
    APP.trail.push({ t: nowStamp(), m: msg });
    var box = $("#trailbox");
    if (box) {
      var row = el("div", "", "<b>" + esc(nowStamp()) + "</b> — " + esc(msg));
      box.appendChild(row);
      box.scrollTop = box.scrollHeight;
    }
  }
  function markDone(id) { APP.done[id] = true; }

  /* ---------------- RDKit loader ---------------- */
  function bootRDKit() {
    var factory = window.initRDKitModule || window.initRDKit;
    if (!factory) return Promise.reject(new Error("RDKit script not found — check your connection"));
    return factory().then(function (m) { APP.RDKit = m; window.RDKIT = m; APP.rdkitReady = true; window.setTimeout(renderFloatMols, 0); return m; });
  }
  function loadRDKit() {
    return new Promise(function (resolve, reject) {
      var factory = window.initRDKitModule || window.initRDKit;
      if (factory) { bootRDKit().then(resolve)["catch"](reject); return; }
      if (window.RDKit && window.RDKit.get_mol) { APP.RDKit = window.RDKit; window.RDKIT = window.RDKit; return resolve(); }
      var s = document.createElement("script");
      s.src = "vendor/RDKit_minimal.js";
      s.onload = function () { bootRDKit().then(resolve)["catch"](reject); };
      s.onerror = function () { reject(new Error("RDKit script not found — check your connection")); };
      document.head.appendChild(s);
    });
  }

  /* ---------------- molecular helpers (RDKit) ---------------- */
  function rdMol(smiles, name) {
    var mol = APP.RDKit.get_mol(smiles);
    if (!mol || !mol.is_valid()) { if (mol) mol.delete(); throw new Error("Invalid SMILES for " + name); }
    return mol;
  }
  var FEAT_LABEL = { mw: "MW", clogp: "cLogP", tpsa: "TPSA", hba: "HBA", hbd: "HBD", rotb: "RotB", rings: "Rings", arom: "Aromatic", heavy: "Heavy", fsp3: "Fsp3" };
  function molFeatures(mol) {
    var d = JSON.parse(mol.get_descriptors());
    return { mw: num(d.exactmw), clogp: num(d.CrippenClogP), tpsa: num(d.tpsa),
      hba: num(d.lipinskiHBA), hbd: num(d.lipinskiHBD), rotb: num(d.NumRotatableBonds),
      rings: num(d.NumRings), arom: num(d.NumAromaticRings), heavy: num(d.NumHeavyAtoms),
      fsp3: num(d.FractionCSP3), atoms: num(d.NumAtoms), het: num(d.NumHeteroatoms) };
  }
  function molSvg(mol) { return mol.get_svg(320, 260) || ""; }
  function fingerprint(mol) {
    var u = mol.get_morgan_fp_as_uint8array ? mol.get_morgan_fp_as_uint8array() : null;
    if (!u) return null;
    var b = new Uint8Array(u.length);
    for (var i = 0; i < u.length; i++) b[i] = u[i] > 0 ? 1 : 0;
    return b;
  }
  function tanimoto(a, b) {
    if (!a || !b) return null;
    var n = Math.min(a.length, b.length), inter = 0, union = 0;
    for (var i = 0; i < n; i++) { if (a[i] && b[i]) inter++; if (a[i] || b[i]) union++; }
    return union ? inter / union : 0;
  }

  /* ---------------- ridge (training-set refit) ---------------- */
  var FEATURES = [];
  function gaussSolve(A, b) {
    var n = A.length, aug = A.map(function (row, i) { return row.concat([b[i]]); });
    for (var c = 0; c < n; c++) {
      var mx = c;
      for (var r = c + 1; r < n; r++) if (Math.abs(aug[r][c]) > Math.abs(aug[mx][c])) mx = r;
      if (mx !== c) { var t = aug[c]; aug[c] = aug[mx]; aug[mx] = t; }
      var pv = aug[c][c];
      if (Math.abs(pv) < 1e-12) continue;
      for (var r2 = c + 1; r2 < n; r2++) {
        var f0 = aug[r2][c] / pv;
        for (var k = c; k <= n; k++) aug[r2][k] -= f0 * aug[c][k];
      }
    }
    var x = new Array(n).fill(0);
    for (var i = n - 1; i >= 0; i--) {
      var s = aug[i][n];
      for (var j = i + 1; j < n; j++) s -= aug[i][j] * x[j];
      x[i] = Math.abs(aug[i][i]) > 1e-12 ? s / aug[i][i] : 0;
    }
    return x;
  }
  function fitRidge(X, y, lam, st) {
    var n = X.length, f = X[0].length;
    var Xt = X.map(function (r) { return [1].concat(r.map(function (v, j) { return st[j].sd > 1e-9 ? (v - st[j].mean) / st[j].sd : 0; })); });
    var A = [], b = new Array(f + 1).fill(0), i, j, k;
    for (i = 0; i <= f; i++) {
      A.push(new Array(f + 1).fill(0));
      for (j = 0; j <= f; j++) {
        var s = 0;
        for (k = 0; k < n; k++) s += Xt[k][i] * Xt[k][j];
        A[i][j] = s + (i > 0 ? lam : 0);
      }
    }
    for (i = 0; i <= f; i++) for (k = 0; k < n; k++) b[i] += Xt[k][i] * y[k];
    return gaussSolve(A, b);
  }
  function predictRidge(w, x, st) {
    var p = w[0];
    for (var j = 0; j < x.length; j++) {
      var z = st[j].sd > 1e-9 ? (x[j] - st[j].mean) / st[j].sd : 0;
      p += w[j + 1] * z;
    }
    return p;
  }
  function featStats(list, feats) {
    return feats.map(function (ft) {
      var vals = list.map(function (r) { return num(r[ft]); });
      var m = mean(vals);
      var sd = Math.sqrt(mean(vals.map(function (v) { var d = v - m; return d * d; })));
      return { name: ft, mean: m, sd: sd || 1 };
    });
  }
  function statsFromMatrix(X) {
    var f = X[0].length, n = X.length, out = [];
    for (var j = 0; j < f; j++) {
      var m = 0, sd;
      for (var i = 0; i < n; i++) m += X[i][j];
      m /= n;
      var s2 = 0;
      for (var k = 0; k < n; k++) { var d = X[k][j] - m; s2 += d * d; }
      sd = Math.sqrt(s2 / n) || 1;
      out.push({ mean: m, sd: sd });
    }
    return out;
  }
  function fitAffinityModel() {
    var tr = M.train;
    var st = featStats(tr, FEATURES);
    var X = tr.map(function (r) { return FEATURES.map(function (ft) { return num(r[ft]); }); });
    var y = tr.map(function (r) { return num(r.dG); });
    return { w: fitRidge(X, y, 1.0, st), stats: st };
  }

  /* ---------------- target resolution ---------------- */
  function resolveTarget(q) {
    q = (q || "").trim();
    if (!q) return null;
    var lq = q.toLowerCase();
    var T = M.targets;
    var t = T.find(function (x) { return x.key === lq; });
    if (!t) t = T.find(function (x) { return x.name.toLowerCase() === lq; });
    if (!t) t = T.find(function (x) { return lq.indexOf(x.name.toLowerCase()) >= 0; });
    if (!t) t = T.find(function (x) { return lq.indexOf(x.protein.toLowerCase()) >= 0; });
    if (!t) t = T.find(function (x) { return x.organism.toLowerCase().indexOf(lq) >= 0; });
    if (!t) t = T.find(function (x) { return (x.focus || "").toLowerCase().indexOf(lq) >= 0; });
    return t || null;
  }
  function calibrantFor(t) {
    if (t && t.anchorKey && t.anchorFeat) return t;
    var g = M.targets.find(function (x) { return x.key === "ache_insect"; });
    return g;
  }
  function affinityForTarget(t, featVec, model) {
    var base = calibrantFor(t);
    var anchorPred = predictRidge(model.w, base.anchorFeat, model.stats);
    var offset = base.anchorDG - anchorPred;
    return predictRidge(model.w, featVec, model.stats) + offset;
  }

  /* ---------------- Molecule Lab ---------------- */
  function doseBand(aff) {
    var potency = Math.exp(-(aff + 5));
    var pmax = Math.exp(5.0), pmin = Math.exp(0.6);
    var p = clamp01((potency - pmin) / (pmax - pmin));
    var idx = 0.85 - 0.55 * p;
    var mg = 5 + idx * 75;
    var band = mg <= 15 ? "Low loading (5–15 mg/mL)" : mg <= 35 ? "Moderate (15–35 mg/mL)" : mg <= 60 ? "Upper band (35–60 mg/mL)" : "Maximum tested (60–80 mg/mL)";
    return { mg: round(mg, 0), band: band };
  }
  function labAnalyze(outSel) {
    var out = $(outSel || "#lab-out");
    if (!out) return null;
    out.innerHTML = "<div class='dim'>Working…</div>";
    var smiA = $("#smi-a").value.trim();
    var smiB = $("#smi-b").value.trim();
    var targetName = $("#target").value.trim();
    if (!smiA || !smiB) { out.innerHTML = "<div class='warn'>Enter two SMILES strings.</div>"; return null; }
    var t = resolveTarget(targetName) || null;
    APP.targetKey = t ? t.key : null;
    if ($("#target-hint")) {
      $("#target-hint").innerHTML = t
        ? "Target resolved: <b>" + esc(t.name) + "</b> — " + esc(t.organism) + (t.accession && t.accession !== "—" ? " (accession " + esc(t.accession) + ")" : "") + (t.anchorDG ? " · anchor " + esc(t.anchor) + " " + t.anchorDG + " kcal/mol" : " · no anchor bundled — generic calibration") : "";
    }
    var molA, molB, svgA = "", svgB = "";
    try {
      molA = rdMol(smiA, "molecule A");
      molB = rdMol(smiB, "molecule B");
      svgA = molSvg(molA); svgB = molSvg(molB);
    } catch (e) { out.innerHTML = "<div class='warn'>" + esc(e.message) + "</div>"; return null; }
    var fA = molFeatures(molA), fB = molFeatures(molB);
    var tani = tanimoto(fingerprint(molA), fingerprint(molB));
    try { molA.delete(); molB.delete(); } catch (e) {}

    var model = fitAffinityModel();
    var xA = FEATURES.map(function (ft) { return fA[ft]; });
    var xB = FEATURES.map(function (ft) { return fB[ft]; });
    var xC = xA.map(function (v, i) { return (v + xB[i]) / 2; });
    var affA = affinityForTarget(t, xA, model);
    var affB = affinityForTarget(t, xB, model);
    var affC = affinityForTarget(t, xC, model);
    var lowConf = !(t && t.anchorKey && t.anchorFeat);

    var complement = tani == null ? 0.5 : 1 - tani;
    var logpMatch = 1 - Math.abs(fA.clogp - fB.clogp) / (Math.max(1, Math.abs(fA.clogp)) + Math.max(1, Math.abs(fB.clogp)));
    var sizeMix = 1 - Math.abs(fA.mw - fB.mw) / (fA.mw + fB.mw || 1);
    var affGain = clamp01((Math.min(affA, affB) - affC) / 2.0);
    var synergy = 0.35 * complement + 0.2 * logpMatch + 0.2 * sizeMix + 0.25 * affGain;
    function synLabel(s) {
      if (s >= 0.72) return { t: "Very high", c: "good" };
      if (s >= 0.55) return { t: "High", c: "good" };
      if (s >= 0.35) return { t: "Moderate", c: "warn" };
      return { t: "Low — likely additive at best", c: "bad" };
    }
    var sl = synLabel(synergy);
    var dbA = doseBand(affA), dbB = doseBand(affB), dbC = doseBand(affC);
    var tLabel = t ? esc(t.name) + " · " + esc(t.organism) : esc(targetName || "custom receptor") + " (generic AChE calibration)";

    function featRows(f) {
      return [
        ["Molecular weight (g/mol)", round(f.mw, 1)],
        ["cLogP (lipophilicity)", round(f.clogp, 2)],
        ["TPSA (Å²)", round(f.tpsa, 1)],
        ["H-bond acceptors / donors", f.hba + " / " + f.hbd],
        ["Rotatable bonds / rings (aromatic)", f.rotb + " / " + f.rings + " (" + f.arom + ")"],
        ["Heavy atoms / heteroatoms", f.heavy + " / " + f.het],
        ["Fraction Csp3", round(f.fsp3, 2)]
      ].map(function (r) { return "<tr><td class='dim'>" + r[0] + "</td><td>" + r[1] + "</td></tr>"; }).join("");
    }
    function molPanel(tag, svg, f) {
      return "<div class='panel'><div class='panel-h'>Molecule " + tag + "</div>" +
        (svg ? "<div class='svgbox'>" + svg + "</div>" : "") +
        "<table class='kv'><tbody>" + featRows(f) + "</tbody></table></div>";
    }

    var html = "";
    html += "<div class='grid c2'>" + molPanel("A", svgA, fA) + molPanel("B", svgB, fB) + "</div>";
    html += "<div class='grid c2'>";
    html += "<div class='panel'><div class='panel-h'>Estimated binding affinity — " + tLabel + "</div>" +
      "<p class='sml dim'>Surrogate ridge QSAR on 20 reference ligands, calibrated to the target anchor (" +
      (t && t.anchor ? esc(t.anchor) + ", " + t.anchorDG + " kcal/mol" : "generic insect AChE scaffold") + ")." +
      (lowConf ? " <span class='warn-bc'>Low confidence — this target has no bundled anchor.</span>" : "") + "</p>" +
      "<table class='kv'><tbody>" +
      "<tr><td>Molecule A</td><td class='num big'>" + round(affA, 2) + " kcal/mol</td></tr>" +
      "<tr><td>Molecule B</td><td class='num big'>" + round(affB, 2) + " kcal/mol</td></tr>" +
      "<tr><td>1:1 combination (mean descriptor)</td><td class='num big'>" + round(affC, 2) + " kcal/mol</td></tr>" +
      "</tbody></table></div>";
    html += "<div class='panel'><div class='panel-h'>Synergy index" + "</div>" +
      "<p class='sml dim'>Loewe-inspired heuristic: structural complementarity (1 − Tanimoto), LogP balance, size mixing, predicted affinity gain of the combination.</p>" +
      "<div class='syn'><span style='width:" + Math.round(synergy * 100) + "%' class='syn-'></span></div>" +
      "<div class='syn-label " + sl.c + "'>" + sl.t + " · index " + round(synergy, 2) + "</div>" +
      "<table class='kv sml'><tbody>" +
      "<tr><td>Structural complementarity (1 − Tanimoto)</td><td class='num'>" + (tani == null ? "n/a" : round(complement, 2)) + "</td></tr>" +
      "<tr><td>Tanimoto similarity</td><td class='num'>" + (tani == null ? "n/a" : round(tani, 2)) + "</td></tr>" +
      "<tr><td>LogP complementarity</td><td class='num'>" + round(logpMatch, 2) + "</td></tr>" +
      "<tr><td>Size mixing balance</td><td class='num'>" + round(sizeMix, 2) + "</td></tr>" +
      "<tr><td>Affinity gain of combination</td><td class='num'>" + round(affGain, 2) + "</td></tr></tbody></table></div></div>";
    html += "<div class='panel'><div class='panel-h'>Dose guidance (" + esc(targetName || "custom receptor") + ")</div>" +
      "<p class='sml dim'>Mapped to the study's tested topical window <b>5–80 mg/mL</b>.</p>" +
      "<table class='kv'><tbody>" +
      "<tr><td>Molecule A</td><td class='num'>" + dbA.mg + " mg/mL — " + dbA.band + "</td></tr>" +
      "<tr><td>Molecule B</td><td class='num'>" + dbB.mg + " mg/mL — " + dbB.band + "</td></tr>" +
      "<tr><td>Combination</td><td class='num'>" + dbC.mg + " mg/mL — " + dbC.band + "</td></tr></tbody></table></div>";
    html += "<div class='panel warn-panel'><b>Caveat.</b> " + esc(M.caveat) + "</div>";

    out.innerHTML = html;
    APP.lastLab = { smiA: smiA, smiB: smiB, target: t ? (t.name + " (" + t.organism + ")") : targetName, targetKey: (t || {}).key,
      affA: affA, affB: affB, affC: affC, synergy: synergy, dose: dbC.mg, fA: fA, fB: fB, lowConf: lowConf };
    APP.ledger.unshift({
      t: nowStamp(), smiA: smiA, smiB: smiB,
      target: t ? (t.name + " (" + t.organism + ")") : targetName, targetKey: (t || {}).key,
      affA: round(affA, 2), affB: round(affB, 2), affC: round(affC, 2), synergy: round(synergy, 2),
      syn: sl.t, dose: dbC.mg + " mg/mL", lowConf: lowConf, mwA: round(fA.mw, 1), mwB: round(fB.mw, 1)
    });
    if (APP.ledger.length > 25) APP.ledger.length = 25;
    logTrail("AgroDockX/AgroSynergyX: analyzed " + (t ? t.name : targetName) + " — affinity " + round(affC, 2) + " kcal/mol, synergy " + round(synergy, 2));
    return APP.lastLab;
  }
  var PHYTO_SMILES = [
    /* Benchmark phytochemicals (Commiphora swynnertonii GC-MS set).
       Ordered most-specific-first: "7alpha-methylcholesterol" must hit
       "methylcholesterol" before the generic "cholest" scaffold.
       Steroid rows marked scaffold carry a representative scaffold —
       curate exact stereochemistry (PubChem / Suppl. File 1) before
       production docking. All SMILES RDKit-validated. */
    ["methylcholesterol", "CC(C)CCCC(C)C1CCC2C3CCC4CC(O)CCC4(C)C3CCC12C"],
    ["chlorfenvinphos", "CCOP(=O)(OCC)OC(=CCl)C1=CC=C(Cl)C=C1Cl"],
    ["hexadecanoic", "CCCCCCCCCCCCCCCC(=O)O"],
    ["palmitic", "CCCCCCCCCCCCCCCC(=O)O"],
    ["pentadecanoate", "CCCCCCCCCCCCCCC(=O)OCC"],
    ["hentriacontanol", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCO"],
    ["phenylaniline", "NC1=CC=CC=C1C1=CC=CC=C1"],
    ["aminobiphenyl", "NC1=CC=CC=C1C1=CC=CC=C1"],
    ["caryophyllene", "CC1=CCCC(=CCCC1)C"],
    ["humulene", "CC1=CCCC(=CCCC1)C"],
    ["cedrol", "CC1CCC2C1C1CCCC1(C)C2O"],
    ["farnesol", "CC(=CCCC(=CCCC(=CCO)C)C)C"],
    ["valencene", "CC1=CCC2CCC(C)C2C1=C(C)C"],
    ["copaene", "CC1=CCC2C3CC(C2C1)C3(C)C"],
    ["curzerene", "C=CC1=C(C)C2=C(O1)C(C)CC2"],
    ["dendrolasin", "CC(=CCCC(=CCOC1=CC=CO1)C)C"],
    ["torreyol", "CC1CCC(C)C2CCC(C)(O)C2C1"],
    ["phytol", "CC(C)CCCC(C)CCCC(C)CCCC(=CCO)C"],
    ["lathosterol", "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C"],
    ["ergost", "CC(C)C(C)CCC(C)C1CCC2C3CC(O)C(O)CC3CCC12C"],
    ["androst", "CC12CCC3C4CCCC4CCC3C1CCC2O"],
    ["spirost", "CC(C)CCCC(C)C1CCC2C3CCC4CC(O)CCC4(C)C3CCC12C"],
    ["furanmethanol", "CC1CCC(C)(C2CCC(C)=CC2)OC1"],
    ["cholest", "CC(C)CCCC(C)C1CCC2C3CCC4CC(O)CCC4(C)C3CCC12C"],
    ["pinene", "CC1=CCC2CC1C2(C)C"],
    ["limonene", "CC1=CCC(CC1)=C(C)C"],
    ["linalool", "CC(C)=CCCC(C)(O)C=C"],
    ["geraniol", "CC(C)=CCCC(C)=CCO"],
    ["thymol", "CC(C)c1cc(C)ccc1O"],
    ["carvacrol", "CC(C)c1ccc(O)cc1C"],
    ["eugenol", "COc1cc(CC=C)ccc1O"],
    ["menthol", "CC(C)C1CCC(C)CC1O"],
    ["camphor", "CC1(C)C2CCC1(C)C(=O)C2"],
    ["squalene", "CC(=CCCC(=CCCC(=CCCC=C(C)C)C)C)C"],
    ["quercetin", "O=C1c2c(O)cc(O)cc2OC(c2cc(O)c(O)c(O)c2)=C1O"],
    ["kaempferol", "O=C1c2c(O)cc(O)cc2OC(c2ccc(O)cc2)=C1O"]
  ];
  /* Resolve a phytochemical *name* to its SMILES entry.
     Longest-key-first substring match, so "7alpha-methylcholesterol"
     resolves to "methylcholesterol" rather than the generic "cholest"
     scaffold. Returns {key, smi} or null (→ manual SMILES entry). */
  function phytoSmilesFor(name) {
    var lname = String(name || "").toLowerCase();
    if (!lname) return null;
    var best = null;
    PHYTO_SMILES.forEach(function (p) {
      if (lname.indexOf(p[0]) >= 0 && (!best || p[0].length > best[0].length)) best = p;
    });
    return best ? { key: best[0], smi: best[1] } : null;
  }
  /* Default pair prefers compounds WITH a curated SMILES, so handoffs work
     even if AgroDockX is opened before AgroPhytoX. */
  function phytoDefaultPair() {
    var names = (D && D.phytochemicals ? D.phytochemicals : []).map(function (r) { return r.compound; });
    var mapped = names.filter(function (n) { return phytoSmilesFor(n); });
    var dA = mapped[0] || names[0] || "";
    var rest = mapped.filter(function (n) { return n !== dA; });
    var dB = rest[0] || mapped[0] || names.filter(function (n) { return n !== dA; })[0] || "";
    return [dA, dB];
  }
  function useTopPhytochemicals() {
    var filled = [];
    if (D && D.phytochemicals) {
      var top = D.phytochemicals.slice().sort(function (a, b) { return b.area - a.area; });
      for (var i = 0; i < top.length && filled.length < 2; i++) {
        var lname = String(top[i].compound || "").toLowerCase();
        var hit = PHYTO_SMILES.find(function (p) { return lname.indexOf(p[0]) >= 0; });
        if (hit) filled.push({ name: top[i].compound, smi: hit[1] });
      }
    }
    var mintop = { name: "thymol (monoterpene)", smi: "CC(C)c1cc(C)ccc1O" };
    var mineug = { name: "eugenol (phenylpropanoid)", smi: "COc1cc(CC=C)ccc1O" };
    var pair = [];
    if (filled.length >= 2) pair = [filled[0], filled[1]];
    else if (filled.length === 1) pair = [filled[0], mineug];
    else pair = [mintop, mineug];
    $("#smi-a").value = pair[0].smi;
    $("#smi-b").value = pair[1].smi;
    $("#target").value = "Acetylcholinesterase (Rhipicephalus microplus)";
    $("#target").dispatchEvent(new Event("input"));
    var hint = el("div", "sml dim");
    hint.style.cssText = "margin-top:8px";
    var src = filled.length >= 2 ? "top two SMILES-mapped phytochemicals from the library (" + filled[0].name + " + " + filled[1].name + ")" :
      (filled.length === 1 ? "top SMILES-mapped phytochemical from the library (" + filled[0].name + ") + eugenol" : "curated library analogues for the two top compound classes");
    hint.textContent = "Filled with " + src + ". Edit freely — any valid SMILES pair works.";
    $("#smi-a").focus();
    return hint;
  }

  /* ---------------- Home render ---------------- */
  function homeRender() {
    $("#b-records").textContent = D ? D.stats.records : "…";
    $("#b-phyt").textContent = D ? D.stats.phytochemicals : "…";
    $("#b-dock").textContent = D ? D.stats.docking_hits : "…";
    $("#b-val").textContent = D ? D.stats.validation_complexes : "…";
    $("#b-lc50").textContent = D ? D.stats.lc50_entries : "…";
    renderModgrid();
    renderPstrip();
    renderRefTable();
    renderPhases();
  }
  function renderModgrid() {
    var wrap = $("#modgrid");
    if (!wrap) return;
    wrap.innerHTML = "";
    M.modules.forEach(function (mod, i) {
      var c = el("div", "modcard");
      c.innerHTML = "<div class='mchip'>" + String(i + 1).padStart(2, "0") + "</div>" +
        "<div class='mname'>" + esc(mod.name) + "</div>" +
        "<div class='mfunc'>" + esc(mod.func) + "</div>" +
        "<div class='mdesc'>" + esc(mod.desc) + "</div>" +
        "<div class='mlaunch'>Open module →</div>";
      c.addEventListener("click", function () { openPipelineAt(i); });
      wrap.appendChild(c);
    });
  }
  function renderPstrip() {
    var wrap = $("#pstrip");
    if (!wrap) return;
    wrap.innerHTML = "";
    M.modules.forEach(function (mod, i) {
      var row = el("div", "prow");
      row.innerHTML = "<div class='pdot'>" + (i + 1) + "</div><div>" + esc(mod.name.replace("™", "")) + "</div>" +
        (i < M.modules.length - 1 ? "<div class='pbar'></div>" : "");
      wrap.appendChild(row);
    });
  }
  function renderRefTable() {
    var t = $("#ref-table");
    if (!t || !M.refparams) return;
    var rows = "";
    Object.keys(M.refparams).forEach(function (k) {
      var r = M.refparams[k];
      rows += "<tr><td><b>" + esc(k) + "</b></td><td>accession " + esc(r.accession) + " · template " + esc(r.template) +
        "</td><td class='num'>centre " + esc(r.center) + "<br>box " + esc(r.box) + " · exhaust. " + esc(r.exhaust) + "</td></tr>";
    });
    t.innerHTML = "<tbody>" + rows + "</tbody>";
  }
  function renderPhases() {
    var w = $("#phases");
    if (!w) return;
    w.innerHTML = "";
    M.phases.forEach(function (p) {
      w.appendChild(el("div", "sml", "• " + esc(p)));
    });
  }

  /* ---------------- Pipeline ---------------- */
  function openPipelineAt(i) {
    showTab("pipeline");
    APP.pipe = Math.max(0, Math.min(M.modules.length - 1, i));
    renderStepper();
    renderStage();
    window.scrollTo({ top: 0, behavior: "smooth" });
  }
  function renderStepper() {
    var st = $("#stepper");
    if (!st) return;
    st.innerHTML = "";
    M.modules.forEach(function (mod, i) {
      var s = el("div", "step" + (i === APP.pipe ? " active" : "") + (APP.done[mod.id] ? " done" : ""));
      s.innerHTML = "<div class='sn'>" + (APP.done[mod.id] ? "✓" : (i + 1)) + "</div><div>" + esc(mod.name.replace("™", "")) + "</div>";
      s.addEventListener("click", function () { APP.pipe = i; renderStepper(); renderStage(); });
      st.appendChild(s);
    });
  }
  function stageNav() {
    var nxt = (APP.pipe < M.modules.length - 1)
      ? '<button class="primary" id="pipe-next">Next: ' + esc(M.modules[APP.pipe + 1].name.replace("™", "")) + " →</button>" : "";
    var prev = (APP.pipe > 0) ? '<button class="ghost" id="pipe-prev">← Previous</button>' : "";
    return "<div class='navbtns'>" + prev + nxt + "</div>";
  }
  function wireNav() {
    var nxt = $("#pipe-next"), prv = $("#pipe-prev");
    if (nxt) nxt.addEventListener("click", function () { APP.pipe++; renderStepper(); renderStage(); });
    if (prv) prv.addEventListener("click", function () { APP.pipe--; renderStepper(); renderStage(); });
  }

  function renderStage() {
    var body = $("#pipe-body");
    if (!body) return;
    var mod = M.modules[APP.pipe];
    body.innerHTML = "";
    var card = el("div", "panel");
    card.innerHTML = "<div class='panel-h'>" + esc(mod.name) + " <span class='tag mid'>" + esc(mod.func) + "</span></div>" +
      "<p class='sml dim'>" + esc(mod.desc) + "</p>";
    body.appendChild(card);
    var content = el("div", "");
    body.appendChild(content);

    var fns = {};
    fns.datahub = function () {
      logTrail("AgroDataHub: benchmark dataset loaded");
      markDone("datahub");
      var inv = el("div", "panel");
      inv.innerHTML = "<div class='panel-h'>Dataset inventory (benchmark bundle)</div>";
      var rows = (D ? D.datasets : []).map(function (ds) {
        return "<tr><td>" + esc(ds.file) + "</td><td class='num'>" + ds.rows + " rows</td><td>" + esc(ds.type) + "</td></tr>";
      }).join("");
      inv.innerHTML += "<table class='data'><thead><tr><th>File</th><th>Records</th><th>Contents</th></tr></thead><tbody>" + rows + "</tbody></table>";
      inv.innerHTML += "<p class='sml dim' style='margin-top:10px'>Import standards: GC-MS peak lists, SDF/SMILES, FASTA (targets), VCF (populations), and bioassay/docking tables are accepted. Every ingestion writes a provenance record for AgroReport.</p>";
      content.appendChild(inv);
    };
    fns.phytox = function () {
      logTrail("AgroPhytoX: phytochemical library browsed");
      markDone("phytox");
      var p = el("div", "panel");
      p.innerHTML = "<div class='panel-h'>Phytochemical library (" + (D.stats ? D.stats.phytochemicals : 0) + " profiled)</div>" +
        "<p class='sml dim'>Pick <b>two phytochemicals by name</b> — their SMILES travel with the session into AgroDockX.</p>" +
        "<div class='grid c2 mlgens'>" +
        "<div><label>Phytochemical A — name</label><select id='phy-a'></select></div>" +
        "<div><label>Phytochemical B — name</label><select id='phy-b'></select></div></div>" +
        "<div style='display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-top:4px'>" +
        "<button class='primary' id='phy-send'>Send pair to AgroDockX →</button>" +
        "<span class='sml dim' id='phy-pair-note'></span></div>" +
        "<div style='margin-top:14px'><input id='phy-q' type='text' placeholder='Filter name / class / part / solvent' style='width:100%;max-width:460px;background:var(--panel2);color:var(--txt);border:1px solid var(--line2);border-radius:8px;padding:9px'>" +
        "<div id='phy-rows' style='margin-top:10px'></div></div>";
      content.appendChild(p);
      var names = (D.phytochemicals || []).map(function (r) { return r.compound; });
      var selA = $("#phy-a"), selB = $("#phy-b");
      selA.innerHTML = names.map(function (n) { return "<option>" + esc(n) + "</option>"; }).join("");
      selB.innerHTML = names.map(function (n) { return "<option>" + esc(n) + "</option>"; }).join("");
      /* Defaults prefer compounds that HAVE a curated SMILES, so the
         AgroDockX handoff visibly works out of the box. */
      var dp = phytoDefaultPair();
      if (!(APP.phytoA && names.indexOf(APP.phytoA) >= 0)) APP.phytoA = dp[0];
      if (!(APP.phytoB && names.indexOf(APP.phytoB) >= 0) || APP.phytoB === APP.phytoA) {
        APP.phytoB = (dp[1] && dp[1] !== APP.phytoA) ? dp[1] : dp[0];
      }
      selA.value = APP.phytoA; selB.value = APP.phytoB;
      function pairNote() {
        APP.phytoA = selA.value; APP.phytoB = selB.value;
        var hA = phytoSmilesFor(APP.phytoA), hB = phytoSmilesFor(APP.phytoB);
        var note = $("#phy-pair-note");
        note.innerHTML = "A: " + (hA ? "SMILES ✓ (" + esc(hA.key) + ")" : "no curated SMILES — paste manually in AgroDockX") +
          " · B: " + (hB ? "SMILES ✓ (" + esc(hB.key) + ")" : "no curated SMILES — paste manually in AgroDockX") +
          ((!hA || !hB) ? " · tip: Phytol, Copaene, Curzerene, Valencene, Lathosterol… carry curated SMILES" : "");
        note.style.color = (hA && hB) ? "var(--em)" : "var(--amber)";
        logTrail("AgroPhytoX: pair selected — " + APP.phytoA + " + " + APP.phytoB);
      }
      selA.addEventListener("change", pairNote);
      selB.addEventListener("change", pairNote);
      pairNote();
      $("#phy-send").addEventListener("click", function () {
        pairNote();
        var i = M.modules.findIndex(function (m) { return m.id === "dockx"; });
        openPipelineAt(i >= 0 ? i : APP.pipe);
      });
      var draw = function () {
        var q = (($("#phy-q") || {}).value || "").trim().toLowerCase();
        var list = (D.phytochemicals || []).slice(0, 40).filter(function (r) {
          if (!q) return true;
          return (r.compound + r.cls + r.plant_part + r.solvent).toLowerCase().indexOf(q) >= 0;
        });
        var rows = list.map(function (r) {
          return "<tr><td>" + esc(r.compound) + "</td><td>" + esc(r.cls) + "</td><td>" + esc(r.plant_part) + "</td><td>" + esc(r.solvent) + "</td><td class='num'>" + round(r.area, 1) + "%</td></tr>";
        }).join("") || "<tr><td colspan='5' class='dim'>No matches.</td></tr>";
        $("#phy-rows").innerHTML = "<table class='data'><thead><tr><th>Compound</th><th>Class</th><th>Plant part</th><th>Solvent</th><th>Relative area</th></tr></thead><tbody>" + rows + "</tbody></table>";
      };
      draw();
      setTimeout(function () { var qi = $("#phy-q"); if (qi) qi.addEventListener("input", draw); }, 0);
    };
    fns.syntheticx = function () {
      markDone("syntheticx");
      logTrail("AgroSyntheticX: drawing workbench opened");
      var p = el("div", "panel");
      p.innerHTML = "<div class='panel-h'>AgroSyntheticX — draw two synthetic molecules → SMILES</div>" +
        "<p class='sml dim'>Sketch two candidates on the in-dashboard editor. Each is stored <b>separately</b> as <b>Synthetic A</b> / <b>Synthetic B</b>: inside the drawer press <b>Export to SMILES</b>, then <b>Use as A</b> / <b>Use as B</b>. AgroDockX can then take either one.</p>" +
        "<div class='grid c2 mlgens'>" +
        "<div><div class='panel-h sml'>Synthetic A</div><div class='kv sml'>SMILES: <span class='kbd' id='syn-smi-a'>—</span></div>" +
        "<div style='margin-top:6px'><button class='primary' id='syn-draw-a'>✏ Draw A</button></div></div>" +
        "<div><div class='panel-h sml'>Synthetic B</div><div class='kv sml'>SMILES: <span class='kbd' id='syn-smi-b'>—</span></div>" +
        "<div style='margin-top:6px'><button class='primary' id='syn-draw-b'>✏ Draw B</button></div></div></div>" +
        "<div style='display:flex;gap:10px;flex-wrap:wrap;margin:10px 0'>" +
        "<button class='ghost' id='syn-to-dock'>Send both to AgroDockX →</button></div>" +
        "<p class='sml dim' id='syn-note'></p>";
      content.appendChild(p);
      function showSlot(w) { var n = $("#syn-smi-" + w.toLowerCase()); if (n) n.textContent = APP["synth" + w] || "—"; }
      showSlot("A"); showSlot("B");
      function openDrawer() {
        if (window.QXSketch) QXSketch.open($("#smi-a").value, $("#smi-b").value);
        else $("#syn-note").textContent = "Drawer unavailable (sketcher.js not loaded).";
      }
      if (window.QXSketch) {
        /* Chain onto the drawer callback (same A/B → Molecule Lab effect
           as the lab's own handler) and store each slot separately. */
        QXSketch.onUse(function (which, smi) {
          var feat = null;
          try {
            if (APP.RDKit) {
              var mm = rdMol(smi, "synthetic");
              feat = molFeatures(mm);
              try { mm.delete(); } catch (e2) {}
            }
          } catch (e) { feat = null; }
          APP.synthSmi = smi; APP.synthFeat = feat;
          if (which === "A") { APP.synthA = smi; APP.synthFeatA = feat; $("#smi-a").value = smi; }
          else { APP.synthB = smi; APP.synthFeatB = feat; $("#smi-b").value = smi; }
          showSlot("A"); showSlot("B");
          $("#syn-note").textContent = "Captured from drawer ✓ — stored as Synthetic " + which + " (+ Molecule Lab " + which + ").";
          logTrail("AgroSyntheticX: molecule drawn → " + smi + " (Synthetic " + which + ")");
        });
      }
      $("#syn-draw-a").addEventListener("click", openDrawer);
      $("#syn-draw-b").addEventListener("click", openDrawer);
      $("#syn-to-dock").addEventListener("click", function () {
        var i = M.modules.findIndex(function (m) { return m.id === "dockx"; });
        openPipelineAt(i >= 0 ? i : APP.pipe);
      });
    };
    fns.targetx = function () {
      logTrail("AgroTargetX: target library browsed");
      markDone("targetx");
      var p = el("div", "panel");
      p.innerHTML = "<div class='panel-h'>Target discovery — pick a pest/pathogen class and a receptor</div>" +
        "<div class='chips' id='tc-chips'></div><div id='tc-cards' class='tcards' style='margin-top:12px'></div>" +
        "<p class='sml dim' id='tc-sel'></p>";
      content.appendChild(p);
      var chipsBox = $("#tc-chips"), cards = $("#tc-cards"), sel = $("#tc-sel");
      var active = "mite";
      M.targetClasses.forEach(function (c) {
        var chip = el("span", "chip" + (c.id === active ? " chip-active" : ""), c.label);
        chip.style.cssText = c.id === active ? "border-color:var(--accent);color:var(--txt)" : "";
        chip.addEventListener("click", function () { active = c.id; draw(); });
        chipsBox.appendChild(chip);
      });
      function draw() {
        cards.innerHTML = "";
        M.targets.filter(function (t) { return t.cls === active; }).forEach(function (t) {
          var selT = APP.targetKey === t.key;
          var c = el("div", "tcard" + (selT ? " sel" : ""));
          c.innerHTML = "<div class='tn'>" + esc(t.name) + "</div><div class='to'>" + esc(t.organism) + "</div>" +
            "<div class='sml dim'>" + esc(t.focus) + (t.accession && t.accession !== "—" ? " · " + esc(t.accession) : "") + "</div>" +
            "<div class='sml' style='margin-top:6px'>" + (t.anchor ? "anchor " + esc(t.anchor) + " (" + t.anchorDG + " kcal/mol)" : "no anchor — generic calibration") + "</div>";
          c.addEventListener("click", function () {
            APP.targetKey = t.key;
            $("#target").value = t.name + " (" + t.organism + ")";
            $("#target").dispatchEvent(new Event("input"));
            sel.innerHTML = "Receptor selected: <b>" + esc(t.name) + "</b> (" + esc(t.organism) + "). Continue to AgroDockX to test molecules against it.";
            logTrail("AgroTargetX: target selected — " + t.name + " (" + t.organism + ")");
            draw();
          });
          cards.appendChild(c);
        });
      }
      draw();
    };
    fns.sitemap = function () {
      logTrail("AgroSiteMap: binding-site context loaded");
      markDone("sitemap");
      var p = el("div", "panel");
      p.innerHTML = "<div class='panel-h'>Binding-site &amp; druggability context</div>";
      var rt = "";
      Object.keys(M.refparams).forEach(function (k) {
        var r = M.refparams[k];
        rt += "<tr><td><b>" + esc(k) + "</b></td><td>" + esc(r.template) + "</td><td class='num'>" + esc(r.center) + "</td><td class='num'>" + esc(r.box) + "</td></tr>";
      });
      p.innerHTML += "<table class='data'><thead><tr><th>Model</th><th>Template</th><th>Docking centre</th><th>Box size</th></tr></thead><tbody>" + rt + "</tbody></table>";
      var panels = M.sitepanels.map(function (s) {
        return "<div class='vrow'><b>" + esc(s.site) + "</b><div class='sml dim'>" + esc(s.note) + "</div></div>";
      }).join("");
      p.innerHTML += "<div style='margin-top:14px'><div class='panel-h'>Pocket notes</div>" + panels +
        "<p class='sml dim'>Cryptic/induced-fit pockets require MD sampling (AgroMD) — static maps can miss them.</p></div>";
      content.appendChild(p);
    };
    fns.dockx = function () {
      markDone("dockx");
      logTrail("AgroDockX: docking-estimate workbench opened");
      var lab = el("div", "panel");
      lab.innerHTML = "<div class='panel-h'>AgroDockX — protein–ligand affinity estimate</div>" +
        "<p class='sml dim'>Enter any SMILES pair and the receptor (set in AgroTargetX or type free text) and run the surrogate. Stronger (more negative) affinity = better predicted pose. " +
        "To bring the SMILES of the two AgroPhytoX phytochemicals, press <b>Use AgroPhytoX pair</b>.</p>" +
        "<div class='grid c2 mlgens'>" +
        "<div><label>Molecule A</label><textarea id='dock-a' rows='2'></textarea></div>" +
        "<div><label>Molecule B</label><textarea id='dock-b' rows='2'></textarea></div></div>" +
        "<div style='display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-top:4px'>" +
        "<button class='primary' id='dock-run'>Estimate affinity</button>" +
        "<button class='ghost' id='dock-phyto'>Use AgroPhytoX pair</button>" +
        "<button class='ghost' id='dock-syn-a'>Take Synthetic A</button>" +
        "<button class='ghost' id='dock-syn-b'>Take Synthetic B</button>" +
        "<button class='ghost' id='dock-sync'>Use Molecule Lab pair</button>" +
        "<span class='sml dim' id='dock-note'>Output below.</span></div>" +
        "<div style='margin-top:10px;padding-top:10px;border-top:1px dashed var(--line)'>" +
        "<div class='panel-h sml'>Mixed pair — e.g. phytochemical as A + synthetic as B</div>" +
        "<div class='grid c2 mlgens'>" +
        "<div><label>Molecule A from</label><select id='mix-a'><option value='keep'>Keep current</option><option value='phyA'>AgroPhytoX A</option><option value='synA'>Synthetic A</option><option value='labA'>Molecule Lab A</option></select></div>" +
        "<div><label>Molecule B from</label><select id='mix-b'><option value='keep'>Keep current</option><option value='phyB'>AgroPhytoX B</option><option value='synB'>Synthetic B</option><option value='labB'>Molecule Lab B</option></select></div></div>" +
        "<div style='margin-top:6px'><button class='ghost' id='mix-apply'>Build mixed pair →</button></div></div>" +
        "<div id='dock-out' style='margin-top:12px'></div>";
      content.appendChild(lab);
      var a = $("#dock-a"), b = $("#dock-b");
      a.value = $("#smi-a").value; b.value = $("#smi-b").value;
      function pullPhytoPair() {
        var dp0 = phytoDefaultPair();
        if (!APP.phytoA) APP.phytoA = dp0[0];
        if (!APP.phytoB) APP.phytoB = dp0[1];
        var notes = [], filled = 0;
        [["phytoA", "A", a], ["phytoB", "B", b]].forEach(function (t3) {
          var nm = APP[t3[0]];
          if (!nm) { notes.push(t3[1] + ": no phytochemical chosen in AgroPhytoX"); return; }
          var hit = phytoSmilesFor(nm);
          if (hit) { t3[2].value = hit.smi; filled++; notes.push(t3[1] + ": " + nm + " ✓"); }
          else { notes.push(t3[1] + ": " + nm + " — no curated SMILES, paste manually (tip: Phytol, Copaene, Curzerene…)"); }
        });
        var noteEl = $("#dock-note");
        noteEl.textContent = filled + "/2 SMILES filled — " + notes.join(" · ");
        noteEl.style.color = filled === 2 ? "var(--em)" : "var(--amber)";
        logTrail("AgroDockX: SMILES pulled from AgroPhytoX pair (" + notes.join("; ") + ")");
      }
      if (APP.phytoA || APP.phytoB) pullPhytoPair();
      $("#dock-phyto").addEventListener("click", pullPhytoPair);
      function takeSynthetic(which) {
        var smi = which === "A" ? APP.synthA : APP.synthB;
        var noteEl = $("#dock-note");
        if (smi) {
          if (which === "A") a.value = smi; else b.value = smi;
          noteEl.textContent = "Synthetic " + which + " taken into Molecule " + which + " ✓";
          noteEl.style.color = "var(--em)";
          logTrail("AgroDockX: took Synthetic " + which + " from AgroSyntheticX");
        } else {
          noteEl.textContent = "No Synthetic " + which + " stored yet — draw it in AgroSyntheticX first.";
          noteEl.style.color = "var(--amber)";
        }
      }
      $("#dock-syn-a").addEventListener("click", function () { takeSynthetic("A"); });
      $("#dock-syn-b").addEventListener("click", function () { takeSynthetic("B"); });
      function sourceSmiles(kind) {
        var slot = kind.slice(-1);
        if (kind.indexOf("phy") === 0) {
          var dp1 = phytoDefaultPair();
          var nm = slot === "A" ? (APP.phytoA || dp1[0]) : (APP.phytoB || dp1[1]);
          var hit = nm && phytoSmilesFor(nm);
          return hit ? { smi: hit.smi, label: nm } : null;
        }
        if (kind.indexOf("syn") === 0) {
          var s = slot === "A" ? APP.synthA : APP.synthB;
          return s ? { smi: s, label: "Synthetic " + slot } : null;
        }
        var v = (slot === "A" ? $("#smi-a") : $("#smi-b")).value || "";
        v = v.trim();
        return v ? { smi: v, label: "Molecule Lab " + slot } : null;
      }
      $("#mix-apply").addEventListener("click", function () {
        var ka = $("#mix-a").value, kb = $("#mix-b").value;
        var notes = [], n = 0;
        [["A", ka, a], ["B", kb, b]].forEach(function (t3) {
          if (t3[1] === "keep") { notes.push(t3[0] + ": kept"); return; }
          var r = sourceSmiles(t3[1]);
          if (r) { t3[2].value = r.smi; n++; notes.push(t3[0] + ": " + r.label + " ✓"); }
          else { notes.push(t3[0] + ": source empty (pick it first in its module)"); }
        });
        var noteEl = $("#dock-note");
        noteEl.textContent = "Mixed pair — " + notes.join(" · ");
        noteEl.style.color = "var(--em)";
        logTrail("AgroDockX: mixed pair built (" + notes.join("; ") + ")");
      });
      $("#dock-run").addEventListener("click", function () {
        $("#smi-a").value = a.value; $("#smi-b").value = b.value;
        labAnalyze("#dock-out");
        if ($("#dock-out").querySelector) {
          content.classList = "";
        }
      });
      $("#dock-sync").addEventListener("click", function () {
        a.value = $("#smi-a").value; b.value = $("#smi-b").value;
        var noteEl = $("#dock-note");
        noteEl.textContent = "Molecule Lab pair copied into A/B ✓";
        noteEl.style.color = "var(--em)";
        logTrail("AgroDockX: Molecule Lab pair copied into A/B");
      });
    };
    fns.dosex = function () {
      markDone("dosex");
      logTrail("AgroDoseX: dose-response workbench opened");
      var p = el("div", "panel");
      p.innerHTML = "<div class='panel-h'>AgroDoseX — dose–response &amp; LC50/LC90</div>" +
        "<p class='sml dim'>A linear surrogate (mortality ∼ log concentration + exposure + species + extract) is fitted live on the 240 benchmark records. Predict mean mortality at any point; compare with reported probit LC50/LC90 below.</p>" +
        "<div class='grid c2' style='margin-top:8px'>" +
        "<div class='kv'>Concentration (mg/mL): <input id='dd-conc' type='number' value='20' min='1' step='1' style='width:120px;background:var(--panel2);color:var(--txt);border:1px solid var(--line2);border-radius:7px;padding:7px'></div>" +
        "<div class='kv'>Exposure (h): <input id='dd-exp' type='number' value='24' min='1' step='1' style='width:120px;background:var(--panel2);color:var(--txt);border:1px solid var(--line2);border-radius:7px;padding:7px'></div>" +
        "<div class='kv'>Extract: <select id='dd-ex'></select></div>" +
        "<div class='kv'>Species: <select id='dd-sp'><option>Rhipicephalus microplus</option><option>Rhipicephalus decoloratus</option></select></div>" +
        "</div><div style='margin-top:8px'><button class='primary' id='dd-run'>Predict mortality</button></div>" +
        "<div id='dd-out' style='margin-top:12px'></div>";
      content.appendChild(p);
      var exSel = $("#dd-ex");
      exSel.innerHTML = D.extracts.map(function (e) { return "<option>" + esc(e.extract) + "</option>"; }).join("");
      function concVector(r) { return [Math.log10(Math.max(r.conc, 0.1)), r.exposure, r.sp_microp, r.sp_decol, r.exudate, r.stembark, r.dcm, r.meoh, r.adult, r.larval]; }
      function pred() {
        var conc = num($("#dd-conc").value), exp = num($("#dd-exp").value);
        var ex = $("#dd-ex").value, sp = $("#dd-sp").value;
        var X = D.mortality.map(function (r) { return concVector(r); });
        var y = D.mortality.map(function (r) { return r.mor; });
        var st = statsFromMatrix(X);
        var w = fitRidge(X, y, 1.0, st);
        var probe = concVector({ conc: conc, exposure: exp, sp_microp: /microplus/.test(sp) ? 1 : 0, sp_decol: /decol/.test(sp) ? 1 : 0, exudate: /Exudate/.test(ex) ? 1 : 0, stembark: /Stem/.test(ex) ? 1 : 0, dcm: /DCM/.test(ex) ? 1 : 0, meoh: /MeOH/.test(ex) ? 1 : 0, adult: 1, larval: 0 });
        var pred = clamp01(predictRidge(w, probe, st) / 100);
        $("#dd-out").innerHTML = "<div class='panel ok-panel'><b>Predicted mean mortality ≈ " + Math.round(pred * 100) + "%</b> at " + conc + " mg/mL for " + esc(ex) + " on " + esc(sp) + " (" + exp + " h). (Linear surrogate — validate with probit.)</div>";
        logTrail("AgroDoseX: mortality predicted " + Math.round(pred * 100) + "% at " + conc + " mg/mL, " + exp + " h, " + ex + " / " + sp);
      }
      $("#dd-run").addEventListener("click", pred);
      var lc = el("div", "panel");
      lc.innerHTML = "<div class='panel-h'>Reported LC50 / LC90 (probit)</div>" +
        "<table class='data'><thead><tr><th>Assay</th><th>Species</th><th>Extract</th><th>LC50 (mg/mL)</th><th>95% CI</th><th>LC90</th><th>Slope</th><th>χ²/df</th></tr></thead><tbody>" +
        (D.lc50 || []).map(function (r) { return "<tr><td>" + esc(r.assay) + "</td><td>" + esc(r.species) + "</td><td>" + esc(r.extract) + "</td><td>" + round(r.lc50, 2) + "</td><td class='sml dim'>" + esc(r.lc50_ci) + "</td><td>" + round(r.lc90, 2) + "</td><td>" + round(r.slope, 2) + "</td><td>" + esc(r.chi2) + "</td></tr>"; }).join("") + "</tbody></table>";
      content.appendChild(lc);
    };
    fns.qml = function () {
      markDone("qml");
      logTrail("AgroQML: methodology layer reviewed");
      var q = el("div", "panel");
      q.innerHTML = "<div class='panel-h'>AgroQML — quantum machine learning layer</div>" +
        "<p>In the Quantum_AgroX architecture, QML is <b>one computational layer</b>, not a claim of advantage. Features are normalized and angle-encoded into qubits; a feature-map circuit builds a Hilbert-space representation; then a quantum fidelity kernel, QSVM-style classifier, variational circuit (VQC) or hybrid regression learns the mapping. Classical baselines remain the reference and every comparison uses identical splits, preprocessing and leakage controls.</p>" +
        "<table class='kv sml'><tbody>" +
        "<tr><td>Qubits</td><td class='num'>8</td></tr>" +
        "<tr><td>Feature-map depth (re-uploading)</td><td class='num'>2</td></tr>" +
        "<tr><td>Variational layers</td><td class='num'>3</td></tr>" +
        "<tr><td>Shots per circuit</td><td class='num'>1024</td></tr>" +
        "<tr><td>Models in scope</td><td class='num'>QKernel · QSVM-style · VQC · QNN (next) · QGNN (future) · QBayesOpt (future)</td></tr>" +
        "<tr><td>Model comparison / audit</td><td class='num'>offline backend audit trail (not shown in this UI)</td></tr>" +
        "</tbody></table>";
      var ph = M.phases.map(function (p) { return "<div class='sml'>• " + esc(p) + "</div>"; }).join("");
      q.innerHTML += "<div style='margin-top:12px'><div class='panel-h sml'>Roadmap</div>" + ph + "</div>";
      content.appendChild(q);
    };
    fns.md = function () {
      markDone("md");
      logTrail("AgroMD: dynamics pipeline described");
      var m = el("div", "panel");
      m.innerHTML = "<div class='panel-h'>AgroMD — molecular dynamics</div>" +
        "<p>Workflow: parametrize complex (GROMACS), solvate, equilibrate, run production (tens to hundreds of ns), then compute <b>RMSD</b> (complex stability), <b>RMSF</b> (residue flexibility), <b>interaction persistence</b> (H-bonds/salt bridges/π-stacking lifetimes) and <b>MM/PBSA</b> binding free energies. These confirm or interrogate the static docking poses and surface cryptic pockets.</p>" +
        "<table class='kv sml'><tbody>" +
        "<tr><td>Engine</td><td class='num'>GROMACS (backend, not bundled in this demo UI)</td></tr>" +
        "<tr><td>Inputs</td><td class='num'>Docked complex · force field · protonation state</td></tr>" +
        "<tr><td>Outputs</td><td class='num'>RMSD/RMSF traces · interaction occupancy · ΔG<sub>bind</sub> (MM/PBSA)</td></tr>" +
        "<tr><td>Status</td><td class='num'>schema + parameters ready — backend to attach in Phase 3/4</td></tr></tbody></table>";
      content.appendChild(m);
    };
    fns.resistance = function () {
      var p = el("div", "panel");
      p.innerHTML = "<div class='panel-h'>AgroResistanceScan — mutation → structure → affinity</div>" +
        "<p class='sml dim'>Pick a known resistance mutation. The scan reports the structural hypothesis, the expected affinity change on the selected receptor, and the literature cue. <b>Hypothesis-generation only.</b></p>" +
        "<div class='mlgens'><label>Mutation</label><select id='mut-sel'></select></div>" +
        "<div id='mut-out' style='margin-top:10px'></div>";
      content.appendChild(p);
      var sel = $("#mut-sel");
      var applicable = M.mutations.filter(function (m) { return !m.target || m.target.indexOf(APP.targetKey) >= 0; });
      sel.innerHTML = applicable.map(function (m) { return "<option>" + esc(m.name) + "</option>"; }).join("");
      function draw() {
        var name = sel.value;
        var m = M.mutations.find(function (x) { return x.name === name; });
        if (!m) return;
        var t = M.targets.find(function (x) { return x.key === APP.targetKey; });
        var aff = APP.lastLab && APP.lastLab.targetKey === APP.targetKey ? APP.lastLab.affC : null;
        markDone("resistance");
        logTrail("AgroResistanceScan: mutation " + name + " evaluated" + (t ? " on " + t.name : ""));
        var d = null;
        if (m.delta != null) d = m.delta;
        var before = aff;
        var after = (aff != null && d != null) ? aff + d : null;
        $("#mut-out").innerHTML =
          "<table class='kv sml'><tbody>" +
          "<tr><td>Mutation / event</td><td><b>" + esc(m.name) + "</b></td></tr>" +
          "<tr><td>Gene / target</td><td>" + esc(m.gene) + (t ? " · " + esc(t.name) : "") + "</td></tr>" +
          "<tr><td>Expected effect</td><td>" + esc(m.effect) + "</td></tr>" +
          "<tr><td>Predicted affinity change</td><td class='num'>" + (d != null ? "+" + d.toFixed(1) + " kcal/mol (weaker binding → resistance)" : "metabolic — not structural") + "</td></tr>" +
          (before != null && after != null ? "<tr><td>Affinity before → after (last lab combo)</td><td class='num'>" + round(before, 2) + " → " + round(after, 2) + " kcal/mol</td></tr>" : "") +
          "<tr><td>Context</td><td>" + esc(m.note) + "</td></tr>" +
          "<tr><td>Lit. cue</td><td>" + esc(m.lit) + "</td></tr></tbody></table>" +
          "<p class='sml warn-bc'>Structural hypothesis only — confirm with molecular dynamics and bioassay under selection pressure.</p>";
      }
      draw();
      sel.addEventListener("change", draw);
    };
    fns.select = function () {
      markDone("select");
      logTrail("AgroSelect: selectivity rules reviewed");
      var s = el("div", "panel");
      s.innerHTML = "<div class='panel-h'>AgroSelect — pest vs beneficial selectivity</div>";
      var rows = M.selectivity.map(function (r) {
        var flag = r.flag === "HIGH" ? "tag bad" : r.flag === "MOD" ? "tag mid" : "tag ok";
        return "<tr><td>" + esc(r.family) + "</td><td class='sml dim'>" + esc(r.chemistry) + "</td><td>" + esc(r.beneficial) + "</td><td><span class='" + flag + "'>" + esc(r.flag) + "</span></td><td class='sml dim'>" + esc(r.note) + "</td></tr>";
      }).join("");
      s.innerHTML += "<table class='data'><thead><tr><th>Family</th><th>Chemistry cue</th><th>Beneficial impact</th><th>Flag</th><th>Note</th></tr></thead><tbody>" + rows + "</tbody></table>" +
        "<p class='sml dim' style='margin-top:10px'>Pair this with the Molecule Lab hit: compute chemistry cues (logP/HBA/volatility) and read the appropriate row above. A positive pest-vs-beneficial margin means the candidate controls the pest at doses that spare key pollinators/parasitoids.</p>";
      content.appendChild(s);
    };
    fns.ecorisk = function () {
      markDone("ecorisk");
      logTrail("AgroEcoRisk: environmental flags reviewed");
      var e = el("div", "panel");
      e.innerHTML = "<div class='panel-h'>AgroEcoRisk — non-target risk prioritization</div>";
      var rows = M.ecorisk.map(function (r) {
        var lv = r.level === "HIGH" ? "tag bad" : r.level === "MOD" ? "tag mid" : "tag ok";
        return "<tr><td>" + esc(r.rule) + "</td><td>" + esc(r.flag) + "</td><td><span class='" + lv + "'>" + esc(r.level) + "</span></td></tr>";
      }).join("");
      e.innerHTML += "<table class='data'><thead><tr><th>Property cue</th><th>Risk flag</th><th>Level</th></tr></thead><tbody>" + rows + "</tbody></table>";
      var lab = APP.lastLab;
      if (lab) {
        var flags = [];
        [lab.fA, lab.fB].forEach(function (f) {
          if (f.clogp >= 4) flags.push("logP ≥ 4 (“" + round(f.clogp, 1) + "”) — bioconcentration");
          if (f.clogp <= 2) flags.push("soluble low-logP (“" + round(f.clogp, 1) + "”) — aquatic mobility");
        });
        e.innerHTML += "<div style='margin-top:12px'><div class='panel-h sml'>Applied to last Molecule Lab pair</div>" +
          (flags.length ? "<div>" + flags.map(function (f0) { return "<div class='sml'>• " + esc(f0) + "</div>"; }).join("") + "</div>" : "<div class='sml dim'>No strong physical-chemical flags from that pair.</div>") + "</div>";
      }
      e.innerHTML += "<p class='sml dim' style='margin-top:10px'>Risk flags prioritize follow-up (residue, bee, aquatic) studies — they are not regulatory assessments.</p>";
      content.appendChild(e);
    };
    fns.synergyx = function () {
      markDone("synergyx");
      logTrail("AgroSynergyX: synergy workbench opened");
      var s = el("div", "panel");
      s.innerHTML = "<div class='panel-h'>AgroSynergyX — combination / synergy</div>" +
        "<p class='sml dim'>Uses the same pair plus receptor as the Molecule Lab. Run to refresh the synergy index into this stage.</p>" +
        "<div style='display:flex;gap:10px;flex-wrap:wrap;margin:8px 0'>" +
        "<button class='primary' id='syn-run'>Run synergy now</button>" +
        "<button class='ghost' id='syn-to-lab'>Open in Molecule Lab</button></div><div id='syn-out'></div>";
      content.appendChild(s);
      $("#syn-run").addEventListener("click", function () { labAnalyze("#syn-out"); });
      $("#syn-to-lab").addEventListener("click", function () { showTab("lab"); });
    };
    fns.optimize = function () {
      markDone("optimize");
      logTrail("AgroOptimize: multi-objective ranking computed");
      var o = el("div", "panel");
      o.innerHTML = "<div class='panel-h'>AgroOptimize — multi-objective ranking</div>" +
        "<p class='sml dim'>Composite = 0.55 × normalized efficacy + 0.30 × inverse LC50 + 0.15 × evidence completeness (demo weights — tunable per project).</p>" +
        "<table class='data' id='opt-t'><thead><tr><th>Extract</th><th>Records</th><th>Mean mortality</th><th>Best LC50 (mg/mL)</th><th>Composite</th></tr></thead><tbody id='opt-b'></tbody></table>" +
        "<p class='sml dim' style='margin-top:10px'>Lab hits, selectivity flags and eco-risk flags join this matrix in the full engine (Phase 7 AgroX decision score).</p>";
      content.appendChild(o);
      var effs = D.extracts.map(function (e) { return e.mean_mortality / 100; });
      var effMin = mean(effs) - 0 * mean(effs); 
      var emax = Math.max.apply(null, effs), emin = Math.min.apply(null, effs);
      var lcRows = D.lc50 || [];
      var invLc = D.extracts.map(function (e) {
        var matches = lcRows.filter(function (r) { return r.extract === e.extract; });
        if (!matches.length) return 0.5;
        var vals = matches.map(function (r) { return r.lc50 > 0 ? r.lc50 : 40; });
        var best = Math.min.apply(null, vals);
        return 1 / (1 + best / 10);
      });
      var lmax = Math.max.apply(null, invLc), lmin = Math.min.apply(null, invLc);
      var body = D.extracts.map(function (e, i) {
        var effN = (emax > emin) ? (effs[i] - emin) / (emax - emin) : 0.5;
        var lcN = (lmax > lmin) ? (invLc[i] - lmin) / (lmax - lmin) : 0.5;
        var comp = 0.55 * effN + 0.30 * lcN + 0.15 * 0.5;
        var match = lcRows.filter(function (r) { return r.extract === e.extract; });
        var bestLc = match.length ? Math.min.apply(null, match.map(function (r) { return r.lc50 > 0 ? r.lc50 : 40; })) : null;
        return { e: e, comp: comp, bestLc: bestLc };
      }).sort(function (a, b) { return b.comp - a.comp; });
      $("#opt-b").innerHTML = body.map(function (x) {
        return "<tr><td><b>" + esc(x.e.extract) + "</b></td><td>" + x.e.records + "</td><td>" + round(x.e.mean_mortality, 1) + "%</td><td>" + (x.bestLc != null ? round(x.bestLc, 2) : "—") + "</td><td><b>" + round(x.comp, 3) + "</b></td></tr>";
      }).join("");
    };
    fns.report = function () {
      markDone("report");
      logTrail("AgroReport: reproducible report generated");
      var r = el("div", "panel report");
      r.innerHTML = "<div class='panel-h'>AgroReport — reproducible summary</div>" +
        "<p class='sml dim'>Parameters grouped under their sub-module. Quantum_AgroX™ · <i>Commiphora swynnertonii</i> vs <i>R. microplus</i> &amp; <i>R. decoloratus</i>.</p>" +
        "<div id='rep-body'></div>";
      content.appendChild(r);
      var b = $("#rep-body");
      function txt() {
        var html = "";
        var runs = APP.ledger.length;
        var best = runs ? APP.ledger.slice().sort(function (a, b) { return a.affC - b.affC; })[0] : null;
        var phytoUsed = !!(APP.phytoA || APP.phytoB);
        var synthUsed = !!(APP.synthA || APP.synthB);
        var n = 0;
        function h(title) { n++; return "<h3>" + n + " · " + title + "</h3>"; }

        /* Brief summary */
        html += "<h3>Summary</h3><table class='kv sml'><tbody>" +
          "<tr><td>Pair</td><td class='num'>" + (best ? esc(best.smiA) + " + " + esc(best.smiB) : esc((APP.phytoA || "—") + " + " + (APP.phytoB || "—"))) + "</td></tr>" +
          (best ? "<tr><td>Best affinity / synergy / dose</td><td class='num'>" + best.affC + " kcal/mol · " + best.synergy + " (" + esc(best.syn) + ") · " + esc(best.dose) + " · " + esc(best.target) + "</td></tr>" : "") +
          "<tr><td>Caveat</td><td class='num'>Pre-screening hypotheses — docking scores are not proof of inhibition.</td></tr>" +
          "</tbody></table>";

        if (!phytoUsed && !synthUsed && !runs) {
          html += "<div class='dim'>Nothing used yet — start at AgroPhytoX, AgroSyntheticX or the Molecule Lab.</div>";
        }

        if (phytoUsed) {
          html += h("AgroPhytoX — phytochemical pair") + "<table class='kv sml'><tbody>" +
            [["Phytochemical A", APP.phytoA], ["Phytochemical B", APP.phytoB]].map(function (row) {
              var hit = phytoSmilesFor(row[1]);
              return "<tr><td>" + row[0] + "</td><td class='num'>" + esc(row[1] || "—") +
                (hit ? "<br><span class='kbd'>" + esc(hit.smi) + "</span>" : "<br>no curated SMILES") + "</td></tr>";
            }).join("") + "</tbody></table>";
        }

        if (synthUsed) {
          function featRow(f) {
            if (!f) return "";
            return "<tr><td>MW / cLogP</td><td class='num'>" + round(f.mw, 1) + " / " + round(f.clogp, 2) + "</td></tr>";
          }
          html += h("AgroSyntheticX — drawn molecules") + "<table class='kv sml'><tbody>" +
            (APP.synthA ? "<tr><td>Synthetic A</td><td class='num kbd'>" + esc(APP.synthA) + "</td></tr>" + featRow(APP.synthFeatA) : "") +
            (APP.synthB ? "<tr><td>Synthetic B</td><td class='num kbd'>" + esc(APP.synthB) + "</td></tr>" + featRow(APP.synthFeatB) : "") +
            "</tbody></table>";
        }

        if (runs) {
          var anchored = !best.lowConf;
          html += h("AgroTargetX / AgroDockX — target & best run") + "<table class='kv sml'><tbody>" +
            "<tr><td>Receptor target</td><td class='num'>" + esc(best.target) + " (" + (anchored ? "bundled anchor" : "generic calibration") + ")</td></tr>" +
            "<tr><td>Best combo</td><td class='num kbd'>" + esc(best.smiA) + " + " + esc(best.smiB) + "</td></tr>" +
            "<tr><td>Affinity / synergy / dose</td><td class='num'>" + best.affC + " kcal/mol · " + best.synergy + " (" + esc(best.syn) + ") · " + esc(best.dose) + "</td></tr></tbody></table>";
          html += h("AgroDoseX / AgroQML — model parameters") + "<table class='kv sml'><tbody>" +
            "<tr><td>Affinity surrogate</td><td class='num'>Ridge QSAR, " + (M.trainligands ? M.trainligands.length : 20) + " reference ligands</td></tr>" +
            "<tr><td>Synergy heuristic</td><td class='num'>0.35 · (1 − Tanimoto) + 0.20 · LogP balance + 0.20 · size mixing + 0.25 · affinity gain</td></tr>" +
            "<tr><td>Dose window</td><td class='num'>5–80 mg/mL topical (study-tested range)</td></tr>" +
            "</tbody></table>";
        }

        html += "<div style='margin-top:8px'><details><summary class='sml dim'>Audit trail (" + APP.trail.length + " events — click to expand)</summary><div class='trailbox' style='margin-top:8px'>" +
          (APP.trail.length ? APP.trail.map(function (t0) { return "<div>" + esc(t0.t) + " — " + esc(t0.m) + "</div>"; }).join("") : "<div class='dim'>No session actions yet.</div>") +
          "</div></details></div>";

        html += "<p class='sml dim'>Generated " + esc(nowStamp()) + " · Quantum_AgroX interactive layer (demo).</p>";
        html += "<div style='margin-top:10px;display:flex;gap:10px;flex-wrap:wrap'><button class='ghost' id='rep-print'>Print / Save as PDF</button>" +
          "<a class='ghost' href='assets/AgroReport.docx' download='AgroReport.docx' style='display:inline-block;padding:13px 18px'>Download full AgroReport (.docx)</a></div>";
        return html;
      }
      b.innerHTML = txt();
      $("#rep-print").addEventListener("click", function () { window.print(); });
    };

    var fn = fns[mod.id];
    if (fn) { fn(); } else {
      content.innerHTML = "<div class='panel'><p class='dim'>Module staged — schema ready for backend.</p></div>";
    }
    body.appendChild(el("div", "panel", stageNav()));
    wireNav();
  }

  /* ---------------- Evidence ---------------- */
  function litRender() {
    var wrap = $("#lit-tables");
    var q = (($("#lit-q") || {}).value || "").trim().toLowerCase();
    if (!D) return;
    function renderTable(title, header, rows) {
      var flt = rows.filter(function (r) {
        if (!q) return true;
        return header.map(function (h, i) { return String(r[i] != null ? r[i] : r[h] || "").toLowerCase(); })
          .some(function (v) { return v.indexOf(q) >= 0; });
      });
      var box = el("div", "panel");
      box.appendChild(el("div", "panel-h", title + ' <span class="dim">(' + flt.length + " shown)</span>"));
      var tbl = el("table", "data");
      tbl.innerHTML = "<thead><tr>" + header.map(function (h) { return "<th>" + esc(h) + "</th>"; }).join("") + "</tr></thead><tbody>" +
        flt.map(function (r) {
          return "<tr>" + header.map(function (h, i) {
            var v = (r[i] != null ? r[i] : r[h]);
            return "<td>" + esc(v) + "</td>";
          }).join("") + "</tr>";
        }).join("") + "</tbody>";
      if (!flt.length) tbl.querySelector("tbody").innerHTML = "<tr><td colspan='" + header.length + "' class='dim'>No rows match.</td></tr>";
      box.appendChild(tbl);
      wrap.appendChild(box);
    }
    var dock = (D.docking || []).map(function (r) {
      return [r.species, r.extract, r.compound, round(r.affinity, 2) + " kcal/mol", pctLink(r.compound + " AChE binding", "verify"), r.interactions];
    });
    var lc = (D.lc50 || []).map(function (r) {
      return [r.assay, r.species, r.extract, round(r.lc50, 2), r.lc50_ci, round(r.lc90, 2), round(r.slope, 2), r.chi2];
    });
    var val = (D.validation || []).map(function (r) {
      return [r.pdb, round(r.calc, 1), round(r.exp, 1), round(r.abs_err, 1), r.rmsd_pass ? "PASS" : "review", round(r.rmsd, 2) + " Å RMSD"];
    });
    var ach = (D.ach_models || []).map(function (r) {
      return [r.species, r.generator, round(r.verify3d, 1) + "%", round(r.errrat, 1), round(r.prosa, 2), round(r.rama, 1) + "%", round(r.disallowed, 1) + "%"];
    });
    var phyt = (D.phytochemicals || []).slice(0, 30).map(function (r) {
      return [r.compound, r.cls, r.plant_part, r.solvent, round(r.area, 2) + "%", pctLink(r.compound + " acaricide", "verify")];
    });
    wrap.innerHTML = "";
    renderTable("Top docking candidates — reported re-docking", ["Species", "Extract", "Compound", "Affinity", "External check", "Interactions"], dock);
    renderTable("LC50 / LC90 probit — reported", ["Assay", "Species", "Extract", "LC50 (mg/mL)", "95% CI", "LC90", "Slope", "χ²/df"], lc);
    renderTable("Docking validation (PDB complexes)", ["PDB", "Calc ΔG kJ/mol", "Exp ΔG kJ/mol", "Abs. error", "RMSD check", "Notes"], val);
    renderTable("AChE homology model QA — reported", ["Species", "Generator", "Verify3D", "ERRAT", "ProSA Z", "Ramachandran", "Disallowed"], ach);
    renderTable("Major phytochemicals — 30 shown", ["Compound", "Class", "Plant part", "Solvent", "Area %", "External check"], phyt);
  }
  function litVerificationPanel() {
    var wrap = $("#lit-verify");
    var rows = (window.AGROX.LIT || []).map(function (l) {
      var k = l.status && l.status.indexOf("Reported") === 0;
      return "<div class='panel vrow'><div><b>" + esc(l.focus) + "</b> <span class='" + (k ? "tag ok" : "tag bad") + "'>" + esc(l.status) + "</span></div>" +
        "<div class='sml dim'>" + esc(l.detail) + "</div><div class='sml'><i>" + esc(l.verify) + "</i></div></div>";
    }).join("");
    if (wrap) wrap.innerHTML = rows;
  }

  /* ---------------- Chat ---------------- */
  var CHAT = {
    replies: [
      { keys: ["hello", "hi ", "hey", "whats up"], ans: function () { return "Hello — I'm the Quantum_AgroX assistant. I cover the 15 modules, the benchmark statistics, target library and safeguards. Try a pilot question below."; } },
      { keys: ["start pipeline", "start the pipeline", "run pipeline", "begin pipeline"], ans: function () { return "Click 'Start Agro Chemical Discovery Pipeline' on Home — or I can open it for you: it walks AgroDataHub → AgroPhytoX → AgroSyntheticX → AgroTargetX → AgroSiteMap → AgroDockX → AgroDoseX → AgroQML → AgroMD → AgroResistanceScan → AgroSelect → AgroEcoRisk → AgroSynergyX → AgroOptimize → AgroReport."; } },
      { keys: ["best extract", "most effective", "top extract"], ans: function () { var top = (D.extracts || [])[0]; return top ? "<b>" + esc(top.extract) + "</b> leads on reported mean mortality at " + top.mean_mortality + "% (peak " + top.max_mortality + "%, n=" + top.records + ")." : "No extract summary loaded."; } },
      { keys: ["top docking", "best docking", "binding affinity", "docking hit"], ans: function () { var d0 = (D.docking || [])[0]; if (!d0) return "No docking rows loaded."; var best = D.docking.slice().sort(function (a, b) { return num(a.affinity) - num(b.affinity); })[0]; return "Strongest reported re-docking: <b>" + esc(best.compound) + "</b> in " + esc(best.extract) + " vs " + esc(best.species) + " at " + round(best.affinity, 2) + " kcal/mol. Interactions: " + esc(best.interactions || "n/a") + "."; } },
      { keys: ["lc50", "lethal concentration"], ans: function () { var lc = D.lc50 || []; var min = lc.slice().sort(function (a, b) { return num(a.lc50) - num(b.lc50); })[0]; return min ? "Lowest reported LC50: " + round(min.lc50, 2) + " mg/mL for " + esc(min.extract) + " (" + esc(min.assay) + ", " + esc(min.species) + "); LC90 " + round(min.lc90, 2) + " mg/mL." : "No LC50 rows loaded."; } },
      { keys: ["average mortality", "mean mortality"], ans: function () { var tot = 0, cnt = 0; (D.extracts || []).forEach(function (e) { tot += e.mean_mortality * e.records; cnt += e.records; }); return cnt ? "Across " + D.stats.records + " records the overall mean mortality is about " + round(tot / cnt, 1) + "%. See the dose-response stage for per-point predictions." : "No data loaded."; } },
      { keys: ["synergy", "combination", "additive", "loewe"], ans: function () { return "AgroSynergyX uses a Loewe-inspired heuristic — structural complementarity (1 − Tanimoto), LogP balance, size mixing and affinity gain of a 1:1 mix. High index = pair worth an isobologram study. Run it in the Molecule Lab or Pipeline → AgroSynergyX."; } },
      { keys: ["dosage", "dose", "concentration"], ans: function () { return "Dosing is guided by the study's tested topical window (5–80 mg/mL). The Molecule Lab returns banded guidance per candidate; AgroDoseX predicts mortality at chosen concentration × exposure."; } },
      { keys: ["caveat", "limitation", "reliable", "risk", "hazard", "safeguard"], ans: function () { return "Safeguards: (1) paper-derived aggregate means pending replicate release; (2) docking/surrogates are hypotheses, not AChE-inhibition proof; (3) no field efficacy or environmental safety from docking alone; (4) QML is one layer, never a superiority claim; (5) wear PPE with botanical extracts; (6) check local pesticide law."; } },
      { keys: ["source", "where does the data", "provenance", "dataset"], ans: function () { return "The Evidence tab exposes the provenance ledger. Data: 12 benchmark CSVs + the Quantum_AgroX architecture document. External-verification links (Scholar/PubMed) are embedded per row."; } },
      { keys: ["qml", "quantum", "qubit", "kernel"], ans: function () { return "AgroQML angle-encodes features into 8 qubits (re-uploading depth 2, variational layers 3, 1024 shots) and learns with quantum fidelity kernels, QSVM-style classifiers, VQC and (next) QNN/QGNN. It is one computational layer with classical baselines as reference; comparisons live in the offline audit trail — deliberately not marketed as accuracy in this UI."; } },
      { keys: ["accuracy", "performance", "score", "f1", "roc", "rmse", "mae", "r2"], ans: function () { return "This UI deliberately does not parade model accuracy: the architecture guidance says a quantum model is not superior just because it scores well. Reproducible comparisons (same splits/preprocessing/leakage controls) are tracked in the backend audit trail and AgroReport."; } },
      { keys: ["resistance", "mutation"], ans: function () { return "AgroResistanceScan works mutation → structure → affinity. Try Pipeline → AgroResistanceScan: pick e.g. AChE G119S (OP/carbamate) or VGSC kdr T929I+L1014F (pyrethroid) and see the predicted binding change on the selected receptor."; } },
      { keys: ["md ", "dynamics", "mmpbsa", "rmsd", "rmsf"], ans: function () { return "AgroMD runs the GROMACS backend (not bundled in the demo UI): solvate, equilibrate, produce RMSD/RMSF traces, interaction persistence and MM/PBSA ΔG for docked poses — this is what interrogates cryptic pockets."; } },
      { keys: ["selectivity", "bee", "pollinator", "beneficial"], ans: function () { return "AgroSelect scores pest-vs-beneficial margin. High-flag families: OP/carbamate- and neonicotinoid-type actives (pollinator caution). Lower: terpenoid and phenolic botanicals (short residual, repellent). See Pipeline → AgroSelect."; } },
      { keys: ["ecorisk", "ecotoxic", "ecology", "environment"], ans: function () { return "AgroEcoRisk prioritizes non-target flags: logP≥4 (bioconcentration), soluble-low-logP (aquatic mobility), neuroactive pharmacophores (non-target neurotoxicity), residue persistence. Flags drive study priorities, not verdicts."; } },
      { keys: ["optimize", "multi-objective", "ranking", "candidate"], ans: function () { return "AgroOptimize composes efficacy, inverse-LC50, selectivity and eco-risk into a weighted candidate ranking. Demo weights in Pipeline → AgroOptimize; the production engine (Phase 7) makes the AgroX decision score fully auditable."; } },
      { keys: ["report", "audit", "print"], ans: function () { return "AgroReport gives a brief summary plus parameters for only the sub-modules used in the session — then Print → Save as PDF. Run it at Pipeline → AgroReport."; } },
      { keys: ["target", "receptor"], ans: function () { return "The target library now includes mites/ticks, insects, nematodes, fungi, weeds and detox enzymes — e.g. AChE (CAA06980/6ARX or CAA11702/5YDH), VGSC, nAChR, CYP51, SDH, β-tubulin, PSII-D1, EPSPS, levamisole-site AChR. Type any name in the Molecule Lab; if unknown, it falls back to generic AChE calibration with a warning."; } },
      { keys: ["reproduce", "how to run", "install", "pipeline locally"], ans: function () { return "Clone git@github.com:SheenaPravin/Quantum_AgroX.git (branch Feat/Sheena), pip install -r Quantum_AgroX_QML_Project/requirements.txt, python main.py. Scripts 01–05 mirror the pipeline's classical + QML stages."; } },
      { keys: ["medxai", "who", "project", "team"], ans: function () { return "Quantum_AgroX is a MEDxAI innovation project: an evidence-integrating agrochemical intelligence platform, first benchmarked on Commiphora swynnertonii botanical acaricide discovery against cattle ticks."; } }
    ],
    fallback: function (w) {
      if (/\d/.test(w)) return "I don't see a matching question — try 'best extract', 'top docking hit', 'LC50', 'synergy', 'resistance', 'optimize', or 'how to start the pipeline'.";
      for (var i = 0; i < M.modules.length; i++) {
        var m = M.modules[i];
        var stem = m.name.replace("™", "").replace(/Agro/, "").toLowerCase();
        if (w.indexOf(stem) >= 0) return m.name + " — " + m.func + ". " + m.desc + " Open it in the Pipeline to run the live demo of this stage.";
      }
      return "I can help with module guidance, benchmark statistics, the target library and safeguards. Try: 'best extract', 'top docking hit', 'LC50', 'QML', 'resistance', or the pilots above.";
    }
  };
  function chatReply(text) {
    var w = text.toLowerCase();
    var hits = CHAT.replies.filter(function (r) {
      return r.keys.some(function (k) { return w.indexOf(k) >= 0; });
    });
    hits.sort(function (a, b) {
      return b.keys.filter(function (k) { return w.indexOf(k) >= 0; }).length - a.keys.filter(function (k) { return w.indexOf(k) >= 0; }).length;
    });
    return hits.length ? hits[0].ans() : CHAT.fallback(w);
  }
  function chatSend() {
    var inp = $("#chat-in");
    var text = inp.value.trim();
    if (!text) return;
    inp.value = "";
    var box = $("#chat-msgs");
    box.appendChild(el("div", "msg user", esc(text)));
    box.appendChild(el("div", "msg bot", chatReply(text)));
    box.scrollTop = box.scrollHeight;
  }
  function buildPilots() {
    var wrap = $("#pilots");
    var list = ["Start the pipeline", "Best extract", "Top docking hit", "What is AgroQML?", "Resistance", "AgroSelect", "Caveats"];
    list.forEach(function (p) {
      var c = el("span", "pilot", p);
      c.addEventListener("click", function () {
        $("#chat-in").value = p;
        chatSend();
      });
      wrap.appendChild(c);
    });
    var welcome = el("div", "msg bot", "Hi — Quantum_AgroX assistant here. Ask about the 15 modules, benchmark stats, targets, or safeguards.");
    $("#chat-msgs").appendChild(welcome);
  }

  /* ---------------- Tabs / init ---------------- */
  var TAB_HASH = { home: "#home", pipeline: "#pipeline", lab: "#lab", evid: "#evid", chat: "#chat" };
  function showTab(name) {
    APP.tab = name;
    $$(".tablink").forEach(function (b) { b.classList.toggle("active", b.dataset.tab === name); });
    $$(".tabpane").forEach(function (p) { p.classList.toggle("active", p.id === "pane-" + name); });
    if (name === "evid") { litRender(); litVerificationPanel(); }
    if (name === "pipeline" && M) { renderStepper(); renderStage(); }
    try { history.replaceState(null, "", TAB_HASH[name] || "#home"); } catch (e) {}
  }
  function initTabs() {
    $$(".tablink").forEach(function (b) {
      b.addEventListener("click", function () { showTab(b.dataset.tab); });
    });
  }
  function deepLinkInit() {
    if (typeof location === "undefined") return;
    var h = (location.hash || "#home").replace(/^#/, "");
    var names = { pipeline: "pipeline", lab: "lab", evid: "evid", chat: "chat" };
    if (names[h]) showTab(names[h]); else showTab("home");
    var pipeParam = (typeof URLSearchParams !== "undefined" && location.search) ? new URLSearchParams(location.search).get("pipe") : null;
    if (pipeParam) {
      APP.pendingPipe = pipeParam;
      if (M) {
        var mi = M.modules.findIndex(function (m) { return m.id === pipeParam; });
        if (mi >= 0) openPipelineAt(mi);
      }
    }
    var s1 = (typeof URLSearchParams !== "undefined" && location.search) ? new URLSearchParams(location.search).get("s1") : null;
    var s2 = (typeof URLSearchParams !== "undefined" && location.search) ? new URLSearchParams(location.search).get("s2") : null;
    if (s1) $("#smi-a").value = s1;
    if (s2) $("#smi-b").value = s2;
    if (s1 && s2) {
      loadRDKit().then(function () { labAnalyze(); })["catch"](function (e) { $("#lab-out").innerHTML = "<div class='warn'>" + esc(e.message) + "</div>"; });
    }
  }
  function register() {
    initTabs();
    buildPilots();
    $("#start-pipe").addEventListener("click", function () { openPipelineAt(0); logTrail("Pipeline started from Home"); });
    $("#browse-modules").addEventListener("click", function () { showTab("home"); var g = $("#modgrid"); if (g) g.scrollIntoView({ behavior: "smooth", block: "start" }); });
    $("#analyze-btn").addEventListener("click", function () { labAnalyze("#lab-out"); });
    $("#lab-use-phyto").addEventListener("click", function () {
      var hint = useTopPhytochemicals();
      if (hint) {
        var hc = $("#target-hint");
        if (hc && hc.parentNode) hc.parentNode.insertBefore(hint, hc.parentNode.children[0]);
      }
    });
    $("#smi-a").addEventListener("keydown", function (e) { if ((e.metaKey || e.ctrlKey) && e.key === "Enter") labAnalyze("#lab-out"); });
    $("#smi-b").addEventListener("keydown", function (e) { if ((e.metaKey || e.ctrlKey) && e.key === "Enter") labAnalyze("#lab-out"); });
    $("#target").addEventListener("input", function () { updateTargetHint(); });
    $("#lit-q").addEventListener("input", function () { litRender(); });
    $("#chat-send").addEventListener("click", chatSend);
    $("#chat-in").addEventListener("keydown", function (e) { if (e.key === "Enter") chatSend(); });
    if ($("#ex-a")) $("#ex-a").addEventListener("change", function () { var v = $("#ex-a").value; if (v) { $("#smi-a").value = v; } });
    if ($("#ex-b")) $("#ex-b").addEventListener("change", function () { var v = $("#ex-b").value; if (v) { $("#smi-b").value = v; } });
    if ($("#ex-target")) $("#ex-target").addEventListener("change", function () { var v = $("#ex-target").value; if (v) { $("#target").value = v; updateTargetHint(); } });
    if ($("#draw-open")) $("#draw-open").addEventListener("click", function () {
      if (!window.QXSketch) return;
      QXSketch.open($("#smi-a").value, $("#smi-b").value);
      QXSketch.onUse(function (which, smi) {
        if (which === "A") $("#smi-a").value = smi; else $("#smi-b").value = smi;
      });
    });
    /* defaults */
    $("#smi-a").value = "CC1=CCC2CC1C2(C)C";   /* alpha-pinene */
    $("#smi-b").value = "COc1cc(CC=C)ccc1O";   /* eugenol */
    $("#target").value = "Acetylcholinesterase (Rhipicephalus microplus)";
    renderTargetDatalist();
    updateTargetHint();
  }
  function updateTargetHint() {
    var v = $("#target").value.trim();
    var h = $("#target-hint");
    if (!M || !h) return;
    var t = resolveTarget(v);
    if (!v) h.innerHTML = "";
    else if (t) h.innerHTML = "Matched: <b>" + esc(t.name) + "</b> — " + esc(t.organism) + (t.anchorDG ? " · anchor " + esc(t.anchor) + " (" + t.anchorDG + " kcal/mol)" : "");
    else h.innerHTML = "<span class='warn-bc'>Not in the bundled library — will use generic insect-AChE calibration.</span>";
  }
  function renderExampleSelects() {
    var ta = $("#ex-a"), tb = $("#ex-b"), exT = $("#ex-target");
    if (!M || (!ta && !tb && !exT)) return;
    var smiOps = "<option value=''>— pick a structure —</option>" + M.smilesExamples.map(function (e) {
      return "<option value='" + esc(e.smi) + "'>" + esc(e.name) + " [" + esc(e.cls) + "]</option>";
    }).join("");
    if (ta) ta.innerHTML = smiOps;
    if (tb) tb.innerHTML = smiOps;
    if (exT) {
      exT.innerHTML = "<option value=''>— pick a receptor —</option>" + M.targets.map(function (t) {
        var label = t.name + " (" + t.organism + ")";
        return "<option value='" + esc(label) + "'>" + esc(label) + " [" + esc(t.cls) + (t.anchorKey ? " · anchored" : "") + "]</option>";
      }).join("");
    }
  }
  function renderTargetDatalist() {
    var dl = $("#targets");
    if (!dl || !M) return;
    var seen = {};
    (M.targets || []).forEach(function (t) {
      [t.name + " (" + t.organism + ")", t.name, t.organism, t.protein].forEach(function (opt) {
        if (!seen[opt] && opt && opt !== "—") { seen[opt] = 1; dl.appendChild(el("option", "", esc(""))); }
      });
    });
    /* datalist option text via attribute */
    dl.innerHTML = "";
    Object.keys(seen).filter(function (k) { return k; }).forEach(function (k) {
      var o = el("option", "", "");
      o.value = k;
      dl.appendChild(o);
    });
  }

  function fillKpis() {
    function set(id, v) { var e = $(id); if (e) e.textContent = v; }
    set("#kp-modules", M.modules.length);
    set("#kp-targets", M.targets.length);
    set("#kp-examples", M.smilesExamples ? M.smilesExamples.length : "–");
    set("#kp-datasets", D.datasets ? D.datasets.length : "–");
    set("#kp-bench", D.bench ? D.bench.length : "–");
  }
  function renderFloatMols() {
    if (!window.RDKIT || !M || !M.smilesExamples) return;
    var map = {};
    M.smilesExamples.forEach(function (e) { map[e.name.toLowerCase()] = e.smi; });
    $$(".fl-svg[data-ext]").forEach(function (el) {
      var nm = (el.getAttribute("data-ext") || "").toLowerCase();
      var smi = map[nm];
      if (!smi) return;
      try {
        var m = RDKIT.get_mol(smi);
        if (m && m.is_valid()) el.innerHTML = m.get_svg(160, 110);
        if (m && m.delete) m.delete();
      } catch (ex) {}
    });
  }
  function loadData() {
    return fetch("assets/dashboard_data.json").then(function (r) { return r.json(); }).then(function (d) {
      D = d;
      M = window.AGROX.MODEL;
      FEATURES.length = 0; M.features.forEach(function (f0) { FEATURES.push(f0); });
      homeRender();
      renderExampleSelects();
      renderTargetDatalist();
      updateTargetHint();
      fillKpis();
      renderFloatMols();
      if (APP.tab === "evid") { litRender(); litVerificationPanel(); }
      if (APP.pendingPipe && M) {
        var miP = M.modules.findIndex(function (mo) { return mo.id === APP.pendingPipe; });
        APP.pendingPipe = null;
        if (miP >= 0) openPipelineAt(miP);
      } else if (APP.tab === "pipeline") { renderStepper(); renderStage(); }
    });
  }

  function init() {
    loadData()["catch"](function (e) { console.error("data load failed", e); });
    register();
    deepLinkInit();
    var rdEl = $("#rdkit-status");
    loadRDKit().then(function () {
      sd(rdEl).textContent = "RDKit core ready — structures render in-browser.";
    })["catch"](function (e) {
      rdEl.textContent = "RDKit core failed (" + e.message + ") — descriptors/affinity disabled.";
    });
  }
  function sd(x){ return x; }

  document.addEventListener("DOMContentLoaded", init);
})();