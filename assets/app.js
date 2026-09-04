/* Quantum_AgroX — interactive dashboard app */
(function () {
  "use strict";

  /* ---------- tiny utils ---------- */
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
  function meanOr(arr) { return arr.length ? mean(arr) : 0; }
  function mean(arr) { return arr.reduce(function (a, b) { return a + num(b); }, 0) / arr.length; }
  function maxOr(arr) { return arr.length ? Math.max.apply(null, arr.map(num)) : 0; }
  function pctLink(q, label) {
    var s = "https://scholar.google.com/scholar?q=" + encodeURIComponent(q);
    var p = "https://pubmed.ncbi.nlm.nih.gov/?term=" + encodeURIComponent(q);
    return '<a class="lnk" target="_blank" rel="noopener" href="' + s + '">' + (label || "Scholar") + "</a> &nbsp;·&nbsp; " +
           '<a class="lnk" target="_blank" rel="noopener" href="' + p + '">PubMed</a>';
  }

  var APP = {
    data: null,
    RDKit: null,
    rdkitReady: false,
    tab: "home",
    chat: []
  };
  window.APP = APP;

  /* ---------- RDKit.js loader (wasm auto-locates beside the script) ---------- */
  function bootRDKit() {
    var factory = window.initRDKitModule || window.initRDKit;
    if (!factory) return Promise.reject(new Error("RDKit script not found — check your connection"));
    return factory().then(function (m) {
      APP.RDKit = m;
      APP.rdkitReady = true;
      return m;
    });
  }
  function loadRDKit() {
    return new Promise(function (resolve, reject) {
      var factory = window.initRDKitModule || window.initRDKit;
      if (factory) { bootRDKit().then(resolve)["catch"](reject); return; }
      if (window.RDKit && window.RDKit.get_mol) { APP.RDKit = window.RDKit; return resolve(); }
      /* fallback: CDN may be blocked (offline/corporate); load the vendored copy */
      var s = document.createElement("script");
      s.src = "vendor/RDKit_minimal.js";
      s.onload = function () { bootRDKit().then(resolve)["catch"](reject); };
      s.onerror = function () { reject(new Error("RDKit script not found — check your connection")); };
      document.head.appendChild(s);
    });
  }

  function rdMolFromSmiles(smiles, name) {
    var mol = APP.RDKit.get_mol(smiles);
    if (!mol || !mol.is_valid()) { if (mol) mol.delete(); throw new Error("Invalid SMILES for " + name); }
    return mol;
  }
  function molFeatures(mol) {
    var d = JSON.parse(mol.get_descriptors());
    return {
      mw: +num(d.exactmw), clogp: +num(d.CrippenClogP), tpsa: +num(d.tpsa),
      hba: num(d.lipinskiHBA), hbd: num(d.lipinskiHBD), rotb: num(d.NumRotatableBonds),
      rings: num(d.NumRings), arom: num(d.NumAromaticRings),
      heavy: num(d.NumHeavyAtoms), fsp3: +num(d.FractionCSP3),
      atoms: num(d.NumAtoms), het: num(d.NumHeteroatoms)
    };
  }
  function molSvg(mol) {
    var h = mol.get_svg(320, 260);
    return h || "<div class='dim'>structure render unavailable</div>";
  }
  function fingerprint(mol) {
    var u = mol.get_morgan_fp_as_uint8array ? mol.get_morgan_fp_as_uint8array() : null;
    if (!u) return null;
    var bits = new Uint8Array(u.length);
    for (var i = 0; i < u.length; i++) bits[i] = u[i] > 0 ? 1 : 0;
    return bits;
  }
  function tanimoto(a, b) {
    if (!a || !b) return null;
    var n = Math.min(a.length, b.length), inter = 0, union = 0;
    for (var i = 0; i < n; i++) {
      if (a[i] && b[i]) inter++;
      if (a[i] || b[i]) union++;
    }
    return union ? inter / union : 0;
  }

  /* ---------- Ridge QSAR (training-set refit, standardised features) ---------- */
  var FEATURES = window.AGROX.MODEL.features;
  function gaussSolve(A, b) {
    var n = A.length, aug = A.map(function (row, i) { return row.concat([b[i]]); });
    for (var c = 0; c < n; c++) {
      var mx = c;
      for (var r = c + 1; r < n; r++) if (Math.abs(aug[r][c]) > Math.abs(aug[mx][c])) mx = r;
      if (mx !== c) { var t = aug[c]; aug[c] = aug[mx]; aug[mx] = t; }
      var pv = aug[c][c];
      if (Math.abs(pv) < 1e-12) continue;
      for (var r2 = c + 1; r2 < n; r2++) {
        var f = aug[r2][c] / pv;
        for (var k = c; k <= n; k++) aug[r2][k] -= f * aug[c][k];
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
  function fitRidge(X, y, lambda, featStats) {
    var n = X.length, f = X[0].length;
    var Xt = X.map(function (row) { return [1].concat(row.map(function (v, j) {
      var st = featStats[j];
      return st.sd > 1e-9 ? (v - st.mean) / st.sd : 0;
    })); });
    var A = [], i, j, k;
    for (i = 0; i <= f; i++) {
      A.push(new Array(f + 1).fill(0));
      for (j = 0; j <= f; j++) {
        var s = 0;
        for (k = 0; k < n; k++) s += Xt[k][i] * Xt[k][j];
        A[i][j] = s + (i > 0 ? lambda : 0);
      }
    }
    var b = new Array(f + 1).fill(0);
    for (i = 0; i <= f; i++) for (k = 0; k < n; k++) b[i] += Xt[k][i] * y[k];
    return gaussSolve(A, b);
  }
  function predictRidge(w, x, featStats) {
    var p = w[0];
    for (var j = 0; j < x.length; j++) {
      var st = featStats[j];
      var z = st.sd > 1e-9 ? (x[j] - st.mean) / st.sd : 0;
      p += w[j + 1] * z;
    }
    return p;
  }
  function fitAffinityModel() {
    var train = window.AGROX.MODEL.train;
    var featStats = FEATURES.map(function (ft) {
      var vals = train.map(function (r) { return num(r[ft]); });
      var m = mean(vals);
      var sd = Math.sqrt(mean(vals.map(function (v) { var d = v - m; return d * d; })));
      return { name: ft, mean: m, sd: sd || 1, vals: vals };
    });
    var X = train.map(function (r) { return FEATURES.map(function (ft) { return num(r[ft]); }); });
    var y = train.map(function (r) { return num(r.dG); });
    var w = fitRidge(X, y, 1.0, featStats);
    return { w: w, stats: featStats, train: train };
  }

  /* ---------- Molecule Lab ---------- */
  function labAnalyze() {
    var out = $("#lab-out");
    out.innerHTML = "<div class='dim'>Working …</div>";
    var smilesA = $("#smi-a").value.trim();
    var smilesB = $("#smi-b").value.trim();
    var target = $("#target").value;
    if (!smilesA || !smilesB) { out.innerHTML = "<div class='warn'>Please enter two SMILES strings.</div>"; return; }
    var A, B, molA, molB;
    try {
      molA = rdMolFromSmiles(smilesA, "molecule A");
      molB = rdMolFromSmiles(smilesB, "molecule B");
    } catch (e) {
      out.innerHTML = "<div class='warn'>" + esc(e.message) + "</div>";
      return;
    }
    var fA = molFeatures(molA), fB = molFeatures(molB);
    var fpA = fingerprint(molA), fpB = fingerprint(molB);
    var tani = tanimoto(fpA, fpB);
    var model = fitAffinityModel();
    var xA = FEATURES.map(function (ft) { return fA[ft]; });
    var xB = FEATURES.map(function (ft) { return fB[ft]; });
    var xC = xA.map(function (v, i) { return (v + xB[i]) / 2; });
    var pA = predictRidge(model.w, xA, model.stats);
    var pB = predictRidge(model.w, xB, model.stats);
    var pC = predictRidge(model.w, xC, model.stats);
    var rec = window.AGROX.MODEL.receptor[target] || window.AGROX.MODEL.receptor["R. microplus"];
    var controlFeat = FEATURES.map(function (ft) {
      return num(window.AGROX.MODEL.train.filter(function (r) { return r.name === rec.controlName; })[0][ft]);
    });
    var pControl = predictRidge(model.w, controlFeat, model.stats);
    var calib = rec.control - pControl;
    var affA = pA + calib, affB = pB + calib, affC = pC + calib;

    /* synergy decomposition */
    var complement = tani == null ? 0.5 : 1 - tani;
    var logpMatch = 1 - Math.abs(fA.clogp - fB.clogp) / (Math.max(1, Math.abs(fA.clogp)) + Math.max(1, Math.abs(fB.clogp)));
    var sizeMix = 1 - Math.abs(fA.mw - fB.mw) / (fA.mw + fB.mw || 1);
    var affGain = clamp01((Math.min(affA, affB) - affC) / 2.0);
    var synergy = 0.35 * complement + 0.2 * logpMatch + 0.2 * sizeMix + 0.25 * affGain;
    function synergyLabel(s) {
      if (s >= 0.72) return { t: "Very high", c: "good" };
      if (s >= 0.55) return { t: "High", c: "good" };
      if (s >= 0.35) return { t: "Moderate", c: "warn" };
      return { t: "Low — likely additive at best", c: "bad" };
    }
    var sl = synergyLabel(synergy);

    function doseBand(aff) {
      var potency = Math.exp(-(aff + 5));            /* relative activity proxy */
      var pmax = Math.exp(5.0), pmin = Math.exp(0.6);
      var p = clamp01((potency - pmin) / (pmax - pmin));
      var doseIdx = 0.85 - 0.55 * p;                  /* map into 5–80 mg/mL window */
      var mg = 5 + doseIdx * 75;
      var band = mg <= 15 ? "Low loading (5–15 mg/mL)" : mg <= 35 ? "Moderate (15–35 mg/mL)" : mg <= 60 ? "Upper band (35–60 mg/mL)" : "Maximum tested (60–80 mg/mL)";
      return { mg: round(mg, 0), band: band };
    }
    var dbA = doseBand(affA), dbB = doseBand(affB), dbC = doseBand(affC);

    function featRows(f, other) {
      var rows = [
        ["Molecular weight (g/mol)", round(f.mw, 1)],
        ["cLogP (lipophilicity)", round(f.clogp, 2)],
        ["TPSA (Å²)", round(f.tpsa, 1)],
        ["H-bond acceptors / donors", f.hba + " / " + f.hbd],
        ["Rotatable bonds", f.rotb],
        ["Rings (aromatic)", f.rings + " (" + f.arom + ")"],
        ["Heavy atoms / heteroatoms", f.heavy + " / " + f.het],
        ["Fraction Csp3", round(f.fsp3, 2)]
      ];
      return rows.map(function (r) {
        return "<tr><td class='dim'>" + r[0] + "</td><td>" + r[1] + "</td></tr>";
      }).join("");
    }

    var html = "";
    html += "<div class='grid lab-grid'>";
    ["A", "B"].forEach(function (tag) {
      var mol = tag === "A" ? molA : molB, f = tag === "A" ? fA : fB, svg = tag === "A" ? molSvg(molA) : molSvg(molB);
      html += "<div class='panel'><div class='panel-h'>Molecule " + tag + "</div>" +
        "<div class='svgbox'>" + svg + "</div>" +
        "<table class='kv'><tbody>" + featRows(f) + "</tbody></table></div>";
    });
    html += "</div>";

    html += "<div class='grid'>";
    html += "<div class='panel'><div class='panel-h'>Predicted AChE binding affinity — " + esc(target) + "</div>" +
      "<p class='dim sml'>Surrogate ridge QSAR fitted live on 20 reference ligands, calibrated to " + esc(rec.controlDisplay) +
      " control (" + rec.control + " kcal/mol as reported in the study).</p>" +
      "<table class='kv'><tbody>" +
      "<tr><td>Molecule A affinity</td><td class='num big'>" + round(affA, 2) + " kcal/mol</td></tr>" +
      "<tr><td>Molecule B affinity</td><td class='num big'>" + round(affB, 2) + " kcal/mol</td></tr>" +
      "<tr><td>1:1 combination (mean descriptor)</td><td class='num big'>" + round(affC, 2) + " kcal/mol</td></tr>" +
      "</tbody></table></div>";

    html += "<div class='panel'><div class='panel-h'>Synergy index (structural complementarity)</div>" +
      "<p class='sml dim'>Combines Tucker-like structural overlap (Tanimoto), LogP complementarity, size mixing and predicted affinity gain in a Loewe-inspired heuristic — a pre-screen, not a formal Loewe CI.</p>" +
      "<div class='syn'><span style='width:" + Math.round(synergy * 100) + "%' class='syn-'></span></div>" +
      "<div class='syn-label " + sl.c + "'>" + esc(sl.t) + " &middot; index " + round(synergy, 2) + "</div>" +
      "<table class='kv sml'><tbody>" +
      "<tr><td>Structural complementarity (1 − Tanimoto)</td><td class='num'>" + (tani == null ? "n/a" : round(complement, 2)) + "</td></tr>" +
      "<tr><td>Tanimoto similarity (Morgan-like)</td><td class='num'>" + (tani == null ? "n/a" : round(tani, 2)) + "</td></tr>" +
      "<tr><td>LogP complementarity</td><td class='num'>" + round(logpMatch, 2) + "</td></tr>" +
      "<tr><td>Size mixing balance</td><td class='num'>" + round(sizeMix, 2) + "</td></tr>" +
      "<tr><td>Affinity gain of combination</td><td class='num'>" + round(affGain, 2) + "</td></tr>" +
      "</tbody></table></div></div>";

    html += "<div class='panel'><div class='panel-h'>Dose guidance (" + esc(target) + ")</div>" +
      "<p class='sml dim'>Mapped to the study's tested topical window <b>5–80 mg/mL</b> (bovine tick bioassays):</p>" +
      "<table class='kv'><tbody>" +
      "<tr><td>Molecule A</td><td class='num'>" + dbA.mg + " mg/mL — " + dbA.band + "</td></tr>" +
      "<tr><td>Molecule B</td><td class='num'>" + dbB.mg + " mg/mL — " + dbB.band + "</td></tr>" +
      "<tr><td>Combination</td><td class='num'>" + dbC.mg + " mg/mL — " + dbC.band + "</td></tr>" +
      "</tbody></table><p class='sml dim'>Low-to-mid concentrations follow the extract dilution ladder (5/10/20/40/80 mg/mL). Always allow for formulation, skin irritation and off-target effects.</p></div>";

    html += "<div class='panel warn-panel'><b>Caveat.</b> " + esc(window.AGROX.MODEL.caveat) + " This is a <i>pre-screening surrogate</i> — validate via molecular docking (AutoDock), then wet-lab repellency/mortality assays before any recommendation.</div>";

    var hook = window.AGROX.MODEL.train.map(function (r) { return r.name; }).join(", ");
    html += "<div class='panel dim sml'>Reference ligands used by the live fit: " + esc(hook) + ". Features standardised and ridge λ=1.0. Chlorfenvinphos is the study's positive control acaricide.</div>";

    out.innerHTML = html;
    ["A", "B"].forEach(function (tag) {
      var mol = tag === "A" ? molA : molB;
      try { mol.delete(); } catch (e) {}
    });
  }

  /* ---------- Literature tab ---------- */
  function litRender() {
    var D = APP.data;
    var wrap = $("#lit-tables");
    var q = ($("#lit-q") && $("#lit-q").value || "").trim().toLowerCase();

    function renderTable(id, title, header, rows) {
      var flt = rows.filter(function (r) {
        if (!q) return true;
        return header.map(function (h, i) { return String(r[i] || r[h] || "").toLowerCase(); })
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
      if (!flt.length) tbl.querySelector("tbody").innerHTML = "<tr><td colspan='" + header.length + "' class='dim'>No rows match your filter.</td></tr>";
      box.appendChild(tbl);
      wrap.appendChild(box);
    }

    wrap.innerHTML = "<div class='dim'>Loading tables…</div>";

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
    renderTable("dock", "Top docking candidates — reported re-docking (per extract × species)", ["Species", "Extract", "Compound", "Affinity", "External check", "Reported interactions"], dock);
    renderTable("lc50", "LC50 / LC90 probit output — reported", ["Assay", "Species", "Extract", "LC50 (mg/mL)", "95% CI", "LC90", "Slope", "χ²/df"], lc);
    renderTable("val", "Docking validation (PDB complexes)", ["PDB", "Calc. ΔG (kJ/mol)", "Exp. ΔG (kJ/mol)", "Abs. error", "RMSD check", "Notes"], val);
    renderTable("ach", "AChE homology model QA — reported", ["Species", "Generator", "Verify3D", "ERRAT", "ProSA Z", "Ramachandran", "Disallowed"], ach);
    renderTable("phyt", "Major phytochemicals (GC-MS) — 30 shown", ["Compound", "Class", "Plant part", "Solvent", "Area %", "External check"], phyt);
  }

  function litVerificationPanel() {
    var wrap = $("#lit-verify");
    var rows = (window.AGROX.LIT || []).map(function (l) {
      var k = l.detail.indexOf(" ") >= 0 ? l.detail.slice(0, 60) : "";
      return "<div class='panel vrow'><div><b>" + esc(l.focus) + "</b> <span class='" + (l.status.indexOf("Reported") === 0 ? "tag ok" : "tag bad") + "'>" + esc(l.status) + "</span></div>" +
        "<div class='sml dim'>" + esc(l.detail) + "</div><div class='sml'><i>" + esc(l.verify) + "</i></div></div>";
    }).join("");
    wrap.innerHTML = rows;
  }

  /* ---------- Home ---------- */
  function homeRender() {
    var D = APP.data;
    $("#stat-records").textContent = D.stats.records;
    $("#stat-phyt").textContent = D.stats.phytochemicals;
    $("#stat-dock").textContent = D.stats.docking_hits;
    $("#stat-val").textContent = D.stats.validation_complexes;

    var ext = el("div", "bars");
    (D.extracts || []).forEach(function (e) {
      var b = el("div", "bar-row");
      b.innerHTML = "<div class='bar-lbl'>" + esc(e.extract) + "</div>" +
        "<div class='bar-track'><div class='bar' style='width:" + Math.round(e.mean_mortality) + "%'></div></div>" +
        "<div class='bar-val'>" + round(e.mean_mortality, 1) + "%</div>";
      ext.appendChild(b);
    });
    $("#ext-chart").innerHTML = "";
    $("#ext-chart").appendChild(ext);

    var spec = el("div", "grid c3");
    (D.species || []).forEach(function (s) {
      var c = el("div", "panel card");
      c.innerHTML = "<div class='num-lg'>" + s.count + "</div><div class='dim'>" + esc(s.species) + " bioassay records</div>";
      spec.appendChild(c);
    });
    $("#spec-cards").innerHTML = "";
    $("#spec-cards").appendChild(spec);

    var cls = el("div", "chips");
    (D.phytochemical_classes || []).forEach(function (c) {
      cls.appendChild(el("span", "chip", esc(c.cls) + " × " + c.count));
    });
    $("#cls-chips").innerHTML = "";
    $("#cls-chips").appendChild(cls);
  }

  /* ---------- Chatbot ---------- */
  var CHAT = {
    replies: [
      {
        keys: ["hello", "hi ", "hey", "whats up", "hiya"],
        ans: function () {
          return "Hello! I can explain the quantum methodology, or answer statistics about the dataset. Try: “what is the best extract?”, “average mortality”, “best docking hit”, “what is QML?”, or “how do I run the pipeline?”.";
        }
      },
      {
        keys: ["average mortality", "mean mortality", "average kill", "mean kill"],
        ans: function () {
          var D = APP.data, ex = D.extracts[0];
          var tot = 0, cnt = 0;
          (D.extracts || []).forEach(function (e) { tot += e.mean_mortality * e.records; cnt += e.records; });
          var over = cnt ? (tot / cnt) : 0;
          return "Across " + D.stats.records + " aggregate bioassay records the mean mortality is about " +
            round(over, 1) + "%. The strongest extract by mean mortality is <b>" + esc(ex.extract) + "</b> at " + ex.mean_mortality + "%.";
        }
      },
      {
        keys: ["best extract", "most effective extract", "top extract", "best performing"],
        ans: function () {
          if (!APP.data.extracts) return "No extract summary loaded yet.";
          var top = APP.data.extracts[0];
          return "By reported mean mortality, <b>" + esc(top.extract) + "</b> leads with " + top.mean_mortality +
            "% mean (peak " + top.max_mortality + "%, n=" + top.records + "). Exudate-DCM and Stem bark-DCM are the front-running fractions in the paper.";
        }
      },
      {
        keys: ["best docking", "top docking", "docking hit", "affinity", "binding"],
        ans: function () {
          var dock = APP.data.docking || [];
          if (!dock.length) return "No docking rows loaded.";
          var best = dock.slice().sort(function (a, b) { return num(a.affinity) - num(b.affinity); })[0];
          return "Strongest reported re-docking hit: <b>" + esc(best.compound) + "</b> in " + esc(best.extract) +
            " against " + esc(best.species) + " with " + round(best.affinity, 2) + " kcal/mol. Reported interactions: " + esc(best.interactions || "n/a") + ".";
        }
      },
      {
        keys: ["lc50", "lethal concentration", "potency"],
        ans: function () {
          var lc = APP.data.lc50 || [];
          if (!lc.length) return "LC50 table not loaded.";
          var min = lc.slice().sort(function (a, b) { return num(a.lc50) - num(b.lc50); })[0];
          return "Lowest reported LC50 is " + round(min.lc50, 2) + " mg/mL for " + esc(min.extract) + " (" + esc(min.assay) + ", " + esc(min.species) + "); LC90 " + round(min.lc90, 2) + " mg/mL.";
        }
      },
      {
        keys: ["qml", "quantum machine learning", "what is quantum", "quantum method", "circuit", "qubit"],
        ans: function () {
          return "Quantum Machine Learning (QML) here = encoding molecular/assay features into qubits (feature map) and learning with a parameterised quantum circuit (VQC or quantum kernel). We use 8 qubits, depth 2 data re-uploading and 3 variational layers. It is a research-grade proof of concept — classical baselines still beat it for small tabular data, and the dashboard runs transparent classical surrogates in-browser.";
        }
      },
      {
        keys: ["run pipeline", "run the pipeline", "reproduce", "how to run", "install", "execute"],
        ans: function () {
          return "Clone git@github.com:SheenaPravin/Quantum_AgroX.git (branch Feat/Sheena), then pip install -r Quantum_AgroX_QML_Project/requirements.txt and run 'python main.py' — it prepares data (01), classical baselines (02), then the quantum kernel regressor/classifier (03–04) and VQC (05). Pennylane is required.";
        }
      },
      {
        keys: ["baseline", "classical", "accuracy", "r2", "mae", "rmse", "score"],
        ans: function () {
          var b = APP.data.baseline_ridge || {};
          return "Classical ridge-baseline on this benchmark: MAE " + (b.MAE_mean || "—") + " ± " + (b.MAE_std || 0) +
            ", RMSE " + (b.RMSE_mean || "—") + ", R² " + (b.R2_mean || "—") + " (5-fold shuffle CV on " + (b.n_samples || 240) + " records). Quantum results are the research goal.";
        }
      },
      {
        keys: ["synergy", "combination", "mix two", "additive", "loewe"],
        ans: function () {
          return "Synergy here is a Loewe-inspired heuristic: structural complementarity (1 − Tanimoto), LogP balance, size mixing and predicted affinity gain of the 1:1 combination vs its components. High index = candidates worth a proper isobologram / Combination Index study. See the Molecule Lab tab — it computes this live for any two SMILES.";
        }
      },
      {
        keys: ["dosage", "dose", "concentration", "mg/ml", "apply how much"],
        ans: function () {
          return "Dosing is mapped to the study's tested topical ladder 5–80 mg/mL. The Molecule Lab estimates a working band from predicted potency; low bands suit repellency screens, upper bands suit mortality assays. Confirm with probit before field use.";
        }
      },
      {
        keys: ["caveat", "limitation", "reliable", "trust", "risk", "hazard", "toxic"],
        ans: function () {
          return "Honest caveats: (1) reported values are paper-derived aggregate means — replicate counts pending; (2) docking/affinity here are lightweight surrogates, not validated AutoDock runs; (3) QML is proof-of-concept; (4) botanical extracts can irritate skin — wear PPE; (5) always cross-check pesticide legality in your jurisdiction.";
        }
      },
      {
        keys: ["source", "dataset", "where does the data", "provenance", "cake"],
        ans: function () {
          return "Data come from the study's benchmark CSVs (12 tables: phytochemistry, larval/adult mortality, LC50/LC90, docking, AChE QA, validation) + the Quantum_AgroX architecture doc. See the Literature tab for the verification ledger and external citation links.";
        }
      },
      {
        keys: ["who", "team", "author", "medxai", "project"],
        ans: function () {
          return "Quantum_AgroX is a MEDxAI innovation project applying quantum machine learning to botanical acaricide discovery, benchmarked on Commiphora swynnertonii against Rhipicephalus microplus and R. decoloratus (cattle ticks).";
        }
      }
    ],
    fallback: function (text) {
      var words = text.toLowerCase();
      if (/\d/.test(words)) return "I don't see a matching question — try asking about mortality, docking, LC50, synergy, dosage, or how to run the pipeline. (I reason over the loaded benchmark tables, not free text from the web.)";
      return "I can help with statistics on this dataset and explain the QML methodology. Try: “best extract”, “top docking hit”, “average mortality”, “synergy meaning”, or “caveats”.";
    }
  };
  function chatReply(text) {
    var words = text.toLowerCase();
    var hits = CHAT.replies.filter(function (r) {
      var score = 0;
      r.keys.forEach(function (k) { if (words.indexOf(k) >= 0) score += 1; });
      return score > 0 && (!r.when || r.when(words));
    });
    hits.sort(function (a, b) { return b.keys.filter(function (k) { return words.indexOf(k) >= 0; }).length -
      a.keys.filter(function (k) { return words.indexOf(k) >= 0; }).length; });
    return hits.length ? hits[0].ans() : CHAT.fallback(words);
  }
  function chatSend() {
    var inp = $("#chat-in");
    var text = inp.value.trim();
    if (!text) return;
    inp.value = "";
    var box = $("#chat-msgs");
    var u = el("div", "msg user", esc(text));
    box.appendChild(u);
    var r = el("div", "msg bot", chatReply(text));
    box.appendChild(r);
    box.scrollTop = box.scrollHeight;
  }

  /* ---------- Navigation & init ---------- */
  function showTab(name) {
    APP.tab = name;
    $$(".tablink").forEach(function (b) { b.classList.toggle("active", b.dataset.tab === name); });
    $$(".tabpane").forEach(function (p) { p.classList.toggle("active", p.id === "pane-" + name); });
    if (name === "lit") { litRender(); litVerificationPanel(); }
    if (name === "doc") { goQml(); }
  }
  function goQml() {
    var b = APP.data && APP.data.baseline_ridge;
    if (b) {
      $("#b-mae").textContent = round(b.MAE_mean, 3) + " ± " + round(b.MAE_std, 3);
      $("#b-rmse").textContent = round(b.RMSE_mean, 3) + " ± " + round(b.RMSE_std, 3);
      $("#b-r2").textContent = round(b.R2_mean, 3) + " ± " + round(b.R2_std, 3);
      $("#b-n").textContent = b.n_samples;
    } else {
      ["#b-mae", "#b-rmse", "#b-r2", "#b-n"].forEach(function (s) { $(s).textContent = "—"; });
    }
  }
  function buildPilots() {
    var list = ["Best extract", "Top docking hit", "Average mortality", "LC50", "What is QML?", "Synergy meaning", "How to run pipeline", "Caveats"];
    var wrap = $("#pilots");
    list.forEach(function (p) {
      var c = el("span", "pilot", p);
      c.addEventListener("click", function () {
        $("#chat-in").value = p;
        chatSend();
      });
      wrap.appendChild(c);
    });
    var box = $("#chat-msgs");
    var welcome = el("div", "msg bot", "Hi — I'm the Quantum_AgroX reasoning assistant. Ask me for statistics on the benchmark tables or about the QML methodology. Try a pilot question above.");
    box.appendChild(welcome);
  }
  function initTabs() {
    $$(".tablink").forEach(function (b) {
      b.addEventListener("click", function () {
        showTab(b.dataset.tab);
        var h = { home: "#home", lab: "#lab", lit: "#lit", doc: "#doc", chat: "#chat" }[b.dataset.tab];
        try { history.replaceState(null, "", h); } catch (e) {}
      });
    });
  }
  function deepLinkInit() {
    if (typeof location === "undefined" || typeof URLSearchParams === "undefined") return;
    /* share/deep-link: #lab + ?s1=&s2= auto-runs a pair */
    var h = (location.hash || "#home").replace(/^#/, "");
    if (h === "lab" || h === "lit" || h === "doc" || h === "chat") { showTab(h); } else { showTab("home"); }
    var params = new URLSearchParams(location.search);
    var s1 = params.get("s1"), s2 = params.get("s2");
    if (s1) $("#smi-a").value = s1;
    if (s2) $("#smi-b").value = s2;
    if (s1 && s2) {
      loadRDKit().then(function () {
        labAnalyze();
      })["catch"](function (e) { $("#lab-out").innerHTML = "<div class='warn'>" + esc(e.message) + "</div>"; });
    }
  }

  function register() {
    initTabs();
    buildPilots();
    $("#analyze-btn").addEventListener("click", function () { labAnalyze(); });
    $("#lit-q").addEventListener("input", function () { litRender(); });
    $("#smi-a").addEventListener("keydown", function (e) { if ((e.metaKey || e.ctrlKey) && e.key === "Enter") labAnalyze(); });
    $("#smi-b").addEventListener("keydown", function (e) { if ((e.metaKey || e.ctrlKey) && e.key === "Enter") labAnalyze(); });
    $("#chat-send").addEventListener("click", chatSend);
    $("#chat-in").addEventListener("keydown", function (e) { if (e.key === "Enter") chatSend(); });
    $("#smi-a").value = "CC1=CCC2CC1C2(C)C";   /* alpha-pinene */
    $("#smi-b").value = "COc1cc(CC=C)ccc1O";   /* eugenol */
  }

  function loadData() {
    return fetch("assets/dashboard_data.json").then(function (r) { return r.json(); }).then(function (d) {
      APP.data = d;
      homeRender();
      if (APP.tab === "lit") { litRender(); litVerificationPanel(); }
      if (APP.tab === "doc") { goQml(); }
    });
  }

  function init() {
    loadData().catch(function (e) { console.error("data load failed", e); $("#stat-records").textContent = "—"; });
    register();
    deepLinkInit();
    var rdEl = document.getElementById("rdkit-status");
    loadRDKit().then(function () {
      rdEl.textContent = "RDKit core ready — molecules will render as 2D structures.";
      rdEl.classList.add("ok");
    })["catch"](function (e) {
      rdEl.textContent = "RDKit core failed to load (" + e.message + ") — structure rendering disabled.";
    });
  }

  document.addEventListener("DOMContentLoaded", init);
})();