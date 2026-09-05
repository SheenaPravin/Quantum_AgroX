/* QXSketch — draw a molecule on a hexagonal lattice in the browser,
   export it as SMILES via a spanning-tree SMILES writer, then validate
   and canonicalize with the vendored RDKit (WebAssembly). */

(function (global) {
  "use strict";

  var CELL = 34, COLS = 10, ROWS = 8;
  var S = {};

  var ELEMENTS = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Se"];
  var VALENCE = { C: 4, N: 3, O: 2, S: 2, P: 3, F: 1, Cl: 1, Br: 1, I: 1, B: 3, Se: 2 };

  function K(r, c) { return r + "," + c; }
  function PK(k) { var p = k.split(","); return { r: +p[0], c: +p[1] }; }

  function center(r, c) {
    var x = CELL + c * (CELL * 1.72) + (r % 2 ? CELL * 0.86 : 0);
    var y = CELL + r * CELL * 1.55;
    return { x: x, y: y };
  }

  function ringAt(r, c) {
    var o = [];
    var cc = r % 2 ? [c, c + 1] : [c - 1, c];
    for (var i = 0; i < 6; i++) {
      var dr = Math.round(Math.sin(i * Math.PI / 3));
      var dc = Math.round(Math.cos(i * Math.PI / 3));
      o.push([r + dr, c + dc]);
    }
    return o;
  }
  function hexNbrs(r, c) {
    var n = [];
    for (var i = 0; i < 6; i++) {
      var dr = Math.round(Math.sin(i * Math.PI / 3));
      var dc = Math.round(Math.cos(i * Math.PI / 3));
      n.push([r + dr, c + dc]);
    }
    return n;
  }
  function adjacent(a, b) {
    var n = hexNbrs(a.r, a.c);
    for (var i = 0; i < n.length; i++) if (n[i][0] === b.r && n[i][1] === b.c) return true;
    return false;
  }

  var state = {
    atoms: {}, bonds: {}, tool: "atom", el: "C", order: 1,
    sel1: null, seed: 100
  };

  function atomHas(r, c) { return state.atoms[K(r, c)] !== undefined; }
  function bondKey(a, b) { return a < b ? a + "~" + b : b + "~" + a; }

  function addAtom(r, c, el) {
    el = el || "C";
    if (atomHas(r, c)) { state.atoms[K(r, c)].el = el; return; }
    state.atoms[K(r, c)] = { r: r, c: c, el: el, id: state.seed++ };
  }
  function rmAtom(r, c) {
    var k = K(r, c);
    delete state.atoms[k];
    Object.keys(state.bonds).forEach(function (bk) {
      var p = bk.split("~");
      if (p[0] === k || p[1] === k) delete state.bonds[bk];
    });
  }
  function setBond(a, b, order) {
    if (order === 0) { delete state.bonds[bondKey(a, b)]; return; }
    state.bonds[bondKey(a, b)] = { order: order };
  }
  function clearState() { state.atoms = {}; state.bonds = {}; }

  /* ---------------- SMILES writer ---------------- */

  function collectGraph() {
    var ids = [];
    var el = {}, nbr = {}, edge = {}, edgeOf = {};
    var eid = 0;
    Object.keys(state.atoms).sort().forEach(function (k) {
      var id = state.atoms[k].id;
      ids.push(id);
      el[id] = state.atoms[k].el;
      nbr[id] = nbr[id] || [];
    });
    var idOfKey = {};
    Object.keys(state.atoms).forEach(function (k) { idOfKey[k] = state.atoms[k].id; });
    Object.keys(state.bonds).forEach(function (bk) {
      var p = bk.split("~");
      if (!idOfKey[p[0]] || !idOfKey[p[1]]) return;
      var a = idOfKey[p[0]], b = idOfKey[p[1]];
      var o = state.bonds[bk].order;
      nbr[a].push({ to: b, o: o, id: eid });
      nbr[b].push({ to: a, o: o, id: eid });
      edge[eid] = [a, b, o];
      edgeOf[eid] = bk;
      eid++;
    });
    return { ids: ids, el: el, nbr: nbr, edge: edge, edgeOf: edgeOf };
  }

  function componentOf(g, start, visit) {
    var stack = [start];
    while (stack.length) {
      var v = stack.pop();
      if (visit[v]) continue;
      visit[v] = true;
      g.nbr[v].forEach(function (e) { if (!visit[e.to]) stack.push(e.to); });
    }
  }

  function dfsAssign(g, root) {
    var pre = [], idx = {}, parent = {}, ringD = {}, seen = {};
    var child = {};
    var nextDigit = 1;
    function dfs(v, parEdge) {
      seen[v] = true;
      idx[v] = pre.length;
      pre.push(v);
      parent[v] = parEdge;
      child[v] = [];
      g.nbr[v].forEach(function (e) {
        if (e.id === parEdge) return;
        if (!seen[e.to]) {
          child[v].push(e);
          dfs(e.to, e.id);
        } else {
          if (ringD[e.id] !== undefined) return;
          if (idx[e.to] < idx[v]) ringD[e.id] = nextDigit++;
        }
      });
    }
    dfs(root, -1);
    return { pre: pre, idx: idx, parent: parent, ringD: ringD, child: child };
  }

  var BONDCH = { 1: "", 2: "=", 3: "#" };

  function fmtDigit(d) { return d < 10 ? String(d) : "%" + d; }

  function emit(g, t, v) {
    var s = g.el[v];
    g.nbr[v].forEach(function (e) {
      if (t.ringD[e.id] === undefined) return;
      s += (BONDCH[e.o] || "") + fmtDigit(t.ringD[e.id]);
    });
    var children = t.child[v] || [];
    children.forEach(function (e, i) {
      if (i === 0) { s += BONDCH[e.o]; s += emit(g, t, e.to); }
      else { s += "(" + BONDCH[e.o] + emit(g, t, e.to) + ")"; }
    });
    return s;
  }

  function graphToSmiles() {
    var g = collectGraph();
    if (!g.ids.length) return "";
    var visit = {};
    var out = [];
    g.ids.slice().sort(function (a, b) { return a - b; }).forEach(function (id) {
      if (visit[id]) return;
      var comp = [];
      var stack = [id];
      visit[id] = true;
      while (stack.length) {
        var v = stack.pop();
        comp.push(v);
        g.nbr[v].forEach(function (e) { if (!visit[e.to]) { visit[e.to] = true; stack.push(e.to); } });
      }
      var nbrs = g.nbr;
      var root = comp.reduce(function (best, v) { return (nbrs[v].length > nbrs[best].length) ? v : best; }, comp[0]);
      var t = dfsAssign(g, root);
      out.push(emit(g, t, root));
    });
    return out.join(".");
  }

  /* ---------------- RDKit bridge ---------------- */

  var rdkitReady = function () { return !!(global.RDKIT); };

  function validate(smiles, altr) {
    if (!rdkitReady()) return { ok: false, msg: "RDKit core not ready", smiles: smiles };
    try {
      var m = global.RDKIT.get_mol(smiles);
      if (m && m.is_valid()) {
        var canon = m.get_smiles();
        var svg = m.get_svg(360, 260);
        m.delete();
        return { ok: true, canon: canon, svg: svg, smiles: smiles };
      }
      var hint = "";
      var ats = Object.keys(state.atoms).length;
      var bds = Object.keys(state.bonds).length;
      if (ats && bds === 0 && ats > 1) hint = "Atoms are placed but no bonds connect them yet.";
      else hint = "Valence check failed — check double-bond positions in rings, or reduce atom valency.";
      return { ok: false, msg: hint, smiles: smiles };
    } catch (e) {
      return { ok: false, msg: "RDKit rejected the structure: " + e.message, smiles: smiles };
    }
  }

  /* ---------------- rendering ---------------- */

  function draw() {
    var wrap = document.getElementById("qxs-grid");
    if (!wrap) return;
    var W = CELL * 2 + (COLS - 1) * CELL * 1.72 + CELL * 0.86 + 24;
    var H = CELL * 2 + (ROWS - 1) * CELL * 1.55 + 24;
    var svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ' + W + ' ' + H + '" width="100%" height="' + H + '">';
    svg += '<defs><radialGradient id="qxs-glow" cx="50%" cy="50%" r="50%">'
      + '<stop offset="0%" stop-color="#2ee6a8" stop-opacity="0.16"/>'
      + '<stop offset="100%" stop-color="#2ee6a8" stop-opacity="0"/>'
      + '</radialGradient></defs>';
    svg += '<rect x="0" y="0" width="' + W + '" height="' + H + '" fill="#071210" rx="14"/>';
    for (var r = 0; r < ROWS; r++) {
      for (var c = 0; c < COLS; c++) {
        var p = center(r, c);
        svg += '<circle cx="' + p.x.toFixed(1) + '" cy="' + p.y.toFixed(1) + '" r="2" fill="#123c2e" opacity="0.55"/>';
      }
    }
    Object.keys(state.bonds).forEach(function (bk) {
      var p = bk.split("~");
      var a = PK(p[0]), b = PK(p[1]);
      var pa = center(a.r, a.c), pb = center(b.r, b.c);
      var o = state.bonds[bk].order;
      svg += bondSvg(pa, pb, o);
    });
    Object.keys(state.atoms).forEach(function (k) {
      var a = state.atoms[k];
      var p = center(a.r, a.c);
      svg += '<circle cx="' + p.x.toFixed(1) + '" cy="' + p.y.toFixed(1) + '" r="14" fill="url(#qxs-glow)"/>';
      var s = a.el;
      svg += '<text x="' + p.x.toFixed(1) + '" y="' + (p.y + 5).toFixed(1) + '" text-anchor="middle" font-family="ui-monospace,Menlo,Consolas,monospace" font-size="14" font-weight="650" fill="#d7f8e9">' + s + '</text>';
    });
    Object.keys(state.atoms).forEach(function (k) {
      var a = state.atoms[k];
      var p = center(a.r, a.c);
      svg += '<circle cx="' + p.x.toFixed(1) + '" cy="' + p.y.toFixed(1) + '" r="17" fill="transparent" stroke="transparent" data-key="atom:' + k + '"/>';
    });
    Object.keys(state.bonds).forEach(function (bk) {
      var p = bk.split("~");
      var a = PK(p[0]), b = PK(p[1]);
      var pa = center(a.r, a.c), pb = center(b.r, b.c);
      var mx = (pa.x + pb.x) / 2, my = (pa.y + pb.y) / 2;
      svg += '<circle cx="' + mx.toFixed(1) + '" cy="' + my.toFixed(1) + '" r="10" fill="transparent" stroke="transparent" data-key="bond:' + bk + '"/>';
    });
    svg += "</svg>";
    wrap.innerHTML = svg;
  }

  function bondSvg(a, b, o) {
    var dx = b.x - a.x, dy = b.y - a.y;
    var len = Math.sqrt(dx * dx + dy * dy) || 1;
    var px = -dy / len, py = dx / len;
    var off = 4;
    var s = "";
    if (o === 1) {
      s += '<line x1="' + a.x.toFixed(1) + '" y1="' + a.y.toFixed(1) + '" x2="' + b.x.toFixed(1) + '" y2="' + b.y.toFixed(1) + '" stroke="#3ff0b4" stroke-width="2.6" stroke-linecap="round"/>';
    } else if (o === 2) {
      s += '<line x1="' + (a.x + px * off).toFixed(1) + '" y1="' + (a.y + py * off).toFixed(1) + '" x2="' + (b.x + px * off).toFixed(1) + '" y2="' + (b.y + py * off).toFixed(1) + '" stroke="#3ff0b4" stroke-width="2.4" stroke-linecap="round"/>';
      s += '<line x1="' + (a.x - px * off).toFixed(1) + '" y1="' + (a.y - py * off).toFixed(1) + '" x2="' + (b.x - px * off).toFixed(1) + '" y2="' + (b.y - py * off).toFixed(1) + '" stroke="#3ff0b4" stroke-width="2.4" stroke-linecap="round"/>';
    } else {
      s += '<line x1="' + (a.x + px * 2 * off).toFixed(1) + '" y1="' + (a.y + py * 2 * off).toFixed(1) + '" x2="' + (b.x + px * 2 * off).toFixed(1) + '" y2="' + (b.y + py * 2 * off).toFixed(1) + '" stroke="#3ff0b4" stroke-width="2.2" stroke-linecap="round"/>';
      s += '<line x1="' + a.x.toFixed(1) + '" y1="' + a.y.toFixed(1) + '" x2="' + b.x.toFixed(1) + '" y2="' + b.y.toFixed(1) + '" stroke="#3ff0b4" stroke-width="2.2"/>';
      s += '<line x1="' + (a.x - px * 2 * off).toFixed(1) + '" y1="' + (a.y - py * 2 * off).toFixed(1) + '" x2="' + (b.x - px * 2 * off).toFixed(1) + '" y2="' + (b.y - py * 2 * off).toFixed(1) + '" stroke="#3ff0b4" stroke-width="2.2" stroke-linecap="round"/>';
    }
    return s;
  }

  function svgToPt(e, wrap) {
    var svgEl = wrap.querySelector("svg");
    var vb = svgEl.viewBox.baseVal;
    var r = svgEl.getBoundingClientRect();
    var x = (e.clientX - r.left) / r.width * vb.width;
    var y = (e.clientY - r.top) / r.height * vb.height;
    var best = null, bestD = 22;
    for (var cr = 0; cr < ROWS; cr++) {
      for (var cc = 0; cc < COLS; cc++) {
        var p = center(cr, cc);
        var d = Math.hypot(x - p.x, y - p.y);
        if (d < bestD) { bestD = d; best = [cr, cc]; }
      }
    }
    return best;
  }

  function onClick(e) {
    var wrap = document.getElementById("qxs-grid");
    if (!wrap) return;
    var t = e.target;
    var dk = t && t.getAttribute && t.getAttribute("data-key");
    var mode = state.tool;
    if (dk && dk.indexOf("bond:") === 0) {
      var bkv = dk.slice(5);
      var parts = bkv.split("~");
      if (mode === "erase") { delete state.bonds[bkv]; }
      else {
        var cur = state.bonds[bkv] ? state.bonds[bkv].order : 1;
        var nxt = (cur % 3) + 1;
        if (cur === 1 && state.bonds[bkv] === undefined) nxt = 1;
        setBond(parts[0], parts[1], nxt);
      }
      setStatus("Bond order toggled — click again for single → double → triple → single.");
      draw();
      return;
    }
    if (dk && dk.indexOf("atom:") === 0) {
      var kk = dk.slice(5);
      var aa = PK(kk);
      if (mode === "erase") { rmAtom(aa.r, aa.c); setStatus("Atom removed."); }
      else if (mode === "atom") {
        addAtom(aa.r, aa.c, state.el);
        setStatus("Atom set to " + state.el + ".");
      } else if (mode === "bond") {
        if (!state.sel1) { state.sel1 = kk; setStatus("First atom selected — pick a neighbour to bond."); }
        else {
          var kn = buildBondIfAdjacent(state.sel1, kk, state.order);
          state.sel1 = null;
          if (kn) setStatus("Bond drawn (order " + state.order + ").");
          else setStatus("Choose an adjacent atom — atoms in different chains must touch first.");
        }
      }
      draw();
      return;
    }
    var cell = svgToPt(e, wrap);
    if (!cell) return;
    var cr = cell[0], cc = cell[1];
    if (atomHas(cr, cc)) {
      if (mode === "atom") { addAtom(cr, cc, state.el); setStatus("Atom changed to " + state.el + "."); draw(); }
      return;
    }
    if (mode === "atom") {
      addAtom(cr, cc, state.el); setStatus("Added " + state.el + "."); draw();
    } else if (mode === "bond") {
      setStatus("Click an existing atom first, then an empty neighbour.");
    } else if (mode === "erase") {
      setStatus("Nothing to erase there.");
    }
  }

  function buildBondIfAdjacent(aKey, bKey, order) {
    var a = PK(aKey), b = PK(bKey);
    if (!adjacent(a, b)) return false;
    setBond(aKey, bKey, order);
    return true;
  }

  function placeTemplate(name, centerRC) {
    var c = centerRC || [3, 4];
    var cr = c[0], cc = c[1];
    var axis;
    if (name === "benzene" || name === "pyridine" || name === "cyclohexane") {
      var k = "C";
      var elRing = [];
      var start = cr % 2 ? cc + 1 : cc;
      var pos = [
        [cr, start], [cr, start + 1],
        [cr + 1, start + 1], [cr + 1, start],
        [cr + 2, start], [cr + 2, start + 1]
      ];
      for (var i = 0; i < pos.length; i++) addAtom(pos[i][0], pos[i][1], k);
      var ring = [K(pos[0][0], pos[0][1]), K(pos[1][0], pos[1][1]), K(pos[2][0], pos[2][1]), K(pos[3][0], pos[3][1]), K(pos[4][0], pos[4][1]), K(pos[5][0], pos[5][1])];
      for (var i2 = 0; i2 < 6; i2++) {
        var e1 = ring[i2].split(","), e2 = ring[(i2 + 1) % 6].split(",");
        setBond(e1[0] + "," + e1[1], e2[0] + "," + e2[1], (name === "cyclohexane" ? 1 : (i2 % 2 === 0 ? 2 : 1)));
      }
      if (name === "pyridine") {
        var tmp = state.atoms[ring[0]]; if (tmp) tmp.el = "N";
      }
    } else if (name === "furan") {
      var pos5 = [[cr, cc], [cr, cc + 1], [cr + 1, cc + 1], [cr + 1, cc], [cr - 1, cc]];
      var fos = [K(pos5[0][0], pos5[0][1]), K(pos5[1][0], pos5[1][1]), K(pos5[2][0], pos5[2][1]), K(pos5[3][0], pos5[3][1]), K(pos5[4][0], pos5[4][1])];
      for (var i3 = 0; i3 < 5; i3++) addAtom(pos5[i3][0], pos5[i3][1], "C");
      for (var i4 = 0; i4 < 5; i4++) {
        var p1 = fos[i4].split(","), p2 = fos[(i4 + 1) % 5].split(",");
        setBond(p1[0] + "," + p1[1], p2[0] + "," + p2[1], (i4 === 0 || i4 === 4 ? 1 : 2));
      }
      var otmp = state.atoms[fos[4]]; if (otmp) otmp.el = "O";
    }
    draw();
  }

  function setStatus(m) {
    var el = document.getElementById("qxs-status");
    if (el) el.textContent = m;
  }

  /* ---------------- public API ---------------- */

  function open(aMolA, aMolB) {
    state.sel1 = null;
    var bg = document.getElementById("qxs-backdrop");
    if (bg) bg.hidden = false;
    var m = document.getElementById("qxs-modal");
    if (m) m.hidden = false;
    if (aMolA !== undefined) state.aMolA = aMolA;
    if (aMolB !== undefined) state.aMolB = aMolB;
    draw();
    updateStatus("Select an atom element or ring template, then click cells to assemble your molecule.");
  }

  function close() {
    var bg = document.getElementById("qxs-backdrop");
    if (bg) bg.hidden = true;
    var m = document.getElementById("qxs-modal");
    if (m) m.hidden = true;
  }

  function exportSmiles() {
    if (!Object.keys(state.atoms).length) return Promise.resolve({ ok: false, msg: "Canvas is empty — place atoms first." });
    var raw = graphToSmiles();
    var res = validate(raw);
    if (res.ok) {
      var out = document.getElementById("qxs-out-smiles");
      var box = document.getElementById("qxs-preview");
      if (out && box) {
        out.textContent = res.canon;
        box.innerHTML = res.svg;
        box.style.display = "";
      }
      setStatus("Valid structure — canonical SMILES produced.");
    } else {
      var out2 = document.getElementById("qxs-out-smiles");
      var box2 = document.getElementById("qxs-preview");
      if (out2 && box2) {
        out2.textContent = raw;
        box2.innerHTML = '<p class="qxs-err">' + res.msg + "</p>";
        box2.style.display = "";
      }
      setStatus(res.msg);
    }
    return Promise.resolve({ ok: res.ok, smiles: raw, canon: res.canon, msg: res.msg });
  }

  function useAs(which) {
    var out = document.getElementById("qxs-out-smiles");
    var smi = out ? out.textContent : "";
    smi = (smi || "").trim();
    if (!smi || /^Invalid/.test(smi)) { setStatus("Export a valid structure before using it."); return; }
    if (state.onUse) state.onUse(which, smi);
    close();
  }

  /* exposed test helpers */
  function mockDraw(spec) {
    clearState();
    var atoms = Array.isArray(spec) ? spec : spec.atoms || [];
    atoms.forEach(function (x) { addAtom(x.r, x.c, x.el || "C"); });
    var bonds = (Array.isArray(spec) ? spec.bonds : spec.bonds) || [];
    bonds.forEach(function (b) { setBond(b[0], b[1], b[2] || 1); });
    draw();
  }

  function ready(cb) { if (cb) cb(); }

  S.open = open; S.close = close; S.exportSmiles = exportSmiles; S.useAs = useAs;
  S.onUse = function (f) { state.onUse = f; };
  S.mockDraw = mockDraw; S.ready = ready;

  /* init DOM wiring */
  function init() {
    var wrap = document.getElementById("qxs-grid");
    if (!wrap) return;
    wrap.addEventListener("click", onClick);
    ["qxs-atom-C", "qxs-atom-N", "qxs-atom-O", "qxs-atom-S", "qxs-atom-P", "qxs-atom-F", "qxs-atom-Cl", "qxs-atom-Br", "qxs-atom-I", "qxs-atom-B", "qxs-atom-Se"].forEach(function (id) {
      var el = document.getElementById(id);
      if (el) el.addEventListener("click", function () {
        state.tool = "atom"; state.el = id.slice(9);
        setActive("qxs-tools", "bk"); refreshToolUI();
      });
    });
    ["qxs-bond1", "qxs-bond2", "qxs-bond3"].forEach(function (id) {
      var el = document.getElementById(id);
      if (el) el.addEventListener("click", function () {
        state.tool = "bond"; state.order = +id.slice(8); refreshToolUI();
      });
    });
    var erase = document.getElementById("qxs-erase");
    if (erase) erase.addEventListener("click", function () { state.tool = "erase"; refreshToolUI(); });
    var clear = document.getElementById("qxs-clear");
    if (clear) clear.addEventListener("click", function () { clearState(); setStatus("Canvas cleared."); draw(); });
    ["qxs-tpl-benzene", "qxs-tpl-pyridine", "qxs-tpl-furan", "qxs-tpl-cyclohexane"].forEach(function (id) {
      var el = document.getElementById(id);
      if (el) el.addEventListener("click", function () { placeTemplate(id.slice(8)); });
    });
    var exp = document.getElementById("qxs-export");
    if (exp) exp.addEventListener("click", function () { exportSmiles(); });
    var useA = document.getElementById("qxs-use-a"), useB = document.getElementById("qxs-use-b");
    if (useA) useA.addEventListener("click", function () { useAs("A"); });
    if (useB) useB.addEventListener("click", function () { useAs("B"); });
    var closeX = document.getElementById("qxs-close");
    if (closeX) closeX.addEventListener("click", close);
    var bg = document.getElementById("qxs-backdrop");
    if (bg) bg.addEventListener("click", function (e) { if (e.target === bg) close(); });
    document.addEventListener("keydown", function (e) {
      if ((e.key === "Escape" || e.key === "Esc") && document.getElementById("qxs-modal") && !document.getElementById("qxs-modal").hidden) close();
    });
    refreshToolUI();
  }

  function refreshToolUI() {
    var atomEls = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Se"];
    atomEls.forEach(function (el) {
      var b = document.getElementById("qxs-atom-" + el);
      if (b) b.classList.toggle("active-tool", state.tool === "atom" && state.el === el);
    });
    [1, 2, 3].forEach(function (o) {
      var b = document.getElementById("qxs-bond" + o);
      if (b) b.classList.toggle("active-tool", state.tool === "bond" && state.order === o);
    });
    var eb = document.getElementById("qxs-erase");
    if (eb) eb.classList.toggle("active-tool", state.tool === "erase");
  }
  function setActive() {}
  function updateStatus(m) { setStatus(m); }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }

  global.QXSketch = S;

  /* headless self-test entry: ?drawtest=1 */
  (function () {
    try {
      var q = (typeof location !== "undefined" && location.search) || "";
      if (q.indexOf("drawtest") < 0) return;
      var tries = 0;
      var t = setInterval(function () {
        tries++;
        if (!global.RDKIT) { if (tries > 60) { clearInterval(t); setStatus("DRAWTEST|TIMEOUT|"); } return; }
        clearInterval(t);
        mockDraw({
          atoms: [
            { r: 3, c: 3, el: "C" }, { r: 3, c: 4, el: "C" }, { r: 4, c: 4, el: "C" },
            { r: 4, c: 3, el: "C" }, { r: 5, c: 3, el: "C" }, { r: 5, c: 4, el: "C" }
          ],
          bonds: [
            ["3,3", "3,4", 2], ["3,4", "4,4", 1], ["4,4", "4,3", 2],
            ["4,3", "5,3", 1], ["5,3", "5,4", 2], ["5,4", "3,3", 1]
          ]
        });
        exportSmiles().then(function (r) {
          var out = document.getElementById("qxs-out-smiles");
          var txt = "DRAWTEST|" + (r.ok ? "OK" : "ERR") + "|" + (r.canon || r.smiles || "");
          if (out) out.textContent = txt;
          setStatus(txt);
          if (q.indexOf("drawdemo") >= 0) open();
        });
      }, 250);
    } catch (e) {
      try { setStatus("DRAWTEST|ERR|" + e.message); } catch (e2) {}
    }
  })();
})(window);