/* Rough-path-driven adjoint systems and gradient-based data assimilation.

   Forward (rough pendulum, Hamiltonian H0 = p²/2 − cos q, H1 = −σ cos q):
     dq = p dt,   dp = −sin q (dt + σ dZ_t),
   where Z is a sample path of fractional Brownian motion with Hurst index H.
   Write one step as x' = x + F(x) with F(x) = (p h, −sin q (h + σ ΔZ)).

   Loss (4D-Var): J(x₀) = ½ Σ_k |x_{n_k} − y_k|² over noisy observations y_k.
   Adjoint: λ_N = r_N,  λ_n = (adjoint step)(λ_{n+1}) + r_n [if n is an observation time],
   gradient ∇J = λ₀.

   Two forward/adjoint pairings:
   • RSPRK (Gauss, one stage): implicit midpoint forward; its symplectic partner backward,
       (I − Aᵀ/2) λ_n = (I + Aᵀ/2) λ_{n+1},  A = DF(x̄_n).
     This is exactly the transpose of the forward Jacobian, so the adjoint
     quasi-conservation law holds to round-off and λ₀ is the exact gradient of the discrete J.
   • Heun / Heun: explicit Heun forward; the continuous adjoint equation integrated backward
     with Heun. Not a symplectic pair, so the quasi-conservation law fails and the gradient is biased. */
(function () {
  'use strict';
  var D = window.Demo;
  var TF = 6, NOBS = 10, MAX_ITER = 200; // gradient descent stops after MAX_ITER steps
  var state = {
    H: 0.4, sigma: 0.8, N: 60, lr: 0.3, seed: 7, running: true,
    truth0: [1.3, 0.2], guess0: [0.3, 1.0], obsNoise: 0.05
  };

  // ---------- fractional Gaussian noise by Cholesky ----------
  function fgn(N, H, h, rng) {
    var g = function (k) { k = Math.abs(k); return 0.5 * Math.pow(h, 2 * H) * (Math.pow(k + 1, 2 * H) - 2 * Math.pow(k, 2 * H) + Math.pow(Math.abs(k - 1), 2 * H)); };
    var Lm = [];
    for (var i = 0; i < N; i++) {
      Lm.push(new Float64Array(N));
      for (var j = 0; j <= i; j++) {
        var s = g(i - j);
        for (var k = 0; k < j; k++) s -= Lm[i][k] * Lm[j][k];
        Lm[i][j] = i === j ? Math.sqrt(Math.max(s, 1e-300)) : s / Lm[j][j];
      }
    }
    var z = [], out = new Float64Array(N);
    for (var a = 0; a < N; a++) z.push(rng.normal());
    for (var b = 0; b < N; b++) { var v = 0; for (var c = 0; c <= b; c++) v += Lm[b][c] * z[c]; out[b] = v; }
    return out;
  }

  // ---------- vector field and its Jacobian for one step ----------
  function F(x, h, dz) { return [x[1] * h, -Math.sin(x[0]) * (h + state.sigma * dz)]; }
  function DF(x, h, dz) { return [0, h, -Math.cos(x[0]) * (h + state.sigma * dz), 0]; } // row-major 2×2
  function mv(A, v) { return [A[0] * v[0] + A[1] * v[1], A[2] * v[0] + A[3] * v[1]]; }
  function mtv(A, v) { return [A[0] * v[0] + A[2] * v[1], A[1] * v[0] + A[3] * v[1]]; }
  function solve2(M, b) { var det = M[0] * M[3] - M[1] * M[2]; return [(M[3] * b[0] - M[1] * b[1]) / det, (M[0] * b[1] - M[2] * b[0]) / det]; }
  function dot(a, b) { return a[0] * b[0] + a[1] * b[1]; }

  // implicit midpoint step by Newton's method
  function midStep(x, h, dz) {
    var y = [x[0] + x[1] * h, x[1] - Math.sin(x[0]) * (h + state.sigma * dz)];
    for (var it = 0; it < 30; it++) {
      var m = [(x[0] + y[0]) / 2, (x[1] + y[1]) / 2], f = F(m, h, dz), A = DF(m, h, dz);
      var G = [y[0] - x[0] - f[0], y[1] - x[1] - f[1]];
      var J = [1 - A[0] / 2, -A[1] / 2, -A[2] / 2, 1 - A[3] / 2];
      var d = solve2(J, G);
      y = [y[0] - d[0], y[1] - d[1]];
      if (Math.abs(d[0]) + Math.abs(d[1]) < 1e-16) break;
    }
    return y;
  }

  var SCHEMES = {
    rsprk: {
      name: 'Symplectic: RSPRK', color: D.colors.sp,
      forward: function (xs, dZ, h) {
        for (var n = 0; n < dZ.length; n++) xs.push(midStep(xs[n], h, dZ[n]));
      },
      A: function (xs, n, h, dz) { return DF([(xs[n][0] + xs[n + 1][0]) / 2, (xs[n][1] + xs[n + 1][1]) / 2], h, dz); },
      tangent: function (xs, n, h, dz, v) { // (I − A/2) v' = (I + A/2) v
        var A = this.A(xs, n, h, dz), r = mv(A, v);
        return solve2([1 - A[0] / 2, -A[1] / 2, -A[2] / 2, 1 - A[3] / 2], [v[0] + r[0] / 2, v[1] + r[1] / 2]);
      },
      adjoint: function (xs, n, h, dz, lam) { // (I − Aᵀ/2) λ_n = (I + Aᵀ/2) λ_{n+1}
        var A = this.A(xs, n, h, dz), r = mtv(A, lam);
        return solve2([1 - A[0] / 2, -A[2] / 2, -A[1] / 2, 1 - A[3] / 2], [lam[0] + r[0] / 2, lam[1] + r[1] / 2]);
      }
    },
    heun: {
      name: 'Not symplectic: Heun / Heun', color: D.colors.other,
      forward: function (xs, dZ, h) {
        for (var n = 0; n < dZ.length; n++) {
          var x = xs[n], f1 = F(x, h, dZ[n]), xt = [x[0] + f1[0], x[1] + f1[1]], f2 = F(xt, h, dZ[n]);
          xs.push([x[0] + 0.5 * (f1[0] + f2[0]), x[1] + 0.5 * (f1[1] + f2[1])]);
        }
      },
      tangent: function (xs, n, h, dz, v) { // exact derivative of the Heun map
        var x = xs[n], A0 = DF(x, h, dz), f1 = F(x, h, dz), A1 = DF([x[0] + f1[0], x[1] + f1[1]], h, dz);
        var a0v = mv(A0, v), a1 = mv(A1, [v[0] + a0v[0], v[1] + a0v[1]]);
        return [v[0] + 0.5 * (a0v[0] + a1[0]), v[1] + 0.5 * (a0v[1] + a1[1])];
      },
      adjoint: function (xs, n, h, dz, lam) { // Heun on dλ = −DFᵀλ, backward from t_{n+1} to t_n
        var g1 = mtv(DF(xs[n + 1], h, dz), lam), lt = [lam[0] + g1[0], lam[1] + g1[1]];
        var g2 = mtv(DF(xs[n], h, dz), lt);
        return [lam[0] + 0.5 * (g1[0] + g2[0]), lam[1] + 0.5 * (g1[1] + g2[1])];
      }
    }
  };
  var KEYS = ['rsprk', 'heun'];

  function trajectory(scheme, x0) { var xs = [x0.slice()]; scheme.forward(xs, state.dZ, state.h); return xs; }
  function loss(xs) {
    var J = 0;
    state.obsIdx.forEach(function (n, k) { var r0 = xs[n][0] - state.obs[k][0], r1 = xs[n][1] - state.obs[k][1]; J += 0.5 * (r0 * r0 + r1 * r1); });
    return J;
  }

  // adjoint sweep: gradient λ₀ and the quasi-conservation residual Q_n
  function adjointSweep(scheme, xs) {
    var N = state.N, h = state.h, dZ = state.dZ;
    var resid = {};
    state.obsIdx.forEach(function (n, k) { resid[n] = [xs[n][0] - state.obs[k][0], xs[n][1] - state.obs[k][1]]; });
    // tangent δx_n from δx_0
    var dv = [[1 / Math.sqrt(1.25), 0.5 / Math.sqrt(1.25)]];
    for (var n = 0; n < N; n++) dv.push(scheme.tangent(xs, n, h, dZ[n], dv[n]));
    var lam = new Array(N + 1);
    lam[N] = resid[N] ? resid[N].slice() : [0, 0];
    for (var m = N - 1; m >= 0; m--) {
      lam[m] = scheme.adjoint(xs, m, h, dZ[m], lam[m + 1]);
      if (resid[m] && m > 0) { lam[m][0] += resid[m][0]; lam[m][1] += resid[m][1]; }
    }
    // Q_n = <λ_n, δx_n> − Σ_{k ≥ n, obs} <r_k, δx_k>  (≡ 0 for the exact discrete adjoint)
    var Q = new Array(N + 1), tail = 0, scale = 0;
    for (var j = N; j >= 0; j--) {
      if (resid[j] && j > 0) { var c = dot(resid[j], dv[j]); tail += c; scale += Math.abs(c); }
      Q[j] = dot(lam[j], dv[j]) - tail;
    }
    scale = Math.max(scale, 1e-300);
    var series = [];
    for (var i = 0; i < N; i++) series.push([i * h, Q[i] / scale]); // Q_N = 0 by construction for both schemes
    return { grad: lam[0], Q: series };
  }

  function fdGradient(scheme, x0) {
    var e = 1e-6, g = [0, 0];
    for (var i = 0; i < 2; i++) {
      var a = x0.slice(), b = x0.slice(); a[i] += e; b[i] -= e;
      g[i] = (loss(trajectory(scheme, a)) - loss(trajectory(scheme, b))) / (2 * e);
    }
    return g;
  }

  // ---------- problem set-up ----------
  function setup() {
    var rng = D.makeRng(state.seed);
    state.h = TF / state.N;
    state.dZ = fgn(state.N, state.H, state.h, rng);
    state.Z = [0];
    for (var n = 0; n < state.N; n++) state.Z.push(state.Z[n] + state.dZ[n]);
    // reference truth: midpoint on a 20× finer grid along the piecewise-linear path
    var sub = 20, hs = state.h / sub, x = state.truth0.slice(), fine = [x.slice()];
    state.truthAt = [x.slice()];
    for (var m = 0; m < state.N; m++) {
      for (var s = 0; s < sub; s++) { x = midStep(x, hs, state.dZ[m] / sub); fine.push(x.slice()); }
      state.truthAt.push(x.slice());
    }
    state.truthFine = fine;
    state.obsIdx = []; state.obs = [];
    for (var k = 1; k <= NOBS; k++) {
      var idx = Math.round(k * state.N / NOBS);
      state.obsIdx.push(idx);
      state.obs.push([state.truthAt[idx][0] + state.obsNoise * rng.normal(), state.truthAt[idx][1] + state.obsNoise * rng.normal()]);
    }
    restart();
  }

  function restart() {
    state.iter = 0;
    state.runs = {};
    KEYS.forEach(function (k) {
      state.runs[k] = { alpha: state.lr, stalled: false, x0: state.guess0.slice(), loss: [], gnorm: [], gerr: [], Q: [], xs: null, path: [state.guess0.slice()] };
      evaluate(k);
    });
  }

  function evaluate(k) {
    var sc = SCHEMES[k], run = state.runs[k];
    run.xs = trajectory(sc, run.x0);
    var J = loss(run.xs), adj = adjointSweep(sc, run.xs), fd = fdGradient(sc, run.x0);
    var gn = Math.hypot(fd[0], fd[1]);
    run.grad = adj.grad; run.Q = adj.Q;
    run.loss.push([state.iter, J]);
    run.gnorm.push([state.iter, gn]);
    run.gerr.push([state.iter, Math.hypot(adj.grad[0] - fd[0], adj.grad[1] - fd[1])]);
  }

  // gradient descent with Armijo backtracking on each scheme's own discrete loss,
  // using that scheme's adjoint gradient as the search direction
  function iterate() {
    state.iter++;
    KEYS.forEach(function (k) {
      var run = state.runs[k], sc = SCHEMES[k], g = run.grad, g2 = g[0] * g[0] + g[1] * g[1];
      var J0 = run.loss[run.loss.length - 1][1], alpha = run.alpha;
      var accepted = false;
      for (var tries = 0; tries < 40 && g2 > 0; tries++) {
        var s = alpha, len = s * Math.sqrt(g2);
        if (len > 0.4) s *= 0.4 / len;
        var cand = [run.x0[0] - s * g[0], run.x0[1] - s * g[1]];
        if (loss(trajectory(sc, cand)) <= J0 - 1e-4 * s * g2) { run.x0 = cand; accepted = true; break; }
        alpha *= 0.5;
      }
      run.alpha = accepted ? Math.min(alpha * 1.5, 50) : alpha;
      run.stalled = !accepted;
      run.path.push(run.x0.slice());
      evaluate(k);
    });
  }

  // ---------- drawing ----------
  var phase = document.getElementById('phase');
  var pv = null;
  var plots = {
    loss: new D.TimePlot(document.getElementById('plot-loss'), { log: true, aspect: 0.7, legend: false, endLabels: true, xlabel: 'iteration', xint: true,
      title: 'Data misfit J(x₀)', subtitle: 'each scheme\'s own discrete loss, per descent step' }),
    gnorm: new D.TimePlot(document.getElementById('plot-gnorm'), { log: true, aspect: 0.7, legend: false, endLabels: true, xlabel: 'iteration', xint: true,
      title: 'True gradient |∇J| at the current x₀', subtitle: 'reaches 0 only if the adjoint gradient is exact' }),
    gerr: new D.TimePlot(document.getElementById('plot-gerr'), { log: true, aspect: 0.7, legend: false, endLabels: true, xlabel: 'iteration', xint: true,
      title: 'Adjoint gradient error |λ₀ − ∇J|', subtitle: 'vs central finite differences (floor ≈ 1e-10)' }),
    Q: new D.TimePlot(document.getElementById('plot-q'), { log: true, aspect: 0.7, legend: false, endLabels: false,
      refLines: [{ y: 1e-16, label: 'machine precision' }],
      title: 'Adjoint conservation law |Qₙ|', subtitle: 'along the trajectory; exactly 0 for a symplectic pairing' })
  };

  function resize() { pv = D.fitCanvas(phase, 0.68); Object.keys(plots).forEach(function (k) { plots[k].resize(); }); }

  function bounds() {
    var pts = state.truthFine.concat(state.obs, [state.guess0]);
    KEYS.forEach(function (k) { pts = pts.concat(state.runs[k].xs); });
    var b = [Infinity, -Infinity, Infinity, -Infinity];
    pts.forEach(function (p) { if (!isFinite(p[0])) return; b[0] = Math.min(b[0], p[0]); b[1] = Math.max(b[1], p[0]); b[2] = Math.min(b[2], p[1]); b[3] = Math.max(b[3], p[1]); });
    var px = 0.1 * (b[1] - b[0] || 1), py = 0.1 * (b[3] - b[2] || 1);
    return [b[0] - px, b[1] + px, b[2] - py, b[3] + py];
  }

  function drawPhase() {
    var ctx = pv.ctx, W = pv.w, Hh = pv.h, P = D.ink();
    var items = [
      { c: SCHEMES.rsprk.color, t: 'Symplectic: RSPRK fit', k: 'line' },
      { c: SCHEMES.heun.color, t: 'Not symplectic: Heun / Heun fit', k: 'line' },
      { c: D.colors.truth, t: 'true trajectory, true x₀ (×)', k: 'dash' },
      { c: D.colors.navy, t: 'observations', k: 'dot' },
      { c: P.muted, t: 'descent path of x₀', k: 'dots' }
    ];
    ctx.font = '600 10px Inter, -apple-system, "Segoe UI", Roboto, sans-serif';
    var rows = [], row = [], lx = 40;
    items.forEach(function (it) {
      var w = 26 + ctx.measureText(it.t).width + 16;
      if (row.length && lx + w > W - 8) { rows.push(row); row = []; lx = 40; }
      row.push({ it: it, x: lx }); lx += w;
    });
    rows.push(row);
    var LH = 14 * rows.length + 8;
    var stripH = Math.round(Hh * 0.2), PH = Hh - stripH - 8;
    ctx.clearRect(0, 0, W, Hh);
    ctx.fillStyle = P.bg; ctx.fillRect(0, 0, W, Hh);
    var b = bounds();
    var X = function (q) { return 40 + (W - 50) * (q - b[0]) / (b[1] - b[0]); };
    var Y = function (p) { return LH + 10 + (PH - LH - 20) * (1 - (p - b[2]) / (b[3] - b[2])); };
    ctx.strokeStyle = P.grid; ctx.strokeRect(40, LH + 10, W - 50, PH - LH - 20);
    ctx.fillStyle = P.axis; ctx.font = '11px Inter, -apple-system, "Segoe UI", Roboto, sans-serif';
    ctx.fillText('q (angle)', W - 62, PH - 14); ctx.save(); ctx.translate(14, (PH + LH) / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'center'; ctx.fillText('p (momentum)', 0, 0); ctx.restore();

    var line = function (pts, color, width, dash) {
      ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash || []);
      ctx.beginPath();
      pts.forEach(function (p, i) { if (i === 0) ctx.moveTo(X(p[0]), Y(p[1])); else ctx.lineTo(X(p[0]), Y(p[1])); });
      ctx.stroke(); ctx.setLineDash([]);
    };
    line(state.truthFine, D.colors.truth, 1.5, [5, 4]);
    KEYS.forEach(function (k) {
      var run = state.runs[k];
      line(run.xs, SCHEMES[k].color, 2);
      line(run.path, SCHEMES[k].color, 1, [2, 3]);
      ctx.fillStyle = SCHEMES[k].color;
      ctx.beginPath(); ctx.arc(X(run.x0[0]), Y(run.x0[1]), 4.5, 0, 2 * Math.PI); ctx.fill();
    });
    ctx.fillStyle = D.colors.navy;
    state.obs.forEach(function (o) { ctx.beginPath(); ctx.arc(X(o[0]), Y(o[1]), 3.5, 0, 2 * Math.PI); ctx.fill(); });
    ctx.strokeStyle = D.colors.truth; ctx.lineWidth = 2;
    var t0 = state.truth0, cx = X(t0[0]), cy = Y(t0[1]);
    ctx.beginPath(); ctx.moveTo(cx - 6, cy - 6); ctx.lineTo(cx + 6, cy + 6); ctx.moveTo(cx + 6, cy - 6); ctx.lineTo(cx - 6, cy + 6); ctx.stroke();

    // driving path strip
    var top = PH, zmin = Math.min.apply(null, state.Z), zmax = Math.max.apply(null, state.Z);
    ctx.strokeStyle = P.grid; ctx.strokeRect(40, top, W - 50, stripH);
    ctx.strokeStyle = D.colors.navy; ctx.lineWidth = 1.2; ctx.beginPath();
    state.Z.forEach(function (z, i) {
      var px = 40 + (W - 50) * i / state.N, py = top + 4 + (stripH - 8) * (1 - (z - zmin) / ((zmax - zmin) || 1));
      if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    });
    ctx.stroke();
    ctx.fillStyle = P.axis;
    ctx.fillText('driving path Z(t), fBm with H = ' + state.H.toFixed(2), 46, top + 14);

    // legend band above the phase plane
    ctx.textBaseline = 'middle';
    rows.forEach(function (r, ri) {
      var yy = 10 + 14 * ri;
      r.forEach(function (e) {
        var it = e.it, x0 = e.x;
        ctx.strokeStyle = it.c; ctx.fillStyle = it.c; ctx.lineWidth = 2;
        if (it.k === 'dot') { ctx.beginPath(); ctx.arc(x0 + 9, yy, 3, 0, 2 * Math.PI); ctx.fill(); }
        else {
          ctx.setLineDash(it.k === 'dash' ? [5, 4] : it.k === 'dots' ? [2, 3] : []);
          ctx.beginPath(); ctx.moveTo(x0, yy); ctx.lineTo(x0 + 18, yy); ctx.stroke(); ctx.setLineDash([]);
        }
        ctx.fillStyle = P.text; ctx.fillText(it.t, x0 + 24, yy);
      });
    });
    ctx.textBaseline = 'alphabetic'; ctx.lineWidth = 1;
  }

  function drawPlots() {
    ['loss', 'gnorm', 'gerr', 'Q'].forEach(function (k) {
      plots[k].setSeries(KEYS.map(function (key) { return { name: SCHEMES[key].name, color: SCHEMES[key].color, data: state.runs[key][k] }; }));
      plots[k].draw();
    });
  }

  function readout() {
    document.getElementById('readout').innerHTML = 'iteration ' + state.iter + (state.iter >= MAX_ITER ? ' (finished; press Restart descent to run again)' : '') + ' · ' + KEYS.map(function (k) {
      var r = state.runs[k];
      return '<b style="color:' + SCHEMES[k].color + '">' + SCHEMES[k].name + '</b>: x₀ = (' + r.x0[0].toFixed(4) + ', ' + r.x0[1].toFixed(4) + ')';
    }).join(' · ') + ' · truth x₀ = (' + state.truth0[0].toFixed(2) + ', ' + state.truth0[1].toFixed(2) + ')';
  }

  var tick = 0;
  function frame() {
    if (state.running && state.iter < MAX_ITER && (tick++ % 3 === 0)) iterate();
    drawPhase(); drawPlots(); readout();
    requestAnimationFrame(frame);
  }

  // ---------- controls ----------
  D.bindRange('hurst', function (v) { return v.toFixed(2); }, function (v) { state.H = v; });
  D.bindRange('sigma', function (v) { return v.toFixed(2); }, function (v) { state.sigma = v; });
  D.bindRange('nsteps', function (v) { return v + ' steps'; }, function (v) { state.N = v; });
  D.bindRange('lr', function (v) { return v.toFixed(2); }, function (v) { state.lr = v; });
  ['hurst', 'sigma', 'nsteps'].forEach(function (id) { document.getElementById(id).addEventListener('change', setup); });
  document.getElementById('btn-restart').addEventListener('click', restart);
  document.getElementById('btn-path').addEventListener('click', function () { state.seed = (state.seed * 7919 + 13) % 100000; setup(); });
  var pauseBtn = document.getElementById('btn-pause');
  pauseBtn.addEventListener('click', function () { state.running = !state.running; pauseBtn.textContent = state.running ? 'Pause' : 'Play'; });

  window.addEventListener('resize', resize);
  resize();
  setup();
  // handle for debugging from the console
  window.DemoDebug = { state: state, step: iterate, reset: setup };
  requestAnimationFrame(frame);
})();
