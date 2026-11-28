/* Stochastic multisymplectic NLS with the Preissmann box scheme.
     dψ = i (ψ_xx + 2γ|ψ|²ψ) dt + i ψ_xx ∘ ε dW_D(t) − σ T ψ ∘ dW_T(t),
     T u = ½ (ξ u_x + (ξ u)_x)   (skew-adjoint stochastic transport)
   on x ∈ [−L, L) with periodic boundary conditions.

   Box scheme: midpoint in time AND midpoint in space. On each cell [x_j, x_{j+1}] × [t_n, t_{n+1}]
   the first-order (multisymplectic) system for z = (Re ψ, Im ψ, Re ψ_x, Im ψ_x) is discretised with
     δ_t M_x  and  δ_x M_t,   where  δ_x u = (u_{j+1} − u_j)/Δx,  M_x u = (u_{j+1} + u_j)/2
   (and likewise in t), with every right-hand side evaluated at the cell centre. Eliminating the
   auxiliary variable ψ_x shows that the cell averages φ = M_x ψ evolve by the implicit midpoint rule
   in time with the box derivative D = δ_x M_x⁻¹ in space. D and D² are circulant, with Fourier symbols
     D  ↦ (2i/Δx) tan(πk/N),     D² ↦ −(4/Δx²) tan²(πk/N),
   so the implicit stages are solved with discrete Fourier transforms. N is odd so that M_x is
   invertible (an even N leaves a spurious checkerboard mode).

   The scheme satisfies a discrete 2-form conservation law, so the global symplectic form
   ω(δ₁, δ₂) = Σ Im(conj(δ₁) δ₂) Δx of any two solutions of the linearised scheme is conserved,
   together with the discrete mass Σ |φ|² Δx. For linear problems (γ = 0) with translation-invariant
   noise it also conserves the discrete momentum.

   The comparison scheme is linearly implicit: the same box operators and midpoint rule for every
   linear term (dispersion and both noises) but the cubic nonlinearity taken explicitly. */
(function () {
  'use strict';
  var D = window.Demo;
  var N = 135, L = 16, dx = 2 * L / N, ROWS = 240; // N = 3³·5, odd

  // cell centres and the Fourier symbols of the box derivative D = δ_x M_x⁻¹
  var xs = new Float64Array(N), k1 = new Float64Array(N), k2 = new Float64Array(N);
  for (var j = 0; j < N; j++) {
    xs[j] = -L + (j + 0.5) * dx;
    k1[j] = 2 / dx * Math.tan(Math.PI * j / N); // D ↦ i k1
    k2[j] = k1[j] * k1[j];                      // D² ↦ −k2
  }

  // ---------- complex arrays and a mixed-radix FFT for any N ----------
  function C() { return { re: new Float64Array(N), im: new Float64Array(N) }; }
  function copyInto(a, b) { b.re.set(a.re); b.im.set(a.im); }

  var TWR = new Float64Array(N), TWI = new Float64Array(N);
  for (var t = 0; t < N; t++) { TWR[t] = Math.cos(2 * Math.PI * t / N); TWI[t] = Math.sin(2 * Math.PI * t / N); }
  function smallestFactor(n) { for (var f = 2; f * f <= n; f++) if (n % f === 0) return f; return n; }

  // Recursive decimation in time by the smallest prime factor p of n, out of place and without
  // allocation: sub-transform r is written to out[r·m … r·m + m), then combined in place.
  var tmpR = new Float64Array(16), tmpI = new Float64Array(16);
  function rec(iR, iI, iOff, st, oR, oI, oOff, n, sign) {
    if (n === 1) { oR[oOff] = iR[iOff]; oI[oOff] = iI[iOff]; return; }
    var p = smallestFactor(n), m = n / p, sc = N / n, r, k, q;
    for (r = 0; r < p; r++) rec(iR, iI, iOff + r * st, st * p, oR, oI, oOff + r * m, m, sign);
    for (k = 0; k < m; k++) {
      for (r = 0; r < p; r++) { tmpR[r] = oR[oOff + r * m + k]; tmpI[r] = oI[oOff + r * m + k]; }
      for (q = 0; q < p; q++) {
        var kk = k + q * m, ar = 0, ai = 0;
        for (r = 0; r < p; r++) {
          var idx = ((r * kk) % n) * sc, wr = TWR[idx], wi = sign * TWI[idx];
          ar += tmpR[r] * wr - tmpI[r] * wi; ai += tmpR[r] * wi + tmpI[r] * wr;
        }
        oR[oOff + q * m + k] = ar; oI[oOff + q * m + k] = ai;
      }
    }
  }
  var inR = new Float64Array(N), inI = new Float64Array(N);
  function fft(re, im, inverse) {
    inR.set(re); inI.set(im);
    rec(inR, inI, 0, 1, re, im, 0, N, inverse ? 1 : -1);
    if (inverse) for (var q = 0; q < N; q++) { re[q] /= N; im[q] /= N; }
  }

  var tA = C(), tB = C();
  function deriv(u, out) { // out = D u, the box derivative δ_x M_x⁻¹
    copyInto(u, out);
    fft(out.re, out.im, false);
    for (var j = 0; j < N; j++) { var r = out.re[j], i = out.im[j]; out.re[j] = -k1[j] * i; out.im[j] = k1[j] * r; }
    fft(out.re, out.im, true);
  }

  // spatially varying part of the transport field
  var xiT = new Float64Array(N);
  for (var j2 = 0; j2 < N; j2++) xiT[j2] = 0.8 * Math.cos(Math.PI * xs[j2] / L) + 0.4 * Math.sin(2 * Math.PI * xs[j2] / L);

  // out = T̃ u = ½ (ξ̃ u_x + (ξ̃ u)_x)
  function transportTilde(u, out) {
    deriv(u, tA);
    for (var j = 0; j < N; j++) { tB.re[j] = xiT[j] * u.re[j]; tB.im[j] = xiT[j] * u.im[j]; }
    deriv(tB, out);
    for (var q = 0; q < N; q++) {
      out.re[q] = 0.5 * (out.re[q] + xiT[q] * tA.re[q]);
      out.im[q] = 0.5 * (out.im[q] + xiT[q] * tA.im[q]);
    }
  }

  var state = {
    noise: 'transport-x', amp: 0.3, dt: 0.005, speed: 3, gamma: 1, scenario: 'collision',
    seed: 5, running: true, t: 0
  };

  /* One step of a linear-implicit problem
       u' − u = diag(ū) + R(ū),   ū = (u + u')/2,
     with diag having Fourier symbol −i a_k (purely imaginary), and R handled by fixed-point
     iteration. Returns the number of iterations (≥ 200 means it did not converge). */
  var Bh = C(), Rv = C(), Mid = C(), Nw = C();
  function implicitStep(u, a, remainder, out) {
    copyInto(u, Bh);
    fft(Bh.re, Bh.im, false);
    for (var k = 0; k < N; k++) { // (1 − i a/2) û
      var c = 0.5 * a[k], x = Bh.re[k], y = Bh.im[k];
      Bh.re[k] = x + y * c; Bh.im[k] = y - x * c;
    }
    copyInto(u, out);
    var it, prev = Infinity;
    for (it = 0; it < 200; it++) {
      for (var j = 0; j < N; j++) { Mid.re[j] = 0.5 * (u.re[j] + out.re[j]); Mid.im[j] = 0.5 * (u.im[j] + out.im[j]); }
      remainder(Mid, Rv);
      fft(Rv.re, Rv.im, false);
      for (var q = 0; q < N; q++) { // (B̂ + R̂) / (1 + i a/2)
        var cc = 0.5 * a[q], xr = Bh.re[q] + Rv.re[q], yi = Bh.im[q] + Rv.im[q], den = 1 + cc * cc;
        Nw.re[q] = (xr + yi * cc) / den; Nw.im[q] = (yi - xr * cc) / den;
      }
      fft(Nw.re, Nw.im, true);
      var err = 0, nrm = 0;
      for (var p = 0; p < N; p++) {
        err = Math.max(err, Math.abs(Nw.re[p] - out.re[p]) + Math.abs(Nw.im[p] - out.im[p]));
        nrm = Math.max(nrm, Math.abs(Nw.re[p]) + Math.abs(Nw.im[p]));
      }
      copyInto(Nw, out);
      if (err <= 1e-15 * (1 + nrm) || (it > 6 && err >= prev)) break;
      prev = err;
    }
    for (var r = 0; r < N; r++) { Mid.re[r] = 0.5 * (u.re[r] + out.re[r]); Mid.im[r] = 0.5 * (u.im[r] + out.im[r]); }
    return it;
  }

  // ---------- step data ----------
  var symb = new Float64Array(N);
  var step = { dWT: 0, dWD: 0, hasTilde: false };
  function prepareStep() {
    var dt = state.dt;
    var dW = D.truncatedIncrement(state.rng, dt);
    var cD = dt, xiBar = 0;
    step.dWT = 0; step.hasTilde = false;
    if (state.noise === 'dispersion') cD = dt + state.amp * dW;
    else if (state.noise === 'transport-const') { xiBar = 1; step.dWT = state.amp * dW; }
    else if (state.noise === 'transport-x') { step.dWT = state.amp * dW; step.hasTilde = true; }
    // diag part: i cD ∂xx − dWT ξ̄ ∂x  →  symbol −i (cD k² + dWT ξ̄ k)
    for (var k = 0; k < N; k++) symb[k] = cD * k2[k] + step.dWT * xiBar * k1[k];
  }

  var T1 = C();
  function addTransport(u, out) { // out −= dWT T̃ u
    if (!step.hasTilde || step.dWT === 0) return;
    transportTilde(u, T1);
    for (var j = 0; j < N; j++) { out.re[j] -= step.dWT * T1.re[j]; out.im[j] -= step.dWT * T1.im[j]; }
  }
  function zero(out) { out.re.fill(0); out.im.fill(0); }

  // i dt γ 2|v|² w   added to out
  function addCubic(v, out, scale) {
    for (var j = 0; j < N; j++) {
      var a2 = v.re[j] * v.re[j] + v.im[j] * v.im[j];
      out.re[j] -= scale * 2 * a2 * v.im[j];
      out.im[j] += scale * 2 * a2 * v.re[j];
    }
  }
  // i dt γ (4|v|² w + 2 v² conj(w))  added to out
  function addCubicLinearised(v, w, out, scale) {
    for (var j = 0; j < N; j++) {
      var pr = v.re[j], pi = v.im[j], a2 = pr * pr + pi * pi, A = pr * pr - pi * pi, B = 2 * pr * pi;
      var zr = 4 * a2 * w.re[j] + 2 * (A * w.re[j] + B * w.im[j]);
      var zi = 4 * a2 * w.im[j] + 2 * (B * w.re[j] - A * w.im[j]);
      out.re[j] -= scale * zi;
      out.im[j] += scale * zr;
    }
  }

  // ---------- the two schemes, each advancing ψ and two tangent vectors ----------
  var PsiOld = C(), PsiMid = C(), Const = C();
  function advanceMidpoint(s) {
    var g = state.gamma * state.dt;
    var its = implicitStep(s.psi, symb, function (mid, out) { zero(out); addTransport(mid, out); addCubic(mid, out, g); }, s.next);
    copyInto(Mid, PsiMid);
    var tmp = s.psi; s.psi = s.next; s.next = tmp;
    s.tan.forEach(function (d) {
      implicitStep(d.v, symb, function (mid, out) { zero(out); addTransport(mid, out); addCubicLinearised(PsiMid, mid, out, g); }, d.w);
      var t2 = d.v; d.v = d.w; d.w = t2;
    });
    return its;
  }
  function advanceLinImplicit(s) {
    var g = state.gamma * state.dt;
    copyInto(s.psi, PsiOld);
    zero(Const); addCubic(PsiOld, Const, g);
    var its = implicitStep(s.psi, symb, function (mid, out) { copyInto(Const, out); addTransport(mid, out); }, s.next);
    var tmp = s.psi; s.psi = s.next; s.next = tmp;
    s.tan.forEach(function (d) {
      zero(Const); addCubicLinearised(PsiOld, d.v, Const, g);
      implicitStep(d.v, symb, function (mid, out) { copyInto(Const, out); addTransport(mid, out); }, d.w);
      var t2 = d.v; d.v = d.w; d.w = t2;
    });
    return its;
  }

  var META = [
    { key: 'ms', name: 'Multisymplectic box scheme', color: D.colors.sp, advance: advanceMidpoint, canvas: 'st-ms' },
    { key: 'li', name: 'Linearly implicit (explicit nonlinearity)', color: D.colors.other, advance: advanceLinImplicit, canvas: 'st-li' }
  ];

  // ---------- invariants ----------
  function mass(u) { var s = 0; for (var j = 0; j < N; j++) s += u.re[j] * u.re[j] + u.im[j] * u.im[j]; return s * dx; }
  var Dx = C();
  function momentum(u) {
    deriv(u, Dx);
    var s = 0;
    for (var j = 0; j < N; j++) s += u.re[j] * Dx.im[j] - u.im[j] * Dx.re[j];
    return s * dx;
  }
  function omega(a, b) { var s = 0; for (var j = 0; j < N; j++) s += a.re[j] * b.im[j] - a.im[j] * b.re[j]; return s * dx; }

  // ---------- initial data ----------
  var SCENARIOS = {
    collision: function (u) {
      var sol = [[1.0, 2.0, -6], [0.8, -2.0, 6]];
      for (var j = 0; j < N; j++) {
        u.re[j] = 0; u.im[j] = 0;
        sol.forEach(function (s) {
          var amp = s[0] / Math.cosh(s[0] * (xs[j] - s[2])), ph = s[1] * xs[j] / 2;
          u.re[j] += amp * Math.cos(ph); u.im[j] += amp * Math.sin(ph);
        });
      }
    },
    breather: function (u) { for (var j = 0; j < N; j++) { u.re[j] = 2 / Math.cosh(xs[j]); u.im[j] = 0; } },
    single: function (u) {
      for (var j = 0; j < N; j++) {
        var amp = 1.2 / Math.cosh(1.2 * (xs[j] + 4)), ph = 1.5 * xs[j] / 2;
        u.re[j] = amp * Math.cos(ph); u.im[j] = amp * Math.sin(ph);
      }
    }
  };

  function reset() {
    state.t = 0; state.steps = 0;
    state.rng = D.makeRng(state.seed);
    var u0 = C();
    SCENARIOS[state.scenario](u0);
    var d1 = C(), d2 = C();
    for (var j = 0; j < N; j++) {
      d1.re[j] = Math.exp(-xs[j] * xs[j] / 4);
      d2.im[j] = Math.exp(-(xs[j] - 1) * (xs[j] - 1) / 4) * (1 + 0.3 * Math.sin(xs[j]));
      d2.re[j] = 0.2 * Math.exp(-(xs[j] + 2) * (xs[j] + 2) / 4);
    }
    state.M0 = mass(u0); state.P0 = momentum(u0); state.W0 = omega(d1, d2);
    var amax = 0;
    for (var q = 0; q < N; q++) amax = Math.max(amax, Math.hypot(u0.re[q], u0.im[q]));
    state.amax = amax * 1.5;
    state.methods = {};
    META.forEach(function (m) {
      var psi = C(); copyInto(u0, psi);
      var a = C(), b = C(); copyInto(d1, a); copyInto(d2, b);
      state.methods[m.key] = {
        psi: psi, next: C(), tan: [{ v: a, w: C() }, { v: b, w: C() }],
        rows: [], mass: [], omega: [], mom: [], dead: false, iters: 0
      };
    });
    recordRows();
  }

  function recordRows() {
    META.forEach(function (m) {
      var s = state.methods[m.key], row = new Float32Array(N);
      for (var j = 0; j < N; j++) row[j] = Math.hypot(s.psi.re[j], s.psi.im[j]);
      s.rows.push(row);
      if (s.rows.length > ROWS) s.rows.shift();
    });
  }

  function advance() {
    prepareStep();
    state.t += state.dt; state.steps++;
    META.forEach(function (m) {
      var s = state.methods[m.key];
      if (s.dead) return;
      s.iters = Math.max(s.iters, m.advance(s));
      var M = mass(s.psi);
      if (!isFinite(M) || M > 1e3 * state.M0) { s.dead = true; return; }
      var W = omega(s.tan[0].v, s.tan[1].v), P = momentum(s.psi);
      D.pushThinned(s.mass, [state.t, (M - state.M0) / state.M0], 2000);
      D.pushThinned(s.omega, [state.t, (W - state.W0) / Math.abs(state.W0)], 2000);
      D.pushThinned(s.mom, [state.t, (P - state.P0) / state.M0], 2000);
    });
    if (state.steps % Math.max(1, Math.round(0.04 / state.dt)) === 0) recordRows();
  }

  // ---------- drawing ----------
  var plots = {
    mass: new D.TimePlot(document.getElementById('plot-mass'), { log: true, title: 'Mass error |M − M₀|/M₀' }),
    omega: new D.TimePlot(document.getElementById('plot-omega'), { log: true, title: '2-form error |ω − ω₀|/|ω₀|' }),
    mom: new D.TimePlot(document.getElementById('plot-mom'), { log: true, title: 'Momentum error |P − P₀|/M₀' })
  };
  var views = {};
  var off = document.createElement('canvas'); off.width = N; off.height = ROWS;
  var offCtx = off.getContext('2d'), img = offCtx.createImageData(N, ROWS);
  // |u| colour map: zero blends into the plot background, large amplitude stands out from it
  var STOPS_LIGHT = [[255, 255, 255], [243, 227, 184], [198, 146, 20], [24, 43, 73]];
  var STOPS_DARK = [[23, 33, 46], [104, 76, 18], [198, 146, 20], [246, 232, 196]];
  var STOPS = STOPS_LIGHT;
  function cmap(v) {
    v = Math.max(0, Math.min(1, v)) * (STOPS.length - 1);
    var i = Math.min(Math.floor(v), STOPS.length - 2), f = v - i, a = STOPS[i], b = STOPS[i + 1];
    return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])];
  }

  function resize() {
    META.forEach(function (m) { views[m.key] = D.fitCanvas(document.getElementById(m.canvas), 0.36); });
    Object.keys(plots).forEach(function (k) { plots[k].resize(); });
  }

  function drawSpaceTime(m) {
    var v = views[m.key], ctx = v.ctx, s = state.methods[m.key], data = img.data, P = D.ink();
    STOPS = D.isDark() ? STOPS_DARK : STOPS_LIGHT;
    for (var r = 0; r < ROWS; r++) {
      var row = s.rows[s.rows.length - ROWS + r];
      for (var j = 0; j < N; j++) {
        var c = row ? cmap(row[j] / state.amax) : STOPS[0];
        var p = 4 * ((ROWS - 1 - r) * N + j);
        data[p] = c[0]; data[p + 1] = c[1]; data[p + 2] = c[2]; data[p + 3] = 255;
      }
    }
    offCtx.putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(off, 0, 0, v.w, v.h);
    ctx.font = '600 12px Inter, -apple-system, "Segoe UI", Roboto, sans-serif';
    var label = m.name + (s.dead ? ' — blew up' : '');
    var w = ctx.measureText(label).width;
    ctx.fillStyle = P.label;
    ctx.fillRect(6, 6, w + 12, 20);
    ctx.fillStyle = m.color === D.colors.sp ? D.colors.navy : m.color;
    ctx.fillText(label, 12, 20);
  }

  function frame() {
    if (state.running) for (var i = 0; i < state.speed; i++) advance();
    META.forEach(drawSpaceTime);
    ['mass', 'omega', 'mom'].forEach(function (k) {
      plots[k].setSeries(META.map(function (m) { return { color: m.color, data: state.methods[m.key][k] }; }));
      plots[k].draw();
    });
    var notes = [];
    if (state.noise === 'transport-x' && state.amp > 0) notes.push('ξ depends on x, so momentum is not conserved by the exact flow');
    if (state.gamma && state.noise !== 'transport-x') notes.push('with the cubic term, the box scheme conserves discrete momentum only approximately');
    var bad = META.filter(function (m) { return state.methods[m.key].iters >= 200; });
    if (bad.length) notes.push('nonlinear solve did not fully converge: reduce Δt or the noise');
    document.getElementById('readout').innerHTML = 't = ' + state.t.toFixed(2) + (notes.length ? ' · ' + notes.join(' · ') : '');
    requestAnimationFrame(frame);
  }

  // ---------- controls ----------
  D.bindRange('amp', function (v) { return v.toFixed(2); }, function (v) { state.amp = v; });
  D.bindRange('dt', function (v) { return v.toFixed(3); }, function (v) { state.dt = v; });
  D.bindRange('speed', function (v) { return v + ' steps/frame'; }, function (v) { state.speed = v; });
  document.getElementById('noise').addEventListener('change', function (e) { state.noise = e.target.value; reset(); });
  document.getElementById('scenario').addEventListener('change', function (e) { state.scenario = e.target.value; reset(); });
  document.getElementById('nonlinear').addEventListener('change', function (e) { state.gamma = e.target.checked ? 1 : 0; reset(); });
  document.getElementById('btn-reset').addEventListener('click', reset);
  document.getElementById('btn-noise').addEventListener('click', function () { state.seed = (state.seed * 7919 + 13) % 100000; reset(); });
  var pauseBtn = document.getElementById('btn-pause');
  pauseBtn.addEventListener('click', function () { state.running = !state.running; pauseBtn.textContent = state.running ? 'Pause' : 'Play'; });

  window.addEventListener('resize', resize);
  resize();
  reset();
  // handle for debugging from the console
  window.DemoDebug = { state: state, step: advance, reset: reset };
  requestAnimationFrame(frame);
})();
