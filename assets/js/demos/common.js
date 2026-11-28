/* Shared helpers for the interactive demos: seeded noise, HiDPI canvases,
   and a small time-series plot with optional log scale. */
(function (global) {
  'use strict';

  // Seeded uniform RNG (mulberry32) and standard normals (Box–Muller).
  function makeRng(seed) {
    var s = seed >>> 0;
    function uniform() {
      s = (s + 0x6D2B79F5) >>> 0;
      var t = s;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    }
    var spare = null;
    function normal() {
      if (spare !== null) { var v = spare; spare = null; return v; }
      var u = 0, w = 0;
      while (u === 0) u = uniform();
      w = uniform();
      var r = Math.sqrt(-2 * Math.log(u));
      spare = r * Math.sin(2 * Math.PI * w);
      return r * Math.cos(2 * Math.PI * w);
    }
    return { uniform: uniform, normal: normal };
  }

  // Brownian increment sqrt(h)*N(0,1), truncated at A = sqrt(2|ln h|) * sqrt(h)
  // (Milstein–Repin–Tretyakov) so that implicit schemes always have a solution.
  function truncatedIncrement(rng, h) {
    var sh = Math.sqrt(h);
    var A = Math.sqrt(2 * Math.abs(Math.log(h)));
    var z = rng.normal();
    if (z > A) z = A; else if (z < -A) z = -A;
    return sh * z;
  }

  // Resize a canvas to its CSS width at device pixel ratio; returns a 2D context
  // whose units are CSS pixels.
  function fitCanvas(canvas, aspect) {
    var dpr = Math.min(window.devicePixelRatio || 1, 2);
    var w = canvas.clientWidth || canvas.parentNode.clientWidth;
    var h = Math.round(w * aspect);
    canvas.style.height = h + 'px';
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
    var ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx: ctx, w: w, h: h };
  }

  // Time-series plot. series: [{name, color, data: [[t, y], ...]}]
  function TimePlot(canvas, opts) {
    this.canvas = canvas;
    this.opts = opts || {};
    this.series = [];
    this.resize();
  }
  TimePlot.prototype.resize = function () {
    var f = fitCanvas(this.canvas, this.opts.aspect || 0.62);
    this.ctx = f.ctx; this.w = f.w; this.h = f.h;
  };
  TimePlot.prototype.setSeries = function (series) { this.series = series; };
  TimePlot.prototype.draw = function () {
    var ctx = this.ctx, W = this.w, H = this.h, o = this.opts, P = ink();
    var log = !!o.log;
    var padL = 46, padR = o.endLabels ? 52 : 10, padB = 22;
    // legend row(s) for named series, laid out under the title so they never cover the curves
    var named = this.series.filter(function (s) { return s.name && !s.hidden; });
    var legendY = o.subtitle ? 44 : 30, legendRows = [];
    if (named.length && o.legend !== false) {
      ctx.font = '600 10px Inter, -apple-system, "Segoe UI", Roboto, sans-serif';
      var row = [], x = padL;
      named.forEach(function (s) {
        var w = 18 + ctx.measureText(s.name).width + 14;
        if (row.length && x + w > W - 6) { legendRows.push(row); row = []; x = padL; }
        row.push({ s: s, x: x }); x += w;
      });
      legendRows.push(row);
    }
    var padT = legendY - 8 + 14 * legendRows.length + (legendRows.length ? 6 : 0);
    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = P.bg;
    ctx.fillRect(0, 0, W, H);

    var tmax = 1e-9, ymin = Infinity, ymax = -Infinity;
    var tr = function (y) { return log ? Math.log10(Math.max(Math.abs(y), o.floor || 1e-16)) : y; };
    this.series.forEach(function (s) {
      if (s.hidden) return;
      for (var i = 0; i < s.data.length; i++) {
        var p = s.data[i];
        if (p[0] > tmax) tmax = p[0];
        var y = tr(p[1]);
        if (isFinite(y)) { if (y < ymin) ymin = y; if (y > ymax) ymax = y; }
      }
    });
    if (!isFinite(ymin)) { ymin = log ? -16 : 0; ymax = log ? 0 : 1; }
    (o.refLines || []).forEach(function (r) { var y = tr(r.y); if (y < ymin) ymin = y; if (y > ymax) ymax = y; });
    if (log) { ymin = Math.floor(ymin); ymax = Math.max(Math.ceil(ymax), ymin + 1); }
    else if (ymax - ymin < 1e-12) { ymin -= 0.5; ymax += 0.5; }
    else { var pad = 0.08 * (ymax - ymin); ymin -= pad; ymax += pad; }

    var X = function (t) { return padL + (W - padL - padR) * t / tmax; };
    var Y = function (y) { return padT + (H - padT - padB) * (1 - (y - ymin) / (ymax - ymin)); };

    // grid + labels
    ctx.font = '11px Inter, -apple-system, "Segoe UI", Roboto, sans-serif';
    ctx.fillStyle = P.axis;
    ctx.strokeStyle = P.grid;
    ctx.lineWidth = 1;
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    var ticks = [];
    if (log) {
      var stepL = Math.max(1, Math.ceil((ymax - ymin) / 5));
      for (var k = ymin; k <= ymax; k += stepL) ticks.push(k);
    } else {
      for (var j = 0; j <= 4; j++) ticks.push(ymin + (ymax - ymin) * j / 4);
    }
    ticks.forEach(function (v) {
      var yy = Y(v);
      ctx.beginPath(); ctx.moveTo(padL, yy); ctx.lineTo(W - padR, yy); ctx.stroke();
      var label = log ? '1e' + v : (Math.abs(v) >= 100 || (Math.abs(v) < 0.01 && v !== 0) ? v.toExponential(1) : v.toFixed(3));
      ctx.fillText(label, padL - 5, yy);
    });
    ctx.textAlign = 'right';
    ctx.textBaseline = 'alphabetic';
    ctx.fillText((o.xlabel || 't') + ' = ' + (o.xint ? Math.round(tmax) : tmax.toFixed(1)), W - padR, H - 6);

    // title
    ctx.textAlign = 'left';
    ctx.fillStyle = P.title;
    ctx.font = '600 12px Inter, -apple-system, "Segoe UI", Roboto, sans-serif';
    ctx.fillText(o.title || '', padL, 15);
    if (o.subtitle) {
      ctx.fillStyle = P.axis;
      ctx.font = '11px Inter, -apple-system, "Segoe UI", Roboto, sans-serif';
      ctx.fillText(o.subtitle, padL, 29);
    }

    // dashed reference lines, e.g. machine precision
    ctx.font = '10px Inter, -apple-system, "Segoe UI", Roboto, sans-serif';
    (o.refLines || []).forEach(function (r) {
      var yy = Y(tr(r.y));
      ctx.strokeStyle = P.ref; ctx.setLineDash([4, 3]);
      ctx.beginPath(); ctx.moveTo(padL, yy); ctx.lineTo(W - padR, yy); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = P.ref; ctx.textAlign = 'left';
      ctx.fillText(r.label, padL + 4, yy - 3);
    });

    // series
    this.series.forEach(function (s) {
      if (s.hidden || !s.data.length) return;
      ctx.strokeStyle = s.color;
      ctx.lineWidth = 1.8;
      ctx.beginPath();
      var started = false;
      var stride = Math.max(1, Math.floor(s.data.length / 600));
      for (var i = 0; i < s.data.length; i += stride) {
        var y = tr(s.data[i][1]);
        if (!isFinite(y)) continue;
        var px = X(s.data[i][0]), py = Y(Math.min(Math.max(y, ymin), ymax));
        if (!started) { ctx.moveTo(px, py); started = true; } else ctx.lineTo(px, py);
      }
      ctx.stroke();
    });

    // value of each series at its last point, at the right edge
    if (o.endLabels) {
      var ends = [];
      this.series.forEach(function (s) {
        if (s.hidden || !s.data.length) return;
        var v = s.data[s.data.length - 1][1], y = tr(v);
        if (!isFinite(y)) return;
        ends.push({ y: Y(Math.min(Math.max(y, ymin), ymax)), text: Math.abs(v) < 1e-300 ? '0' : Math.abs(v).toExponential(0), color: s.color });
      });
      ends.sort(function (a, b) { return a.y - b.y; });
      for (var e = 1; e < ends.length; e++) if (ends[e].y - ends[e - 1].y < 12) ends[e].y = ends[e - 1].y + 12;
      ctx.font = '600 10px Inter, -apple-system, "Segoe UI", Roboto, sans-serif';
      ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
      ends.forEach(function (en) { ctx.fillStyle = en.color; ctx.fillText(en.text, W - padR + 4, en.y); });
      ctx.textBaseline = 'alphabetic';
    }

    // legend rows
    ctx.font = '600 10px Inter, -apple-system, "Segoe UI", Roboto, sans-serif';
    ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    legendRows.forEach(function (row, ri) {
      var yy = legendY + 14 * ri;
      row.forEach(function (it) {
        ctx.fillStyle = it.s.color; ctx.fillRect(it.x, yy - 2, 12, 4);
        ctx.fillStyle = P.text; ctx.fillText(it.s.name, it.x + 17, yy);
      });
    });
    ctx.textBaseline = 'alphabetic';
  };

  // Plot colours for the current colour mode (see _includes/head/custom.html); read on every draw,
  // so the plots follow the light/dark switch without a reload.
  var PALETTES = {
    light: { bg: '#ffffff', grid: '#EAEAEA', axis: '#666666', title: '#182B49', text: '#333333', ref: '#9aa3ad', navy: '#182B49', muted: '#999999', label: 'rgba(255,255,255,0.85)' },
    dark: { bg: '#17212E', grid: '#2A3646', axis: '#9AA7B6', title: '#E8EEF6', text: '#D5DCE5', ref: '#6F7C8C', navy: '#C9D6E8', muted: '#7D8896', label: 'rgba(23,33,46,0.85)' }
  };
  function isDark() { return document.documentElement.getAttribute('data-theme') === 'dark'; }
  function ink() { return PALETTES[isDark() ? 'dark' : 'light']; }

  // Keep a series bounded in length by thinning every other point.
  function pushThinned(arr, point, maxLen) {
    arr.push(point);
    if (arr.length > maxLen) {
      var out = [];
      for (var i = 0; i < arr.length; i += 2) out.push(arr[i]);
      arr.length = 0;
      Array.prototype.push.apply(arr, out);
    }
  }

  function bindRange(id, fmt, onChange) {
    var el = document.getElementById(id);
    var out = document.querySelector('output[for="' + id + '"]');
    function upd() { if (out) out.textContent = fmt(parseFloat(el.value)); if (onChange) onChange(parseFloat(el.value)); }
    el.addEventListener('input', upd);
    upd();
    return el;
  }

  global.Demo = {
    makeRng: makeRng,
    truncatedIncrement: truncatedIncrement,
    fitCanvas: fitCanvas,
    TimePlot: TimePlot,
    pushThinned: pushThinned,
    bindRange: bindRange,
    ink: ink,
    isDark: isDark,
    colors: { sp: '#C69214', other: '#D1495B', grey: '#666666', truth: '#2BB673' }
  };
  // dark series colour (observations, driving path): navy on light, pale slate on dark
  Object.defineProperty(global.Demo.colors, 'navy', { get: function () { return ink().navy; }, enumerable: true });
})(window);
