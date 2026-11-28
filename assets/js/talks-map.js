/* World map of selected talks (Talks page). Countries are shaded by the number of talks given there,
   and cities are marked with dots sized by their number of talks; a second panel enlarges Europe, where
   most of the cities are. Data comes from _data/talks.yml via the #talks-data JSON block; country
   outlines from the world-atlas package (Natural Earth, 1:110m and 1:50m). */
(function () {
  var worldHost = document.getElementById('talks-map-world');
  var europeHost = document.getElementById('talks-map-europe');
  var dataEl = document.getElementById('talks-data');
  if (!worldHost || !dataEl || !window.d3 || !window.topojson) return;

  var ATLAS = 'https://cdn.jsdelivr.net/npm/world-atlas@2.0.2/';
  // display names in _data/talks.yml that differ from the Natural Earth names in the atlas
  var ATLAS_NAMES = { 'United States': 'United States of America' };
  // area shown in the Europe panel (lon/lat corners)
  var WORLD_ASPECT = 0.45;
  var EUROPE = { type: 'MultiPoint', coordinates: [[-6.5, 41.5], [13, 41.5], [-6.5, 57], [13, 57]] };

  var talks = JSON.parse(dataEl.textContent);

  // group by country (keyed by atlas name) and by city
  var countries = {}, cityMap = {};
  talks.forEach(function (t) {
    var key = ATLAS_NAMES[t.country] || t.country;
    (countries[key] = countries[key] || { name: t.country, talks: [] }).talks.push(t);
    var ck = t.city + '|' + t.country;
    (cityMap[ck] = cityMap[ck] || { city: t.city, country: t.country, lat: +t.lat, lon: +t.lon, talks: [] }).talks.push(t);
  });
  var cities = Object.keys(cityMap).map(function (k) { return cityMap[k]; })
    .sort(function (a, b) { return b.talks.length - a.talks.length; });
  var maxCountry = d3.max(Object.keys(countries), function (k) { return countries[k].talks.length; });
  var maxCity = d3.max(cities, function (c) { return c.talks.length; });

  // colours follow the site's theme variables, so the map matches the colour scheme and the
  // light/dark mode. One hue for the number of talks: on light, light accent -> accent -> navy;
  // on dark, dim accent -> accent -> light accent (more talks always stand out more).
  var ramp, rampLo, rampHi, LAND, DOT;
  function plural(n) { return n + ' talk' + (n === 1 ? '' : 's'); }
  function fill(n) { return ramp(maxCountry > 1 ? rampLo + (rampHi - rampLo) * (n - 1) / (maxCountry - 1) : (rampLo + rampHi) / 2); }
  function setColours() {
    var css = getComputedStyle(document.documentElement);
    var v = function (name, fallback) { return (css.getPropertyValue(name) || '').trim() || fallback; };
    var dark = document.documentElement.getAttribute('data-theme') === 'dark';
    var accent = v('--gold', '#C69214'), accentLight = v('--gold-light', '#F3E3B8'), navy = v('--navy', '#182B49');
    if (dark) {
      ramp = d3.piecewise(d3.interpolateLab, [d3.interpolateLab(v('--panel', '#121A25'), accent)(0.4), accent, accentLight]);
      rampLo = 0.2; rampHi = 0.85;
      DOT = v('--ink', '#E8EEF6');
    } else {
      ramp = d3.piecewise(d3.interpolateLab, [accentLight, accent, navy]);
      rampLo = 0.22; rampHi = 0.72;
      DOT = navy;
    }
    LAND = v('--map-land', '#E4E8EE');
    // legend: one swatch per count, and the dot meaning
    var legend = document.getElementById('talks-map-legend');
    if (legend) {
      var html = '';
      for (var n = 1; n <= maxCountry; n++) html += '<span><i style="background:' + fill(n) + '"></i>' + plural(n) + '</span>';
      html += '<span><i class="dot" style="background:' + DOT + '"></i>city</span>';
      legend.innerHTML = html;
    }
  }
  setColours();

  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); }

  function makePanel(host, opts) {
    var tip = document.createElement('div');
    tip.className = 'talks-map-tip';
    tip.hidden = true;
    host.appendChild(tip);
    var svg;

    function showTip(evt, title, list) {
      tip.innerHTML = '<b>' + esc(title) + '</b><ul>' + list.map(function (t) {
        return '<li>' + esc(t.event) + ' <span>' + esc(t.date) + (t.online ? ', online' : '') + '</span></li>';
      }).join('') + '</ul>';
      tip.hidden = false;
      var r = host.getBoundingClientRect(), x = evt.clientX - r.left, y = evt.clientY - r.top;
      var w = tip.offsetWidth, h = tip.offsetHeight;
      tip.style.left = Math.max(4, Math.min(r.width - w - 4, x - w / 2)) + 'px';
      tip.style.top = (y - h - 14 >= 0 ? y - h - 14 : y + 18) + 'px';
    }
    function hideTip() { tip.hidden = true; }

    return {
      hide: hideTip,
      draw: function (features) {
        var width = host.clientWidth, height = opts.height ? opts.height(width) : Math.round(width * opts.aspect);
        var projection = opts.projection().fitExtent([[6, 6], [width - 6, height - 6]], opts.fit(features));
        var path = d3.geoPath(projection);
        hideTip();
        if (svg) svg.remove();
        svg = d3.select(host).insert('svg', ':first-child')
          .attr('viewBox', '0 0 ' + width + ' ' + height).attr('width', width).attr('height', height);

        svg.append('g').selectAll('path').data(features).join('path')
          .attr('class', function (f) { return countries[f.properties.name] ? 'country has-talks' : 'country'; })
          .attr('d', path)
          .attr('fill', function (f) { var c = countries[f.properties.name]; return c ? fill(c.talks.length) : LAND; })
          .on('pointermove', function (evt, f) {
            var c = countries[f.properties.name];
            if (c) showTip(evt, c.name + ': ' + plural(c.talks.length), c.talks); else hideTip();
          })
          .on('pointerleave', hideTip);

        var r = d3.scaleSqrt().domain([1, maxCity]).range(opts.radius);
        svg.append('g').selectAll('circle').data(cities).join('circle')
          .attr('class', 'city')
          .attr('transform', function (c) { return 'translate(' + projection([c.lon, c.lat]) + ')'; })
          .attr('r', function (c) { return r(c.talks.length); })
          .attr('fill', DOT)
          .on('pointermove', function (evt, c) { showTip(evt, c.city + ', ' + c.country + ': ' + plural(c.talks.length), c.talks); })
          .on('pointerleave', hideTip);

        if (opts.frame) {
          var b = path.bounds(opts.frame);
          svg.append('rect').attr('class', 'inset-frame')
            .attr('x', b[0][0]).attr('y', b[0][1]).attr('width', b[1][0] - b[0][0]).attr('height', b[1][1] - b[0][1]);
        }
      }
    };
  }

  var panels = [];
  var worldPanel = makePanel(worldHost, {
    projection: d3.geoNaturalEarth1,
    aspect: WORLD_ASPECT,
    fit: function (features) { return { type: 'FeatureCollection', features: features }; },
    radius: [2.5, 4],
    frame: EUROPE
  });
  var europePanel = europeHost && makePanel(europeHost, {
    projection: function () { return d3.geoConicConformal().rotate([-3, 0]).parallels([42, 56]); },
    // side by side, match the world map's height; stacked (narrow screens), use its own shape
    height: function (width) {
      var a = document.getElementById('talks-map-world-panel'), b = document.getElementById('talks-map-europe-panel');
      var sideBySide = a && b && Math.abs(a.offsetTop - b.offsetTop) < 2;
      return sideBySide ? Math.round(worldHost.clientWidth * WORLD_ASPECT) : Math.round(width * 0.62);
    },
    fit: function () { return EUROPE; },
    radius: [3.5, 6]
  });

  function load(file) {
    return fetch(ATLAS + file).then(function (r) { return r.json(); }).then(function (topo) {
      return topojson.feature(topo, topo.objects.countries).features
        .filter(function (f) { return f.properties.name !== 'Antarctica'; });
    });
  }

  function failed(host) {
    host.classList.add('talks-map-failed');
    host.textContent = 'The map could not be loaded; the list of countries below gives the same information.';
  }

  load('countries-110m.json').then(function (world) {
    panels.push([worldPanel, world]);
    worldPanel.draw(world);
  }).catch(function () { failed(worldHost); if (europeHost) europeHost.hidden = true; });

  // the Europe panel uses the finer 1:50m outlines, so coastlines hold up when enlarged
  if (europePanel) {
    load('countries-50m.json').then(function (world) {
      panels.push([europePanel, world]);
      europePanel.draw(world);
    }).catch(function () { europeHost.parentNode.hidden = true; });
  }

  var last = worldHost.clientWidth;
  window.addEventListener('resize', function () {
    if (worldHost.clientWidth === last) return;
    last = worldHost.clientWidth;
    panels.forEach(function (p) { p[0].draw(p[1]); });
  });
  // light/dark switch: recolour and redraw
  window.addEventListener('themechange', function () {
    setColours();
    panels.forEach(function (p) { p[0].draw(p[1]); });
  });
})();
