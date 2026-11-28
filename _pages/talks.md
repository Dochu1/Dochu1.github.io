---
layout: archive
title: "Talks"
permalink: /talks/
author_profile: true
---

{% comment %} Talks are listed in _data/talks.yml, which also drives the map at the end of the page. {% endcomment %}

## Selected talks

{% assign talks_by_year = site.data.talks | group_by: "year" %}
{% for y in talks_by_year %}
### {{ y.name }}

{% for t in y.items %}- **{{ t.event }}**, {{ t.kind }}, {{ t.venue }}, {{ t.country }}{% if t.online %} (online){% endif %} ({{ t.date }})
{% endfor %}
{% endfor %}

I have also given numerous invited internal talks at Imperial College London and at STUOD workshops.

## Poster presentations

- **Foundations of Computational Mathematics (FoCM) 2026**, University of Vienna, Austria (Jul 2026)
- **AGU Ocean Sciences Meeting 2024**, two posters, New Orleans, United States (Feb 2024)

{% assign by_country = site.data.talks | group_by: "country" | sort: "size" | reverse %}
## Talks map

<figure class="talks-map-fig">
  <div class="talks-map-row">
    <div class="talks-map-panel" id="talks-map-world-panel">
      <div class="talks-map-label">World</div>
      <div id="talks-map-world" class="talks-map" role="img" aria-label="World map of the selected talks above, shaded by the number of talks given in each country"></div>
    </div>
    <div class="talks-map-panel" id="talks-map-europe-panel">
      <div class="talks-map-label">Europe</div>
      <div id="talks-map-europe" class="talks-map" role="img" aria-label="Enlarged map of Western Europe showing the cities where talks were given"></div>
    </div>
  </div>
  <div class="talks-map-legend" id="talks-map-legend" aria-hidden="true"></div>
</figure>

<ul class="talks-by-country">
{% for c in by_country %}<li><b>{{ c.name }}</b> {{ c.size }} talk{% if c.size != 1 %}s{% endif %}</li>
{% endfor %}</ul>

<script id="talks-data" type="application/json">{{ site.data.talks | jsonify }}</script>
<script src="https://cdn.jsdelivr.net/npm/d3@7.9.0/dist/d3.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/topojson-client@3.1.0/dist/topojson-client.min.js"></script>
<script src="{{ '/assets/js/talks-map.js' | relative_url }}"></script>
