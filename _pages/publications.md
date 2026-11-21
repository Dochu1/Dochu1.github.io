---
layout: archive
title: "Publications"
permalink: /publications/
author_profile: true
---

You can also find my publications on [my Google Scholar profile](https://scholar.google.co.uk/citations?user=6_olQZoAAAAJ&hl=en). The preprint versions of these publications can also be found on [my arXiv profile](https://arxiv.org/a/hu_r_3).

{% comment %} The lists are generated from _data/publications.yml. Edit that file to add papers; type "preprint" puts a paper under Preprints. {% endcomment %}

## Preprints

{% assign preprints = site.data.publications | where: "type", "preprint" %}
{% include pub-list.html pubs=preprints %}

## Published

{% assign published = site.data.publications | where_exp: "p", "p.type != 'preprint'" | where_exp: "p", "p.type != 'abstract'" %}
{% include pub-list.html pubs=published by_year=true %}

## Conference abstracts

{% assign abstracts = site.data.publications | where: "type", "abstract" %}
{% include pub-list.html pubs=abstracts %}
