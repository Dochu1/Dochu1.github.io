---
permalink: /
title: "About me"
excerpt: "About me"
author_profile: true
redirect_from: 
  - /about/
  - /about.html
---
I am a **Stefan E. Warschawski (S.E.W) Visiting Assistant Professor** in the Department of Mathematics at the University of California San Diego. I am also affiliated with Imperial College London as a Visiting Researcher.

Previously, I was a research associate in the department of Mathematics at Imperial College London working on the **Stochastic Parameterization of Ocean Turbulence for Observational Networks** (SPOT-ON) project.

Research interests
------
<div class="box-grid one">
<div class="box-item">
<h3>Stochastic &amp; rough geometric mechanics</h3>
<p>Variational principles with transport noise and rough paths, Euler–Poincaré and Lie–Poisson reduction, and the geometric structure of stochastic fluid models such as SALT and SFLT.</p>
</div>
<div class="box-item">
<h3>Structure-preserving numerical methods</h3>
<p>Symplectic, multisymplectic and variational integrators for ODEs, PDEs and SPDEs that respect conservation laws, Casimirs and symmetries, together with geometric methods for adjoint systems. Some of this work can be explored on the <a href="/science/#demos">Science page</a>.</p>
</div>
<div class="box-item">
<h3>Geophysical fluid dynamics</h3>
<p>Wave–mean flow and wave–current interaction, vertical-slice and thin-domain models, stochastic parameterisation of ocean turbulence, and nonlinear waves and solitons.</p>
</div>
<div class="box-item">
<h3>Data assimilation &amp; machine learning</h3>
<p>Geometric approaches to uncertainty quantification and data assimilation, and machine learning algorithms and architectures informed by geometric mechanics.</p>
</div>
</div>

Selected publications & preprints
------
The full list is on the [publications page](/publications/), and also on [Google Scholar](https://scholar.google.co.uk/citations?user=6_olQZoAAAAJ&hl=en) and [arXiv](https://arxiv.org/a/hu_r_3).

{% assign selected = site.data.publications | where: "selected", true %}
<details class="pubs">
<summary>View Selected Publications &amp; Preprints</summary>
<div class="pubs-content">
{% include pub-list.html pubs=selected %}
</div>
</details>

Contact
------
<div class="contact-info" markdown="1">
{% if site.author.email %}**Email:** [{{ site.author.email }}](mailto:{{ site.author.email }})
{% endif %}
**Mailing address:**  
Department of Mathematics  
University of California San Diego  
9500 Gilman Drive  
La Jolla, CA 92093, USA
</div>
