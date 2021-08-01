---
layout: front_page
title: Single cell omics technologies
---


## Welcome!

The following pages include my full code for processing of various projects while working in the <a href="https://www.ohsu.edu/school-of-medicine/oroak-lab">O'Roak</a> and <a href="https://adeylab.org">Adey labs</a>. More recently, I joined <a href="https://www.ohsu.edu/knight-cancer-institute/cedar">CEDAR</a>, with actively developed code available as well. I also added some background information about the state of the field for various single-cell assays (as of 2020-ish) which can be found in the "Background" section.

Some of the basic processing (bcl2fastq to generation of a counts matrix) was performed in our bash environment or with the Adey lab developed <a href="https://github.com/adeylab/scitools-dev">scitools helper functions</a>. The rest was performed in either R or python environments. All of this is a work of progress and will likely change quite a bit as projects progress.

<h3>Background</h3>
<ul>
    {% for doc in site.pages | sort:"order" %}
      {% if doc.category == "background" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
        {% endif %}
    {% endfor %}
</ul>

<h3>sciATAC</h3>
<ul>
    {% for doc in site.pages %}
      {% if doc.category == "sciATAC" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>

<h3>s3 Assays</h3>
<ul>
    {% for doc in site.pages %}
      {% if doc.category == "s3processing" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>

<h3>CEDAR Projects</h3>
<ul>
    {% for doc in site.pages %}
      {% if doc.category == "CEDAR" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>

<h3>Alternative Assays</h3>
<ul>
    {% for doc in site.pages %}
      {% if doc.category == "alternative" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>

<h3>Side Projects</h3>
<ul>
    {% for doc in site.pages %}
      {% if doc.category == "side_projects" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>



