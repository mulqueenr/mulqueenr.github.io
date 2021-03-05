---
layout: default
title: Single cell omics technologies
---
![omics overview](/assets/images/omics_overviews.png){: .center-image width="60%"}

## Welcome!

The following pages include my full code for processing of various projects while working in the <a href="https://www.ohsu.edu/school-of-medicine/oroak-lab">O'Roak</a> and <a href="https://adeylab.org">Adey labs</a>.

Some of the basic processing (bcl2fastq to generation of a counts matrix) was performed in our bash environment or with the Adey lab developed <a href="https://github.com/adeylab/scitools-dev">scitools helper functions</a>. The rest was performed in either R or python environments. All of this is a work of progress and will likely change quite a bit.

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
      {% if doc.category == "s3" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>

<h3>s3 Assays</h3>
<ul>
    {% for doc in site.pages %}
      {% if doc.category == "s3" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>
