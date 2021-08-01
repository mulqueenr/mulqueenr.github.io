---
layout: front_page
title: Single cell omics technologies
---


## Welcome!

My field of study is in the development of single cell technologies and their applications to cancer early detection and neurodevelopment. A summary of single cell methods for various regulatory features can be found here. 

<h3>Background</h3>
<ul>
    {% for doc in site.pages | sort:"order" %}
      {% if doc.category == "background" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
        {% endif %}
    {% endfor %}
</ul>

My full code for processing of various projects while working in the <a href="https://www.ohsu.edu/school-of-medicine/oroak-lab">O'Roak</a> and <a href="https://adeylab.org">Adey labs</a> for my graduate work is available. Some of the basic processing (bcl2fastq to generation of a counts matrix) was performed in our bash environment or with the Adey lab developed <a href="https://github.com/adeylab/scitools-dev">scitools helper functions</a>. The rest was performed in either R or python environments. All of this is a work of progress and will likely change quite a bit as projects progress. 

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

<h3>Alternative Assays</h3>
<ul>
    {% for doc in site.pages %}
      {% if doc.category == "alternative" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>

More recently, I joined <a href="https://www.ohsu.edu/knight-cancer-institute/cedar">CEDAR</a>, with actively developed code available as well.

<h3>CEDAR Projects</h3>
<ul>
    {% for doc in site.pages %}
      {% if doc.category == "CEDAR" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>

Finally, I include some side projects from hobbies including 3D printing and rendering. Have a look around!

<h3>Side Projects</h3>
<ul>
    {% for doc in site.pages %}
      {% if doc.category == "side_projects" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>



