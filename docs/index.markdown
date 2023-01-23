---
layout: front_page
title: Single cell omics technologies
---


## Welcome!

My field of study is in the development of single cell technologies and their applications to cancer early detection and neurodevelopment. A summary of the motivation and state of single cell methods for various regulatory features can be found here:

<h3>Background</h3>
<ul>
    {% for doc in site.pages | sort:"order" %}
      {% if doc.category == "background" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
        {% endif %}
    {% endfor %}
</ul>
  
My full code for processing of various projects while working towards my graduate thesis in the <a href="https://www.ohsu.edu/school-of-medicine/oroak-lab">O'Roak</a> and <a href="https://adeylab.org">Adey labs</a> is available below. Some of the basic processing (bcl2fastq to generation of a counts matrix) was performed in our bash environment or with the Adey lab developed <a href="https://github.com/adeylab/scitools-dev">scitools helper functions</a>. The rest was performed in either R or python environments. 

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

I joined <a href="https://www.ohsu.edu/knight-cancer-institute/cedar">CEDAR</a>, at OHSU for roughly one year to begin my transistion into a cancer researcher. My actively developed code for projects there are available as well.

<h3>CEDAR Projects</h3>
<ul>
    {% for doc in site.pages %}
      {% if doc.category == "CEDAR" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>

More recently, I joined up with Nick Navin at MD Anderson to explore the effect of copy number variation on the epigenome.

<h3>CEDAR Projects</h3>
<ul>
    {% for doc in site.pages %}
      {% if doc.category == "mdanderson" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>

Finally, I include some side projects from hobbies including 3D printing and rendering. 

<h3>Side Projects</h3>
<ul>
    {% for doc in site.pages %}
      {% if doc.category == "side_projects" %}
        <li><a href="{{ doc.url }}">{{ doc.title }}</a></li>
      {% endif %}
    {% endfor %}
</ul>

Have a look around!


