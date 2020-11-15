---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
title:"Single cell omics technologies"
---


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
