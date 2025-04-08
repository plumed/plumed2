# Summary of manual coverage

{% assign failed = ''  | split: ',' %}
{% assign noexamples = ''  | split: ',' %}
{% assign undockey = ''  | split: ',' %}
{% assign date = "now" | date: "%Y-%m-%d %H:%M" %}

{% for item in site.data.actionlist %}
   {% if item.nfail > 0 %}
     {% assign failed = failed | push: item %}
   {% endif %}
   {% if item.ninp == 0 %}
     {% assign noexamples = noexamples | push: item %}
   {% endif %}
   {% if item.nukey > 0 %}
      {% assign undockey = undockey | push: item %}
   {% endif %}
{% endfor %}

Data on this page was generated on {{ date }}.

__List of actions that have broken examples__

There are {{ failed.size }} action pages with failing inputs.

{:#browse-table .display}
| Name | Module | # fails |
|:----:|:-------:|:-------:|
{% for item in failed %} | [{{ item.name }}]( {{ item.path }}) | {{ item.module }} | {{ item.nfail }} |
{% endfor %}

__List of actions that have undocumented keywords__

There are {{ undockey.size }} actions with undocumented keywords

{:#browse-table .display}
| Name | Module | # undocumented |
|:----:|:-------:|:-------:|
{% for item in undockey %} | [{{ item.name }}]( {{ item.path }}) | {{ item.module }} | {{ item.nukey }} |
{% endfor %}

__List of actions that have no examples in manual__

There are {{ noexamples.size }} action pages with no examples.

{:#browse-table .display}
| Name | Module | # fails |
|:----:|:-------:|:-------:|
{% for item in noexamples %} | [{{ item.name }}]( {{ item.path }}) | {{ item.module }} | {{ item.nfail }} |
{% endfor %}
