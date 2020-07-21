<div align="center">
  <img src="https://raw.githubusercontent.com/Sander-w/breakwater/master//doc/_figures/bw-icon-html_icon.png" width="500"><br>
</div>

-----------------

# breakwater: a Python package for the conceptual design of breakwaters
[![Documentation Status](https://readthedocs.org/projects/breakwater/badge/?version=latest)](https://breakwater.readthedocs.io/en/latest/?badge=latest)
![GitHub issues](https://img.shields.io/github/issues-raw/Sander-w/breakwater)
![PyPI](https://img.shields.io/pypi/v/breakwater)

<table class="tg">
  <tr>
    <td align="left" class="tg-0pky">Date</td>
    <td align="left" class="tg-0pky">July 2020</td>
  </tr>
  <tr>
    <td class="tg-0pky">Author</td>
    <td class="tg-0pky">Sander Winkel</td>
  </tr>
  <tr>
    <td class="tg-0pky">Contact</td>
    <td class="tg-0pky"><a href="mailto:pybreakwater@gmail.com">pybreakwater@gmail.com</a></td>
  </tr>  
  <tr>
    <td class="tg-0pky">Master Thesis</td>
    <td class="tg-0pky"><a href="https://repository.tudelft.nl/">Developing a design automation tool for the conceptual design of breakwaters</a></td>
  </tr>
  <tr>
    <td class="tg-0pky">Documentation</td>
    <td class="tg-0pky"><a href="https://breakwater.readthedocs.io/en/latest/index.html">https://breakwater.readthedocs.io/</a></td>
  </tr>
  <tr>
    <td class="tg-0pky">License</td>
    <td class="tg-0pky">CC-BY-NC-SA 4.0</td>
  </tr>
</table>

## Purpose

The main goal of **breakwater** is to support the designer in exploring different
types of breakwaters, see Figure 1.1 for a cross section of the breakwater
types defined by CIRIA, CUR, CETMEF (2007). The main use case is thus to make a
conceptual design of a breakwater using one of the design classes, see [Chapter 7]
of the documentation. However, because all functions and classes are also available
from `breakwater.core` it is also possible to develop your own design automation
script.

   [Chapter 7]: https://breakwater.readthedocs.io/en/latest/types.html

<div align="center">
  <img src="https://raw.githubusercontent.com/Sander-w/breakwater/master/doc/_figures/breakwater-types.png" width="800"><br>
</div>   

###### Note
> Due to the limited time not all breakwater types defined by CIRIA, CUR, CETMEF (2007)
  could be implemented. Currently only the following structures have been implemented:
  conventional rubble mound breakwaters with rock and armour units as armour layer,
  caisson breakwaters and vertically composite breakwaters.

## Main Features

**breakwater** offers the following features to support the designer:

- Design a rubble mound breakwater with rock or concrete armour units as armour
  layer, with `bw.RockRubbleMound` or `bw.ConcreteRubbleMound`.
- Design a vertical or vertically composite breakwater with `bw.Caisson`.
- Design with an interactive design application by using `bw.interactive_design`.
  With this application several parameters can be changed with sliders, to assess
  the influence of certain parameters on the design and cost.
- Design multiple breakwaters at ones by using `bw.Configurations`. With this class
  multiple breakwater types can be designed at ones, these concepts can than
  be assess by using a multi-criteria analysis or with the [DesignExplorer].
- Use the functions and classes from `breakwater.core` to create your own
  design automation approach. `breakwater.core` consist of all functions and
  classes defined in Chapters [8] to [11].

   [DesignExplorer]: http://tt-acm.github.io/DesignExplorer
   [8]: https://breakwater.readthedocs.io/en/latest/stability.html
   [11]: https://breakwater.readthedocs.io/en/latest/geo.html

## Getting Started

The latest released version is available at the [Python package index]

```sh
pip install breakwater
```

Alternatively, the source code can be downloaded, or cloned, from the GitHub repository.  

```sh
python setup.py install
```

A tutorial of **breakwater** can be found in [Chapter 3] of the documentation. This
[documentation] also provides a full overview of the implemented failure mechanisms,
and is thus a full overview of all features of **breakwater**.

   [Python package index]: https://pypi.org/project/breakwater/
   [documentation]: https://breakwater.readthedocs.io/en/latest/index.html
   [Chapter 3]: https://breakwater.readthedocs.io/en/latest/tutorial.html
