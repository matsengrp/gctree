Open source code repository
===========================

All code is freely available at `<https://github.com/matsengrp/gctree>`_

Developer tools
===============

Developer install::

  make install

.. warning:: 
  
    If you're working on an ARM Mac, you may run into trouble installing the PyQt5 dependency
    (which is also a dependency of the ``historydag`` package) via pip. As a workaround,
    we recommend commenting out the the ``PyQt5``, ``ete3``, and ``historydag`` lines from
    the ``setup.py`` and installing those via Conda.

Run tests::

  make test

Format code::

  make format

Lint::

  make lint

Build docs locally (you can then see the generated documentation in ``docs/_build/html/index.html``)::

  make docs

Docs are automatically deployed to github pages via a workflow on push to the main branch.

Todo list
=========

.. todolist::
