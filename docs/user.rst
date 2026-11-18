Command-line tools
==================

The web interface is the normal way to search results, configure exceptional
fits, and choose the public candidate. Command-line tools remain available for
service operation, automation, and recovery. ``COMMAND --help`` is the
authoritative option reference.

Current catalog and package tools
---------------------------------

``sdf-fits``
  Serve the complete exported catalog, or inspect, fit, render, reprioritize,
  and select retained candidates for one package. See :doc:`server-operations`
  and :doc:`joint-fitting`.

``sdf-fit``
  Fit one IPAC photometry table without running the catalog service. The input
  is staged into a self-contained package and fitted through the same retained-
  candidate core as server-managed systems. Use ``--stellar`` and repeatable
  ``--disk`` options for an exact model combination, or omit them for automatic
  model selection.

``sdf-index``
  Rebuild or verify the disposable SQLite projection of package-local results.

``sdf-jobs``
  Inspect the durable fit queue or run a headless queue consumer. The normal
  ``sdf-fits serve`` process supervises its own workers.

Compatibility tools
-------------------

``sdf-sample``
  Generate sample products. The current hierarchy-aware mode uses
  ``--index``; the form without an index uses the legacy result database.

``sdf-cleanup``
  Legacy MultiNest cleanup for the earlier sample/static-site workflow.
