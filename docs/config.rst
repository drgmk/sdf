Configuration
=============

The complete configuration is ``sdf/sdf.toml`` in the source tree and is
installed with the package. It supplies the defaults used by
``sdf-fit``, ``sdf-fits``, and the Python API. It includes the fitting and model
settings used by current code as well as settings retained for specialist
model-generation and historical database/sample utilities.

Most installations should not copy or replace it. Put only machine-specific
changes in ``~/.sdf.conf`` using INI syntax. A full server commonly needs just::

   [file]
   model_root = /path/to/full/model/library/

Other tables and setting names are the same as those in ``sdf.toml``. For
example, a local MultiNest override could be::

   [fitting]
   n_live = 300
   verb = false

Relative paths in ``sdf.toml`` are interpreted relative to that file; paths in
``~/.sdf.conf`` are interpreted relative to the home directory. Unknown TOML
tables and settings are rejected so spelling mistakes cannot silently change a
fit.

For an unusual deployment, ``SDF_CONFIG`` may select an explicit partial TOML
or INI override instead of ``~/.sdf.conf``::

   export SDF_CONFIG=/absolute/path/to/override.toml

Per-system model choices do not belong in the process-wide configuration. They
remain in the web UI or package-local ``sdf-fit.yml``.

Module reference
----------------

.. automodule:: sdf.config
   :members: load_config
