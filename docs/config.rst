Configuration
=============

Most users do not need an SDF configuration file. The package contains a small
model set and the command-line defaults used by both ``sdf-fit`` and
``sdf-fits``. A full server installation normally configures only the external
model library.

Copy ``sdf.example.toml`` to ``sdf.toml`` in the working directory, or install
it as ``~/.config/sdf/config.toml``. An explicit file can be selected with::

   export SDF_CONFIG=/absolute/path/to/sdf.toml

Relative paths are interpreted relative to the TOML file, rather than the
process working directory. This makes one configuration usable from a shell,
launchd, or another service manager.

The supported tables are:

``[models]``
   ``root`` selects the model-library directory. Optional role patterns and
   wavelength-grid settings support custom model libraries.

``[fitting]``
   Process-wide numerical defaults: the evidence threshold, MultiNest live
   points and update interval, retained posterior sample count, normalisation
   ranges, and global filter policies. Per-system model choices belong in the
   web UI or package-local ``sdf-fit.yml``.

``[paths]``
   ``spectra`` supports historical rawphot files whose spectra are relative to
   a shared directory. Modern SDB exports carry a package-local spectra
   manifest and do not need it.

Unknown settings are rejected so spelling mistakes cannot silently change a
fit. See ``sdf.example.toml`` for the complete supported surface.

Legacy configuration
--------------------

When ``SDF_CONFIG`` is unset, an existing ``~/.sdf.conf`` is still read before
the TOML discovery locations. It exists for older database, calibration, and
model-generation utilities; new installations should use TOML. The complete
internal defaults remain in the package's ``sdf.conf`` but should not be
edited.

Module reference
----------------

.. automodule:: sdf.config
   :members: load_config
