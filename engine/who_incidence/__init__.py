"""Offline WHO TB incidence snapshots: schema, validation, adapters and trend interfaces.

Nothing in this package performs network access. Maintainers create snapshots with
``scripts/import_who_incidence.py`` from local source files; ordinary application
sessions read the validated, checksummed snapshot bundled with the repository.
"""
