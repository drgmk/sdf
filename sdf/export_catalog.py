"""Read the current sample views carried by SDB export manifests."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path


@dataclass(frozen=True)
class ExportSample:
    name: str
    identifiers: tuple[str, ...]
    run_id: int
    manifest_path: Path
    completed_at: str | None
    sample_date: str | None = None
    note: str | None = None


class ExportCatalog:
    """Small read-only catalog projection over an SDB export root."""

    def __init__(self, root):
        self.root = Path(root).expanduser().resolve()

    def samples(self):
        """Return sample views represented by the current export root."""

        manifests = []
        for path in sorted(self.root.glob('export-*-manifest.json')):
            try:
                payload = json.loads(path.read_text(encoding='utf-8'))
            except (OSError, json.JSONDecodeError):
                continue
            if payload.get('schema') != 'sdb-fit-package-export':
                continue
            if payload.get('status') != 'completed':
                continue
            manifests.append((int(payload.get('run_id', 0)), path, payload))

        snapshots = [
            value for value in manifests
            if isinstance(value[2].get('samples'), list)
        ]
        if snapshots:
            run_id, path, payload = max(snapshots, key=lambda value: value[0])
            samples = []
            for row in payload['samples']:
                if not isinstance(row, dict):
                    continue
                name = row.get('name')
                identifiers = row.get('sdbids')
                if (
                    not isinstance(name, str) or not name.strip()
                    or not isinstance(identifiers, list)
                ):
                    continue
                samples.append(ExportSample(
                    name=name.strip(),
                    identifiers=tuple(dict.fromkeys(
                        str(value) for value in identifiers
                    )),
                    run_id=run_id,
                    manifest_path=path,
                    completed_at=payload.get('completed_at'),
                    sample_date=row.get('sample_date'),
                    note=row.get('note'),
                ))
            return tuple(sorted(samples, key=lambda value: value.name.lower()))

        samples = {}
        for run_id, path, payload in manifests:
            selection = payload.get('selection', {})
            if selection.get('kind') != 'sample':
                continue
            name = selection.get('sample') or selection.get('value')
            if not isinstance(name, str) or not name.strip():
                continue
            identifiers = selection.get('selected_sdbids', ())
            if not isinstance(identifiers, list):
                continue
            sample = ExportSample(
                name=name.strip(),
                identifiers=tuple(dict.fromkeys(
                    str(value) for value in identifiers
                )),
                run_id=run_id,
                manifest_path=path,
                completed_at=payload.get('completed_at'),
            )
            previous = samples.get(sample.name)
            if previous is None or sample.run_id > previous.run_id:
                samples[sample.name] = sample
        return tuple(samples[name] for name in sorted(samples, key=str.lower))

    def sample(self, name):
        clean = str(name).strip()
        for sample in self.samples():
            if sample.name == clean:
                return sample
        raise KeyError('exported sample not found: {}'.format(clean))
