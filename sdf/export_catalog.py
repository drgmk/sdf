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


class ExportCatalog:
    """Small read-only catalog projection over an SDB export root."""

    def __init__(self, root):
        self.root = Path(root).expanduser().resolve()

    def samples(self):
        """Return the newest completed manifest for every exported sample."""

        samples = {}
        for path in sorted(self.root.glob('export-*-manifest.json')):
            try:
                payload = json.loads(path.read_text(encoding='utf-8'))
            except (OSError, json.JSONDecodeError):
                continue
            if payload.get('schema') != 'sdb-fit-package-export':
                continue
            if payload.get('status') != 'completed':
                continue
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
                run_id=int(payload.get('run_id', 0)),
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
