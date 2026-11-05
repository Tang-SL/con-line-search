"""The run folder: state, results, plots, notes and a step log."""
from pathlib import Path
from datetime import datetime, timezone
import json
import math
import numpy as np


def _clean(value):
    """JSON-safe copy (numpy scalars/arrays, NaN -> None, Paths -> str)."""
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return _clean(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        value = float(value)
        return value if math.isfinite(value) else None
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


class Run:
    def __init__(self, path):
        self.path = Path(path).expanduser().resolve()
        for sub in ['', 'plots', 'products', 'notes', 'cache']:
            (self.path/sub).mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------- json
    def save(self, name, obj):
        target = self.path/f'{name}.json'
        target.write_text(json.dumps(_clean(obj), indent=1))
        return target

    def load(self, name, default=None):
        target = self.path/f'{name}.json'
        if not target.exists():
            if default is not None:
                return default
            raise FileNotFoundError(f'{target.name} is missing: run the step that makes it first.')
        return json.loads(target.read_text())

    def has(self, name):
        return (self.path/f'{name}.json').exists()

    # -------------------------------------------------------------- state
    @property
    def state(self):
        return self.load('run', {})

    def update(self, **values):
        state = self.state
        state.update(_clean(values))
        self.save('run', state)
        return state

    # -------------------------------------------------------------- misc
    def plot(self, name):
        return self.path/'plots'/f'{name}.png'

    def product(self, name):
        return self.path/'products'/name

    def log(self, step, message, **data):
        entry = dict(utc=datetime.now(timezone.utc).isoformat(timespec='seconds'), step=step, message=message, **_clean(data))
        with (self.path/'log.jsonl').open('a') as stream:
            stream.write(json.dumps(entry)+'\n')

    def cubes(self):
        """SPW key -> cube record (path, pb, mask, ...), from the inventory or discovery."""
        state = self.state
        if 'cubes' not in state:
            raise ValueError('No cubes registered: run `data` first.')
        return {c['key']: c for c in state['cubes']}
