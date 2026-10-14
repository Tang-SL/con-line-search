"""Link existing local pilot inputs; no large data are copied or committed."""
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('workspace', type=Path, help='Existing CONs workspace root')
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
relative = Path('targets/ic860/alma/espada_b5_all')
items = ['inputs', 'qa/spectra.npz']
for name in items:
    source = (args.workspace / relative / name).resolve()
    if not source.exists():
        parser.error(f'Missing prepared input: {source}')
for name in items:
    source = (args.workspace / relative / name).resolve()
    destination = root / relative / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        if destination.resolve() != source:
            parser.error(f'Refusing to replace existing path: {destination}')
    else:
        destination.symlink_to(source, target_is_directory=source.is_dir())
    print(f'Linked {name}')
