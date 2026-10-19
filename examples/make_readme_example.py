"""
Regenerate the README example (figure + text output) with the packaged default model.
Run from the repository root:  python examples/make_readme_example.py
Writes docs/observe_and_etc_point.png, docs/observe_and_etc_extended.png and docs/example_output_*.txt
"""
import contextlib
import io
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
from sevendt_etc import observe_and_etc

DOCS = Path('docs')
DOCS.mkdir(exist_ok=True)
TARGET = dict(ra=30.968, dec=-4.245, date='2026-09-14', seeing=2.0)
CASES = {'point': dict(value=19.5, snr=5),
         'extended': dict(value=24.0, snr=5, mode='extended', area=100)}
for name, kw in CASES.items():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        observe_and_etc(**TARGET, **kw, save_path=DOCS / f'observe_and_etc_{name}.png')
    (DOCS / f'example_output_{name}.txt').write_text(buf.getvalue())
    print(f'--- {name} ---\n{buf.getvalue()}')
