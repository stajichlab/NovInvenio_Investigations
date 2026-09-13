import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))
sys.path.insert(0, '/bigdata/stajichlab/jstajich/projects/NovInvenio/lib')
import generate_busco_map as gbm  # noqa: E402


def test_filter_wanted_keeps_only_requested_ids():
    rows = [
        ('100036at4751', 'Ncra', 'protA', 500),
        ('999999at4751', 'Ncra', 'protZ', 200),  # not wanted -- must be dropped
        ('100149at4751', 'Ncra', 'protB', 300),
    ]
    wanted = {'100036at4751', '100149at4751'}
    mapping = gbm.filter_wanted(rows, wanted)
    assert mapping == {'100036at4751': 'protA', '100149at4751': 'protB'}


def test_main_writes_two_column_tsv(tmp_path):
    full_table = tmp_path / 'full_table.tsv'
    full_table.write_text(
        "# BUSCO\n"
        "100036at4751\tComplete\tprotA:1-100\t0\t500\n"
        "999999at4751\tComplete\tprotZ:1-50\t0\t200\n"
    )
    wanted_csv = tmp_path / 'controls.csv'
    wanted_csv.write_text(
        "control_id,class,expected_call,anchor_type,anchor,proteome_short,gene_name,expected_origin,source,notes\n"
        "NEG_BUSCO01,negative,core,busco,100036at4751,,x,y,z,note\n"
    )
    out = tmp_path / 'out.tsv'
    import subprocess
    subprocess.run([
        sys.executable, str(REPO / 'bin' / 'generate_busco_map.py'),
        '--full-table', str(full_table), '--species', 'Ncra',
        '--controls', str(wanted_csv), '--output', str(out),
    ], check=True)
    assert out.read_text() == '100036at4751\tprotA\n'
