import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'bin'))
import compare_cluster_tiers as cct  # noqa: E402


def test_parse_summary_tsv(tmp_path):
    p = tmp_path / 'summary.tsv'
    p.write_text('metric\tvalue\nrecall\t0.4\nfp_rate\t0.0\n')
    assert cct.parse_summary_tsv(p) == {'recall': '0.4', 'fp_rate': '0.0'}


def test_build_tier_comparison_tags_rows_with_tier():
    tier_summaries = {
        'P': {'recall': '1.0', 'fp_rate': '0.0'},
        'C+H': {'recall': '0.4', 'fp_rate': '0.0'},
    }
    tier_per_control = {
        'P': [{'control_id': 'POS_HEX1', 'outcome': 'hit'}],
        'C+H': [{'control_id': 'POS_HEX1', 'outcome': 'miss'}],
    }
    per_control, summaries = cct.build_tier_comparison(tier_summaries, tier_per_control)
    assert {'control_id': 'POS_HEX1', 'outcome': 'hit', 'tier': 'P'} in per_control
    assert {'control_id': 'POS_HEX1', 'outcome': 'miss', 'tier': 'C+H'} in per_control
    assert {'tier': 'P', 'recall': '1.0', 'fp_rate': '0.0'} in summaries
    assert {'tier': 'C+H', 'recall': '0.4', 'fp_rate': '0.0'} in summaries
