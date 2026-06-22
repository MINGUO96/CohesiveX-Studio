from pathlib import Path

from cohesivex_studio import __version__
from cohesivex_studio.kernel import compare_backends, generate_cohesive_inp


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


def test_version_is_release_100():
    assert __version__ == "1.0.0"


def test_two_grain_example_generates_one_grain_boundary_element(tmp_path: Path):
    summary = generate_cohesive_inp(
        EXAMPLES / "two_grain" / "two_grain.inp",
        tmp_path / "two_grain_coh.inp",
        interface_scope="grain_boundary",
        fast_mode=True,
        verbose=False,
    )
    assert summary["solid_elements"] == 2
    assert summary["cohesive_elements"] == 1
    assert summary["duplicated_nodes"] == 2
    assert summary["cohesive_family_counts"] == {"GB_COH": 1}
    assert (tmp_path / "two_grain_coh.inp").is_file()
    assert (tmp_path / "two_grain_coh_cae_preview.inp").is_file()
    assert (tmp_path / "two_grain_coh_summary.json").is_file()


def test_two_domain_example_generates_one_intragranular_element(tmp_path: Path):
    summary = generate_cohesive_inp(
        EXAMPLES / "two_domain" / "two_domain.inp",
        tmp_path / "two_domain_coh.inp",
        interface_scope="intragranular",
        fast_mode=False,
        verbose=False,
    )
    assert summary["solid_elements"] == 2
    assert summary["cohesive_elements"] == 1
    assert summary["duplicated_nodes"] == 2
    assert summary["cohesive_family_counts"] == {"INTRA_COH": 1}
    assert (tmp_path / "two_domain_coh.inp").is_file()
    assert (tmp_path / "two_domain_coh_cae_preview.inp").is_file()
    assert (tmp_path / "two_domain_coh_summary.json").is_file()


def test_backends_are_topology_consistent_for_example(tmp_path: Path):
    comparison = compare_backends(
        EXAMPLES / "two_grain" / "two_grain.inp",
        tmp_path / "two_grain_compare",
        interface_scope="grain_boundary",
    )
    assert comparison["topology_consistent"] is True
    assert Path(comparison["comparison_report"]).is_file()
