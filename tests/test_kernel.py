from pathlib import Path

from cohesivex_studio.kernel import compare_backends, generate_cohesive_inp, run_self_tests


def test_kernel_regression_suite():
    run_self_tests()


def test_timing_and_backend_comparison(tmp_path: Path):
    inp = tmp_path / "two_domain.inp"
    inp.write_text(
        """
*Heading
*Node
1, 0.0, 0.0
2, 1.0, 0.0
3, 2.0, 0.0
4, 0.0, 1.0
5, 1.0, 1.0
6, 2.0, 1.0
*Element, type=CPE4
1, 1, 2, 5, 4
2, 2, 3, 6, 5
*Solid Section, elset=ALL, material=M1
,
*Material, name=M1
*Elastic
1., 0.3
*Step
*Static
0.1, 1.
*End Step
""".strip()
        + "\n",
        encoding="utf-8",
    )
    summary = generate_cohesive_inp(inp, tmp_path / "single.inp", interface_scope="intragranular", fast_mode=False, verbose=False)
    assert summary["backend"] == "pure_python"
    assert summary["total_preprocessing_time_seconds"] > 0
    assert "face_detection" in summary["timings_seconds"]

    comparison = compare_backends(inp, tmp_path / "bench", interface_scope="intragranular")
    assert comparison["topology_consistent"] is True
    assert Path(comparison["comparison_report"]).is_file()
