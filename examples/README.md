# CohesiveX Studio examples

These examples are intentionally compact so that installation and CI-style tests can run without Abaqus. They are not intended to represent calibrated material models. Their purpose is to verify parsing, interface detection, node duplication, cohesive element writing and report generation.

## two_grain

`two_grain/two_grain.inp` is a two-element CPE4 model with two element sets, `GRAIN-1` and `GRAIN-2`.

Run:

```bash
python -m cohesivex_studio examples/two_grain/two_grain.inp output/two_grain_coh.inp --interface-scope grain_boundary
```

Expected summary:

```text
solid_elements: 2
num_domains: 2
generated_cohesive_elements: 1
duplicated_nodes: 2
cohesive_family_counts: {"GB_COH": 1}
```

## two_domain

`two_domain/two_domain.inp` is a two-element CPE4 model without grain sets. It is treated as one domain and is useful for intragranular insertion tests.

Run:

```bash
python -m cohesivex_studio examples/two_domain/two_domain.inp output/two_domain_coh.inp --interface-scope intragranular
```

Expected summary:

```text
solid_elements: 2
num_domains: 1
generated_cohesive_elements: 1
duplicated_nodes: 2
cohesive_family_counts: {"INTRA_COH": 1}
```

## Backend comparison

Run:

```bash
python -m cohesivex_studio examples/two_grain/two_grain.inp output/two_grain_compare.inp --interface-scope grain_boundary --compare-backends
```

The comparison report should state that the NumPy and pure-Python backends are topology-consistent. On very small examples, the pure-Python backend may appear faster because NumPy setup overhead dominates; acceleration is expected to become useful for larger meshes.
