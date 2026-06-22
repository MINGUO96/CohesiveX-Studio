# CohesiveX Studio

**Version:** v1.0.0  
**Author:** guomin  
**Repository:** <https://github.com/MINGUO96/CohesiveX-Studio>  
**License:** BSD 3-Clause License

CohesiveX Studio is an open-source Python platform for generating zero-thickness cohesive interface elements in existing Abaqus INP finite-element models. It is designed for inter-domain interfaces such as grain boundaries or phase boundaries, intra-domain interfaces such as intragranular crack paths, and combined inter-/intra-domain cohesive-zone simulations.

The software provides both a graphical user interface and command-line workflow. It reads an existing Abaqus input file, detects conformal shared element edges or faces, classifies the interface type, duplicates only the nodes required by the selected cohesive topology, assigns cohesive elements to independent property families, and writes solver-ready Abaqus input files together with CAE-preview files and machine-readable reports.

CohesiveX Studio does **not** replace Abaqus or the user's Fortran UEL. It is a preprocessing and verification tool. Abaqus, Fortran compilers, UMAT/VUMAT files and UEL source files must be supplied by the user when running the final simulation.

## Why this software is useful

Manual cohesive-zone preprocessing in Abaqus is time-consuming and easy to get wrong. A user must identify interface faces, split nodes, update solid-element connectivity, define cohesive user elements, assign cohesive properties, preserve boundary conditions on duplicated nodes, and verify face orientation. These operations become especially difficult in large 2D or 3D meshes.

CohesiveX Studio addresses this problem through a conservative and auditable transformation:

- It preserves the original material cards, solid sections, steps, amplitudes, loads, boundary conditions, solver controls and output requests.
- It inserts only the nodes, updated solid connectivity and cohesive UEL keyword blocks required by the selected interface scope.
- It writes separate solver and CAE-preview INP files.
- It creates CSV, JSON and TXT reports for topology inspection and reproducibility.
- It supports a NumPy-accelerated backend and a pure-Python reference backend for verification.

## Main features

- GUI workflow based on Tkinter.
- CLI workflow for batch processing and regression tests.
- Inter-domain cohesive insertion for grain-boundary or phase-boundary interfaces.
- Intra-domain cohesive insertion for intragranular or conventional elastic-plastic cracking paths.
- Combined inter-/intra-domain insertion with independent cohesive property families.
- Coordinate-based conformal face hashing with adjustable tolerance.
- NumPy backend for large meshes and pure-Python backend for reference comparison.
- Deterministic intragranular sampling with a user-provided random seed.
- Solver-ready Abaqus INP output with `*User element`, cohesive `*Element` and `*Uel Property` blocks.
- Abaqus/CAE preview INP output without unsupported UEL keyword blocks.
- Reports for generated cohesive elements, duplicated nodes, interface measures, normal vectors, warnings and timing data.
- Built-in self-tests that run without Abaqus.

## Repository structure

```text
CohesiveX-Studio/
├── docs/
│   ├── quick_start.md
│   ├── softwarex_notes.md
│   └── theory.md
├── examples/
│   ├── README.md
│   ├── two_domain/
│   │   └── two_domain.inp
│   └── two_grain/
│       └── two_grain.inp
├── src/
│   └── cohesivex_studio/
│       ├── __init__.py
│       ├── __main__.py
│       ├── app.py
│       └── kernel.py
├── tests/
│   ├── test_examples.py
│   └── test_kernel.py
├── .gitignore
├── CITATION.cff
├── LICENSE
├── MANIFEST.in
├── README.md
└── pyproject.toml
```

## Requirements

CohesiveX Studio requires Python 3.9 or later.

Required Python packages:

```text
numpy
matplotlib
```

Tkinter is required for the graphical interface. It is included in most standard Python distributions. On some Linux systems, install it through the system package manager, for example:

```bash
sudo apt-get install python3-tk
```

Optional tools for development and testing:

```text
pytest
build
```

## Installation from GitHub

Clone the repository:

```bash
git clone https://github.com/MINGUO96/CohesiveX-Studio.git
cd CohesiveX-Studio
```

Install the package in editable mode:

```bash
python -m pip install -e .
```

For development and tests:

```bash
python -m pip install -e ".[dev]"
```

Check the installed version:

```bash
python -m cohesivex_studio --version
```

Expected output:

```text
CohesiveX Studio 1.0.0
```

## Launching the graphical interface

Use:

```bash
python -m cohesivex_studio
```

or, after installation:

```bash
cohesivex-studio
```

The interface includes project setup, cohesive-law settings, generation options, report inspection, visualization and Abaqus job preparation pages.

## Command-line usage

Generate grain-boundary cohesive elements for a two-grain example:

```bash
python -m cohesivex_studio examples/two_grain/two_grain.inp output/two_grain_coh.inp --interface-scope grain_boundary
```

Generate intragranular cohesive elements for a single-domain example:

```bash
python -m cohesivex_studio examples/two_domain/two_domain.inp output/two_domain_coh.inp --interface-scope intragranular
```

Generate both grain-boundary and intragranular cohesive families:

```bash
python -m cohesivex_studio model.inp output/model_coh.inp \
  --interface-scope both \
  --uel-elset GB_COH \
  --intra-elset INTRA_COH \
  --props 3 1e7 1e7 100 100 5 5 2 1 \
  --intra-props 3 1e7 1e7 180 180 10 10 2 1
```

Use the pure-Python reference backend:

```bash
python -m cohesivex_studio model.inp output/model_python.inp --interface-scope both --backend python
```

Compare the NumPy and pure-Python backends:

```bash
python -m cohesivex_studio model.inp output/backend_case.inp --interface-scope both --compare-backends
```

The comparison command writes NumPy and pure-Python outputs, then reports whether the cohesive topology is consistent between the two backends.

## Built-in test cases

The repository contains two compact Abaqus INP examples:

### `examples/two_grain/two_grain.inp`

A two-element, two-domain CPE4 model. It is used to test grain-boundary cohesive insertion.

Expected behavior:

- two solid elements are parsed;
- two domain sets are detected;
- one inter-domain interface is found;
- one `GB_COH` cohesive element is generated;
- two nodes are duplicated;
- solver INP, CAE-preview INP and reports are written.

### `examples/two_domain/two_domain.inp`

A two-element, single-domain CPE4 model. It is used to test intragranular cohesive insertion.

Expected behavior:

- two solid elements are parsed;
- the mesh is treated as one domain;
- one intra-domain interface is found;
- one `INTRA_COH` cohesive element is generated;
- two nodes are duplicated;
- solver INP, CAE-preview INP and reports are written.

## Running tests

Run the self-test without Abaqus:

```bash
python -m cohesivex_studio --self-test
```

Run the full pytest suite:

```bash
pytest -q
```

The tests check kernel regression behavior, example generation, backend comparison and version consistency. Abaqus is not required for these tests.

## Output files

For an output path such as:

```text
output/model_coh.inp
```

CohesiveX Studio may generate:

```text
model_coh.inp                         solver-ready Abaqus input file
model_coh_cae_preview.inp             Abaqus/CAE preview input file
model_coh_summary.json                structured generation summary
model_coh_mesh_check.txt              human-readable mesh and warning report
model_coh_grain_boundary_table.csv    cohesive interface table
model_coh_duplicated_nodes.csv        original-to-duplicated node map
```

Generated result files are ignored by `.gitignore` unless intentionally added.

## Cohesive theory summary

For a zero-thickness cohesive interface, the displacement jump is

```text
Delta u = u+ - u-
```

The local normal and shear separations are computed by projecting this jump onto the interface basis. Before damage initiation, the traction response is represented by a penalty stiffness relation:

```text
t = K delta
```

For a bilinear cohesive law, the pure-mode final separations are related to the fracture energies by:

```text
delta_n_f = 2 G_Ic  / S_I
delta_s_f = 2 G_IIc / S_II
```

The mixed-mode fracture energy can be described by the Benzeggagh-Kenane expression:

```text
G_c = G_Ic + (G_IIc - G_Ic) (G_s / G_T)^eta
```

CohesiveX Studio writes a nine-parameter UEL property vector for each cohesive family:

```text
[mode, K_I, K_II, S_I, S_II, G_Ic, G_IIc, eta, HEIGHT]
```

The exact constitutive update is implemented in the user's Fortran UEL. CohesiveX Studio provides the mesh topology, Abaqus keyword blocks and verification reports needed to use such a UEL safely.

## Limitations

The v1.0.0 release targets conformal Abaqus meshes where internal interfaces are shared by exactly two owner elements. Non-manifold faces are reported and skipped. Non-matching interfaces, mixed unsupported element blocks and automatic UEL constitutive-code generation are outside the current release scope.

## SoftwareX submission notes

For a SoftwareX submission, cite the public repository, the release tag and the associated manuscript. Recommended metadata:

```text
Current code version: v1.0.0
Repository: https://github.com/MINGUO96/CohesiveX-Studio
License: BSD-3-Clause
Language: Python
Main dependencies: NumPy, Matplotlib, Tkinter
```

## License

CohesiveX Studio is distributed under the BSD 3-Clause License.

Abaqus and user-provided Fortran subroutines are not part of this repository and are not redistributed under this license.

## Citation

If you use CohesiveX Studio in academic work, please cite the associated SoftwareX article and this repository. A machine-readable citation file is provided in:

```text
CITATION.cff
```

## Contact

Please use the GitHub issue tracker for questions, bug reports or feature requests:

<https://github.com/MINGUO96/CohesiveX-Studio/issues>
