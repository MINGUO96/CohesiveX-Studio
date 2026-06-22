# Quick start

1. Install the package with `python -m pip install -e .`.
2. Run the GUI with `python -m cohesivex_studio`.
3. Select an Abaqus input file and an output path.
4. Choose the interface scope: grain-boundary, intragranular, or both.
5. Choose the interface-detection backend: NumPy accelerated or pure Python reference.
6. Define cohesive properties for each active family.
7. Generate the solver input file and CAE preview input file.
8. Inspect the generated reports before submitting the Abaqus job.

The solver input file contains `*User element` and `*Uel Property` blocks and should be submitted to Abaqus/Standard with the appropriate Fortran source. The CAE preview file removes UEL blocks and is intended only for visual inspection.

## Backend timing and speedup reports

The NumPy backend is the default because it vectorizes face-key construction and shared-face grouping. The pure-Python backend implements the same topology rules and is kept as a transparent reference and fallback. The GUI provides a backend selector in the Generate page, and the command line provides:

```bash
python -m cohesivex_studio model.inp model_coh.inp --backend numpy
python -m cohesivex_studio model.inp model_coh.inp --backend python
```

To obtain a direct benchmark report, run:

```bash
python -m cohesivex_studio model.inp benchmark_case.inp --interface-scope both --compare-backends
```

The comparison report checks that the two backends produce identical cohesive-element counts, duplicated-node counts, interface-type counts and cohesive-family counts before reporting the total preprocessing speedup and the face-detection speedup.

## Installation note

The repository follows a standard `pyproject.toml` / `src` package layout. If editable installation fails on Windows because of an existing broken console-script wrapper, use `python -m cohesivex_studio` after installation, or install in a clean virtual environment.
