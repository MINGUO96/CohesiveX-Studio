# SoftwareX preparation notes

Recommended submission materials:

- Public GitHub repository.
- README, LICENSE, CITATION.cff, and pyproject.toml.
- Reproducible examples with input files and generated reports.
- Regression tests that run without Abaqus.
- A manuscript describing the software scope, algorithm, impact, and limitations.

CohesiveX Studio should be described as a preprocessing platform. It generates cohesive topology and Abaqus keyword blocks, while users provide their own solver-side Fortran subroutine.

## Benchmark reporting

For the SoftwareX manuscript, backend acceleration should be reported using the `--compare-backends` command. The comparison is meaningful because both backends use the same parser, interface-scope rules, node-duplication logic and cohesive-family definitions. Only the shared-face detection backend is changed. The generated comparison JSON and TXT files provide reproducible timing data for the manuscript figure.
