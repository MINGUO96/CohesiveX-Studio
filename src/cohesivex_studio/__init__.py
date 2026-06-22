"""CohesiveX Studio: graphical preprocessing tools for Abaqus cohesive-zone models."""
__version__ = "1.0.0"

from .kernel import compare_backends, generate_cohesive_inp, read_text_auto, run_self_tests

__all__ = ["compare_backends", "generate_cohesive_inp", "read_text_auto", "run_self_tests"]
