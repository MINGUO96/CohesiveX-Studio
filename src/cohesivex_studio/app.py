# -*- coding: utf-8 -*-
"""Graphical interface for CohesiveX Studio.

The GUI is imported lazily so that the command-line kernel and regression tests
can run in headless environments without Tkinter.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple
import csv
import json
import math
import os
import queue
import shlex
import signal
import subprocess
import sys
import threading
import time

from . import __version__
from .kernel import compare_backends, generate_cohesive_inp, read_text_auto, run_self_tests

# =============================================================================
# Tkinter graphical user interface
# =============================================================================


try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    TK_AVAILABLE = True
except Exception:  # pragma: no cover
    class _MissingTk:
        class Tk:
            pass
        TclError = Exception

    class _MissingTtk:
        class Frame:
            pass

    tk = _MissingTk()  # type: ignore
    ttk = _MissingTtk()  # type: ignore
    filedialog = None  # type: ignore
    messagebox = None  # type: ignore
    TK_AVAILABLE = False

try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    MATPLOTLIB_AVAILABLE = True
except Exception:  # pragma: no cover
    FigureCanvasTkAgg = None  # type: ignore
    NavigationToolbar2Tk = None  # type: ignore
    Figure = None  # type: ignore
    MATPLOTLIB_AVAILABLE = False


PLATFORM_NAME = "CohesiveX Studio"
PLATFORM_SUBTITLE = "Abaqus inter-/intra-domain cohesive-zone preprocessing platform"
PREF_PATH = Path.home() / ".cohesivex_studio.json"



COHESIVE_THEORY_FORMULAS = r"""
Cohesive-zone theory used by CohesiveX Studio
=============================================

1. Interface kinematics
-----------------------
A zero-thickness cohesive element connects two coincident surfaces.  Let u+ and
u- be the displacement vectors on the positive and negative sides of the
interface.  The displacement jump is

    Delta u = u+ - u-.

With the local interface basis {n, s, t}, the opening and shear separations are

    delta_n  = Delta u dot n
    delta_s  = Delta u dot s
    delta_t  = Delta u dot t
    delta_sh = sqrt(delta_s^2 + delta_t^2).

In 2D, only delta_n and one shear component are active.  The reported 2D
interface measure is edge_length x HEIGHT.  In 3D, the interface measure is the
polygonal face area.

2. Elastic penalty response
---------------------------
Before damage starts, the local traction vector is computed from the penalty
stiffnesses:

    t_n = K_I  delta_n
    t_s = K_II delta_s
    t_t = K_II delta_t.

Equivalently,

    {t} = K {delta},

where

    K = diag(K_I, K_II)          for 2D,
    K = diag(K_I, K_II, K_II)    for 3D.

A penalty stiffness that is too small gives artificial compliance.  A penalty
stiffness that is too large can make the global stiffness matrix ill-conditioned
and can reduce convergence robustness.

3. Damage initiation
--------------------
For a bilinear law, pure-mode initiation separations are

    delta_n0 = S_I  / K_I,
    delta_s0 = S_II / K_II.

A common mixed-mode nominal stress index is

    F_ini = (<t_n>/S_I)^2 + (t_sh/S_II)^2,

where

    <t_n> = max(t_n, 0),
    t_sh  = sqrt(t_s^2 + t_t^2).

Damage initiation occurs when

    F_ini = 1.

4. Fracture energy and final separation
---------------------------------------
For a triangular bilinear traction-separation law,

    G_Ic  = 1/2 S_I  delta_nf,
    G_IIc = 1/2 S_II delta_sf.

Therefore,

    delta_nf = 2 G_Ic  / S_I,
    delta_sf = 2 G_IIc / S_II.

5. Linear softening damage variable
-----------------------------------
For an effective separation delta, damage initiation separation delta0 and final
failure separation deltaf, a standard linear-softening damage variable is

    D = deltaf (delta - delta0) / [ delta (deltaf - delta0) ],

with

    0 <= D <= 1.

The degraded traction response is

    t_n = (1 - D) K_I  delta_n,
    t_s = (1 - D) K_II delta_s,
    t_t = (1 - D) K_II delta_t.

6. Benzeggagh-Kenane mixed-mode fracture criterion
--------------------------------------------------
The modal work terms are

    G_I   = integral <t_n> d(delta_n),
    G_II  = integral t_s d(delta_s),
    G_III = integral t_t d(delta_t).

Define

    G_T = G_I + G_II + G_III,
    G_s = G_II + G_III.

The B-K mixed-mode critical energy is

    G_c = G_Ic + (G_IIc - G_Ic) (G_s/G_T)^eta.

The mixed-mode failure index is

    F_fail = G_T / G_c.

7. Cohesive property vector
---------------------------
CohesiveX writes a nine-parameter UEL property vector for each cohesive family:

    PROPS = [mode, K_I, K_II, S_I, S_II, G_Ic, G_IIc, eta, HEIGHT].

mode = 1 activates opening mode, mode = 2 activates shear mode, and mode = 3
activates mixed mode.  HEIGHT is required for 2D cohesive elements to convert a
line measure into an area-like measure in the UEL and in the reports.
"""

UEL_MAPPING_NOTES = r"""
Abaqus UEL mapping and cohesive families
========================================

1. User-element definition
--------------------------
For a solid face with N_f face nodes, CohesiveX inserts a zero-thickness UEL with
2 N_f nodes:

    *User element, nodes=2N_f, type=<UEL type>, properties=9,
    coordinates=NDIM, variables=NIP * NSVARS_PER_IP
      1, 2           for 2D
      1, 2, 3        for 3D

For the supplied OXFORD cohesive UEL, NSVARS_PER_IP = 1 because only damage is
stored at each cohesive integration point.  If the Fortran UEL is extended to
store additional SDVs, the GUI value and the generated variables count must be
updated together.

2. Cohesive connectivity
------------------------
For a shared interface between elements e1 and e2, the cohesive connectivity is

    e_coh, n1+, n2+, ..., nN+, n1-, n2-, ..., nN-.

The first group belongs to one side of the split interface, and the second group
belongs to the neighbouring side.  CohesiveX compares the interface normal with
the vector from the centroid of e1 to the centroid of e2.  If the orientation is
opposite, the cohesive node order is reversed according to the element topology.

3. Interface scope selected by the user
---------------------------------------
Let d(e) be the domain ID of a solid element.  In a polycrystal, d(e) is the
grain ID.  In a conventional elastic-plastic model without grain sets, all
parsed solid elements are assigned to one domain.

For two neighbouring solid elements e1 and e2:

    grain-boundary / inter-domain interface:    d(e1) != d(e2),
    intragranular / intra-domain interface:     d(e1) == d(e2).

The GUI offers three insertion scopes:

    Grain-boundary only              insert only when d(e1) != d(e2),
    Intragranular only               insert only when d(e1) == d(e2),
    Grain-boundary + intragranular   insert both interface classes.

4. Independent cohesive property families
-----------------------------------------
The same UEL formulation can be used with multiple elsets and different UEL
property vectors:

    *Element, type=U1, elset=GB_COH
    ...
    *Uel Property, elset=GB_COH
    PROPS_GB

    *Element, type=U1, elset=INTRA_COH
    ...
    *Uel Property, elset=INTRA_COH
    PROPS_INTRA

This allows a model to use weaker grain-boundary properties and stronger
intragranular properties, or any other calibrated pair of cohesive parameters.

5. Solver input and CAE preview input
-------------------------------------
Abaqus/CAE cannot import input files containing *User element and *Uel Property
blocks.  CohesiveX therefore writes two files:

    solver INP       contains UEL cohesive elements and is submitted to Abaqus,
    CAE preview INP  removes only UEL blocks for visual mesh inspection.

Submit the solver input with a compatible Fortran subroutine, for example:

    abaqus job=<job_name> input=<solver.inp> user=<subroutine.f> cpus=<N> interactive

6. Cohesive-only preservation rule
----------------------------------
CohesiveX does not regenerate the whole Abaqus model.  It preserves original
materials, solid sections, steps, amplitudes, loads, boundary conditions,
controls, output requests and time incrementation.  It modifies or inserts only
nodes, the parsed solid element connectivity, UEL cohesive blocks, UEL
properties, cohesive elsets and optional duplicated-node supplements to existing
Nsets.
"""

ALGORITHM_VERIFICATION_NOTES = r"""
Algorithm, acceleration and verification
========================================

1. Coordinate-based face hashing
--------------------------------
For each solid element e and each local face f, CohesiveX constructs a face key
from quantized nodal coordinates:

    key(e, f) = sort( round(x_i / tolerance) for nodes i on face f ).

The tolerance is an absolute geometric tolerance.  It should be selected based
on the coordinate unit and coordinate precision of the Abaqus input file.  This
coordinate key allows conformal interfaces with duplicated node IDs to be
identified, while non-matching interfaces are not silently merged.

2. NumPy accelerated backend
----------------------------
When NumPy is available, face-node arrays are assembled as

    faces: (number_of_elements * number_of_faces, number_of_face_nodes),

and unique face keys are computed using vectorized array operations.  The pure
Python backend uses the same face-key definition.  The self-test checks that the
fast NumPy backend and the pure Python backend generate identical cohesive
counts for the regression meshes.

3. Manifold requirement
-----------------------
The current release supports conformal manifold meshes, where each internal
interface face has exactly two owner elements.  If a face is shared by more than
two elements, CohesiveX skips it and reports a warning instead of attempting
pairwise insertion.  This avoids creating non-physical cohesive topology in
invalid or thickness-collapsed meshes.

4. Intragranular sampling
-------------------------
If the intragranular fraction is f, CohesiveX inserts an intra-domain cohesive
element only when

    r(seed, interface_id) <= f,

where r is a deterministic hash-based pseudo-random number in [0, 1).  When
f < 1 and the user does not provide a random seed, seed = 0 is used and a warning
is written.  This makes the sampling reproducible and independent of Python
iteration order or NumPy availability.

5. Node splitting by selected interfaces
----------------------------------------
The selected cohesive interfaces define which solid faces should be separated.
Solid elements connected through unselected internal faces remain in the same
connectivity component.  Nodes are duplicated per component, so only the chosen
interface network is opened.  This is more general than purely grain-based node
splitting and enables grain-boundary-only, intragranular-only and combined
insertion scopes.

6. Recommended verification checks
----------------------------------
For a publishable generation run, verify that

    fast_mode=True and fast_mode=False give the same cohesive topology,
    the selected interface scope gives the expected cohesive count,
    the solver INP preserves non-cohesive Abaqus keywords,
    the CAE preview INP imports without UEL keywords,
    duplicated boundary nodes are supplemented into boundary Nsets when enabled,
    the mesh-check report contains no unexpected warnings,
    grain-boundary and intragranular families have the intended elsets and properties.

7. Current scope and safe claims
-------------------------------
This release targets one supported homogeneous solid element block per generation
run.  Mixed element-type blocks and non-matching interfaces are treated as future
extensions.  The platform is UEL-agnostic: it creates the topology and Abaqus
keyword structure, while the traction-separation law is defined by the user's
Fortran UEL.
"""

class TextEditor(ttk.Frame):  # type: ignore[misc]
    """Small scrollable text editor used by the GUI."""

    def __init__(self, master: Any, wrap: str = "none", height: int = 12) -> None:
        super().__init__(master)
        self.text = tk.Text(
            self,
            wrap=wrap,
            undo=True,
            height=height,
            font=("Consolas", 10),
            bg="#ffffff",
            fg="#0f172a",
            insertbackground="#0f172a",
        )
        y = ttk.Scrollbar(self, orient="vertical", command=self.text.yview)
        x = ttk.Scrollbar(self, orient="horizontal", command=self.text.xview)
        self.text.configure(yscrollcommand=y.set, xscrollcommand=x.set)
        self.text.grid(row=0, column=0, sticky="nsew")
        y.grid(row=0, column=1, sticky="ns")
        x.grid(row=1, column=0, sticky="ew")
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

    def set(self, value: str) -> None:
        self.text.configure(state="normal")
        self.text.delete("1.0", "end")
        self.text.insert("1.0", value)

    def get(self) -> str:
        return self.text.get("1.0", "end-1c")

    def append(self, value: str) -> None:
        self.text.configure(state="normal")
        self.text.insert("end", value)
        self.text.see("end")


class PlotCanvas(ttk.Frame):  # type: ignore[misc]
    """Embedded Matplotlib canvas used by the Visualization page."""

    def __init__(self, master: Any, figsize: Tuple[float, float] = (7.8, 5.6)) -> None:
        super().__init__(master)
        if not MATPLOTLIB_AVAILABLE:
            self.fig = None
            ttk.Label(
                self,
                text="Matplotlib is not available in this Python environment.\nInstall matplotlib to enable the visualization dashboard.",
                style="Muted.TLabel",
                justify="center",
            ).pack(fill="both", expand=True, padx=20, pady=20)
            return
        self.fig = Figure(figsize=figsize, dpi=100, facecolor="white")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        self.toolbar = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
        self.toolbar.grid(row=1, column=0, sticky="ew")
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

    @property
    def figure(self) -> Any:
        return self.fig

    def draw(self) -> None:
        if MATPLOTLIB_AVAILABLE and self.fig is not None:
            self.canvas.draw_idle()


class CohesiveXStudio(tk.Tk):  # type: ignore[misc]
    """All-in-one GUI for cohesive element insertion.

    The design follows the workflow-oriented interface: a left navigation bar,
    independent pages, structured parameter editors, a log window, file preview
    and report inspection.  It calls the scientific kernel
    ``generate_cohesive_inp`` while keeping the GUI layer separate from the
    headless command-line and testing workflows.
    """

    NAV_NAMES = [
        "Dashboard",
        "Project",
        "Cohesive Law",
        "Generate",
        "Reports",
        "Visualization",
        "Run Abaqus",
        "Files",
        "Theory",
    ]

    PROP_LABELS = [
        ("mode", "Fracture mode: 1 opening, 2 shear, 3 mixed"),
        ("KI", "Penalty stiffness for opening mode"),
        ("KII", "Penalty stiffness for shearing mode"),
        ("SI", "Maximum normal stress before failure"),
        ("SII", "Maximum shear stress before failure"),
        ("GCI", "Fracture toughness for opening mode"),
        ("GCII", "Fracture toughness for shearing mode"),
        ("ETA", "Benzeggagh-Kenane exponent"),
        ("HEIGHT", "Height for 2D elements"),
    ]

    PROP_DEFAULTS = ["3", "1e7", "1e7", "100", "100", "5", "5", "2", "1"]

    def __init__(self) -> None:
        super().__init__()
        self.title(PLATFORM_NAME)
        self.geometry("1380x860")
        self.minsize(1050, 680)

        self.log_queue: queue.Queue[str] = queue.Queue()
        self.worker_thread: Optional[threading.Thread] = None
        self.abaqus_process: Optional[subprocess.Popen[str]] = None
        self.summary: Dict[str, Any] = {}
        self.selected_file: Optional[Path] = None

        self.prefs = self._load_prefs()
        self._setup_style()
        self._build_layout()
        self._load_defaults()
        self.after(150, self._drain_log_queue)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    # ------------------------------------------------------------------ style
    def _setup_style(self) -> None:
        self.colors = {
            "bg": "#f6f8fb",
            "panel": "#ffffff",
            "ink": "#0f172a",
            "muted": "#64748b",
            "teal": "#0f766e",
            "teal_dark": "#115e59",
            "line": "#dbe4ef",
            "warn": "#b45309",
            "bad": "#b91c1c",
        }
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        self.configure(bg=self.colors["bg"])
        style.configure("TFrame", background=self.colors["bg"])
        style.configure("Panel.TFrame", background=self.colors["panel"], relief="solid", borderwidth=1)
        style.configure("Header.TFrame", background=self.colors["teal"])
        style.configure("Header.TLabel", background=self.colors["teal"], foreground="white", font=("Segoe UI", 18, "bold"))
        style.configure("SubHeader.TLabel", background=self.colors["teal"], foreground="#d1fae5", font=("Segoe UI", 9, "bold"))
        style.configure("Title.TLabel", background=self.colors["bg"], foreground=self.colors["ink"], font=("Segoe UI", 16, "bold"))
        style.configure("CardTitle.TLabel", background=self.colors["panel"], foreground=self.colors["ink"], font=("Segoe UI", 11, "bold"))
        style.configure("Muted.TLabel", background=self.colors["panel"], foreground=self.colors["muted"], font=("Segoe UI", 9))
        style.configure("TLabel", background=self.colors["bg"], foreground=self.colors["ink"], font=("Segoe UI", 9))
        style.configure("TButton", font=("Segoe UI", 9), padding=(9, 5))
        style.configure("Accent.TButton", background=self.colors["teal"], foreground="white", font=("Segoe UI", 9, "bold"))
        style.map("Accent.TButton", background=[("active", self.colors["teal_dark"])])
        style.configure("Nav.TButton", anchor="w", padding=(14, 9), font=("Segoe UI", 10, "bold"))
        style.configure("Active.Nav.TButton", background=self.colors["teal"], foreground="white", anchor="w", padding=(14, 9), font=("Segoe UI", 10, "bold"))
        style.configure("Treeview", background="#ffffff", fieldbackground="#ffffff", foreground=self.colors["ink"], rowheight=24, font=("Segoe UI", 9))
        style.configure("Treeview.Heading", font=("Segoe UI", 9, "bold"))

    # ---------------------------------------------------------------- layout
    def _build_layout(self) -> None:
        header = ttk.Frame(self, style="Header.TFrame")
        header.pack(side="top", fill="x")
        ttk.Label(header, text=PLATFORM_NAME, style="Header.TLabel").pack(side="left", padx=18, pady=14)
        ttk.Label(
            header,
            text="GUI-assisted Abaqus platform for inter-/intra-domain cohesive-zone generation, verification and visualization",
            style="SubHeader.TLabel",
        ).pack(side="left", padx=10)
        ttk.Label(header, text=f"v{__version__} | Research software platform", style="SubHeader.TLabel").pack(side="right", padx=20)

        body = ttk.Frame(self)
        body.pack(side="top", fill="both", expand=True)

        self.nav = ttk.Frame(body, style="Panel.TFrame", width=180)
        self.nav.pack(side="left", fill="y", padx=(8, 4), pady=8)
        self.nav.pack_propagate(False)
        ttk.Label(self.nav, text="Workflow", style="CardTitle.TLabel").pack(anchor="w", padx=16, pady=(18, 8))
        self.nav_buttons: Dict[str, ttk.Button] = {}
        for name in self.NAV_NAMES:
            b = ttk.Button(self.nav, text=name, style="Nav.TButton", command=lambda n=name: self.show_page(n))
            b.pack(fill="x", padx=10, pady=2)
            self.nav_buttons[name] = b
        self.status_label = ttk.Label(self.nav, text="Status\nReady", style="Muted.TLabel", justify="left")
        self.status_label.pack(anchor="w", padx=18, pady=(26, 8))

        self.content = ttk.Frame(body)
        self.content.pack(side="left", fill="both", expand=True, padx=(4, 8), pady=8)
        self.pages: Dict[str, ttk.Frame] = {}
        for name in self.NAV_NAMES:
            frame = ttk.Frame(self.content)
            frame.grid(row=0, column=0, sticky="nsew")
            self.pages[name] = frame
        self.content.grid_rowconfigure(0, weight=1)
        self.content.grid_columnconfigure(0, weight=1)

        self._build_dashboard_page()
        self._build_project_page()
        self._build_cohesive_law_page()
        self._build_generate_page()
        self._build_reports_page()
        self._build_visualization_page()
        self._build_run_abaqus_page()
        self._build_files_page()
        self._build_theory_page()
        self.show_page("Project")

    def show_page(self, name: str) -> None:
        for n, b in self.nav_buttons.items():
            b.configure(style="Active.Nav.TButton" if n == name else "Nav.TButton")
        self.pages[name].tkraise()
        if name == "Reports":
            self.refresh_reports()
        if name == "Visualization":
            self.refresh_visualization()
        if name == "Files":
            self.refresh_file_browser()

    def card(self, master: Any, title: str, subtitle: str = "") -> ttk.Frame:
        f = ttk.Frame(master, style="Panel.TFrame", padding=12)
        ttk.Label(f, text=title, style="CardTitle.TLabel").pack(anchor="w")
        if subtitle:
            ttk.Label(f, text=subtitle, style="Muted.TLabel").pack(anchor="w", pady=(2, 8))
        return f

    def _entry_row(self, master: Any, label: str, var: tk.StringVar, width: int = 14) -> ttk.Entry:
        row = ttk.Frame(master, style="Panel.TFrame")
        row.pack(fill="x", pady=2)
        ttk.Label(row, text=label, style="Muted.TLabel", width=22).pack(side="left")
        ent = ttk.Entry(row, textvariable=var, width=width)
        ent.pack(side="left", fill="x", expand=True, padx=(4, 0))
        return ent

    def _path_row(self, master: Any, label: str, var: tk.StringVar, cmd: Callable[[], None]) -> None:
        row = ttk.Frame(master, style="Panel.TFrame")
        row.pack(fill="x", pady=4)
        ttk.Label(row, text=label, style="Muted.TLabel", width=22).pack(side="left")
        ttk.Entry(row, textvariable=var).pack(side="left", fill="x", expand=True, padx=6)
        ttk.Button(row, text="Browse", command=cmd).pack(side="left")

    def _combo_row(self, master: Any, label: str, var: tk.StringVar, values: Sequence[str]) -> ttk.Combobox:
        row = ttk.Frame(master, style="Panel.TFrame")
        row.pack(fill="x", pady=2)
        ttk.Label(row, text=label, style="Muted.TLabel", width=22).pack(side="left")
        cb = ttk.Combobox(row, textvariable=var, values=list(values), state="readonly")
        cb.pack(side="left", fill="x", expand=True, padx=(4, 0))
        return cb

    # --------------------------------------------------------------- defaults
    def _load_prefs(self) -> Dict[str, Any]:
        try:
            return json.loads(PREF_PATH.read_text(encoding="utf-8"))
        except Exception:
            return {}

    def _save_prefs(self) -> None:
        try:
            data = {
                "base_dir": self.var_base.get(),
                "input_inp": self.var_input.get(),
                "output_inp": self.var_output.get(),
                "uel_file": self.var_uel.get(),
                "abaqus_cmd": self.var_abaqus_cmd.get(),
                "job_name": self.var_job.get(),
                "run_inp": self.var_run_inp.get(),
                "run_fortran": self.var_run_fortran.get(),
                "run_cpus": self.var_cpus.get(),
                "uel_type": self.var_uel_type.get() if hasattr(self, "var_uel_type") else "U1",
                "uel_elset": self.var_uel_elset.get() if hasattr(self, "var_uel_elset") else "GB_COH",
                "intra_elset": self.var_intra_elset.get() if hasattr(self, "var_intra_elset") else "INTRA_COH",
                "interface_scope": self.var_interface_scope.get() if hasattr(self, "var_interface_scope") else "Grain-boundary only",
                "backend_mode": self.var_backend_mode.get() if hasattr(self, "var_backend_mode") else "NumPy accelerated",
                "geometry_tolerance": self.var_tolerance.get() if hasattr(self, "var_tolerance") else "1e-10",
            }
            PREF_PATH.write_text(json.dumps(data, indent=2), encoding="utf-8")
        except Exception:
            pass

    def _load_defaults(self) -> None:
        base = Path(self.prefs.get("base_dir", Path.cwd())).expanduser()
        if not base.exists():
            base = Path.cwd()
        self.var_base.set(str(base))
        self.var_input.set(self.prefs.get("input_inp", str(base / "abq_hex.inp")))
        self.var_output.set(self.prefs.get("output_inp", str(base / "Job_coh.inp")))
        self.var_uel.set(self.prefs.get("uel_file", str(base / "UEL.f")))
        self.var_abaqus_cmd.set(self.prefs.get("abaqus_cmd", "abaqus"))
        self.var_job.set(self.prefs.get("job_name", "Job_coh"))
        self.var_run_inp.set(self.prefs.get("run_inp", self.var_output.get()))
        self.var_run_fortran.set(self.prefs.get("run_fortran", self.var_uel.get()))
        self.var_cpus.set(self.prefs.get("run_cpus", "1"))
        if hasattr(self, "var_uel_type"):
            self.var_uel_type.set(self.prefs.get("uel_type", "U1"))
        if hasattr(self, "var_uel_elset"):
            self.var_uel_elset.set(self.prefs.get("uel_elset", "GB_COH"))
        if hasattr(self, "var_intra_elset"):
            self.var_intra_elset.set(self.prefs.get("intra_elset", "INTRA_COH"))
        if hasattr(self, "var_interface_scope"):
            self.var_interface_scope.set(self.prefs.get("interface_scope", "Grain-boundary only"))
        if hasattr(self, "var_backend_mode"):
            self.var_backend_mode.set(self.prefs.get("backend_mode", "NumPy accelerated"))
        if hasattr(self, "var_tolerance"):
            self.var_tolerance.set(self.prefs.get("geometry_tolerance", "1e-10"))
        self.update_dashboard()
        self.preview_input_file()

    def _on_close(self) -> None:
        self._save_prefs()
        self.destroy()

    # -------------------------------------------------------------- dashboard
    def _build_dashboard_page(self) -> None:
        p = self.pages["Dashboard"]
        ttk.Label(p, text="Dashboard", style="Title.TLabel").pack(anchor="w", pady=(4, 12))
        grid = ttk.Frame(p)
        grid.pack(fill="x")
        self.dashboard_cards: Dict[str, ttk.Label] = {}
        for i, title in enumerate(["Input", "Element type", "Cohesive elements", "Warnings"]):
            c = self.card(grid, title)
            c.grid(row=0, column=i, sticky="nsew", padx=5, pady=5)
            val = ttk.Label(c, text="—", style="CardTitle.TLabel", font=("Segoe UI", 16, "bold"))
            val.pack(anchor="w", pady=6)
            self.dashboard_cards[title] = val
            grid.grid_columnconfigure(i, weight=1)

        quick = self.card(p, "Quick actions", "Load an Abaqus INP, edit cohesive parameters, then generate the cohesive model.")
        quick.pack(fill="x", pady=10)
        ttk.Button(quick, text="Choose INP", command=self.browse_input).pack(side="left", padx=4)
        ttk.Button(quick, text="Generate cohesive INP", style="Accent.TButton", command=self.generate_in_thread).pack(side="left", padx=4)
        ttk.Button(quick, text="Open output folder", command=lambda: self.open_folder(Path(self.var_output.get()).parent)).pack(side="left", padx=4)
        ttk.Button(quick, text="Visualization", command=lambda: self.show_page("Visualization")).pack(side="left", padx=4)

        self.dashboard_text = TextEditor(p, height=20)
        self.dashboard_text.pack(fill="both", expand=True, pady=8)

    def update_dashboard(self) -> None:
        s = self.summary or {}
        self.dashboard_cards.get("Input", ttk.Label()).configure(text=Path(self.var_input.get()).name or "—")
        self.dashboard_cards.get("Element type", ttk.Label()).configure(text=str(s.get("element_type", "—")))
        self.dashboard_cards.get("Cohesive elements", ttk.Label()).configure(text=str(s.get("cohesive_elements", "—")))
        warns = s.get("warnings", []) if isinstance(s, dict) else []
        self.dashboard_cards.get("Warnings", ttk.Label()).configure(text=str(len(warns)) if warns else "0")
        lines = [
            "CohesiveX Studio project summary",
            "=================================",
            f"Base directory : {self.var_base.get() if hasattr(self, 'var_base') else ''}",
            f"Input INP      : {self.var_input.get() if hasattr(self, 'var_input') else ''}",
            f"Output INP     : {self.var_output.get() if hasattr(self, 'var_output') else ''}",
            f"UEL file       : {self.var_uel.get() if hasattr(self, 'var_uel') else ''}",
            "",
            "Last generation summary:",
        ]
        if s:
            for key in ["element_type", "dimension", "original_nodes", "duplicated_nodes", "total_nodes", "solid_elements", "cohesive_elements", "grains", "nonmanifold_faces"]:
                lines.append(f"  {key:20s}: {s.get(key)}")
            if warns:
                lines.append("\nWarnings:")
                for w in warns:
                    lines.append(f"  - {w}")
        else:
            lines.append("  No cohesive model has been generated in this session yet.")
        if hasattr(self, "dashboard_text"):
            self.dashboard_text.set("\n".join(lines))
        if hasattr(self, "status_label"):
            self.status_label.configure(text="Status\nReady" if not s else "Status\nGenerated")

    # --------------------------------------------------------------- project
    def _build_project_page(self) -> None:
        p = self.pages["Project"]
        ttk.Label(p, text="Project", style="Title.TLabel").pack(anchor="w", pady=(4, 12))
        c = self.card(p, "Project paths", "Select input/output files. Running the script without arguments now opens this GUI.")
        c.pack(fill="x")
        self.var_base = tk.StringVar()
        self.var_input = tk.StringVar()
        self.var_output = tk.StringVar()
        self.var_uel = tk.StringVar()
        self._path_row(c, "Base directory", self.var_base, self.browse_base)
        self._path_row(c, "Input Abaqus INP", self.var_input, self.browse_input)
        self._path_row(c, "Output cohesive INP", self.var_output, self.browse_output)
        self._path_row(c, "UEL Fortran file", self.var_uel, self.browse_uel)
        actions = ttk.Frame(c, style="Panel.TFrame")
        actions.pack(fill="x", pady=(8, 0))
        ttk.Button(actions, text="Preview input", command=self.preview_input_file).pack(side="left", padx=4)
        ttk.Button(actions, text="Save config JSON", command=self.save_config_json).pack(side="left", padx=4)
        ttk.Button(actions, text="Load config JSON", command=self.load_config_json).pack(side="left", padx=4)
        ttk.Button(actions, text="Open base folder", command=lambda: self.open_folder(Path(self.var_base.get()))).pack(side="left", padx=4)
        self.input_preview = TextEditor(p, height=28)
        self.input_preview.pack(fill="both", expand=True, pady=10)

    def browse_base(self) -> None:
        d = filedialog.askdirectory(initialdir=self.var_base.get() or str(Path.cwd()))
        if d:
            self.var_base.set(d)

    def browse_input(self) -> None:
        f = filedialog.askopenfilename(
            initialdir=self.var_base.get() or str(Path.cwd()),
            filetypes=[("Abaqus input", "*.inp *.INP"), ("All files", "*.*")],
        )
        if f:
            p = Path(f)
            self.var_input.set(str(p))
            self.var_base.set(str(p.parent))
            if not self.var_output.get().strip() or Path(self.var_output.get()).name == "Job_coh.inp":
                self.var_output.set(str(p.with_name(p.stem + "_coh.inp")))
            self.preview_input_file()
            self.update_dashboard()

    def browse_output(self) -> None:
        f = filedialog.asksaveasfilename(
            initialdir=self.var_base.get() or str(Path.cwd()),
            initialfile=Path(self.var_output.get() or "Job_coh.inp").name,
            defaultextension=".inp",
            filetypes=[("Abaqus input", "*.inp *.INP"), ("All files", "*.*")],
        )
        if f:
            self.var_output.set(f)

    def browse_uel(self) -> None:
        f = filedialog.askopenfilename(
            initialdir=self.var_base.get() or str(Path.cwd()),
            filetypes=[("Fortran source", "*.f *.for *.f90 *.F"), ("All files", "*.*")],
        )
        if f:
            self.var_uel.set(f)

    def preview_input_file(self) -> None:
        if not hasattr(self, "input_preview"):
            return
        p = Path(self.var_input.get())
        if not p.is_file():
            self.input_preview.set(f"Input file not found:\n{p}")
            return
        try:
            text = read_text_auto(p)
            self.input_preview.set(text[:30000] + ("\n\n... <truncated>" if len(text) > 30000 else ""))
        except OSError as e:
            self.input_preview.set(f"Could not read input file:\n{e}")

    # ---------------------------------------------------------- cohesive law
    def _build_cohesive_law_page(self) -> None:
        p = self.pages["Cohesive Law"]
        ttk.Label(p, text="Cohesive model parameters", style="Title.TLabel").pack(anchor="w", pady=(4, 10))

        # A notebook layout avoids the clipped right-side panel that occurs on
        # smaller screens with a three-column PanedWindow.  It also separates
        # input parsing, UEL settings and cohesive-law parameters more clearly.
        nb = ttk.Notebook(p)
        nb.pack(fill="both", expand=True)

        tab_setup = ttk.Frame(nb, padding=8)
        tab_props = ttk.Frame(nb, padding=8)
        tab_mapping = ttk.Frame(nb, padding=8)
        nb.add(tab_setup, text="Input and UEL")
        nb.add(tab_props, text="Cohesive properties")
        nb.add(tab_mapping, text="Interface scope")

        setup_paned = ttk.PanedWindow(tab_setup, orient="horizontal")
        setup_paned.pack(fill="both", expand=True)
        left = ttk.Frame(setup_paned, padding=(0, 0, 6, 0))
        right = ttk.Frame(setup_paned, padding=(6, 0, 0, 0))
        setup_paned.add(left, weight=1)
        setup_paned.add(right, weight=1)

        c0 = self.card(
            left,
            "Analysis template",
            "Choose a template to fill recommended UEL settings.  The template does not limit the actual INP parser; it is a safety preset for common Abaqus workflows.",
        )
        c0.pack(fill="x", pady=(0, 8))
        self.var_model_template = tk.StringVar(value="2D plane strain - linear cohesive edge")
        self._combo_row(
            c0,
            "Template",
            self.var_model_template,
            [
                "2D plane strain - linear cohesive edge",
                "2D plane stress - linear cohesive edge",
                "2D quadratic cohesive edge",
                "3D tetrahedral cohesive face",
                "3D hexahedral cohesive face",
                "Custom / keep current settings",
            ],
        )
        ttk.Button(c0, text="Apply template", style="Accent.TButton", command=self.apply_model_template).pack(anchor="w", pady=(8, 0))
        ttk.Label(
            c0,
            text="For the current OXFORD-UMAT.f / UEL.f, INTMTD is hard-coded in the Fortran source.  Keep the GUI integration method consistent with that source, normally INTMTD = 1 unless the Fortran file is edited and recompiled.",
            style="Muted.TLabel",
            wraplength=520,
        ).pack(anchor="w", pady=(8, 0))

        c1 = self.card(
            left,
            "Input parsing and mesh options",
            "These fields control how domain/grain sets and geometric checks are interpreted. They do not rewrite material cards or analysis steps.",
        )
        c1.pack(fill="x", pady=(0, 8))
        self.var_grain_prefix = tk.StringVar(value="GRAIN-")
        self.var_tolerance = tk.StringVar(value="1e-10")
        self.var_supplement_nsets = tk.BooleanVar(value=True)
        self.var_write_cae_preview = tk.BooleanVar(value=True)
        self._entry_row(c1, "Grain set prefix", self.var_grain_prefix)
        self._entry_row(c1, "Geometry tolerance", self.var_tolerance)
        ttk.Checkbutton(c1, text="Supplement existing Nsets with duplicated nodes", variable=self.var_supplement_nsets).pack(anchor="w", pady=(6, 0))
        ttk.Checkbutton(c1, text="Generate CAE preview INP without UEL blocks", variable=self.var_write_cae_preview).pack(anchor="w", pady=(2, 0))
        ttk.Label(
            c1,
            text="Cohesive-only mode preserves original *Material, *Solid Section, *Step, *Static, *Boundary, *Cload, *Amplitude, *Controls and *Output blocks.",
            style="Muted.TLabel",
            wraplength=520,
        ).pack(anchor="w", pady=(8, 0))

        c2 = self.card(right, "User-element block", "Abaqus UEL keyword settings inserted into the solver INP.")
        c2.pack(fill="x", pady=(0, 8))
        self.var_uel_type = tk.StringVar(value="U1")
        self.var_uel_elset = tk.StringVar(value="GB_COH")
        self.var_intra_elset = tk.StringVar(value="INTRA_COH")
        self.var_intmtd = tk.StringVar(value="1")
        self.var_nsvars_per_ip = tk.StringVar(value="1")
        self._entry_row(c2, "UEL type", self.var_uel_type)
        self._entry_row(c2, "GB cohesive elset", self.var_uel_elset)
        self._entry_row(c2, "Intra cohesive elset", self.var_intra_elset)
        self._combo_row(c2, "Integration method", self.var_intmtd, ["1", "2", "3"])
        self._entry_row(c2, "SVARS / IP", self.var_nsvars_per_ip)
        ttk.Label(
            c2,
            text="Current OXFORD cohesive UEL stores one state variable per integration point, i.e. damage only. Therefore SVARS/IP = 1 and variables = number of integration points. If the UEL is extended to output more SDVs, this value and the *User element variable count must be updated together.",
            style="Muted.TLabel",
            wraplength=520,
        ).pack(anchor="w", pady=(8, 0))

        c_warn = self.card(right, "2D plane-strain submission checklist", "Use this when a CPE3/CPE4/CPE4R model stops immediately after job submission.")
        c_warn.pack(fill="x", pady=(0, 8))
        ttk.Label(
            c_warn,
            text=(
                "Recommended for the supplied OXFORD cohesive UEL: 2D linear cohesive element = 4 UEL nodes, coordinates = 2, active DOFs = 1,2, integration method = 1, variables = 2, HEIGHT > 0.  Do not import the solver INP into CAE; use the CAE preview file for inspection."
            ),
            style="Muted.TLabel",
            wraplength=520,
        ).pack(anchor="w")

        # ------------------------------ properties tab
        prop_top = self.card(tab_props, "Cohesive property preset", "Select a preset, then edit the nine values below if needed.")
        prop_top.pack(fill="x", pady=(0, 8))
        self.var_property_preset = tk.StringVar(value="Mixed mode - default")
        self._combo_row(
            prop_top,
            "Preset",
            self.var_property_preset,
            [
                "Mixed mode - default",
                "Opening mode only",
                "Shear mode only",
                "High-strength interface",
                "Weak interface",
                "Custom / keep current values",
            ],
        )
        ttk.Button(prop_top, text="Apply property preset", style="Accent.TButton", command=self.apply_property_preset).pack(anchor="w", pady=(8, 0))

        c4 = self.card(tab_props, "Bilinear cohesive UEL properties", "Nine properties written to *Uel Property for the cohesive element set.")
        c4.pack(fill="both", expand=True)
        self.prop_vars: Dict[str, tk.StringVar] = {}
        unit_map = {
            "mode": "-", "KI": "N/mm", "KII": "N/mm", "SI": "MPa", "SII": "MPa",
            "GCI": "energy", "GCII": "energy", "ETA": "-", "HEIGHT": "mm",
        }
        for (name, desc), default in zip(self.PROP_LABELS, self.PROP_DEFAULTS):
            row = ttk.Frame(c4, style="Panel.TFrame")
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=name, style="Muted.TLabel", width=10).pack(side="left")
            var = tk.StringVar(value=default)
            self.prop_vars[name] = var
            ttk.Entry(row, textvariable=var, width=16).pack(side="left", padx=5)
            ttk.Label(row, text=unit_map.get(name, ""), style="Muted.TLabel", width=10).pack(side="left")
            ttk.Label(row, text=desc, style="Muted.TLabel", wraplength=760).pack(side="left", fill="x", expand=True)
        actions = ttk.Frame(c4, style="Panel.TFrame")
        actions.pack(fill="x", pady=(10, 0))
        ttk.Button(actions, text="Opening mode preset", command=lambda: self.apply_prop_preset("opening")).pack(side="left", padx=3)
        ttk.Button(actions, text="Mixed mode preset", command=lambda: self.apply_prop_preset("mixed")).pack(side="left", padx=3)
        ttk.Button(actions, text="Plot cohesive law", style="Accent.TButton", command=lambda: self.show_page("Visualization") or self.draw_visual_law()).pack(side="left", padx=3)
        ttk.Button(actions, text="Reset defaults", command=self.reset_props).pack(side="left", padx=3)

        # ------------------------------ mapping tab
        c_scope = self.card(tab_mapping, "Interface insertion scope", "Select where cohesive elements are inserted and which property family is used.")
        c_scope.pack(fill="x", pady=(0, 8))
        self.var_interface_scope = tk.StringVar(value="Grain-boundary only")
        self._combo_row(
            c_scope,
            "Insert location",
            self.var_interface_scope,
            ["Grain-boundary only", "Intragranular only", "Grain-boundary + intragranular"],
        )
        self.var_intragranular_fraction = tk.StringVar(value="1.0")
        self.var_random_seed = tk.StringVar(value="")
        self._entry_row(c_scope, "Intra fraction", self.var_intragranular_fraction)
        self._entry_row(c_scope, "Random seed", self.var_random_seed)
        ttk.Label(
            c_scope,
            text=(
                "Selection rule: for two neighbouring solid elements e1 and e2, the interface is treated as grain-boundary when domain(e1) != domain(e2), "
                "and as intragranular/intradomain when domain(e1) == domain(e2).  GB interfaces use the GB property vector and elset; "
                "intragranular interfaces use the intra property vector and elset.  The intragranular fraction f selects a reproducible subset of intra interfaces: "
                "select interface i when r(seed, i) <= f.  If no grain/domain sets are found, the parsed solid block is treated as one domain, "
                "so the intragranular option can also be used for ordinary elastic-plastic models."
            ),
            style="Muted.TLabel",
            wraplength=900,
        ).pack(anchor="w", pady=(8, 0))

        c_intra = self.card(tab_mapping, "Intragranular cohesive properties", "Independent property vector for cohesive elements inserted inside a grain/domain.")
        c_intra.pack(fill="x", pady=(0, 8))
        self.intra_prop_vars: Dict[str, tk.StringVar] = {}
        intra_defaults = ["3", "1e7", "1e7", "180", "180", "10", "10", "2", "1"]
        for (name, desc), default in zip(self.PROP_LABELS, intra_defaults):
            row = ttk.Frame(c_intra, style="Panel.TFrame")
            row.pack(fill="x", pady=1)
            ttk.Label(row, text=name, style="Muted.TLabel", width=10).pack(side="left")
            var = tk.StringVar(value=default)
            self.intra_prop_vars[name] = var
            ttk.Entry(row, textvariable=var, width=16).pack(side="left", padx=5)
            ttk.Label(row, text=desc, style="Muted.TLabel", wraplength=650).pack(side="left", fill="x", expand=True)
        intra_actions = ttk.Frame(c_intra, style="Panel.TFrame")
        intra_actions.pack(fill="x", pady=(8, 0))
        ttk.Button(intra_actions, text="Copy GB properties to intra", command=self.copy_gb_to_intra).pack(side="left", padx=3)
        ttk.Button(intra_actions, text="Use stronger intra preset", command=self.reset_intra_props).pack(side="left", padx=3)


        c_notes = self.card(tab_mapping, "Implementation note for CPE3/CPE4", "Matched linear 2D interfaces are conceptually compatible with 4-node cohesive UELs, but the current release focuses on a single supported solid element block per generation run.")
        c_notes.pack(fill="x", pady=8)
        ttk.Label(
            c_notes,
            text="For publishable reliability, use a homogeneous CPE3, CPE4, CPS3 or CPS4 block for the present release. Mixed CPE3/CPE4 blocks and non-matching interfaces should be treated as a planned extension, not as a claimed feature.",
            style="Muted.TLabel",
            wraplength=900,
        ).pack(anchor="w")

    def apply_model_template(self) -> None:
        template = self.var_model_template.get() if hasattr(self, "var_model_template") else ""
        # The supplied OXFORD cohesive UEL uses INTMTD=1 in the Fortran source.
        # Keep GUI-generated variable counts consistent unless the user has edited and recompiled the UEL.
        if "2D" in template:
            self.var_intmtd.set("1")
            self.var_nsvars_per_ip.set("1")
            self.var_uel_type.set(self.var_uel_type.get() or "U1")
            self.var_uel_elset.set(self.var_uel_elset.get() or "UEL")
            if "HEIGHT" in self.prop_vars:
                try:
                    if float(self.prop_vars["HEIGHT"].get().replace("D", "E")) <= 0:
                        self.prop_vars["HEIGHT"].set("1")
                except Exception:
                    self.prop_vars["HEIGHT"].set("1")
        elif "3D" in template:
            self.var_intmtd.set("1")
            self.var_nsvars_per_ip.set("1")
            if "HEIGHT" in self.prop_vars:
                self.prop_vars["HEIGHT"].set("1")
        messagebox.showinfo("Template applied", f"Recommended settings were applied for:\n{template}")

    def apply_property_preset(self) -> None:
        preset = self.var_property_preset.get() if hasattr(self, "var_property_preset") else "Mixed mode - default"
        if preset == "Opening mode only":
            vals = ["1", "1e7", "0", "100", "0", "5", "0", "0", "1"]
        elif preset == "Shear mode only":
            vals = ["2", "0", "1e7", "0", "100", "0", "5", "0", "1"]
        elif preset == "High-strength interface":
            vals = ["3", "1e7", "1e7", "200", "200", "10", "10", "2", "1"]
        elif preset == "Weak interface":
            vals = ["3", "1e7", "1e7", "50", "50", "2", "2", "2", "1"]
        elif preset == "Custom / keep current values":
            return
        else:
            vals = self.PROP_DEFAULTS
        for (name, _), val in zip(self.PROP_LABELS, vals):
            self.prop_vars[name].set(val)
    def apply_prop_preset(self, kind: str) -> None:
        if kind == "opening":
            vals = ["1", "1e7", "0", "100", "0", "5", "0", "0", "1"]
        else:
            vals = self.PROP_DEFAULTS
        for (name, _), val in zip(self.PROP_LABELS, vals):
            self.prop_vars[name].set(val)

    def reset_props(self) -> None:
        for (name, _), val in zip(self.PROP_LABELS, self.PROP_DEFAULTS):
            self.prop_vars[name].set(val)

    def _float_from_var(self, var: tk.StringVar, name: str) -> float:
        try:
            return float(var.get().replace("D", "E").replace("d", "e"))
        except ValueError as e:
            raise ValueError(f"Invalid numeric value for {name}: {var.get()}") from e

    def current_props(self) -> Tuple[float, ...]:
        return tuple(self._float_from_var(self.prop_vars[name], name) for name, _ in self.PROP_LABELS)

    def current_intra_props(self) -> Tuple[float, ...]:
        if not hasattr(self, "intra_prop_vars"):
            return self.current_props()
        return tuple(self._float_from_var(self.intra_prop_vars[name], f"intra {name}") for name, _ in self.PROP_LABELS)

    def copy_gb_to_intra(self) -> None:
        if not hasattr(self, "intra_prop_vars"):
            return
        for name, _ in self.PROP_LABELS:
            self.intra_prop_vars[name].set(self.prop_vars[name].get())

    def reset_intra_props(self) -> None:
        if not hasattr(self, "intra_prop_vars"):
            return
        vals = ["3", "1e7", "1e7", "180", "180", "10", "10", "2", "1"]
        for (name, _), val in zip(self.PROP_LABELS, vals):
            self.intra_prop_vars[name].set(val)

    def _interface_scope_code(self) -> str:
        value = getattr(self, "var_interface_scope", tk.StringVar(value="Grain-boundary only")).get()
        if value == "Intragranular only":
            return "intragranular"
        if value == "Grain-boundary + intragranular":
            return "both"
        return "grain_boundary"

    def _selected_fast_mode(self) -> bool:
        value = getattr(self, "var_backend_mode", tk.StringVar(value="NumPy accelerated")).get()
        return not value.lower().startswith("pure")

    # --------------------------------------------------------------- generate
    def _build_generate_page(self) -> None:
        p = self.pages["Generate"]
        ttk.Label(p, text="Generate cohesive model", style="Title.TLabel").pack(anchor="w", pady=(4, 12))
        bar = self.card(p, "Run generator", "The GUI calls the cohesive-only INP modifier in a background thread.")
        bar.pack(fill="x")
        top_actions = ttk.Frame(bar, style="Panel.TFrame")
        top_actions.pack(fill="x")
        ttk.Button(top_actions, text="Generate cohesive INP", style="Accent.TButton", command=self.generate_in_thread).pack(side="left", padx=4)
        ttk.Button(top_actions, text="Compare backends", command=self.compare_backends_in_thread).pack(side="left", padx=4)
        ttk.Button(top_actions, text="Open output folder", command=lambda: self.open_folder(Path(self.var_output.get()).parent)).pack(side="left", padx=4)
        ttk.Button(top_actions, text="Open CAE preview", command=self.open_cae_preview_file).pack(side="left", padx=4)
        ttk.Button(top_actions, text="Go to reports", command=lambda: self.show_page("Reports")).pack(side="left", padx=4)
        self.generate_status = tk.StringVar(value="Idle")
        ttk.Label(top_actions, textvariable=self.generate_status, style="Muted.TLabel").pack(side="left", padx=14)

        backend_row = ttk.Frame(bar, style="Panel.TFrame")
        backend_row.pack(fill="x", pady=(10, 0))
        self.var_backend_mode = tk.StringVar(value="NumPy accelerated")
        ttk.Label(backend_row, text="Interface-detection backend", style="Muted.TLabel", width=28).pack(side="left")
        ttk.Combobox(
            backend_row,
            textvariable=self.var_backend_mode,
            values=["NumPy accelerated", "Pure Python reference"],
            state="readonly",
            width=28,
        ).pack(side="left", padx=(4, 12))
        ttk.Label(
            backend_row,
            text="Default: NumPy accelerated face hashing. The pure-Python backend keeps the same topology rules and is useful for fallback runs and speedup benchmarks.",
            style="Muted.TLabel",
            wraplength=760,
        ).pack(side="left", fill="x", expand=True)
        self.run_log = TextEditor(p, height=30)
        self.run_log.pack(fill="both", expand=True, pady=10)

    def generate_in_thread(self) -> None:
        if self.worker_thread is not None and self.worker_thread.is_alive():
            messagebox.showinfo("Busy", "Generation is already running.")
            return
        try:
            input_path = Path(self.var_input.get())
            output_path = Path(self.var_output.get())
            if not input_path.is_file():
                messagebox.showerror("Input missing", f"Input INP file not found:\n{input_path}")
                return
            props = self.current_props()
            intmtd = int(float(self.var_intmtd.get()))
            tolerance = self._float_from_var(self.var_tolerance, "geometry tolerance")
            uel_type = self.var_uel_type.get().strip() or "U1"
            uel_elset = self.var_uel_elset.get().strip() or "GB_COH"
            intra_elset = getattr(self, "var_intra_elset", tk.StringVar(value="INTRA_COH")).get().strip() or "INTRA_COH"
            interface_scope = self._interface_scope_code()
            intra_props = self.current_intra_props()
            try:
                intragranular_fraction = float(getattr(self, "var_intragranular_fraction", tk.StringVar(value="1.0")).get() or "1.0")
            except ValueError:
                raise ValueError("Intra fraction must be a number in [0, 1].")
            seed_text = getattr(self, "var_random_seed", tk.StringVar(value="")).get().strip()
            random_seed = int(seed_text) if seed_text else None
            if intragranular_fraction < 1.0 and random_seed is None:
                random_seed = 0
                self.log_queue.put("[settings] Intra fraction < 1.0 and no random seed was provided; using deterministic seed 0.\n")
            grain_prefix = self.var_grain_prefix.get()
            fast_mode = self._selected_fast_mode()
            supplement_nsets = bool(self.var_supplement_nsets.get())
            write_cae_preview_flag = bool(self.var_write_cae_preview.get())
            nsvars_per_ip = int(float(getattr(self, "var_nsvars_per_ip", tk.StringVar(value="1")).get()))
            if nsvars_per_ip != 1:
                raise ValueError("The current generator is configured for the supplied OXFORD cohesive UEL with NSVPN=1. Set SVARS/IP to 1 unless the UEL and generator are updated together.")
            if props[8] <= 0.0:
                raise ValueError("HEIGHT must be positive for 2D cohesive elements because the UEL multiplies the 2D line Jacobian by HEIGHT.")
            if intmtd != 1:
                ok = messagebox.askyesno(
                    "Integration method check",
                    "The supplied OXFORD-UMAT.f / UEL.f has INTMTD hard-coded as 1.\n\n"
                    "Using another GUI integration method can make *User element variables inconsistent with the compiled UEL for higher-order elements.\n\n"
                    "Continue anyway?"
                )
                if not ok:
                    return
        except Exception as e:
            messagebox.showerror("Invalid settings", str(e))
            return

        self.run_log.set("")
        self.generate_status.set("Running...")
        self.status_label.configure(text="Status\nGenerating")
        self.log_queue.put(f"[{time.strftime('%H:%M:%S')}] Starting cohesive generation...\n")

        def worker() -> None:
            try:
                summary = generate_cohesive_inp(
                    input_path,
                    output_path,
                    grain_set_prefix=grain_prefix,
                    intmtd=intmtd,
                    props=props,
                    gb_props=props,
                    intra_props=intra_props,
                    interface_scope=interface_scope,
                    uel_type=uel_type,
                    uel_elset=uel_elset,
                    gb_elset=uel_elset,
                    intra_elset=intra_elset,
                    nsvars_per_ip=nsvars_per_ip,
                    intragranular_fraction=intragranular_fraction,
                    random_seed=random_seed,
                    fast_mode=fast_mode,
                    supplement_nsets=supplement_nsets,
                    write_cae_preview=write_cae_preview_flag,
                    tolerance=tolerance,
                    verbose=False,
                )
                self.summary = summary
                self.log_queue.put(f"[{time.strftime('%H:%M:%S')}] Done. Solver INP: {output_path}\n")
                self.log_queue.put(f"CAE preview INP: {summary.get('cae_preview_file')}\n")
                self.log_queue.put(f"Generated cohesive elements: {summary.get('cohesive_elements')}\n")
                self.log_queue.put(f"Backend used: {summary.get('backend')}\n")
                self.log_queue.put(f"Total preprocessing time: {summary.get('total_preprocessing_time_seconds', 0.0):.6f} s\n")
                self.log_queue.put("Step timings:\n")
                for _name, _value in summary.get("timings_seconds", {}).items():
                    self.log_queue.put(f"  - {_name:<28}: {_value:.6f} s\n")
                self.log_queue.put("__GEN_DONE__")
            except Exception as e:
                self.log_queue.put(f"[{time.strftime('%H:%M:%S')}] ERROR: {e}\n")
                self.log_queue.put("__GEN_ERROR__")

        self.worker_thread = threading.Thread(target=worker, daemon=True)
        self.worker_thread.start()

    def compare_backends_in_thread(self) -> None:
        if self.worker_thread is not None and self.worker_thread.is_alive():
            messagebox.showinfo("Busy", "Generation is already running.")
            return
        try:
            input_path = Path(self.var_input.get())
            output_path = Path(self.var_output.get())
            if not input_path.is_file():
                messagebox.showerror("Input missing", f"Input INP file not found:\n{input_path}")
                return
            props = self.current_props()
            intra_props = self.current_intra_props()
            intmtd = int(float(self.var_intmtd.get()))
            tolerance = self._float_from_var(self.var_tolerance, "geometry tolerance")
            uel_type = self.var_uel_type.get().strip() or "U1"
            uel_elset = self.var_uel_elset.get().strip() or "GB_COH"
            intra_elset = getattr(self, "var_intra_elset", tk.StringVar(value="INTRA_COH")).get().strip() or "INTRA_COH"
            interface_scope = self._interface_scope_code()
            intragranular_fraction = float(getattr(self, "var_intragranular_fraction", tk.StringVar(value="1.0")).get() or "1.0")
            seed_text = getattr(self, "var_random_seed", tk.StringVar(value="")).get().strip()
            random_seed = int(seed_text) if seed_text else None
            nsvars_per_ip = int(float(getattr(self, "var_nsvars_per_ip", tk.StringVar(value="1")).get()))
        except Exception as e:
            messagebox.showerror("Invalid settings", str(e))
            return

        prefix = output_path.with_suffix("")
        self.run_log.set("")
        self.generate_status.set("Benchmarking...")
        self.status_label.configure(text="Status\nBenchmark")
        self.log_queue.put(f"[{time.strftime('%H:%M:%S')}] Comparing NumPy and pure-Python backends...\n")

        def worker() -> None:
            try:
                comparison = compare_backends(
                    input_path,
                    prefix,
                    grain_set_prefix=self.var_grain_prefix.get(),
                    intmtd=intmtd,
                    props=props,
                    gb_props=props,
                    intra_props=intra_props,
                    interface_scope=interface_scope,
                    uel_type=uel_type,
                    uel_elset=uel_elset,
                    gb_elset=uel_elset,
                    intra_elset=intra_elset,
                    nsvars_per_ip=nsvars_per_ip,
                    intragranular_fraction=intragranular_fraction,
                    random_seed=random_seed,
                    supplement_nsets=bool(self.var_supplement_nsets.get()),
                    write_cae_preview=bool(self.var_write_cae_preview.get()),
                    tolerance=tolerance,
                    verbose=False,
                )
                self.summary = comparison.get("numpy", {})
                self.log_queue.put("Backend comparison completed.\n")
                self.log_queue.put(f"Topology consistent: {comparison.get('topology_consistent')}\n")
                if comparison.get("speedup_total") is not None:
                    self.log_queue.put(f"Total preprocessing speedup: {comparison.get('speedup_total'):.3f} x\n")
                if comparison.get("speedup_face_detection") is not None:
                    self.log_queue.put(f"Face-detection speedup: {comparison.get('speedup_face_detection'):.3f} x\n")
                self.log_queue.put(f"Comparison report: {comparison.get('comparison_report')}\n")
                self.log_queue.put("__GEN_DONE__")
            except Exception as e:
                self.log_queue.put(f"[{time.strftime('%H:%M:%S')}] ERROR: {e}\n")
                self.log_queue.put("__GEN_ERROR__")

        self.worker_thread = threading.Thread(target=worker, daemon=True)
        self.worker_thread.start()

    def _drain_log_queue(self) -> None:
        try:
            while True:
                msg = self.log_queue.get_nowait()
                if msg == "__GEN_DONE__":
                    self.generate_status.set("Done")
                    self.status_label.configure(text="Status\nGenerated")
                    self.var_run_inp.set(self.var_output.get())
                    self.var_run_fortran.set(self.var_uel.get())
                    self.var_job.set(Path(self.var_output.get()).stem or self.var_job.get())
                    self.update_abaqus_preview()
                    self.update_dashboard()
                    self.refresh_reports()
                    self.refresh_visualization()
                elif msg == "__GEN_ERROR__":
                    self.generate_status.set("Failed")
                    self.status_label.configure(text="Status\nError")
                elif msg == "__ABAQUS_DONE__":
                    self.abaqus_status.set("Done")
                else:
                    if hasattr(self, "run_log"):
                        self.run_log.append(msg)
                    if hasattr(self, "abaqus_notes") and (msg.startswith("[abaqus]") or self.abaqus_status.get() == "Running..."):
                        self.abaqus_notes.append(msg)
        except queue.Empty:
            pass
        self.after(150, self._drain_log_queue)

    def open_cae_preview_file(self) -> None:
        preview = Path(self.var_output.get()).with_name(Path(self.var_output.get()).stem + "_cae_preview.inp")
        if not preview.is_file():
            messagebox.showinfo("CAE preview", f"CAE preview file not found:\n{preview}\nGenerate the cohesive INP first.")
            return
        try:
            if sys.platform.startswith("win"):
                os.startfile(str(preview))  # type: ignore[attr-defined]
            elif sys.platform == "darwin":
                subprocess.Popen(["open", str(preview)])
            else:
                subprocess.Popen(["xdg-open", str(preview)])
        except OSError as e:
            messagebox.showerror("Open failed", str(e))

    # ---------------------------------------------------------------- reports
    def _build_reports_page(self) -> None:
        p = self.pages["Reports"]
        ttk.Label(p, text="Reports", style="Title.TLabel").pack(anchor="w", pady=(4, 8))
        bar = ttk.Frame(p)
        bar.pack(fill="x", pady=(0, 6))
        ttk.Button(bar, text="Refresh", command=self.refresh_reports).pack(side="left", padx=2)
        ttk.Button(bar, text="Open output folder", command=lambda: self.open_folder(Path(self.var_output.get()).parent)).pack(side="left", padx=2)
        ttk.Button(bar, text="Go to visualization", command=lambda: self.show_page("Visualization")).pack(side="left", padx=2)

        paned = ttk.PanedWindow(p, orient="horizontal")
        paned.pack(fill="both", expand=True)
        left = ttk.Frame(paned, padding=4)
        paned.add(left, weight=1)
        ttk.Label(left, text="Mesh-check report", style="CardTitle.TLabel").pack(anchor="w")
        self.check_text = TextEditor(left, height=24)
        self.check_text.pack(fill="both", expand=True, pady=4)

        right = ttk.Frame(paned, padding=4)
        paned.add(right, weight=2)
        ttk.Label(right, text="Grain-boundary cohesive table", style="CardTitle.TLabel").pack(anchor="w")
        cols = ("coh", "e1", "g1", "e2", "g2", "measure", "nx", "ny", "nz")
        self.gb_tree = ttk.Treeview(right, columns=cols, show="headings", height=22)
        headings = ["Coh ID", "Elem 1", "Grain 1", "Elem 2", "Grain 2", "Measure", "nx", "ny", "nz"]
        widths = [80, 70, 70, 70, 70, 95, 75, 75, 75]
        for c, h, w in zip(cols, headings, widths):
            self.gb_tree.heading(c, text=h)
            self.gb_tree.column(c, width=w, anchor="e")
        y = ttk.Scrollbar(right, orient="vertical", command=self.gb_tree.yview)
        self.gb_tree.configure(yscrollcommand=y.set)
        self.gb_tree.pack(side="left", fill="both", expand=True, pady=4)
        y.pack(side="left", fill="y", pady=4)

    def _report_prefix(self) -> Path:
        return Path(self.var_output.get()).with_suffix("")

    def refresh_reports(self) -> None:
        if not hasattr(self, "check_text"):
            return
        prefix = self._report_prefix()
        check = prefix.with_name(prefix.name + "_mesh_check.txt")
        gb = prefix.with_name(prefix.name + "_grain_boundary_table.csv")
        if check.is_file():
            self.check_text.set(read_text_auto(check))
        else:
            self.check_text.set(f"Mesh-check report not found:\n{check}\nGenerate a cohesive INP first.")
        self.gb_tree.delete(*self.gb_tree.get_children())
        if gb.is_file():
            try:
                with gb.open("r", encoding="utf-8", newline="") as fh:
                    rdr = csv.DictReader(fh)
                    for i, row in enumerate(rdr):
                        if i >= 5000:
                            break
                        self.gb_tree.insert(
                            "",
                            "end",
                            values=(
                                row.get("cohesive_id", ""),
                                row.get("element_1", ""),
                                row.get("grain_1", ""),
                                row.get("element_2", ""),
                                row.get("grain_2", ""),
                                row.get("measure", ""),
                                row.get("normal_x", ""),
                                row.get("normal_y", ""),
                                row.get("normal_z", ""),
                            ),
                        )
            except OSError as e:
                self.check_text.append(f"\nCould not read GB table: {e}\n")

    # ---------------------------------------------------------- visualization
    def _build_visualization_page(self) -> None:
        p = self.pages["Visualization"]
        ttk.Label(p, text="Visualization dashboard", style="Title.TLabel").pack(anchor="w", pady=(4, 8))
        top = self.card(
            p,
            "Generated-model visualization",
            "Inspect cohesive-element statistics from the CSV/JSON reports without importing the UEL solver INP into Abaqus/CAE.",
        )
        top.pack(fill="x", pady=(0, 6))
        ttk.Button(top, text="Refresh data", command=self.refresh_visualization).pack(side="left", padx=4)
        ttk.Button(top, text="Summary bars", style="Accent.TButton", command=self.draw_visual_summary).pack(side="left", padx=4)
        ttk.Button(top, text="Grain-pair counts", command=self.draw_visual_grain_pairs).pack(side="left", padx=4)
        ttk.Button(top, text="Interface measure histogram", command=self.draw_visual_measure_hist).pack(side="left", padx=4)
        ttk.Button(top, text="Normal direction scatter", command=self.draw_visual_normals).pack(side="left", padx=4)
        ttk.Button(top, text="Cohesive law curves", command=self.draw_visual_law).pack(side="left", padx=4)
        ttk.Button(top, text="B-K mixed-mode envelope", command=self.draw_visual_bk).pack(side="left", padx=4)
        ttk.Button(top, text="Open output folder", command=lambda: self.open_folder(Path(self.var_output.get()).parent)).pack(side="left", padx=4)

        paned = ttk.PanedWindow(p, orient="horizontal")
        paned.pack(fill="both", expand=True)
        left = ttk.Frame(paned, padding=4)
        paned.add(left, weight=1)
        ttk.Label(left, text="Visualization data summary", style="CardTitle.TLabel").pack(anchor="w")
        self.visual_text = TextEditor(left, height=28, wrap="word")
        self.visual_text.pack(fill="both", expand=True, pady=4)

        right = ttk.Frame(paned, padding=4)
        paned.add(right, weight=3)
        ttk.Label(right, text="Figure preview", style="CardTitle.TLabel").pack(anchor="w")
        self.visual_canvas = PlotCanvas(right, figsize=(8.6, 6.2))
        self.visual_canvas.pack(fill="both", expand=True, pady=4)

    def _summary_json_path(self) -> Path:
        prefix = self._report_prefix()
        return prefix.with_name(prefix.name + "_summary.json")

    def _grain_boundary_csv_path(self) -> Path:
        prefix = self._report_prefix()
        return prefix.with_name(prefix.name + "_grain_boundary_table.csv")

    def _read_summary_report(self) -> Dict[str, Any]:
        path = self._summary_json_path()
        if not path.is_file():
            return dict(self.summary or {})
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return dict(self.summary or {})

    def _read_gb_rows(self) -> List[Dict[str, str]]:
        path = self._grain_boundary_csv_path()
        if not path.is_file():
            return []
        rows: List[Dict[str, str]] = []
        try:
            with path.open("r", encoding="utf-8", newline="") as fh:
                for row in csv.DictReader(fh):
                    rows.append(dict(row))
        except OSError:
            return []
        return rows

    @staticmethod
    def _to_float(value: Any, default: float = 0.0) -> float:
        try:
            return float(str(value).replace("D", "E").replace("d", "e"))
        except Exception:
            return default

    def refresh_visualization(self) -> None:
        if not hasattr(self, "visual_text"):
            return
        summary = self._read_summary_report()
        rows = self._read_gb_rows()
        gb_path = self._grain_boundary_csv_path()
        sm_path = self._summary_json_path()
        if not summary and not rows:
            self.visual_text.set(
                "No visualization data found.\n\n"
                f"Expected summary file:\n{sm_path}\n\n"
                f"Expected grain-boundary table:\n{gb_path}\n\n"
                "Generate the cohesive INP first, then return to this page."
            )
            if MATPLOTLIB_AVAILABLE and hasattr(self, "visual_canvas") and self.visual_canvas.figure is not None:
                fig = self.visual_canvas.figure
                fig.clear()
                ax = fig.add_subplot(1, 1, 1)
                ax.axis("off")
                ax.text(0.5, 0.55, "No generated report data", ha="center", va="center", fontsize=13, fontweight="bold")
                ax.text(0.5, 0.45, "Generate a cohesive model to enable visualization.", ha="center", va="center", fontsize=10)
                self.visual_canvas.draw()
            return

        pair_counts: Dict[str, int] = {}
        measures: List[float] = []
        for row in rows:
            g1 = str(row.get("grain_1", "?"))
            g2 = str(row.get("grain_2", "?"))
            try:
                a, b = sorted([int(float(g1)), int(float(g2))])
                key = f"{a}-{b}"
            except Exception:
                key = f"{g1}-{g2}"
            pair_counts[key] = pair_counts.get(key, 0) + 1
            if row.get("measure", "") not in ("", None):
                measures.append(self._to_float(row.get("measure"), 0.0))

        lines = [
            f"{PLATFORM_NAME} visualization summary",
            "=" * (len(PLATFORM_NAME) + 22),
            f"Summary JSON       : {sm_path}",
            f"GB table CSV       : {gb_path}",
            "",
            f"Element type       : {summary.get('element_type', '—')}",
            f"Dimension          : {summary.get('dimension', '—')}",
            f"Original nodes     : {summary.get('original_nodes', '—')}",
            f"Duplicated nodes   : {summary.get('duplicated_nodes', '—')}",
            f"Total nodes        : {summary.get('total_nodes', '—')}",
            f"Solid elements     : {summary.get('solid_elements', '—')}",
            f"Cohesive elements  : {summary.get('cohesive_elements', len(rows))}",
            f"Grains             : {summary.get('grains', '—')}",
            f"Grain-pair groups  : {len(pair_counts)}",
            f"Non-manifold faces : {summary.get('nonmanifold_faces', '—')}",
        ]
        if measures:
            lines.extend([
                "",
                "Interface measure statistics",
                f"  count : {len(measures)}",
                f"  min   : {min(measures):.6g}",
                f"  max   : {max(measures):.6g}",
                f"  mean  : {sum(measures) / max(1, len(measures)):.6g}",
            ])
        warnings = summary.get("warnings", []) if isinstance(summary, dict) else []
        lines.append("")
        lines.append("Warnings")
        if warnings:
            for w in warnings:
                lines.append(f"  - {w}")
        else:
            lines.append("  None")
        self.visual_text.set("\n".join(lines))
        self.draw_visual_summary()

    def _clear_visual_figure(self, title: str = "") -> Any:
        if not MATPLOTLIB_AVAILABLE or not hasattr(self, "visual_canvas") or self.visual_canvas.figure is None:
            return None
        fig = self.visual_canvas.figure
        fig.clear()
        fig.patch.set_facecolor("white")
        try:
            fig.set_layout_engine("tight")
        except Exception:
            pass
        if title:
            fig.suptitle(title, fontsize=12, fontweight="bold")
        return fig

    def draw_visual_summary(self) -> None:
        fig = self._clear_visual_figure("Cohesive model summary")
        if fig is None:
            return
        summary = self._read_summary_report()
        keys = ["original_nodes", "duplicated_nodes", "solid_elements", "cohesive_elements", "grains"]
        labels = ["Original nodes", "Duplicated nodes", "Solid elements", "Cohesive elements", "Grains"]
        values = [self._to_float(summary.get(k, 0), 0.0) for k in keys]
        ax = fig.add_subplot(1, 1, 1)
        ax.bar(labels, values)
        ax.set_ylabel("Count")
        ax.set_title("Generated model size and cohesive insertion result", pad=10)
        ax.tick_params(axis="x", labelrotation=20)
        ax.grid(axis="y", alpha=0.25)
        for i, v in enumerate(values):
            ax.text(i, v, f"{int(v) if abs(v-round(v))<1e-9 else v:.0f}", ha="center", va="bottom", fontsize=9)
        self.visual_canvas.draw()

    def draw_visual_grain_pairs(self) -> None:
        fig = self._clear_visual_figure("Grain-pair cohesive interface counts")
        if fig is None:
            return
        rows = self._read_gb_rows()
        pair_counts: Dict[str, int] = {}
        for row in rows:
            g1 = str(row.get("grain_1", "?"))
            g2 = str(row.get("grain_2", "?"))
            try:
                a, b = sorted([int(float(g1)), int(float(g2))])
                key = f"G{a}-G{b}"
            except Exception:
                key = f"{g1}-{g2}"
            pair_counts[key] = pair_counts.get(key, 0) + 1
        ax = fig.add_subplot(1, 1, 1)
        if not pair_counts:
            ax.axis("off")
            ax.text(0.5, 0.5, "No grain-boundary rows were found.", ha="center", va="center", fontsize=12)
            self.visual_canvas.draw()
            return
        items = sorted(pair_counts.items(), key=lambda kv: kv[1], reverse=True)[:20]
        labels = [k for k, _ in items][::-1]
        values = [v for _, v in items][::-1]
        ax.barh(labels, values)
        ax.set_xlabel("Cohesive element count")
        ax.set_title("Top grain-pair interface groups", pad=10)
        ax.grid(axis="x", alpha=0.25)
        self.visual_canvas.draw()

    def draw_visual_measure_hist(self) -> None:
        fig = self._clear_visual_figure("Interface measure distribution")
        if fig is None:
            return
        rows = self._read_gb_rows()
        measures = [self._to_float(row.get("measure"), 0.0) for row in rows if row.get("measure", "") not in ("", None)]
        ax = fig.add_subplot(1, 1, 1)
        if not measures:
            ax.axis("off")
            ax.text(0.5, 0.5, "No interface measure column was found.", ha="center", va="center", fontsize=12)
            self.visual_canvas.draw()
            return
        bins = min(40, max(8, int(math.sqrt(len(measures)))))
        ax.hist(measures, bins=bins)
        ax.set_xlabel("Interface measure: length × height in 2D, area in 3D")
        ax.set_ylabel("Frequency")
        ax.set_title("Distribution of cohesive interface measure", pad=10)
        ax.grid(axis="y", alpha=0.25)
        self.visual_canvas.draw()

    def draw_visual_normals(self) -> None:
        fig = self._clear_visual_figure("Cohesive interface normal directions")
        if fig is None:
            return
        rows = self._read_gb_rows()
        nx = [self._to_float(row.get("normal_x"), 0.0) for row in rows if row.get("normal_x", "") not in ("", None)]
        ny = [self._to_float(row.get("normal_y"), 0.0) for row in rows if row.get("normal_y", "") not in ("", None)]
        nz = [self._to_float(row.get("normal_z"), 0.0) for row in rows if row.get("normal_z", "") not in ("", None)]
        ax = fig.add_subplot(1, 1, 1)
        if not nx or len(nx) != len(ny):
            ax.axis("off")
            ax.text(0.5, 0.5, "No normal-vector columns were found.", ha="center", va="center", fontsize=12)
            self.visual_canvas.draw()
            return
        if nz and len(nz) == len(nx):
            sc = ax.scatter(nx, ny, c=nz, s=18, alpha=0.75)
            cb = fig.colorbar(sc, ax=ax, shrink=0.82)
            cb.set_label("normal_z")
        else:
            ax.scatter(nx, ny, s=18, alpha=0.75)
        ax.set_xlabel("normal_x")
        ax.set_ylabel("normal_y")
        ax.set_aspect("equal", adjustable="box")
        ax.set_title("Projected interface normal distribution", pad=10)
        ax.grid(alpha=0.25)
        self.visual_canvas.draw()


    def draw_visual_law(self) -> None:
        """Plot bilinear normal/shear traction-separation curves from the current property panel."""
        fig = self._clear_visual_figure("Bilinear cohesive traction-separation law")
        if fig is None:
            return
        try:
            mode, kn, ks, sn, ss, gic, giic, eta, height = self.current_props()
        except Exception as e:
            ax = fig.add_subplot(1, 1, 1)
            ax.axis("off")
            ax.text(0.5, 0.5, f"Invalid cohesive properties:\n{e}", ha="center", va="center", fontsize=12)
            self.visual_canvas.draw()
            return

        def bilinear(k: float, strength: float, gc: float) -> Tuple[List[float], List[float], float, float]:
            if k <= 0 or strength <= 0 or gc <= 0:
                return [0.0, 1.0], [0.0, 0.0], 0.0, 0.0
            d0 = strength / k
            df = max(2.0 * gc / strength, d0 * 1.001)
            x = [0.0, d0, df]
            y = [0.0, strength, 0.0]
            return x, y, d0, df

        ax = fig.add_subplot(1, 1, 1)
        xn, yn, dn0, dnf = bilinear(kn, sn, gic)
        xs, ys, ds0, dsf = bilinear(ks, ss, giic)
        ax.plot(xn, yn, marker="o", label=f"Mode I: δ0={dn0:.3g}, δf={dnf:.3g}")
        if mode in (2, 3):
            ax.plot(xs, ys, marker="s", label=f"Mode II: δ0={ds0:.3g}, δf={dsf:.3g}")
        ax.set_xlabel("Separation")
        ax.set_ylabel("Traction")
        ax.set_title("Penalty response, damage initiation and linear softening", pad=10)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc="best")
        note = (
            "δ0 = σmax/K,    δf = 2Gc/σmax\n"
            "Before damage: t = Kδ; after initiation: t = (1-D)Kδ"
        )
        ax.text(0.02, 0.98, note, transform=ax.transAxes, va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cbd5e1", alpha=0.95))
        self.visual_canvas.draw()

    def draw_visual_bk(self) -> None:
        """Plot the Benzeggagh-Kenane mixed-mode fracture-energy envelope."""
        fig = self._clear_visual_figure("Benzeggagh-Kenane mixed-mode fracture envelope")
        if fig is None:
            return
        try:
            _mode, _kn, _ks, _sn, _ss, gic, giic, eta, _height = self.current_props()
        except Exception as e:
            ax = fig.add_subplot(1, 1, 1)
            ax.axis("off")
            ax.text(0.5, 0.5, f"Invalid cohesive properties:\n{e}", ha="center", va="center", fontsize=12)
            self.visual_canvas.draw()
            return
        ax = fig.add_subplot(1, 1, 1)
        x = [i / 100.0 for i in range(101)]
        y = [gic + (giic - gic) * (r ** eta) for r in x]
        ax.plot(x, y, lw=2)
        ax.set_xlabel("Shear energy ratio, Gs / GT")
        ax.set_ylabel("Critical fracture energy, Gc")
        ax.set_title("Gc = GIc + (GIIc - GIc)(Gs/GT)^η", pad=10)
        ax.grid(alpha=0.25)
        ax.text(0.02, 0.98, f"GIc={gic:g}, GIIc={giic:g}, η={eta:g}", transform=ax.transAxes,
                va="top", fontsize=9, bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cbd5e1", alpha=0.95))
        self.visual_canvas.draw()

    # ------------------------------------------------------------- Abaqus run
    def _build_run_abaqus_page(self) -> None:
        p = self.pages["Run Abaqus"]
        ttk.Label(p, text="Run Abaqus", style="Title.TLabel").pack(anchor="w", pady=(4, 12))
        c = self.card(
            p,
            "Abaqus command",
            "Optional helper for submitting a selected solver INP with a selected Fortran UEL/VUMAT/UMAT source.",
        )
        c.pack(fill="x")
        self.var_abaqus_cmd = tk.StringVar(value="abaqus")
        self.var_job = tk.StringVar(value="Job_coh")
        self.var_run_inp = tk.StringVar()
        self.var_run_fortran = tk.StringVar()
        self.var_cpus = tk.StringVar(value="1")
        self._entry_row(c, "Abaqus command", self.var_abaqus_cmd)
        self._entry_row(c, "Job name", self.var_job)
        self._path_row(c, "Solver INP for run", self.var_run_inp, self.browse_run_inp)
        self._path_row(c, "Fortran source", self.var_run_fortran, self.browse_run_fortran)
        self._entry_row(c, "Parallel CPUs", self.var_cpus)
        quick = ttk.Frame(c, style="Panel.TFrame")
        quick.pack(fill="x", pady=(4, 0))
        ttk.Button(quick, text="Use generated INP/UEL", command=self.use_generated_run_files).pack(side="left", padx=4)
        ttk.Button(quick, text="Update command preview", command=self.update_abaqus_preview).pack(side="left", padx=4)
        actions = ttk.Frame(c, style="Panel.TFrame")
        actions.pack(fill="x", pady=(8, 0))
        ttk.Button(actions, text="Run Abaqus", style="Accent.TButton", command=self.run_abaqus).pack(side="left", padx=4)
        ttk.Button(actions, text="Generate run.bat", command=self.save_run_bat).pack(side="left", padx=4)
        ttk.Button(actions, text="Stop", command=self.stop_abaqus).pack(side="left", padx=4)
        ttk.Button(actions, text="Open job folder", command=lambda: self.open_folder(Path(self.var_run_inp.get() or self.var_output.get()).parent)).pack(side="left", padx=4)
        self.abaqus_status = tk.StringVar(value="Idle")
        ttk.Label(actions, textvariable=self.abaqus_status, style="Muted.TLabel").pack(side="left", padx=14)
        self.abaqus_notes = TextEditor(p, height=28)
        self.abaqus_notes.pack(fill="both", expand=True, pady=10)
        self.update_abaqus_preview()

    def browse_run_inp(self) -> None:
        f = filedialog.askopenfilename(
            initialdir=self.var_base.get() or str(Path.cwd()),
            filetypes=[("Abaqus input", "*.inp *.INP"), ("All files", "*.*")],
        )
        if f:
            self.var_run_inp.set(f)
            p = Path(f)
            if not self.var_job.get().strip() or self.var_job.get().strip() == "Job_coh":
                self.var_job.set(p.stem)
            self.update_abaqus_preview()

    def browse_run_fortran(self) -> None:
        f = filedialog.askopenfilename(
            initialdir=self.var_base.get() or str(Path.cwd()),
            filetypes=[("Fortran source", "*.f *.for *.f90 *.F *.FOR *.F90"), ("All files", "*.*")],
        )
        if f:
            self.var_run_fortran.set(f)
            self.update_abaqus_preview()

    def use_generated_run_files(self) -> None:
        self.var_run_inp.set(self.var_output.get())
        self.var_run_fortran.set(self.var_uel.get())
        out = Path(self.var_output.get())
        if out.name:
            self.var_job.set(out.stem)
        self.update_abaqus_preview()

    def _safe_cpu_count(self) -> int:
        try:
            n = int(float(self.var_cpus.get()))
        except (TypeError, ValueError):
            n = 1
        return max(1, n)

    def _quote_cmd_arg(self, value: str) -> str:
        # Windows command shell quoting.  POSIX shlex.quote also works for display on Linux/macOS.
        if sys.platform.startswith("win"):
            v = str(value)
            if not v:
                return '""'
            if any(ch.isspace() for ch in v) or any(ch in v for ch in '()&^%!,;='):
                return '"' + v.replace('\"', '\\"') + '"'
            return v
        return shlex.quote(str(value))

    def _abaqus_command_payload(self) -> Tuple[str, List[str], Path, Path]:
        inp = Path(self.var_run_inp.get() or self.var_output.get()).expanduser()
        user = Path(self.var_run_fortran.get() or self.var_uel.get()).expanduser()
        job = (self.var_job.get() or inp.stem or "Job_coh").strip()
        cpus = self._safe_cpu_count()
        abaqus_cmd = (self.var_abaqus_cmd.get() or "abaqus").strip()
        # Submit from the INP folder and pass the INP as a local file name.
        # This is usually more robust for Abaqus on Windows than absolute input paths.
        inp_arg = inp.name if inp.name else str(inp)
        user_arg = str(user.resolve()) if user.exists() else str(user)
        abaqus_args = [
            abaqus_cmd,
            f"job={job}",
            f"input={inp_arg}",
            f"user={user_arg}",
            f"cpus={cpus}",
            "interactive",
        ]
        display_cmd = " ".join([
            self._quote_cmd_arg(abaqus_cmd),
            f"job={self._quote_cmd_arg(job)}",
            f"input={self._quote_cmd_arg(inp_arg)}",
            f"user={self._quote_cmd_arg(user_arg)}",
            f"cpus={cpus}",
            "interactive",
        ])
        if sys.platform.startswith("win"):
            # Use cmd.exe only as the .bat resolver, while still keeping
            # shell=False.  Passing one quoted command string after /c is more
            # reliable for abaqus.bat and full paths containing spaces than
            # asking CreateProcess to execute a batch file directly.
            run_args = ["cmd.exe", "/d", "/c", display_cmd]
        else:
            run_args = abaqus_args
        return display_cmd, run_args, inp, user

    def _abaqus_command_string(self) -> Tuple[str, Path, Path]:
        display_cmd, _run_args, inp, user = self._abaqus_command_payload()
        return display_cmd, inp, user

    def update_abaqus_preview(self) -> None:
        if not hasattr(self, "abaqus_notes"):
            return
        cmd, inp, user = self._abaqus_command_string()
        cpus = self._safe_cpu_count() if hasattr(self, "var_cpus") else 1
        bat_text = (
            "cd /d " + self._quote_cmd_arg(str(inp.parent.resolve() if inp.parent.exists() else inp.parent)) + "\n"
            + cmd + "\n"
        )
        self.abaqus_notes.set(
            "Abaqus command preview:\n"
            + cmd
            + "\n\nEquivalent Windows .bat content:\n"
            + bat_text
            + "\nSelected run files:\n"
            + f"  Solver INP     : {inp}\n"
            + f"  Fortran source : {user}\n"
            + f"  Parallel CPUs  : {cpus}\n\n"
            + "Important: do NOT import the solver INP with File > Import > Model in Abaqus/CAE,\n"
            + "because CAE cannot import *User element / *Uel Property blocks.\n"
            + "Use the *_cae_preview.inp file only for CAE inspection, and submit the solver INP here or from an Abaqus Command Prompt.\n\n"
            + "If PyCharm reports WinError 2 for 'abaqus', open Abaqus Command Prompt or set Abaqus command to the full path of abaqus.bat. On Windows, direct submission is routed through cmd.exe with shell=False so abaqus.bat can be resolved; the generated .bat file remains the most reproducible workflow."
        )

    def save_run_bat(self) -> None:
        cmd, inp, user = self._abaqus_command_string()
        folder = inp.parent if str(inp.parent) else Path.cwd()
        default = (self.var_job.get().strip() or inp.stem or "Job_coh") + "_run.bat"
        f = filedialog.asksaveasfilename(
            initialdir=str(folder),
            initialfile=default,
            defaultextension=".bat",
            filetypes=[("Windows batch", "*.bat"), ("All files", "*.*")],
        )
        if not f:
            return
        bat = (
            "@echo off\n"
            + "cd /d " + self._quote_cmd_arg(str(folder.resolve() if folder.exists() else folder)) + "\n"
            + "call " + cmd + "\n"
            + "pause\n"
        )
        try:
            Path(f).write_text(bat, encoding="utf-8")
            messagebox.showinfo("run.bat generated", f"Batch file written to:\n{f}\n\nRun it inside an Abaqus Command Prompt if PyCharm cannot find Abaqus.")
        except OSError as e:
            messagebox.showerror("Save failed", str(e))

    def run_abaqus(self) -> None:
        if self.abaqus_process is not None and self.abaqus_process.poll() is None:
            messagebox.showinfo("Busy", "Abaqus is already running.")
            return
        inp = Path(self.var_run_inp.get() or self.var_output.get())
        uel = Path(self.var_run_fortran.get() or self.var_uel.get())
        if not inp.is_file():
            messagebox.showerror("Missing INP", f"Selected solver INP not found:\n{inp}")
            return
        if not uel.is_file():
            messagebox.showerror("Missing Fortran source", f"Selected Fortran source file not found:\n{uel}")
            return
        cpus = self._safe_cpu_count()
        cmd_str, run_args, _inp, _uel = self._abaqus_command_payload()
        self.update_abaqus_preview()
        self.abaqus_notes.append("\n[run] " + cmd_str + "\n")
        self.abaqus_status.set("Running...")

        def worker() -> None:
            try:
                popen_kwargs = dict(
                    cwd=str(inp.parent),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                    shell=False,
                )
                if sys.platform.startswith("win"):
                    popen_kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP  # type: ignore[attr-defined]
                else:
                    popen_kwargs["start_new_session"] = True
                self.abaqus_process = subprocess.Popen(run_args, **popen_kwargs)
                assert self.abaqus_process.stdout is not None
                for line in self.abaqus_process.stdout:
                    self.log_queue.put(line)
                code = self.abaqus_process.wait()
                self.log_queue.put(f"[abaqus] exit code = {code}\n")
            except Exception as e:
                self.log_queue.put(f"[abaqus] ERROR: {e}\n")
            finally:
                self.abaqus_process = None
                self.log_queue.put("__ABAQUS_DONE__")

        threading.Thread(target=worker, daemon=True).start()

    def stop_abaqus(self) -> None:
        """Stop a job launched from the GUI.

        Abaqus on Windows is commonly launched through abaqus.bat, which may
        spawn child processes.  Terminating only the cmd.exe wrapper is often
        insufficient, so the GUI first tries an Abaqus job termination command
        and then kills the process tree for the process it launched.  Jobs
        started outside the GUI, for example by manually running the generated
        .bat file, must be stopped from Abaqus or the operating system.
        """
        proc = self.abaqus_process
        if proc is None or proc.poll() is not None:
            self.abaqus_notes.append("\n[abaqus] no running GUI-launched process to stop\n")
            return
        job = (self.var_job.get() or "").strip()
        abaqus_cmd = (self.var_abaqus_cmd.get() or "abaqus").strip()
        try:
            if job:
                # Ask Abaqus to terminate the analysis cleanly when possible.
                if sys.platform.startswith("win"):
                    term_cmd = ["cmd.exe", "/d", "/c", f"{self._quote_cmd_arg(abaqus_cmd)} terminate job={self._quote_cmd_arg(job)}"]
                else:
                    term_cmd = [abaqus_cmd, "terminate", f"job={job}"]
                subprocess.run(term_cmd, cwd=str(Path(self.var_run_inp.get() or self.var_output.get()).parent),
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False, timeout=5)
        except Exception:
            pass
        try:
            if sys.platform.startswith("win"):
                subprocess.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
            else:
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            self.abaqus_notes.append("\n[abaqus] stop request sent to Abaqus job/process tree\n")
            self.abaqus_status.set("Stopping...")
        except Exception as e:
            try:
                proc.terminate()
                self.abaqus_notes.append("\n[abaqus] terminate signal sent to wrapper process\n")
            except OSError:
                self.abaqus_notes.append(f"\n[abaqus] stop failed: {e}\n")

    # ------------------------------------------------------------------ files
    def _build_files_page(self) -> None:
        p = self.pages["Files"]
        ttk.Label(p, text="Files", style="Title.TLabel").pack(anchor="w", pady=(4, 8))
        bar = ttk.Frame(p)
        bar.pack(fill="x", pady=(0, 4))
        ttk.Button(bar, text="Refresh", command=self.refresh_file_browser).pack(side="left", padx=2)
        ttk.Button(bar, text="Open output folder", command=lambda: self.open_folder(Path(self.var_output.get()).parent)).pack(side="left", padx=2)
        paned = ttk.PanedWindow(p, orient="horizontal")
        paned.pack(fill="both", expand=True)
        left = ttk.Frame(paned, padding=4)
        paned.add(left, weight=1)
        self.file_tree = ttk.Treeview(left, columns=("name", "size"), show="headings", height=24)
        self.file_tree.heading("name", text="Name")
        self.file_tree.heading("size", text="Size")
        self.file_tree.column("name", width=280, anchor="w")
        self.file_tree.column("size", width=90, anchor="e")
        self.file_tree.pack(fill="both", expand=True)
        self.file_tree.bind("<<TreeviewSelect>>", self.preview_selected_file)
        right = ttk.Frame(paned, padding=4)
        paned.add(right, weight=2)
        self.file_editor = TextEditor(right, height=28)
        self.file_editor.pack(fill="both", expand=True)
        ttk.Button(right, text="Save edits", style="Accent.TButton", command=self.save_selected_file).pack(anchor="e", pady=4)

    def refresh_file_browser(self) -> None:
        if not hasattr(self, "file_tree"):
            return
        self.file_tree.delete(*self.file_tree.get_children())
        folders = []
        for d in [Path(self.var_base.get()), Path(self.var_output.get()).parent]:
            if d.exists() and d not in folders:
                folders.append(d)
        for folder in folders:
            for path in sorted(folder.iterdir()):
                if not path.is_file():
                    continue
                if path.suffix.lower() not in {".inp", ".csv", ".txt", ".json", ".f", ".for", ".f90", ".sta", ".msg", ".dat"}:
                    continue
                try:
                    size = path.stat().st_size
                except OSError:
                    continue
                label = path.name if path.parent == Path(self.var_base.get()) else str(path)
                self.file_tree.insert("", "end", iid=str(path), values=(label, f"{size//1024} KB" if size > 1024 else f"{size} B"))

    def preview_selected_file(self, _event: Any = None) -> None:
        sel = self.file_tree.selection()
        if not sel:
            return
        path = Path(sel[0])
        self.selected_file = path
        try:
            text = read_text_auto(path)
            self.file_editor.set(text[:60000] + ("\n\n... <truncated>" if len(text) > 60000 else ""))
        except OSError as e:
            self.file_editor.set(f"Could not read file:\n{e}")

    def save_selected_file(self) -> None:
        if self.selected_file is None:
            messagebox.showinfo("No file", "Select a file first.")
            return
        if not messagebox.askyesno("Overwrite", f"Save edits to\n{self.selected_file}?"):
            return
        try:
            self.selected_file.write_text(self.file_editor.get(), encoding="utf-8")
        except OSError as e:
            messagebox.showerror("Save failed", str(e))

    # ---------------------------------------------------------------- theory
    def _build_theory_page(self) -> None:
        p = self.pages["Theory"]
        ttk.Label(p, text="Cohesive-zone theory and implementation notes", style="Title.TLabel").pack(anchor="w", pady=(4, 8))
        nb = ttk.Notebook(p)
        nb.pack(fill="both", expand=True)

        tab1 = ttk.Frame(nb, padding=6)
        tab2 = ttk.Frame(nb, padding=6)
        tab3 = ttk.Frame(nb, padding=6)
        nb.add(tab1, text="CZM formulas")
        nb.add(tab2, text="Abaqus UEL mapping")
        nb.add(tab3, text="Algorithm and verification")

        ed1 = TextEditor(tab1, height=30, wrap="word")
        ed1.pack(fill="both", expand=True)
        ed1.set(COHESIVE_THEORY_FORMULAS)
        ed1.text.configure(state="disabled")

        ed2 = TextEditor(tab2, height=30, wrap="word")
        ed2.pack(fill="both", expand=True)
        ed2.set(UEL_MAPPING_NOTES)
        ed2.text.configure(state="disabled")

        ed3 = TextEditor(tab3, height=30, wrap="word")
        ed3.pack(fill="both", expand=True)
        ed3.set(ALGORITHM_VERIFICATION_NOTES)
        ed3.text.configure(state="disabled")

    # ------------------------------------------------------------- config I/O
    def current_config(self) -> Dict[str, Any]:
        return {
            "base_dir": self.var_base.get(),
            "input_inp": self.var_input.get(),
            "output_inp": self.var_output.get(),
            "uel_file": self.var_uel.get(),
            "grain_set_prefix": self.var_grain_prefix.get(),
            "intmtd": self.var_intmtd.get(),
            "analysis_template": getattr(self, "var_model_template", tk.StringVar(value="")).get(),
            "nsvars_per_ip": getattr(self, "var_nsvars_per_ip", tk.StringVar(value="1")).get(),
            "property_preset": getattr(self, "var_property_preset", tk.StringVar(value="")).get(),
            "uel_type": self.var_uel_type.get(),
            "uel_elset": self.var_uel_elset.get(),
            "intra_elset": getattr(self, "var_intra_elset", tk.StringVar(value="INTRA_COH")).get(),
            "interface_scope": getattr(self, "var_interface_scope", tk.StringVar(value="Grain-boundary only")).get(),
            "intragranular_fraction": getattr(self, "var_intragranular_fraction", tk.StringVar(value="1.0")).get(),
            "random_seed": getattr(self, "var_random_seed", tk.StringVar(value="")).get(),
            "intra_props": {name: self.intra_prop_vars[name].get() for name, _ in self.PROP_LABELS} if hasattr(self, "intra_prop_vars") else {},
            "geometry_tolerance": self.var_tolerance.get(),
            "supplement_nsets": bool(self.var_supplement_nsets.get()),
            "write_cae_preview": bool(self.var_write_cae_preview.get()),
            "props": {name: self.prop_vars[name].get() for name, _ in self.PROP_LABELS},
            "abaqus_cmd": self.var_abaqus_cmd.get(),
            "job_name": self.var_job.get(),
            "run_inp": self.var_run_inp.get(),
            "run_fortran": self.var_run_fortran.get(),
            "run_cpus": self.var_cpus.get(),
        }

    def apply_config(self, cfg: Dict[str, Any]) -> None:
        self.var_base.set(str(cfg.get("base_dir", self.var_base.get())))
        self.var_input.set(str(cfg.get("input_inp", self.var_input.get())))
        self.var_output.set(str(cfg.get("output_inp", self.var_output.get())))
        self.var_uel.set(str(cfg.get("uel_file", self.var_uel.get())))
        self.var_grain_prefix.set(str(cfg.get("grain_set_prefix", self.var_grain_prefix.get())))
        self.var_intmtd.set(str(cfg.get("intmtd", self.var_intmtd.get())))
        if hasattr(self, "var_model_template"):
            self.var_model_template.set(str(cfg.get("analysis_template", self.var_model_template.get())))
        if hasattr(self, "var_nsvars_per_ip"):
            self.var_nsvars_per_ip.set(str(cfg.get("nsvars_per_ip", self.var_nsvars_per_ip.get())))
        if hasattr(self, "var_property_preset"):
            self.var_property_preset.set(str(cfg.get("property_preset", self.var_property_preset.get())))
        if hasattr(self, "var_uel_type"):
            self.var_uel_type.set(str(cfg.get("uel_type", self.var_uel_type.get())))
        if hasattr(self, "var_uel_elset"):
            self.var_uel_elset.set(str(cfg.get("uel_elset", self.var_uel_elset.get())))
        if hasattr(self, "var_intra_elset"):
            self.var_intra_elset.set(str(cfg.get("intra_elset", self.var_intra_elset.get())))
        if hasattr(self, "var_interface_scope"):
            self.var_interface_scope.set(str(cfg.get("interface_scope", self.var_interface_scope.get())))
        if hasattr(self, "var_intragranular_fraction"):
            self.var_intragranular_fraction.set(str(cfg.get("intragranular_fraction", self.var_intragranular_fraction.get())))
        if hasattr(self, "var_random_seed"):
            self.var_random_seed.set(str(cfg.get("random_seed", self.var_random_seed.get())))
        if hasattr(self, "var_tolerance"):
            self.var_tolerance.set(str(cfg.get("geometry_tolerance", self.var_tolerance.get())))
        if hasattr(self, "var_supplement_nsets"):
            self.var_supplement_nsets.set(bool(cfg.get("supplement_nsets", self.var_supplement_nsets.get())))
        if hasattr(self, "var_write_cae_preview"):
            self.var_write_cae_preview.set(bool(cfg.get("write_cae_preview", self.var_write_cae_preview.get())))
        props = cfg.get("props", {})
        if isinstance(props, dict):
            for name, _ in self.PROP_LABELS:
                if name in props:
                    self.prop_vars[name].set(str(props[name]))
        intra_props_cfg = cfg.get("intra_props", {})
        if isinstance(intra_props_cfg, dict) and hasattr(self, "intra_prop_vars"):
            for name, _ in self.PROP_LABELS:
                if name in intra_props_cfg:
                    self.intra_prop_vars[name].set(str(intra_props_cfg[name]))
        self.var_abaqus_cmd.set(str(cfg.get("abaqus_cmd", self.var_abaqus_cmd.get())))
        self.var_job.set(str(cfg.get("job_name", self.var_job.get())))
        self.var_run_inp.set(str(cfg.get("run_inp", self.var_run_inp.get() or self.var_output.get())))
        self.var_run_fortran.set(str(cfg.get("run_fortran", self.var_run_fortran.get() or self.var_uel.get())))
        self.var_cpus.set(str(cfg.get("run_cpus", self.var_cpus.get())))
        self.update_abaqus_preview()
        self.update_dashboard()
        self.preview_input_file()

    def save_config_json(self) -> None:
        f = filedialog.asksaveasfilename(
            initialdir=self.var_base.get() or str(Path.cwd()),
            initialfile="cohesive_config.json",
            defaultextension=".json",
            filetypes=[("JSON", "*.json"), ("All files", "*.*")],
        )
        if not f:
            return
        Path(f).write_text(json.dumps(self.current_config(), indent=2), encoding="utf-8")
        messagebox.showinfo("Saved", f"Configuration saved to\n{f}")

    def load_config_json(self) -> None:
        f = filedialog.askopenfilename(
            initialdir=self.var_base.get() or str(Path.cwd()),
            filetypes=[("JSON", "*.json"), ("All files", "*.*")],
        )
        if not f:
            return
        try:
            cfg = json.loads(Path(f).read_text(encoding="utf-8"))
            self.apply_config(cfg)
        except Exception as e:
            messagebox.showerror("Load failed", str(e))

    # ------------------------------------------------------------- utilities
    def open_folder(self, path: Path) -> None:
        path = Path(path)
        if not path.exists():
            messagebox.showinfo("Open folder", f"Path does not exist:\n{path}")
            return
        try:
            if sys.platform.startswith("win"):
                os.startfile(str(path))  # type: ignore[attr-defined]
            elif sys.platform == "darwin":
                subprocess.Popen(["open", str(path)])
            else:
                subprocess.Popen(["xdg-open", str(path)])
        except OSError as e:
            messagebox.showerror("Open folder failed", str(e))
def launch_gui() -> None:
    """Launch the graphical user interface."""
    if not TK_AVAILABLE:
        raise RuntimeError(
            "Tkinter is not available in this Python environment. "
            "Install Tkinter or run a command-line generation task instead."
        )
    app = CohesiveXStudio()
    app.mainloop()


def main() -> None:
    """Entry point used by ``python -m cohesivex_studio`` and the console script.

    With no positional input file, the graphical interface is launched.  With an
    input/output pair, the same package can be used as a command-line cohesive
    generator.  The ``--self-test`` option exercises the kernel without opening
    the GUI.
    """
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "CohesiveX Studio: GUI/CLI tool for inserting UEL cohesive elements "
            "at grain-boundary and/or intragranular interfaces."
        )
    )
    parser.add_argument("input", nargs="?", help="Input Abaqus .inp file. If omitted, the GUI opens.")
    parser.add_argument("output", nargs="?", default="Job_coh.inp", help="Output Abaqus .inp file")
    parser.add_argument("--gui", action="store_true", help="Force GUI mode")
    parser.add_argument("--version", action="version", version=f"CohesiveX Studio {__version__}")
    parser.add_argument("--self-test", action="store_true", help="Run kernel regression tests and exit")
    parser.add_argument("--grain-prefix", default="GRAIN-", help="Grain/domain elset prefix")
    parser.add_argument("--intmtd", type=int, default=1, choices=[1, 2, 3], help="Integration method used for UEL variable count")
    parser.add_argument("--uel-type", default="U1", help="Abaqus UEL type name inserted in *User element")
    parser.add_argument("--uel-elset", default="GB_COH", help="Element-set name for grain-boundary cohesive elements")
    parser.add_argument("--intra-elset", default="INTRA_COH", help="Element-set name for intragranular cohesive elements")
    parser.add_argument("--interface-scope", default="grain_boundary", choices=["grain_boundary", "intragranular", "both"], help="Where to insert cohesive elements")
    parser.add_argument("--intra-fraction", type=float, default=1.0, help="Fraction of intragranular interfaces to insert when intragranular scope is active")
    parser.add_argument("--random-seed", type=int, default=None, help="Random seed for optional intragranular sampling")
    parser.add_argument("--tolerance", type=float, default=1.0e-10, help="Absolute geometry tolerance for coordinate matching")
    parser.add_argument("--no-nset-supplement", action="store_true", help="Do not supplement original Nsets with duplicated nodes")
    parser.add_argument("--no-cae-preview", action="store_true", help="Do not write a CAE preview INP")
    parser.add_argument("--props", nargs=9, type=float, metavar=("MODE", "KI", "KII", "SI", "SII", "GCI", "GCII", "ETA", "HEIGHT"), help="Nine grain-boundary cohesive properties")
    parser.add_argument("--intra-props", nargs=9, type=float, metavar=("MODE", "KI", "KII", "SI", "SII", "GCI", "GCII", "ETA", "HEIGHT"), help="Nine intragranular cohesive properties")
    parser.add_argument("--backend", choices=["numpy", "python"], default="numpy", help="Interface-detection backend. Default: numpy")
    parser.add_argument("--no-fast", action="store_true", help="Disable NumPy-accelerated face hashing backend; equivalent to --backend python")
    parser.add_argument("--compare-backends", action="store_true", help="Run NumPy and pure-Python backends and write a backend comparison report")
    args = parser.parse_args()

    if args.self_test:
        run_self_tests()
        return

    if args.gui or not args.input:
        launch_gui()
        return

    fast_mode = (args.backend == "numpy") and not args.no_fast
    if args.compare_backends:
        compare_backends(
            args.input,
            Path(args.output).with_suffix(""),
            grain_set_prefix=args.grain_prefix,
            intmtd=args.intmtd,
            props=tuple(args.props) if args.props else (3.0, 1.0e7, 1.0e7, 100.0, 100.0, 5.0, 5.0, 2.0, 1.0),
            gb_props=tuple(args.props) if args.props else None,
            intra_props=tuple(args.intra_props) if args.intra_props else (3.0, 1.0e7, 1.0e7, 180.0, 180.0, 10.0, 10.0, 2.0, 1.0),
            interface_scope=args.interface_scope,
            uel_type=args.uel_type,
            uel_elset=args.uel_elset,
            gb_elset=args.uel_elset,
            intra_elset=args.intra_elset,
            intragranular_fraction=args.intra_fraction,
            random_seed=args.random_seed,
            supplement_nsets=not args.no_nset_supplement,
            write_cae_preview=not args.no_cae_preview,
            tolerance=args.tolerance,
        )
        return

    generate_cohesive_inp(
        args.input,
        args.output,
        grain_set_prefix=args.grain_prefix,
        intmtd=args.intmtd,
        props=tuple(args.props) if args.props else (3.0, 1.0e7, 1.0e7, 100.0, 100.0, 5.0, 5.0, 2.0, 1.0),
        gb_props=tuple(args.props) if args.props else None,
        intra_props=tuple(args.intra_props) if args.intra_props else (3.0, 1.0e7, 1.0e7, 180.0, 180.0, 10.0, 10.0, 2.0, 1.0),
        interface_scope=args.interface_scope,
        uel_type=args.uel_type,
        uel_elset=args.uel_elset,
        gb_elset=args.uel_elset,
        intra_elset=args.intra_elset,
        intragranular_fraction=args.intra_fraction,
        random_seed=args.random_seed,
        fast_mode=fast_mode,
        supplement_nsets=not args.no_nset_supplement,
        write_cae_preview=not args.no_cae_preview,
        tolerance=args.tolerance,
    )


if __name__ == "__main__":
    main()
