# -*- coding: utf-8 -*-
"""
CohesiveX Studio
================

A GUI-assisted Abaqus preprocessing platform for inter-domain and intra-domain
cohesive-zone model generation.  Design principle
----------------
The generator performs cohesive-only modifications.  It preserves user-defined
materials, sections, steps, amplitudes, loads, controls, output requests and
original time-incrementation settings, while adding only the duplicated nodes,
updated solid connectivity, UEL cohesive blocks and optional node-set
supplements required for the cohesive model.
"""
from __future__ import annotations

from pathlib import Path
from itertools import combinations
import csv
import json
import math
import re
import hashlib
from time import perf_counter
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple


def read_text_auto(path: str | Path, warnings: Optional[List[str]] = None) -> str:
    """Read text without silently replacing undecodable bytes.

    Abaqus input files are usually UTF-8/ASCII, but Windows cases may contain
    GBK comments.  We try common encodings and record a warning when the file is
    not clean UTF-8.  This avoids the old ``errors="replace"`` behaviour,
    which could silently turn non-ASCII bytes into replacement characters.
    """
    path = Path(path)
    data = path.read_bytes()
    for enc in ("utf-8-sig", "utf-8", "gbk"):
        try:
            text = data.decode(enc)
            if enc not in {"utf-8", "utf-8-sig"} and warnings is not None:
                warnings.append(f"Read {path.name} using {enc} encoding; please verify non-ASCII comments if the file is edited externally.")
            return text
        except UnicodeDecodeError:
            continue
    if warnings is not None:
        warnings.append(f"Could not decode {path.name} as UTF-8 or GBK; using latin1 fallback without byte replacement.")
    return data.decode("latin1")


def generate_cohesive_inp(
    input_inp: str | Path,
    output_inp: str | Path = "Job_coh.inp",
    *,
    grain_set_prefix: str = "GRAIN-",
    intmtd: int = 1,
    props: Sequence[float] = (3.0, 1.0e7, 1.0e7, 100.0, 100.0, 5.0, 5.0, 2.0, 1.0),
    uel_type: str = "U1",
    uel_elset: str = "GB_COH",
    interface_scope: str = "grain_boundary",
    gb_props: Optional[Sequence[float]] = None,
    intra_props: Optional[Sequence[float]] = None,
    gb_elset: Optional[str] = None,
    intra_elset: str = "INTRA_COH",
    nsvars_per_ip: int = 1,
    fast_mode: bool = True,
    intragranular_fraction: float = 1.0,
    random_seed: Optional[int] = None,
    supplement_nsets: bool = True,
    write_cae_preview: bool = True,
    tolerance: float = 1.0e-10,
    zero_tolerance: float = 1.0e-14,
    report_prefix: Optional[str | Path] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Insert UEL cohesive elements at selected inter-/intra-domain interfaces.

    This upgraded kernel supports three insertion scopes:

    ``grain_boundary``
        Insert cohesive elements only between neighbouring elements that belong
        to different grain/domain sets.
    ``intragranular``
        Insert cohesive elements only between neighbouring elements that belong
        to the same grain/domain set.  If no grain sets are present, the whole
        parsed solid element block is treated as one domain; this makes the
        tool useful for conventional elastic-plastic models with potential
        intradomain cracking.
    ``both``
        Insert cohesive elements at both grain-boundary and intragranular
        interfaces, assigning different elsets and property vectors to each
        family.

    The operation remains cohesive-only: it preserves the original Abaqus
    material cards, sections, steps, controls, loads, amplitudes, output
    requests and time-incrementation settings.  Only the parsed node block, the
    parsed solid element connectivity, UEL cohesive blocks and optional Nset
    supplements are modified/inserted.
    """

    input_inp = Path(input_inp)
    output_inp = Path(output_inp)
    output_inp.parent.mkdir(parents=True, exist_ok=True)
    if report_prefix is None:
        report_prefix = output_inp.with_suffix("")
    report_prefix = Path(report_prefix)
    report_prefix.parent.mkdir(parents=True, exist_ok=True)

    warnings: List[str] = []
    step_timings: Dict[str, float] = {}
    _t_total = perf_counter()
    _t_last = _t_total

    def mark_step(name: str) -> None:
        nonlocal _t_last
        now = perf_counter()
        step_timings[name] = now - _t_last
        _t_last = now

    def log(msg: str) -> None:
        if verbose:
            print(msg)

    def normalize_scope(value: str) -> str:
        s = str(value or "grain_boundary").strip().lower().replace(" ", "_").replace("-", "_")
        if s in {"gb", "grain", "grain_boundary", "grain_boundary_only", "intergranular", "inter_domain", "interdomain"}:
            return "grain_boundary"
        if s in {"intra", "intragranular", "intragranular_only", "intra_domain", "intradomain", "within_grain"}:
            return "intragranular"
        if s in {"both", "all", "grain_boundary_+_intragranular", "grain_boundary_and_intragranular", "inter_and_intra"}:
            return "both"
        raise ValueError(f"Unknown interface_scope: {value!r}. Use grain_boundary, intragranular, or both.")

    scope = normalize_scope(interface_scope)
    uel_type = str(uel_type or "U1").strip().upper()
    if gb_elset is None:
        gb_elset = str(uel_elset or "GB_COH").strip()
    gb_elset = str(gb_elset or "GB_COH").strip()
    intra_elset = str(intra_elset or "INTRA_COH").strip()
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_\-]*", uel_type):
        raise ValueError(f"Invalid UEL type name: {uel_type!r}")
    for elset_name in (gb_elset, intra_elset):
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_\-]*", elset_name):
            raise ValueError(f"Invalid cohesive elset name: {elset_name!r}")

    def to_props(seq: Optional[Sequence[float]], fallback: Sequence[float]) -> Tuple[float, ...]:
        vals = tuple(float(x) for x in (seq if seq is not None else fallback))
        if len(vals) != 9:
            raise ValueError("Each cohesive property vector must contain exactly 9 values.")
        return vals

    base_props = to_props(props, props)
    gb_props_t = to_props(gb_props, base_props)
    intra_props_t = to_props(intra_props, base_props)
    nsvars_per_ip = max(1, int(nsvars_per_ip))
    intragranular_fraction = float(intragranular_fraction)
    if not (0.0 <= intragranular_fraction <= 1.0):
        raise ValueError("intragranular_fraction must be in [0, 1].")
    if intragranular_fraction < 1.0 and random_seed is None:
        random_seed_effective = 0
        warnings.append(
            "Intragranular fraction is smaller than 1.0 but random_seed was not provided; "
            "using deterministic default seed 0 for reproducible sampling."
        )
    else:
        random_seed_effective = int(random_seed) if random_seed is not None else None

    num_re = re.compile(r"[-+]?(?:\d*\.\d+|\d+\.?)(?:[EeDd][-+]?\d+)?")

    def strip_numeric_payload(line: str) -> str:
        stripped = line.lstrip()
        if not stripped or stripped.startswith("**"):
            return ""
        return re.split(r"[#!]", line, maxsplit=1)[0].strip()

    def numeric_values(line: str) -> List[float]:
        payload = strip_numeric_payload(line)
        if not payload:
            return []
        return [float(x.replace("D", "E").replace("d", "e")) for x in num_re.findall(payload)]

    def numeric_ints(line: str) -> List[int]:
        payload = strip_numeric_payload(line)
        if not payload:
            return []
        vals: List[int] = []
        for token in re.split(r"[,\s]+", payload):
            if not token:
                continue
            try:
                vals.append(int(float(token.replace("D", "E").replace("d", "e"))))
            except ValueError:
                continue
        return vals

    def is_keyword(line: str) -> bool:
        stripped = line.lstrip()
        return stripped.startswith("*") and not stripped.startswith("**")

    def keyword_name(line: str) -> str:
        if not is_keyword(line):
            return ""
        return line.strip().split(",", 1)[0].lower()

    def parse_keyword_options(line: str) -> Dict[str, str]:
        opts: Dict[str, str] = {}
        if not is_keyword(line):
            return opts
        for part in line.strip().split(",")[1:]:
            part = part.strip()
            if not part:
                continue
            if "=" in part:
                k, v = part.split("=", 1)
                opts[k.strip().lower()] = v.strip()
            else:
                opts[part.strip().lower()] = ""
        return opts

    def expand_generate_ids(start: int, stop: int, step: int, context: str) -> List[int]:
        if step == 0:
            warnings.append(f"{context}: ignored generate range with step=0 ({start}, {stop}, {step}).")
            return []
        if (stop - start) * step < 0:
            warnings.append(f"{context}: generate range direction is inconsistent ({start}, {stop}, {step}); no IDs expanded.")
            return []
        return list(range(start, stop + (1 if step > 0 else -1), step))

    def validate_props(p: Sequence[float], label: str) -> None:
        mode = int(round(float(p[0])))
        if mode not in (1, 2, 3):
            raise ValueError(f"{label}: MODE must be 1, 2, or 3.")
        if mode in (1, 3):
            if p[1] <= 0 or p[3] <= 0 or p[5] <= 0:
                raise ValueError(f"{label}: KI, SI and GCI must be positive for opening/mixed mode.")
        if mode in (2, 3):
            if p[2] <= 0 or p[4] <= 0 or p[6] <= 0:
                raise ValueError(f"{label}: KII, SII and GCII must be positive for shear/mixed mode.")
        if p[7] < 0:
            raise ValueError(f"{label}: ETA must be non-negative.")
        if p[8] <= 0:
            raise ValueError(f"{label}: HEIGHT must be positive for 2D reports/UEL property consistency.")

    if scope in ("grain_boundary", "both"):
        validate_props(gb_props_t, "grain-boundary cohesive properties")
    if scope in ("intragranular", "both"):
        validate_props(intra_props_t, "intragranular cohesive properties")

    geom_tol = max(float(tolerance), 1.0e-30)
    zero_tol = max(float(zero_tolerance), 1.0e-30)

    def coord_token(pt: Sequence[float]) -> Tuple[int, ...]:
        return tuple(int(round(float(x) / geom_tol)) for x in pt)

    def vec_sub(a: Sequence[float], b: Sequence[float]) -> List[float]:
        return [float(x) - float(y) for x, y in zip(a, b)]

    def vec_add(a: Sequence[float], b: Sequence[float]) -> List[float]:
        return [float(x) + float(y) for x, y in zip(a, b)]

    def vec_dot(a: Sequence[float], b: Sequence[float]) -> float:
        return float(sum(float(x) * float(y) for x, y in zip(a, b)))

    def vec_cross(a: Sequence[float], b: Sequence[float]) -> List[float]:
        return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]

    def vec_norm(a: Sequence[float]) -> float:
        return math.sqrt(sum(float(x) * float(x) for x in a))

    def vec_unit(a: Sequence[float]) -> List[float]:
        n = vec_norm(a)
        if n <= zero_tol:
            return [0.0 for _ in a]
        return [float(x) / n for x in a]

    def centroid(node_ids: Sequence[int], crd: Dict[int, Tuple[float, ...]]) -> Tuple[float, ...]:
        dim = len(next(iter(crd.values())))
        s = [0.0] * dim
        for nid in node_ids:
            s = vec_add(s, crd[nid])
        return tuple(x / max(1, len(node_ids)) for x in s)

    def dist(a: Sequence[float], b: Sequence[float]) -> float:
        return vec_norm(vec_sub(a, b))

    def element_info(eltype_raw: str, integration_method: int) -> Dict[str, Any]:
        eltype = eltype_raw.upper().strip()
        if eltype in {"CPS3", "CPE3"}:
            return dict(eltype=eltype, eltyp=6, nnpel=3, nsurf=3, nnps=2, dim=2, surfaces=[[0, 1], [1, 2], [2, 0]], numpt=2)
        if eltype in {"CPS4", "CPE4", "CPS4R", "CPE4R"}:
            return dict(eltype=eltype, eltyp=8, nnpel=4, nsurf=4, nnps=2, dim=2, surfaces=[[0, 1], [1, 2], [2, 3], [3, 0]], numpt=2)
        if eltype in {"CPS6", "CPE6", "CPS6R", "CPE6R"}:
            return dict(eltype=eltype, eltyp=12, nnpel=6, nsurf=3, nnps=3, dim=2, surfaces=[[0, 1, 3], [1, 2, 4], [2, 0, 5]], numpt=3 if integration_method in (1, 2) else 2)
        if eltype in {"CPS8", "CPE8", "CPS8R", "CPE8R"}:
            return dict(eltype=eltype, eltyp=16, nnpel=8, nsurf=4, nnps=3, dim=2, surfaces=[[0, 1, 4], [1, 2, 5], [2, 3, 6], [3, 0, 7]], numpt=3 if integration_method in (1, 2) else 2)
        if eltype == "C3D4":
            return dict(eltype=eltype, eltyp=12, nnpel=4, nsurf=4, nnps=3, dim=3, surfaces=[[0, 1, 2], [0, 2, 3], [2, 1, 3], [0, 3, 1]], numpt=1)
        if eltype in {"C3D8", "C3D8R"}:
            return dict(eltype=eltype, eltyp=24, nnpel=8, nsurf=6, nnps=4, dim=3, surfaces=[[3, 2, 1, 0], [4, 5, 6, 7], [0, 4, 7, 3], [1, 2, 6, 5], [0, 1, 5, 4], [2, 3, 7, 6]], numpt=4)
        if eltype == "C3D10":
            return dict(eltype=eltype, eltyp=30, nnpel=10, nsurf=4, nnps=6, dim=3, surfaces=[[0, 1, 2, 4, 5, 6], [0, 2, 3, 6, 9, 7], [2, 1, 3, 5, 8, 9], [0, 3, 1, 7, 8, 4]], numpt=4 if integration_method == 1 else (6 if integration_method == 2 else 3))
        if eltype in {"C3D20", "C3D20R"}:
            return dict(eltype=eltype, eltyp=60, nnpel=20, nsurf=6, nnps=8, dim=3, surfaces=[[3, 2, 1, 0, 10, 9, 8, 11], [4, 5, 6, 7, 12, 13, 14, 15], [0, 4, 7, 3, 16, 15, 19, 11], [1, 2, 6, 5, 9, 18, 13, 17], [0, 1, 5, 4, 8, 17, 12, 16], [2, 3, 7, 6, 10, 19, 14, 18]], numpt=9 if integration_method in (1, 2) else 4)
        raise ValueError(f"Unsupported Abaqus element type: {eltype_raw}")

    def corner_face_nodes(face_nodes: Sequence[int]) -> List[int]:
        n = len(face_nodes)
        if n in (2, 3, 4):
            return list(face_nodes)
        if n == 6:
            return list(face_nodes[:3])
        if n == 8:
            return list(face_nodes[:4])
        return list(face_nodes)

    def polygon_area_3d(node_ids: Sequence[int], crd: Dict[int, Tuple[float, ...]]) -> float:
        pts = [crd[n] for n in node_ids]
        if len(pts) < 3:
            return 0.0
        p0 = pts[0]
        area = 0.0
        for ii in range(1, len(pts) - 1):
            area += 0.5 * vec_norm(vec_cross(vec_sub(pts[ii], p0), vec_sub(pts[ii + 1], p0)))
        return area

    def interface_measure(face_nodes: Sequence[int], crd: Dict[int, Tuple[float, ...]], dim: int, height: float) -> float:
        if dim == 2:
            if len(face_nodes) < 2:
                return 0.0
            return dist(crd[face_nodes[0]], crd[face_nodes[1]]) * height
        return polygon_area_3d(corner_face_nodes(face_nodes), crd)

    def reorder_if_needed(data: List[int], info: Dict[str, Any]) -> List[int]:
        eltyp = info["eltyp"]
        dim = info["dim"]
        if (eltyp == 6 or eltyp == 8) and dim == 2:
            order = [1, 0, 3, 2]
        elif (eltyp == 12 or eltyp == 16) and dim == 2:
            order = [1, 0, 2, 4, 3, 5]
        elif eltyp == 12 and dim == 3:
            order = [0, 2, 1, 3, 5, 4]
        elif eltyp == 30:
            order = [0, 2, 1, 5, 4, 3, 6, 8, 7, 11, 10, 9]
        elif eltyp == 24:
            order = [0, 3, 2, 1, 4, 7, 6, 5]
        elif eltyp == 60:
            order = [0, 3, 2, 1, 7, 6, 5, 4, 8, 11, 10, 9, 15, 14, 13, 12]
        else:
            order = list(range(len(data)))
        if len(order) != len(data):
            warnings.append(f"Internal reordering length mismatch for eltyp={eltyp}; connectivity left unchanged.")
            return data
        return [data[i] for i in order]

    def face_normal(face_nodes: Sequence[int], crd: Dict[int, Tuple[float, ...]], dim: int) -> List[float]:
        corners = corner_face_nodes(face_nodes)
        if dim == 2:
            if len(corners) < 2:
                return [0.0, 0.0]
            t = vec_unit(vec_sub(crd[corners[1]], crd[corners[0]]))
            return [-t[1], t[0]]
        if len(corners) < 3:
            return [0.0, 0.0, 0.0]
        pts = [crd[n] for n in corners]
        for ii, jj in combinations(range(1, len(pts)), 2):
            n = vec_cross(vec_sub(pts[ii], pts[0]), vec_sub(pts[jj], pts[0]))
            if vec_norm(n) > zero_tol:
                return vec_unit(n)
        return [0.0, 0.0, 0.0]

    class DSU:
        def __init__(self, items: Sequence[int]) -> None:
            self.parent = {int(x): int(x) for x in items}
            self.rank = {int(x): 0 for x in items}
        def find(self, x: int) -> int:
            x = int(x)
            p = self.parent[x]
            if p != x:
                self.parent[x] = self.find(p)
            return self.parent[x]
        def union(self, a: int, b: int) -> None:
            ra, rb = self.find(a), self.find(b)
            if ra == rb:
                return
            if self.rank[ra] < self.rank[rb]:
                ra, rb = rb, ra
            self.parent[rb] = ra
            if self.rank[ra] == self.rank[rb]:
                self.rank[ra] += 1

    if not input_inp.is_file():
        raise FileNotFoundError(input_inp)
    lines = read_text_auto(input_inp, warnings).splitlines()

    node_blocks = [i for i, line in enumerate(lines) if keyword_name(line) == "*node"]
    if not node_blocks:
        raise ValueError("Could not locate *Node block.")
    node_start = node_blocks[0]
    if len(node_blocks) > 1:
        warnings.append("Multiple *Node blocks were detected. This release modifies only the first supported mesh block; other blocks are preserved unchanged.")

    coords: Dict[int, Tuple[float, ...]] = {}
    i = node_start + 1
    while i < len(lines) and not is_keyword(lines[i]):
        vals = numeric_values(lines[i])
        if len(vals) >= 3:
            coords[int(vals[0])] = tuple(float(x) for x in vals[1:])
        i += 1
    if not coords:
        raise ValueError("No nodes were parsed from *Node block.")
    numdim = len(next(iter(coords.values())))

    elem_start = None
    elem_type = None
    supported_element_blocks: List[Tuple[int, str]] = []
    for k, line in enumerate(lines):
        if keyword_name(line) == "*element":
            opts = parse_keyword_options(line)
            et = opts.get("type", "").strip().upper()
            if et:
                try:
                    element_info(et, intmtd)
                    supported_element_blocks.append((k, et))
                    if elem_start is None:
                        elem_start = k
                        elem_type = et
                except ValueError:
                    continue
    if elem_start is None or elem_type is None:
        raise ValueError("Could not locate a supported *Element block.")
    if len(supported_element_blocks) > 1:
        block_types = ", ".join(f"line {idx+1}: {typ}" for idx, typ in supported_element_blocks)
        warnings.append("Multiple supported solid *Element blocks were detected. This release processes only the first supported block. " + f"Detected blocks: {block_types}.")
    info = element_info(elem_type, intmtd)
    nnpel = int(info["nnpel"])
    nnps = int(info["nnps"])
    surfaces = info["surfaces"]
    numpt = int(info["numpt"])
    if info["dim"] != numdim:
        warnings.append(f"Element type {elem_type} suggests dim={info['dim']}, but node coordinates have dim={numdim}.")

    solid_conn: Dict[int, List[int]] = {}
    element_order: List[int] = []
    i = elem_start + 1
    while i < len(lines) and not is_keyword(lines[i]):
        vals = numeric_ints(lines[i])
        if not vals:
            i += 1
            continue
        eid = vals[0]
        nodes = vals[1:]
        i += 1
        while len(nodes) < nnpel and i < len(lines) and not is_keyword(lines[i]):
            more = numeric_ints(lines[i])
            if more:
                nodes.extend(more)
            i += 1
        if len(nodes) < nnpel:
            raise ValueError(f"Element {eid} has {len(nodes)} nodes, expected {nnpel}.")
        solid_conn[eid] = nodes[:nnpel]
        element_order.append(eid)
    if not solid_conn:
        raise ValueError("No solid elements were parsed.")

    grain_by_element: Dict[int, int] = {}
    grain_pat = re.compile(re.escape(grain_set_prefix) + r"\s*([+-]?\d+)", re.IGNORECASE)
    for idx, line in enumerate(lines):
        if keyword_name(line) != "*elset":
            continue
        opts = parse_keyword_options(line)
        elset_name = opts.get("elset", "")
        m = grain_pat.search(elset_name)
        if not m:
            continue
        gid = int(m.group(1))
        generate = "generate" in opts or line.lower().strip().endswith("generate")
        ids: List[int] = []
        j = idx + 1
        if generate:
            while j < len(lines) and not is_keyword(lines[j]):
                vals = numeric_ints(lines[j])
                if len(vals) >= 3:
                    ids.extend(expand_generate_ids(vals[0], vals[1], vals[2], f"Elset {elset_name}"))
                    break
                j += 1
        else:
            while j < len(lines) and not is_keyword(lines[j]):
                ids.extend(numeric_ints(lines[j]))
                j += 1
        for eid in ids:
            if eid in solid_conn:
                grain_by_element[eid] = gid

    if not grain_by_element:
        warnings.append("No grain/domain elsets were found with the selected prefix. The parsed solid block is treated as one domain (grain ID 0); use intragranular scope for ordinary elastic-plastic cohesive insertion.")
    missing = [eid for eid in element_order if eid not in grain_by_element]
    if missing:
        if grain_by_element:
            warnings.append(f"{len(missing)} solid elements have no grain/domain ID; they are assigned to domain 0.")
        for eid in missing:
            grain_by_element[eid] = 0
    grain_order = sorted(set(grain_by_element.values()))
    mark_step("parse_input")

    # -------- Fast face hashing backend --------
    face_records: List[Tuple[int, int, Tuple[int, ...], int]] = []
    face_groups: List[List[int]] = []
    used_numpy_backend = False
    try:
        if fast_mode:
            import numpy as np  # type: ignore
            eids_arr = np.asarray(element_order, dtype=np.int64)
            conn_arr = np.asarray([solid_conn[int(eid)] for eid in element_order], dtype=np.int64)
            missing_node_refs = sorted(set(int(x) for x in conn_arr.ravel().tolist()) - set(coords.keys()))
            if missing_node_refs:
                raise ValueError(
                    "Solid connectivity references node IDs that are not present in the parsed *Node block: "
                    + ", ".join(str(x) for x in missing_node_refs[:20])
                    + (" ..." if len(missing_node_refs) > 20 else "")
                )
            surf_arr = np.asarray(surfaces, dtype=np.int64)
            faces_nodes = conn_arr[:, surf_arr]  # (nelem, nsurf, nnps)

            # Map coordinate tokens to compact integer IDs.  The lookup length is
            # based on both coordinate and connectivity node IDs so sparse Abaqus
            # labels remain safe.  A zero token after lookup indicates a broken
            # node reference and is treated as an error rather than silently
            # producing false face matches.
            token_to_id: Dict[Tuple[int, ...], int] = {}
            node_token_id: Dict[int, int] = {}
            for nid, xyz in coords.items():
                tok = coord_token(xyz)
                if tok not in token_to_id:
                    token_to_id[tok] = len(token_to_id) + 1
                node_token_id[nid] = token_to_id[tok]
            max_lookup_id = max(max(coords), int(conn_arr.max()))
            token_lookup = np.zeros(max_lookup_id + 1, dtype=np.int64)
            for nid, tid in node_token_id.items():
                token_lookup[int(nid)] = int(tid)
            face_token = token_lookup[faces_nodes]
            if int(face_token.min()) == 0:
                raise ValueError("One or more solid faces contain an unmapped node token; check the parsed node block and element connectivity.")

            face_keys = np.sort(face_token.reshape(-1, nnps), axis=1)
            _unique_keys, inverse, counts = np.unique(face_keys, axis=0, return_inverse=True, return_counts=True)
            nsurf = len(surfaces)
            flat_faces = faces_nodes.reshape(-1, nnps)
            face_records = []
            for flat_idx, fnodes in enumerate(flat_faces):
                eid = int(eids_arr[flat_idx // nsurf])
                fid = int(flat_idx % nsurf) + 1
                face_records.append((eid, fid, tuple(int(x) for x in fnodes.tolist()), int(grain_by_element[eid])))

            grouped_indices_by_uid: Dict[int, List[int]] = {}
            for flat_idx, uid in enumerate(inverse.tolist()):
                grouped_indices_by_uid.setdefault(int(uid), []).append(int(flat_idx))
            face_groups = [inds for uid, inds in grouped_indices_by_uid.items() if int(counts[int(uid)]) >= 2]
            used_numpy_backend = True
    except Exception as exc:
        warnings.append(f"NumPy fast face hashing was unavailable or failed ({exc}); falling back to pure-Python face hashing.")
        used_numpy_backend = False

    if not used_numpy_backend:
        face_map: Dict[Tuple[Tuple[int, ...], ...], List[int]] = {}
        for eid in element_order:
            nodes = solid_conn[eid]
            gid = grain_by_element[eid]
            for fid, surf in enumerate(surfaces, start=1):
                face_orig = tuple(nodes[ii] for ii in surf)
                key = tuple(sorted(coord_token(coords[nid]) for nid in face_orig))
                face_records.append((eid, fid, face_orig, gid))
                face_map.setdefault(key, []).append(len(face_records) - 1)
        face_groups = [inds for inds in face_map.values() if len(inds) >= 2]

    mark_step("face_detection")

    selected_interfaces: List[Dict[str, Any]] = []
    all_internal_pairs: List[Tuple[int, int]] = []
    selected_pair_keys: set[Tuple[int, int, int, int]] = set()
    nonmanifold_faces = 0

    def stable_sample_unit(*parts: Any) -> float:
        """Return a deterministic pseudo-random number in [0, 1).

        The value depends only on the interface identity and the selected seed,
        not on Python iteration order.  This keeps NumPy fast mode and pure
        Python mode reproducible for the same mesh and settings.
        """
        seed = 0 if random_seed_effective is None else int(random_seed_effective)
        payload = "|".join([str(seed)] + [str(p) for p in parts]).encode("utf-8")
        digest = hashlib.blake2b(payload, digest_size=8).digest()
        return int.from_bytes(digest, "big") / float(1 << 64)

    def should_select(interface_type: str, eid0: int, fid0: int, eid1: int, fid1: int) -> bool:
        if interface_type == "grain_boundary" and scope in ("grain_boundary", "both"):
            return True
        if interface_type == "intragranular" and scope in ("intragranular", "both"):
            if intragranular_fraction >= 1.0:
                return True
            pair_key = (min(eid0, eid1), max(eid0, eid1), min(fid0, fid1), max(fid0, fid1))
            return stable_sample_unit(*pair_key) <= intragranular_fraction
        return False

    for inds in face_groups:
        if len(inds) > 2:
            nonmanifold_faces += 1
            warnings.append(
                f"A face group is shared by {len(inds)} elements; cohesive insertion for this non-manifold face was skipped. "
                "The current release supports conformal manifold meshes where each internal face has at most two owners."
            )
            continue
        for ia, ib in combinations(inds, 2):
            rec0 = face_records[ia]
            rec1 = face_records[ib]
            eid0, fid0, face0, gid0 = rec0
            eid1, fid1, face1, gid1 = rec1
            if eid0 == eid1:
                continue
            pair = tuple(sorted((eid0, eid1)))
            all_internal_pairs.append(pair)
            itype = "grain_boundary" if gid0 != gid1 else "intragranular"
            if should_select(itype, eid0, fid0, eid1, fid1):
                key_pair = (min(eid0, eid1), max(eid0, eid1), min(fid0, fid1), max(fid0, fid1))
                if key_pair in selected_pair_keys:
                    continue
                selected_pair_keys.add(key_pair)
                fam = "GB_COH" if itype == "grain_boundary" else "INTRA_COH"
                selected_interfaces.append(dict(
                    element_1=eid0, face_1=fid0, face_nodes_1=face0, grain_1=gid0,
                    element_2=eid1, face_2=fid1, face_nodes_2=face1, grain_2=gid1,
                    interface_type=itype,
                    family_name=fam,
                    uel_elset=gb_elset if itype == "grain_boundary" else intra_elset,
                    props=gb_props_t if itype == "grain_boundary" else intra_props_t,
                ))

    selected_element_pairs = {tuple(sorted((int(r["element_1"]), int(r["element_2"])))) for r in selected_interfaces}
    mark_step("interface_classification")

    # Build element components connected through unselected internal faces. Nodes are duplicated per component.
    dsu = DSU(element_order)
    for a, b in set(all_internal_pairs):
        if tuple(sorted((a, b))) not in selected_element_pairs:
            dsu.union(a, b)

    conn_dup: Dict[int, List[int]] = {eid: list(nodes) for eid, nodes in solid_conn.items()}
    coords_dup: Dict[int, Tuple[float, ...]] = dict(coords)
    max_node_id = max(coords_dup)
    duplicate_records: List[Dict[str, Any]] = []
    node_occ: Dict[int, List[Tuple[int, int, int, int]]] = {}
    for eid in element_order:
        gid = grain_by_element[eid]
        root = dsu.find(eid)
        for local_idx, nid in enumerate(solid_conn[eid]):
            node_occ.setdefault(nid, []).append((eid, local_idx, gid, root))

    for nid, occs in node_occ.items():
        groups: Dict[int, List[Tuple[int, int, int, int]]] = {}
        for rec in occs:
            groups.setdefault(rec[3], []).append(rec)
        if len(groups) <= 1:
            continue
        roots = sorted(groups)
        for root in roots[1:]:
            max_node_id += 1
            new_nid = max_node_id
            coords_dup[new_nid] = coords[nid]
            gid_set = sorted(set(r[2] for r in groups[root]))
            for eid, local_idx, _gid, _root in groups[root]:
                conn_dup[eid][local_idx] = new_nid
            duplicate_records.append(dict(original_node=nid, new_node=new_nid, component=root, grains=";".join(str(g) for g in gid_set)))

    mark_step("node_duplication")

    family_rows: Dict[str, List[List[int]]] = {gb_elset: [], intra_elset: []}
    family_props: Dict[str, Tuple[float, ...]] = {gb_elset: gb_props_t, intra_elset: intra_props_t}
    family_types: Dict[str, str] = {gb_elset: "grain_boundary", intra_elset: "intragranular"}
    gb_records: List[Dict[str, Any]] = []
    max_elem_id = max(element_order)

    for r in selected_interfaces:
        eid0 = int(r["element_1"]); eid1 = int(r["element_2"])
        fid0 = int(r["face_1"]); fid1 = int(r["face_2"])
        face0_orig = tuple(int(x) for x in r["face_nodes_1"])
        face1_orig = tuple(int(x) for x in r["face_nodes_2"])
        surf0 = surfaces[fid0 - 1]
        surf1 = surfaces[fid1 - 1]
        face0_dup = [conn_dup[eid0][idx] for idx in surf0]
        map_coord_to_dup_1 = {coord_token(coords[orig]): conn_dup[eid1][idx] for orig, idx in zip(face1_orig, surf1)}
        try:
            face1_dup = [map_coord_to_dup_1[coord_token(coords[orig])] for orig in face0_orig]
        except KeyError:
            warnings.append(f"Could not coordinate-match face nodes between elements {eid0} and {eid1}; skipped.")
            continue
        if set(face0_orig) != set(face1_orig):
            warnings.append(f"Interface between elements {eid0} and {eid1} was matched by coordinates, not shared node IDs. Verify this conformal-but-disconnected interface.")
        data = list(face0_dup) + list(face1_dup)
        c0 = centroid(solid_conn[eid0], coords)
        c1 = centroid(solid_conn[eid1], coords)
        ndir = vec_unit(vec_sub(c1, c0))
        nface = face_normal(face0_dup, coords_dup, numdim)
        if vec_dot(ndir, nface) < 0.0:
            data = reorder_if_needed(data, info)
            nface = face_normal(data[:nnps], coords_dup, numdim)
        coh_id = max_elem_id + sum(len(v) for v in family_rows.values()) + 1
        elset = str(r["uel_elset"])
        family_rows.setdefault(elset, []).append([coh_id] + data)
        if elset not in family_props:
            family_props[elset] = tuple(r["props"])
            family_types[elset] = str(r["interface_type"])
        pvec = tuple(float(x) for x in r["props"])
        nrep = list(nface) + [0.0] * max(0, 3 - len(nface))
        gb_records.append(dict(
            cohesive_id=coh_id,
            family_name=r["family_name"],
            interface_type=r["interface_type"],
            uel_type=uel_type,
            uel_elset=elset,
            element_1=eid0,
            face_1=fid0,
            grain_1=r["grain_1"],
            element_2=eid1,
            face_2=fid1,
            grain_2=r["grain_2"],
            measure=interface_measure(face0_orig, coords, numdim, float(pvec[8])),
            normal_x=nrep[0], normal_y=nrep[1], normal_z=nrep[2],
            mode=pvec[0], KI=pvec[1], KII=pvec[2], SI=pvec[3], SII=pvec[4],
            GCI=pvec[5], GCII=pvec[6], ETA=pvec[7], HEIGHT=pvec[8],
        ))

    mark_step("cohesive_connectivity")

    total_coh = sum(len(rows) for rows in family_rows.values())
    log(f"Parsed nodes: {len(coords)}")
    log(f"Parsed solid elements: {len(solid_conn)} of type {elem_type}")
    log(f"Parsed domains/grains: {len(grain_order)}")
    log(f"Selected interface scope: {scope}")
    log(f"Duplicated nodes: {len(duplicate_records)}")
    log(f"Generated cohesive elements: {total_coh}")
    if total_coh == 0:
        warnings.append("No cohesive elements were generated. Check interface scope, domain/grain sets and shared faces.")

    def fmt_float(x: float) -> str:
        return f"{float(x):.10e}"

    def write_id_list_to_lines(ids: Sequence[int], per_line: int = 16) -> List[str]:
        ids = list(ids)
        if not ids:
            return [""]
        return [", ".join(str(int(x)) for x in ids[ii:ii + per_line]) for ii in range(0, len(ids), per_line)]

    def block_data_end(start_idx: int) -> int:
        j = start_idx + 1
        while j < len(lines) and not is_keyword(lines[j]):
            j += 1
        return j

    node_end = block_data_end(node_start)
    elem_end = block_data_end(elem_start)

    duplicate_by_original: Dict[int, List[int]] = {}
    for rec in duplicate_records:
        duplicate_by_original.setdefault(int(rec["original_node"]), []).append(int(rec["new_node"]))

    def parse_nset_ids_from_block(start_idx: int, end_idx: int) -> Tuple[str, List[int], bool, Dict[str, str]]:
        kw = lines[start_idx]
        opts = parse_keyword_options(kw)
        name = opts.get("nset", "")
        generate = "generate" in opts or kw.lower().strip().endswith("generate")
        ids: List[int] = []
        if generate:
            for row in lines[start_idx + 1:end_idx]:
                vals = numeric_ints(row)
                if len(vals) >= 3:
                    ids.extend(expand_generate_ids(vals[0], vals[1], vals[2], f"Nset {name}"))
                    break
        else:
            for row in lines[start_idx + 1:end_idx]:
                ids.extend(numeric_ints(row))
        return name, ids, generate, opts

    def make_nset_supplement_block(name: str, opts: Dict[str, str], add_ids: Sequence[int]) -> List[str]:
        header = f"*Nset, nset={name}"
        if opts.get("instance"):
            header += f", instance={opts['instance']}"
        out = ["** CohesiveX: duplicated nodes added to preserve original node set", header]
        out.extend(write_id_list_to_lines(sorted(set(int(x) for x in add_ids))))
        return out

    def modified_node_lines() -> List[str]:
        out = [lines[node_start]]
        for nid in sorted(coords_dup):
            out.append(str(nid) + ", " + ", ".join(fmt_float(v) for v in coords_dup[nid]))
        return out

    def modified_solid_element_lines() -> List[str]:
        out = [lines[elem_start]]
        for eid in element_order:
            out.append(", ".join(str(int(x)) for x in ([eid] + conn_dup[eid])))
        return out

    def cohesive_header_lines() -> List[str]:
        out: List[str] = []
        out.append("** CohesiveX: user cohesive element definitions inserted below")
        out.append(f"*User element, nodes={nnps * 2}, type={uel_type}, properties=9, coordinates={numdim}, variables={numpt * nsvars_per_ip}")
        out.append("1, 2" if numdim == 2 else "1, 2, 3")
        out.append("**")
        for elset_name, rows in family_rows.items():
            if not rows:
                continue
            out.append(f"** CohesiveX family: {family_types.get(elset_name, 'interface')} -> {elset_name}")
            out.append(f"*Element, type={uel_type}, elset={elset_name}")
            for row in rows:
                out.append(", ".join(str(int(x)) for x in row))
            out.append(f"*Uel Property, elset={elset_name}")
            out.append(", ".join(fmt_float(x) for x in family_props[elset_name]))
            out.append(f"*Elset, elset={elset_name}")
            out.extend(write_id_list_to_lines([row[0] for row in rows]))
            out.append("**")
        out.append("** CohesiveX: end of cohesive element insertion")
        return out

    out_lines: List[str] = [
        "** Generated by CohesiveX Studio in cohesive-only mode",
        f"** Source file: {input_inp.name}",
        f"** Interface scope: {scope}",
        "** Original materials, sections, steps, controls, loads, amplitudes and time increments are preserved.",
        "** Only nodes, parsed solid connectivity, UEL cohesive blocks and duplicated-node set supplements are added/modified.",
        "**",
    ]
    nset_supplement_records: List[Dict[str, Any]] = []
    i = 0
    while i < len(lines):
        if i == node_start:
            out_lines.extend(modified_node_lines())
            i = node_end
            continue
        if i == elem_start:
            out_lines.extend(cohesive_header_lines())
            out_lines.extend(modified_solid_element_lines())
            i = elem_end
            continue
        if supplement_nsets and keyword_name(lines[i]) == "*nset":
            j = block_data_end(i)
            out_lines.extend(lines[i:j])
            nset_name, nset_ids, _generate, opts = parse_nset_ids_from_block(i, j)
            if nset_name:
                add_ids: List[int] = []
                for old_id in nset_ids:
                    add_ids.extend(duplicate_by_original.get(old_id, []))
                add_ids = sorted(set(add_ids))
                if add_ids:
                    out_lines.extend(make_nset_supplement_block(nset_name, opts, add_ids))
                    nset_supplement_records.append({"nset": nset_name, "added_duplicated_nodes": len(add_ids), "instance": opts.get("instance", "")})
            i = j
            continue
        out_lines.append(lines[i])
        i += 1

    output_inp.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    mark_step("write_solver_inp")

    cae_preview_inp = output_inp.with_name(output_inp.stem + "_cae_preview.inp")
    family_elsets_upper = {k.upper() for k, rows in family_rows.items() if rows}
    def write_cae_preview(source: Path, target: Path) -> None:
        src_lines = read_text_auto(source, warnings).splitlines()
        out: List[str] = [
            "** CAE preview file generated by CohesiveX Studio",
            "** This file is for Abaqus/CAE import only.",
            f"** Solver file: {source.name}",
            "** Removed for CAE compatibility: *User element, *Uel Property and inserted cohesive element/elset blocks.",
            "**",
        ]
        i = 0
        while i < len(src_lines):
            line = src_lines[i]
            kname = keyword_name(line)
            opts = parse_keyword_options(line)
            if kname == "*user element":
                out.append("** [CAE preview] skipped *User element block")
                i += 1
                while i < len(src_lines) and not is_keyword(src_lines[i]):
                    i += 1
                continue
            if kname == "*uel property" and opts.get("elset", "").strip().upper() in family_elsets_upper:
                out.append("** [CAE preview] skipped *Uel Property cohesive block")
                i += 1
                while i < len(src_lines) and not is_keyword(src_lines[i]):
                    i += 1
                continue
            if kname == "*element" and opts.get("type", "").strip().upper() == uel_type and opts.get("elset", "").strip().upper() in family_elsets_upper:
                out.append(f"** [CAE preview] skipped *Element, type={uel_type} cohesive block")
                i += 1
                while i < len(src_lines) and not is_keyword(src_lines[i]):
                    i += 1
                continue
            if kname == "*elset" and opts.get("elset", "").strip().upper() in family_elsets_upper:
                out.append("** [CAE preview] skipped cohesive *Elset block")
                i += 1
                while i < len(src_lines) and not is_keyword(src_lines[i]):
                    i += 1
                continue
            out.append(line)
            i += 1
        target.write_text("\n".join(out) + "\n", encoding="utf-8")

    if write_cae_preview:
        write_cae_preview(output_inp, cae_preview_inp)
    mark_step("write_cae_preview")

    gb_csv = report_prefix.with_name(report_prefix.name + "_grain_boundary_table.csv")
    dup_csv = report_prefix.with_name(report_prefix.name + "_duplicated_nodes.csv")
    check_txt = report_prefix.with_name(report_prefix.name + "_mesh_check.txt")
    summary_json = report_prefix.with_name(report_prefix.name + "_summary.json")

    with gb_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = list(gb_records[0].keys()) if gb_records else ["cohesive_id", "family_name", "interface_type", "element_1", "element_2"]
        wr = csv.DictWriter(fh, fieldnames=fieldnames)
        wr.writeheader()
        for rec in gb_records:
            wr.writerow(rec)
    with dup_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = ["original_node", "new_node", "component", "grains"]
        wr = csv.DictWriter(fh, fieldnames=fieldnames)
        wr.writeheader()
        for rec in duplicate_records:
            wr.writerow(rec)

    family_counts = {elset: len(rows) for elset, rows in family_rows.items() if rows}
    interface_type_counts: Dict[str, int] = {}
    for rec in gb_records:
        interface_type_counts[str(rec["interface_type"])] = interface_type_counts.get(str(rec["interface_type"]), 0) + 1

    summary: Dict[str, Any] = dict(
        input_file=str(input_inp), output_file=str(output_inp), cae_preview_file=str(cae_preview_inp),
        element_type=elem_type, dimension=numdim,
        original_nodes=len(coords), duplicated_nodes=len(duplicate_records), total_nodes=len(coords_dup),
        solid_elements=len(solid_conn), cohesive_elements=total_coh,
        interface_scope=scope, intragranular_fraction=intragranular_fraction, random_seed=random_seed_effective,
        interface_type_counts=interface_type_counts, cohesive_family_counts=family_counts,
        grains=len(grain_order), grain_ids=grain_order,
        uel_type=uel_type, nsvars_per_ip=nsvars_per_ip, uel_variable_count=numpt * nsvars_per_ip,
        used_numpy_backend=used_numpy_backend, nonmanifold_faces=nonmanifold_faces,
        nset_supplement_enabled=bool(supplement_nsets), nset_supplemented=len(nset_supplement_records),
        nset_supplements=nset_supplement_records,
        reports=dict(grain_boundary_table=str(gb_csv), duplicated_nodes=str(dup_csv), mesh_check=str(check_txt), summary=str(summary_json), cae_preview=str(cae_preview_inp)),
        backend="numpy" if used_numpy_backend else "pure_python",
        requested_backend="numpy" if fast_mode else "pure_python",
        timings_seconds=dict(step_timings),
        total_preprocessing_time_seconds=0.0,
        warnings=warnings,
    )
    summary["timings_seconds"] = dict(step_timings)
    summary["total_preprocessing_time_seconds"] = perf_counter() - _t_total
    with check_txt.open("w", encoding="utf-8") as fh:
        fh.write("Cohesive generation mesh-check report\n")
        fh.write("====================================\n")
        fh.write(f"Input file             : {input_inp}\n")
        fh.write(f"Solver output file     : {output_inp}\n")
        fh.write(f"CAE preview file       : {cae_preview_inp}\n")
        fh.write(f"Element type           : {elem_type}\n")
        fh.write(f"Dimension              : {numdim}\n")
        fh.write(f"Interface scope        : {scope}\n")
        fh.write(f"Intra fraction         : {intragranular_fraction}\n")
        fh.write(f"Random seed            : {random_seed_effective}\n")
        fh.write(f"NumPy fast backend     : {used_numpy_backend}\n")
        fh.write(f"Requested backend      : {'numpy' if fast_mode else 'pure_python'}\n")
        fh.write(f"Original nodes         : {len(coords)}\n")
        fh.write(f"Duplicated nodes       : {len(duplicate_records)}\n")
        fh.write(f"Total nodes            : {len(coords_dup)}\n")
        fh.write(f"Solid elements         : {len(solid_conn)}\n")
        fh.write(f"Cohesive elements      : {total_coh}\n")
        fh.write(f"Family counts          : {family_counts}\n")
        fh.write(f"Interface type counts  : {interface_type_counts}\n")
        fh.write(f"Domain/grain IDs       : {grain_order}\n")
        fh.write(f"Non-manifold faces     : {nonmanifold_faces}\n")
        fh.write(f"Nset supplements       : {len(nset_supplement_records)} original node sets received duplicated nodes\n")
        fh.write("\nStep timings (seconds)\n")
        for _name, _value in summary.get("timings_seconds", {}).items():
            fh.write(f"  - {_name:<28}: {_value:.6f}\n")
        fh.write(f"  - {'total_preprocessing':<28}: {summary.get('total_preprocessing_time_seconds', 0.0):.6f}\n")
        if nset_supplement_records:
            fh.write("\nSupplemented node sets\n")
            for rec in nset_supplement_records:
                inst = f", instance={rec.get('instance')}" if rec.get("instance") else ""
                fh.write(f"  - {rec.get('nset')}{inst}: +{rec.get('added_duplicated_nodes')} nodes\n")
        fh.write("\nWarnings\n")
        if warnings:
            for w in warnings:
                fh.write(f"  - {w}\n")
        else:
            fh.write("  None\n")
    mark_step("write_reports")
    summary["timings_seconds"] = dict(step_timings)
    summary["total_preprocessing_time_seconds"] = perf_counter() - _t_total
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log(f"Wrote: {output_inp}")
    log(f"Wrote reports with prefix: {report_prefix}")
    return summary


def compare_backends(
    input_inp: str | Path,
    output_prefix: str | Path,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Run the same generation with NumPy and pure-Python backends and compare outputs.

    This function is intended for publishable benchmarks.  It keeps the parser,
    interface-scope rules, cohesive-family definitions and report generation
    unchanged while switching only the interface-detection backend.
    """
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    common = dict(kwargs)
    common.pop("fast_mode", None)
    common.setdefault("verbose", False)

    fast = generate_cohesive_inp(
        input_inp,
        output_prefix.with_name(output_prefix.name + "_numpy.inp"),
        report_prefix=output_prefix.with_name(output_prefix.name + "_numpy"),
        fast_mode=True,
        **common,
    )
    slow = generate_cohesive_inp(
        input_inp,
        output_prefix.with_name(output_prefix.name + "_python.inp"),
        report_prefix=output_prefix.with_name(output_prefix.name + "_python"),
        fast_mode=False,
        **common,
    )

    comparison = {
        "input_file": str(input_inp),
        "numpy_output": fast.get("output_file"),
        "python_output": slow.get("output_file"),
        "topology_consistent": (
            fast.get("cohesive_elements") == slow.get("cohesive_elements")
            and fast.get("duplicated_nodes") == slow.get("duplicated_nodes")
            and fast.get("interface_type_counts") == slow.get("interface_type_counts")
            and fast.get("cohesive_family_counts") == slow.get("cohesive_family_counts")
        ),
        "numpy": fast,
        "pure_python": slow,
        "speedup_total": (
            slow.get("total_preprocessing_time_seconds", 0.0) / fast.get("total_preprocessing_time_seconds", 1.0)
            if fast.get("total_preprocessing_time_seconds", 0.0) > 0 else None
        ),
        "speedup_face_detection": (
            slow.get("timings_seconds", {}).get("face_detection", 0.0) / fast.get("timings_seconds", {}).get("face_detection", 1.0)
            if fast.get("timings_seconds", {}).get("face_detection", 0.0) > 0 else None
        ),
    }
    out_json = output_prefix.with_name(output_prefix.name + "_backend_comparison.json")
    out_txt = output_prefix.with_name(output_prefix.name + "_backend_comparison.txt")
    out_json.write_text(json.dumps(comparison, indent=2), encoding="utf-8")
    with out_txt.open("w", encoding="utf-8") as fh:
        fh.write("CohesiveX backend comparison report\n")
        fh.write("====================================\n")
        fh.write(f"Input file                : {input_inp}\n")
        fh.write(f"Topology consistent       : {comparison['topology_consistent']}\n")
        fh.write(f"NumPy total time          : {fast.get('total_preprocessing_time_seconds', 0.0):.6f} s\n")
        fh.write(f"Pure-Python total time    : {slow.get('total_preprocessing_time_seconds', 0.0):.6f} s\n")
        if comparison["speedup_total"] is not None:
            fh.write(f"Total preprocessing speedup: {comparison['speedup_total']:.3f} x\n")
        if comparison["speedup_face_detection"] is not None:
            fh.write(f"Face-detection speedup    : {comparison['speedup_face_detection']:.3f} x\n")
        fh.write("\nTopology quantities\n")
        for key in ("cohesive_elements", "duplicated_nodes", "interface_type_counts", "cohesive_family_counts"):
            fh.write(f"  {key}: NumPy={fast.get(key)} | Pure-Python={slow.get(key)}\n")
        fh.write("\nNumPy step timings\n")
        for key, value in fast.get("timings_seconds", {}).items():
            fh.write(f"  {key:<28}: {value:.6f} s\n")
        fh.write("\nPure-Python step timings\n")
        for key, value in slow.get("timings_seconds", {}).items():
            fh.write(f"  {key:<28}: {value:.6f} s\n")
    comparison["comparison_json"] = str(out_json)
    comparison["comparison_report"] = str(out_txt)
    return comparison

def run_self_tests() -> None:
    """Run a small regression suite for the cohesive-generation kernel.

    The tests intentionally use tiny Abaqus INP snippets so they can run on a
    normal Python installation without Abaqus.  They check the issues most
    likely to corrupt scientific results: numeric comments inside data blocks,
    fast/slow backend consistency, scope-specific cohesive counts,
    deterministic intragranular sampling, and non-manifold face handling.
    """
    import tempfile

    def write_case(folder: Path, name: str, text: str) -> Path:
        path = folder / name
        path.write_text(text.strip() + "\n", encoding="utf-8")
        return path

    def gen(inp: Path, name: str, **kwargs: Any) -> Dict[str, Any]:
        return generate_cohesive_inp(
            inp,
            inp.with_name(name),
            verbose=False,
            write_cae_preview=True,
            **kwargs,
        )

    two_element_grain_case = """
*Heading
** tiny two-element test with numeric comments inside data blocks
*Node
** seed 0.0 grain 3 size 1.0 -- must be ignored by numeric parser
1, 0.0, 0.0
2, 1.0, 0.0
3, 2.0, 0.0
4, 0.0, 1.0
5, 1.0, 1.0
6, 2.0, 1.0
*Element, type=CPE4
1, 1, 2, 5, 4
** subset 5 created -- must be ignored by connectivity parser
2, 2, 3, 6, 5
*Elset, elset=GRAIN-1
1
*Elset, elset=GRAIN-2
2
*Nset, nset=LEFT
** 123.456 comment should not enter the set
1, 4
*Solid Section, elset=GRAIN-1, material=M1
,
*Solid Section, elset=GRAIN-2, material=M2
,
*Material, name=M1
*Elastic
1., 0.3
*Material, name=M2
*Elastic
1., 0.3
*Step
*Static
0.1, 1.
*End Step
"""

    two_element_one_domain_case = """
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
"""

    three_element_one_domain_case = """
*Heading
*Node
1, 0.0, 0.0
2, 1.0, 0.0
3, 2.0, 0.0
4, 3.0, 0.0
5, 0.0, 1.0
6, 1.0, 1.0
7, 2.0, 1.0
8, 3.0, 1.0
*Element, type=CPE4
1, 1, 2, 6, 5
2, 2, 3, 7, 6
3, 3, 4, 8, 7
*Solid Section, elset=ALL, material=M1
,
*Material, name=M1
*Elastic
1., 0.3
*Step
*Static
0.1, 1.
*End Step
"""

    nonmanifold_case = """
*Heading
*Node
1, 0.0, 0.0
2, 1.0, 0.0
3, 0.0, 1.0
4, 1.0, 1.0
5, 2.0, 0.0
6, 2.0, 1.0
7, 0.0, -1.0
8, 1.0, -1.0
*Element, type=CPE4
1, 1, 2, 4, 3
2, 2, 5, 6, 4
3, 4, 2, 8, 7
*Elset, elset=GRAIN-1
1
*Elset, elset=GRAIN-2
2
*Elset, elset=GRAIN-3
3
*Step
*Static
0.1, 1.
*End Step
"""

    with tempfile.TemporaryDirectory() as td:
        folder = Path(td)
        inp_gb = write_case(folder, "two_gb.inp", two_element_grain_case)
        s_gb = gen(inp_gb, "two_gb_coh.inp", interface_scope="grain_boundary", fast_mode=True)
        assert s_gb["cohesive_elements"] == 1, s_gb
        assert s_gb["interface_type_counts"].get("grain_boundary") == 1, s_gb

        s_intra_empty = gen(inp_gb, "two_intra_empty.inp", interface_scope="intragranular", fast_mode=True)
        assert s_intra_empty["cohesive_elements"] == 0, s_intra_empty

        inp_intra = write_case(folder, "two_domain.inp", two_element_one_domain_case)
        s_fast = gen(inp_intra, "two_domain_fast.inp", interface_scope="intragranular", fast_mode=True)
        s_slow = gen(inp_intra, "two_domain_slow.inp", interface_scope="intragranular", fast_mode=False)
        assert s_fast["cohesive_elements"] == s_slow["cohesive_elements"] == 1, (s_fast, s_slow)
        assert s_fast["duplicated_nodes"] == s_slow["duplicated_nodes"], (s_fast, s_slow)

        inp_three = write_case(folder, "three_domain.inp", three_element_one_domain_case)
        s_sample_fast = gen(inp_three, "sample_fast.inp", interface_scope="intragranular", fast_mode=True, intragranular_fraction=0.5)
        s_sample_slow = gen(inp_three, "sample_slow.inp", interface_scope="intragranular", fast_mode=False, intragranular_fraction=0.5)
        assert s_sample_fast["cohesive_elements"] == s_sample_slow["cohesive_elements"], (s_sample_fast, s_sample_slow)
        assert s_sample_fast["random_seed"] == s_sample_slow["random_seed"] == 0, (s_sample_fast, s_sample_slow)

        inp_non = write_case(folder, "nonmanifold.inp", nonmanifold_case)
        s_non = gen(inp_non, "nonmanifold_coh.inp", interface_scope="both", fast_mode=True)
        assert s_non["nonmanifold_faces"] >= 1, s_non
        assert any("non-manifold" in w.lower() and "skipped" in w.lower() for w in s_non["warnings"]), s_non

    print("CohesiveX Studio self-tests passed.")
