#!/usr/bin/env python3
"""Fail-closed engineering qualification of stored Blender 4.3.2 THICK fibres.

No scene, solver, rendering, or production-array mutation is performed. All
geometry is derived from the admitted float32 keys, promoted to float64. The
continuous geometric inequalities use explicit numerical safety margins; this
is NOT a formal interval-arithmetic certificate or a Cycles ray-precision proof.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
if os.name == 'nt' and __name__ == '__main__' and '--_worker' not in sys.argv:
    from native_cli_watchdog import main as _windows_main
    raise SystemExit(_windows_main(Path(__file__).resolve(),sys.argv[1:]))

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import portable_resources
import signal
import time
import zipfile

import numpy as np
from scipy.spatial import cKDTree

# Historical solver modules use the same bare module name. Load this audited
# standalone kernel under its own name when both systems share one process.
import importlib.util
_spec = importlib.util.spec_from_file_location('_native_clearance_segments',
                                              Path(__file__).with_name('certified_segments.py'))
_kernel = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_kernel)
projection_distance_lower = _kernel.projection_distance_lower

HERE = Path(__file__).resolve().parent
AUDIT = HERE.parents[1]/"carpet-self-clearance-audit"
EPS = np.finfo(np.float64).eps
INTERFIBRE_GAP_M = 1.95e-6


@dataclass(frozen=True)
class GateConfig:
    chord_tolerance_m: float = 2.5e-9
    distance_margin_m: float = 1e-12
    curvature_min_depth: int = 3
    max_depth: int = 16
    max_leaves: int = 1_000_000
    max_native_spans: int = 500_000
    max_input_bytes: int = 64*1024**2
    max_candidate_pairs: int = 30_000_000
    pair_batch_size: int = 8192
    query_batch_size: int = 32
    max_neighbors_per_query: int = 131072
    max_failure_examples: int = 16
    wall_seconds: float = 180.
    memory_limit_mib: int = 900


class Unproven(RuntimeError):
    def __init__(self, code, **details):
        super().__init__(code)
        self.code, self.details = code, details


class Budget:
    def __init__(self, config):
        self.config, self.start = config, time.monotonic()

    def check(self):
        elapsed = time.monotonic()-self.start
        if elapsed > self.config.wall_seconds:
            raise Unproven("wall_budget_exceeded", elapsed_seconds=elapsed)
        # portable_resources normalizes every platform to KiB. A high-water mark,
        # conservative when imported into a process that already did other work.
        peak = portable_resources.getrusage(portable_resources.RUSAGE_SELF).ru_maxrss/1024
        if peak > self.config.memory_limit_mib:
            raise Unproven("memory_budget_exceeded", process_peak_rss_mib=peak)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1<<20), b""):
            digest.update(block)
    return digest.hexdigest()


def array_receipt(value):
    array = np.asarray(value)
    return {"shape": list(array.shape), "dtype": array.dtype.str,
            "sha256_c_order_bytes": hashlib.sha256(array.tobytes(order="C")).hexdigest()}


def native_beziers(keys):
    """Repeated-endpoint uniform Catmull–Rom, one cubic per adjacent key pair."""
    p = np.asarray(keys, dtype=np.float64)
    before = np.concatenate((p[:1], p[:-1]))
    after = np.concatenate((p[1:], p[-1:]))
    tangent = .5*(after-before)
    return np.stack((p[:-1], p[:-1]+tangent[:-1]/3,
                     p[1:]-tangent[1:]/3, p[1:]), axis=1)


def split_beziers(bez):
    first = (bez[:, :-1]+bez[:, 1:])*.5
    second = (first[:, :-1]+first[:, 1:])*.5
    middle = (second[:, 0]+second[:, 1])*.5
    left = np.stack((bez[:, 0], first[:, 0], second[:, 0], middle), axis=1)
    right = np.stack((middle, second[:, 1], first[:, 2], bez[:, 3]), axis=1)
    return np.stack((left, right), axis=1).reshape(-1, 4, 3)


def continuous_bounds(bez, position_pad=0.):
    """Bernstein positive-speed lower and curvature upper bounds on each cubic.

    The quadratic derivative hull projects onto its midpoint direction. Norms
    of the degree-three Bernstein coefficients bound gamma' cross gamma''.
    Derivatives here are with respect to the leaf's own [0,1] parameter.
    """
    derivative = 3*np.diff(bez, axis=1)
    second = 2*np.diff(derivative, axis=1)
    middle = (derivative[:, 0]+2*derivative[:, 1]+derivative[:, 2])*.25
    middle_norm = np.linalg.norm(middle, axis=1)
    direction = np.divide(middle, middle_norm[:, None], out=np.zeros_like(middle),
                          where=middle_norm[:, None]>0)
    direction *= 1-32*EPS
    dscale = np.max(np.linalg.norm(derivative, axis=2), axis=1)
    sscale = np.max(np.linalg.norm(second, axis=2), axis=1)
    derivative_pad = 12*position_pad+128*EPS*dscale
    speed = np.min(np.sum(derivative*direction[:, None], axis=2), axis=1)-derivative_pad
    numerator = np.zeros(len(bez))
    comb2, comb3 = (1, 2, 1), (1, 3, 3, 1)
    for k in range(4):
        coefficient = np.zeros((len(bez), 3))
        for i in range(3):
            j = k-i
            if 0 <= j < 2:
                coefficient += (comb2[i]/comb3[k])*np.cross(derivative[:, i], second[:, j])
        numerator = np.maximum(numerator, np.linalg.norm(coefficient, axis=1))
    # Bounds include subtraction/cross-product roundoff and accumulated control
    # coordinate uncertainty. Deliberately conservative, not interval arithmetic.
    numerator += 256*EPS*dscale*sscale + derivative_pad*(sscale+4*dscale+4*derivative_pad)
    curvature = np.divide(numerator, speed**3, out=np.full(len(bez), np.inf), where=speed>0)
    curvature *= 1+256*EPS
    return speed, curvature


def chord_arc_bounds(bez, position_pad=0.):
    chord = bez[:, 3]-bez[:, 0]
    parameter_error = np.maximum(np.linalg.norm(bez[:, 1]-bez[:, 0]-chord/3, axis=1),
                                 np.linalg.norm(bez[:, 2]-bez[:, 0]-2*chord/3, axis=1))
    # The difference from the affine chord is itself a cubic Bezier curve with
    # zero endpoint controls. Convexity gives a continuous parameterwise bound.
    # Distance to a convex segment is convex too. This second valid bound avoids
    # pointless refinement of a straight span with nonuniform parameter speed.
    relative = bez[:, 1:3]-bez[:, :1]
    chord_square = np.sum(chord*chord, axis=1)
    parameter = np.divide(np.sum(relative*chord[:, None], axis=2), chord_square[:, None],
                          out=np.zeros((len(bez), 2)), where=chord_square[:, None]>0)
    projected = np.clip(parameter, 0., 1.)[:, :, None]*chord[:, None]
    hull_error = np.linalg.norm(relative-projected, axis=2).max(axis=1)
    error = np.minimum(parameter_error, hull_error)*(1+128*EPS)+8*position_pad
    upper = np.sum(np.linalg.norm(np.diff(bez, axis=1), axis=2), axis=1)
    upper = upper*(1+128*EPS)+8*position_pad
    return error, upper


def owner_curvature(bez, radius, config, budget, position_pad, owner):
    pending = bez
    native_span = np.arange(len(bez), dtype=np.int32)
    maximum, speed_lower, accepted = 0., math.inf, 0
    deepest = 0
    for depth in range(config.max_depth+1):
        budget.check()
        if len(pending) > config.max_leaves:
            raise Unproven("curvature_work_budget_exceeded", owner=owner, pending=len(pending))
        speed, curvature = continuous_bounds(pending, position_pad)
        good = np.isfinite(curvature) & np.isfinite(speed) & (speed>0) & (curvature*radius<1)
        if depth < config.curvature_min_depth:
            good[:] = False
        if np.any(good):
            maximum = max(maximum, float(curvature[good].max()))
            speed_lower = min(speed_lower, float((speed[good]*2**depth).min()))
            accepted += int(good.sum())
            deepest = max(deepest, depth)
        if np.all(good):
            break
        if depth == config.max_depth:
            bad = int(np.flatnonzero(~good)[0])
            raise Unproven("curvature_or_speed_unresolved", owner=owner,
                           native_span=int(native_span[bad]), depth=depth,
                           speed_lower=float(speed[bad]),
                           kappa_radius_upper=float(curvature[bad]*radius)
                           if np.isfinite(curvature[bad]) else None)
        if 2*int((~good).sum()) > config.max_leaves:
            raise Unproven("curvature_work_budget_exceeded", owner=owner,
                           prospective_pending=2*int((~good).sum()))
        pending = split_beziers(pending[~good])
        native_span = np.repeat(native_span[~good], 2)
    cutoff = math.inf if maximum == 0 else np.pi/maximum*(1-256*EPS)
    return {"curvature_upper_per_m": maximum, "kappa_radius_upper": maximum*radius,
            "native_parameter_speed_lower_m": speed_lower,
            "local_arc_cutoff_m": cutoff if np.isfinite(cutoff) else None,
            "curvature_bound_leaves": accepted, "maximum_curvature_depth": deepest}, cutoff


def owner_leaves(bez, owner, cutoff, config, budget, position_pad, remaining):
    pending = bez
    spans = np.arange(len(bez), dtype=np.int32)
    starts, ends = np.zeros(len(bez)), np.ones(len(bez))
    accepted = []
    count = 0
    for depth in range(config.max_depth+1):
        budget.check()
        if count+len(pending) > remaining:
            raise Unproven("leaf_budget_exceeded", owner=owner, accepted=count,
                           pending=len(pending), remaining=remaining)
        error, arc = chord_arc_bounds(pending, position_pad)
        good = ((error <= config.chord_tolerance_m) & (arc <= cutoff*.5*(1-256*EPS)) &
                np.isfinite(arc) & (arc>0) & np.isfinite(error))
        if np.any(good):
            accepted.append((pending[good, 0], pending[good, 3], error[good], arc[good],
                             spans[good], starts[good], ends[good]))
            count += int(good.sum())
        if np.all(good):
            break
        if depth == config.max_depth:
            raise Unproven("chord_or_arc_bound_unresolved", owner=owner, depth=depth,
                           native_span=int(spans[np.flatnonzero(~good)[0]]))
        if count+2*int((~good).sum()) > remaining:
            raise Unproven("leaf_budget_exceeded", owner=owner, accepted=count,
                           prospective_pending=2*int((~good).sum()), remaining=remaining)
        pending = split_beziers(pending[~good])
        spans = np.repeat(spans[~good], 2)
        lo, hi = starts[~good], ends[~good]
        mid = (lo+hi)*.5
        starts = np.stack((lo, mid), axis=1).ravel()
        ends = np.stack((mid, hi), axis=1).ravel()
    arrays = [np.concatenate([record[k] for record in accepted]) for k in range(7)]
    order = np.lexsort((arrays[5], arrays[4]))
    arrays = [value[order] for value in arrays]
    a, b, error, arc, span, t0, t1 = arrays
    # Ordering and native span coverage are admission invariants, not an
    # incidental property of breadth-first adaptive subdivision.
    if not (span[0] == 0 and t0[0] == 0 and span[-1] == len(bez)-1 and t1[-1] == 1):
        raise Unproven("leaf_order_or_coverage_error", owner=owner)
    adjacent = ((span[1:]==span[:-1]) & (t0[1:]==t1[:-1])) | (
                (span[1:]==span[:-1]+1) & (t0[1:]==0) & (t1[:-1]==1))
    if not np.all(adjacent):
        raise Unproven("leaf_order_or_coverage_error", owner=owner)
    # Prefix differences are padded for both accumulation and subtraction. A
    # plain subtraction of two rounded-up cumulative sums would not suffice.
    prefix = np.r_[0., np.cumsum(arc)]
    accumulation_pad = 256*EPS*float(prefix[-1])*max(1, len(arc)) + 8*position_pad
    return {"a": a, "b": b, "error": error, "arc_upper": arc,
            "owner": np.full(len(a), owner, dtype=np.int32), "native_span": span,
            "t0": t0, "t1": t1, "arc_before": prefix[:-1], "arc_after": prefix[1:]}, accumulation_pad


def _validate(body, wraps, body_radius, wrap_radius, identity_world, config):
    integer_fields = ("curvature_min_depth", "max_depth", "max_leaves", "max_native_spans",
                      "max_input_bytes", "max_candidate_pairs", "pair_batch_size", "query_batch_size",
                      "max_neighbors_per_query", "max_failure_examples", "memory_limit_mib")
    if any(type(getattr(config, name)) is not int for name in integer_fields):
        raise Unproven("invalid_configuration")
    if identity_world is not True:
        raise Unproven("identity_world_contract_not_asserted")
    for name, data, owners in (("body", body, 187), ("wraps", wraps, 6)):
        if not isinstance(data, np.ndarray) or data.dtype != np.dtype("float32"):
            raise Unproven("stored_float32_required", field=name)
        if data.ndim != 3 or data.shape[0] != owners or data.shape[2] != 3 or data.shape[1] < 2:
            raise Unproven("invalid_curve_shape", field=name, actual_shape=list(data.shape))
        if not np.all(np.isfinite(data)):
            raise Unproven("nonfinite_input", field=name)
        if float(np.max(np.abs(data))) > 1e3:
            raise Unproven("coordinate_range_outside_engineering_contract", field=name)
    for name, value in (("body_radius", body_radius), ("wrap_radius", wrap_radius)):
        data = np.asarray(value)
        if data.dtype != np.dtype("float32") or data.size != 1:
            raise Unproven("constant_float32_radius_required", field=name)
        if not np.isfinite(data).all() or not 1e-12 < float(data.reshape(-1)[0]) <= 1:
            raise Unproven("radius_outside_engineering_contract", field=name)
    native_spans = len(body)*(body.shape[1]-1)+len(wraps)*(wraps.shape[1]-1)
    if native_spans > config.max_native_spans:
        raise Unproven("native_span_budget_exceeded", native_spans=native_spans)
    if body.nbytes+wraps.nbytes > config.max_input_bytes:
        raise Unproven("input_byte_budget_exceeded")
    if not (0 < config.chord_tolerance_m < 1e-3 and 0 < config.distance_margin_m < 1e-3
            and np.isfinite(config.wall_seconds) and config.wall_seconds>0
            and 0 <= config.curvature_min_depth <= config.max_depth <= 24
            and 0 < config.memory_limit_mib <= 1024
            and 0 < config.max_native_spans <= 500_000
            and 0 < config.max_input_bytes <= 64*1024**2
            and 0 < config.max_candidate_pairs <= 100_000_000
            and 0 < config.max_failure_examples <= 128
            and 0 < config.max_leaves <= 1_000_000
            and 0 < config.pair_batch_size <= 65536
            and 0 < config.query_batch_size <= 128
            and 0 < config.max_neighbors_per_query <= 131072):
        raise Unproven("invalid_configuration")


def _leaf_label(leaves, index):
    owner = int(leaves["owner"][index])
    return {"owner": owner, "kind": "body" if owner<187 else "wrap",
            "kind_index": owner if owner<187 else owner-187,
            "native_span": int(leaves["native_span"][index]),
            "parameter_start": float(leaves["t0"][index]),
            "parameter_end": float(leaves["t1"][index]), "ordered_leaf_index": int(index)}


def _contact(leaves, radii, cutoffs, arc_pads, config, budget):
    a, b, error, owner = (leaves[k] for k in ("a", "b", "error", "owner"))
    radius = radii[owner]
    mid = (a+b)*.5
    half = np.linalg.norm(b-a, axis=1)*.5
    extent = radius+error+half
    # Both chord error and physical radius participate in the sphere broad phase.
    numerical_pad = config.distance_margin_m+512*EPS*max(1e-6, float(np.abs(mid).max()))
    maximum_extent = float(extent.max())
    tree = cKDTree(mid)
    budget.check()
    failures, construction_failures = [], []
    count = {"sphere_candidate_pairs": 0, "aabb_pruned_pairs": 0,
             "same_owner_local_pairs": 0, "interfibre_pairs_tested": 0,
             "same_owner_nonlocal_pairs_tested": 0, "unresolved_pairs": 0,
             "interfibre_physical_unresolved_pairs": 0,
             "same_fibre_self_unresolved_pairs": 0,
             "construction_gap_unresolved_pairs": 0}
    minimum_inter, minimum_self = math.inf, math.inf
    for first in range(0, len(a), config.query_batch_size):
        budget.check()
        last = min(len(a), first+config.query_batch_size)
        horizon = extent[first:last]+maximum_extent+INTERFIBRE_GAP_M+numerical_pad
        # Count before materializing neighbor lists; pathological coincident
        # geometry must not allocate an unbounded quadratic pair table.
        sizes = tree.query_ball_point(mid[first:last], horizon, eps=0,
                                      workers=1, return_length=True)
        if int(np.max(sizes)) > config.max_neighbors_per_query:
            raise Unproven("neighbor_budget_exceeded", first_leaf=first,
                           maximum_neighbors=int(np.max(sizes)))
        for index in range(first, last):
            neighbors = tree.query_ball_point(mid[index], horizon[index-first], eps=0, workers=1)
            js = np.asarray(neighbors, dtype=np.int64)
            js = js[js>index]
            count["sphere_candidate_pairs"] += len(js)
            if count["sphere_candidate_pairs"] > config.max_candidate_pairs:
                raise Unproven("candidate_pair_budget_exceeded", **count)
            for begin in range(0, len(js), config.pair_batch_size):
                budget.check()
                j = js[begin:begin+config.pair_batch_size]
                same = owner[j] == owner[index]
                # The inclusive span covers both entire endpoint leaves. For an
                # open curve it upper-bounds every point pair in the rectangle.
                span_upper = leaves["arc_after"][j]-leaves["arc_before"][index]+arc_pads[owner[index]]
                local = same & (span_upper <= cutoffs[owner[index]])
                count["same_owner_local_pairs"] += int(local.sum())
                j, same = j[~local], same[~local]
                if not len(j):
                    continue
                threshold = np.where(same, 0., INTERFIBRE_GAP_M)
                required = radius[index]+radius[j]+error[index]+error[j]+threshold+numerical_pad
                # Distance between AABBs of the two full chords is another
                # continuous lower bound; positive axis separation is sufficient.
                separation = np.maximum(np.maximum(np.minimum(a[j], b[j])-np.maximum(a[index], b[index]),
                                                     np.minimum(a[index], b[index])-np.maximum(a[j], b[j])), 0.)
                far = np.max(separation, axis=1) > required
                count["aabb_pruned_pairs"] += int(far.sum())
                j, same, threshold = j[~far], same[~far], threshold[~far]
                if not len(j):
                    continue
                lower = projection_distance_lower(np.broadcast_to(a[index], (len(j), 3)),
                                                  np.broadcast_to(b[index], (len(j), 3)), a[j], b[j])
                surface = lower-error[index]-error[j]-radius[index]-radius[j]-numerical_pad
                if not np.all(np.isfinite(surface)):
                    raise Unproven("nonfinite_distance_bound", leaf=index)
                count["interfibre_pairs_tested"] += int((~same).sum())
                count["same_owner_nonlocal_pairs_tested"] += int(same.sum())
                if np.any(~same):
                    minimum_inter = min(minimum_inter, float(surface[~same].min()))
                if np.any(same):
                    minimum_self = min(minimum_self, float(surface[same].min()))
                # Prototype qualification requires physical separation STRICTLY
                # >0 for both same and different fibres. The construction gap
                # is independently reported and never silently called preserved.
                bad = surface<=0.
                construction_bad = (~same) & (surface<INTERFIBRE_GAP_M)
                count["unresolved_pairs"] += int(bad.sum())
                count["interfibre_physical_unresolved_pairs"] += int((bad & ~same).sum())
                count["same_fibre_self_unresolved_pairs"] += int((bad & same).sum())
                count["construction_gap_unresolved_pairs"] += int(construction_bad.sum())
                for offset in np.flatnonzero(bad):
                    if len(failures) < config.max_failure_examples:
                        failures.append({"code": "nonlocal_self_clearance_unproven" if same[offset]
                                         else "interfibre_clearance_unproven",
                                         "first": _leaf_label(leaves, index),
                                         "second": _leaf_label(leaves, int(j[offset])),
                                         "surface_clearance_lower_m": float(surface[offset]),
                                         "required_surface_clearance_m": 0.,
                                         "strict_inequality_required": True})
                for offset in np.flatnonzero(construction_bad):
                    if len(construction_failures) < config.max_failure_examples:
                        construction_failures.append({"code": "construction_gap_unproven",
                                         "first": _leaf_label(leaves, index),
                                         "second": _leaf_label(leaves, int(j[offset])),
                                         "surface_clearance_lower_m": float(surface[offset]),
                                         "required_surface_clearance_m": INTERFIBRE_GAP_M})
    return {**count, "interfibre_surface_clearance_lower_m": min(INTERFIBRE_GAP_M, minimum_inter),
            "minimum_tested_interfibre_surface_clearance_lower_m": minimum_inter if np.isfinite(minimum_inter) else None,
            "minimum_tested_nonlocal_self_surface_clearance_lower_m": minimum_self if np.isfinite(minimum_self) else None,
            "complete_radius_sphere_and_aabb_broad_phase": True,
            "same_owner_pairs_dropped_by_owner_only": 0,
            "interfibre_physical_nonpenetration_pass": count["interfibre_physical_unresolved_pairs"] == 0,
            "same_fibre_self_clearance_pass": count["same_fibre_self_unresolved_pairs"] == 0,
            "construction_gap_preserved_pass": count["construction_gap_unresolved_pairs"] == 0,
            "failure_examples": failures, "construction_gap_failure_examples": construction_failures}


def qualify_native(body, wraps, body_radius, wrap_radius, *, identity_world=False, config=None):
    """Qualify one actual per-loop stored float32 array bundle; never mutate it.

    identity_world=True asserts that these are the authoritative stored keys in
    metres with identity object transforms, static THICK representation, constant
    stated radii, and no modifiers/deformation. A later scene builder must verify
    that export assertion against its actual objects before relying on this gate.
    This function does not inspect or qualify an existing Blender scene.
    """
    config = GateConfig() if config is None else config
    budget = Budget(config)
    result = {"qualified": False, "complete": False, "status": "unstarted",
              "schema": "native_curve_qualification_v1", "configuration": asdict(config),
              "formal_outward_rounded_arithmetic": False,
              "scope": "Continuous geometric engineering bounds; no Blender scene or Cycles ray precision qualification.",
              "representation": {"blender_version": "4.3.2", "cycles_shape": "THICK",
                 "basis": "uniform Catmull-Rom", "endpoint_keys": "repeated",
                 "curves": "open; no implicit cyclic seam", "units": "metres",
                 "identity_world_asserted": identity_world is True,
                 "actual_scene_contract_verified": False},
              "interfibre_construction_gap_m": INTERFIBRE_GAP_M,
              "interfibre_physical_target_m": 0., "same_fibre_physical_target_m": 0.,
              "qualification_requires_construction_gap": False,
              "physical_nonpenetration_pass": False, "construction_gap_preserved_pass": False,
              "strict_nonfold_pass": False, "same_fibre_self_clearance_pass": False,
              "same_fibre_nonlocal_test_is_strict": True, "failures": []}
    try:
        _validate(body, wraps, body_radius, wrap_radius, identity_world, config)
        budget.check()
        result["inputs"] = {name: array_receipt(value) for name, value in
                            (("body", body), ("wraps", wraps), ("body_radius", body_radius), ("wrap_radius", wrap_radius))}
        radius_values = (float(np.asarray(body_radius).reshape(-1)[0]), float(np.asarray(wrap_radius).reshape(-1)[0]))
        radii = np.r_[np.full(187, radius_values[0]), np.full(6, radius_values[1])]
        result["constant_radii_m"] = {"body": radius_values[0], "wraps": radius_values[1]}
        owner_receipts, leaf_groups, arc_pads, cutoffs = [], [], [], []
        leaf_count = 0
        for owner in range(193):
            budget.check()
            keys = body[owner] if owner<187 else wraps[owner-187]
            bez = native_beziers(keys)
            position_pad = 256*EPS*max(1e-6, float(np.abs(keys).max()))*(config.max_depth+3)
            curvature_receipt, cutoff = owner_curvature(bez, radii[owner], config, budget, position_pad, owner)
            leaves, arc_pad = owner_leaves(bez, owner, cutoff, config, budget,
                                            position_pad, config.max_leaves-leaf_count)
            leaf_count += len(leaves["a"])
            leaf_groups.append(leaves)
            cutoffs.append(cutoff)
            arc_pads.append(arc_pad)
            owner_receipts.append({"owner": owner, "kind": "body" if owner<187 else "wrap",
                                   "native_spans": len(bez), "leaf_count": len(leaves["a"]),
                                   "arclength_upper_m": float(leaves["arc_after"][-1]+arc_pad),
                                   "maximum_leaf_arc_upper_m": float(leaves["arc_upper"].max()),
                                   **curvature_receipt})
        result["owners"] = owner_receipts
        result["leaf_count"] = leaf_count
        result["maximum_kappa_radius_upper"] = max(x["kappa_radius_upper"] for x in owner_receipts)
        result["strict_nonfold_pass"] = True
        result["minimum_native_parameter_speed_lower_m"] = min(x["native_parameter_speed_lower_m"] for x in owner_receipts)
        leaves = {key: np.concatenate([part[key] for part in leaf_groups]) for key in leaf_groups[0]}
        del leaf_groups
        budget.check()
        result["maximum_continuous_chord_error_m"] = float(leaves["error"].max())
        result["leaf_order"] = "owner, native_span, parameter_start; continuous complete coverage checked"
        result["leaf_data_hashes"] = {key: array_receipt(value) for key, value in leaves.items()}
        contacts = _contact(leaves, radii, np.asarray(cutoffs), np.asarray(arc_pads), config, budget)
        result["contacts"] = contacts
        result["failures"] = contacts["failure_examples"]
        result["qualified"] = contacts["unresolved_pairs"] == 0
        result["physical_nonpenetration_pass"] = contacts["interfibre_physical_nonpenetration_pass"]
        result["construction_gap_preserved_pass"] = contacts["construction_gap_preserved_pass"]
        result["same_fibre_self_clearance_pass"] = contacts["same_fibre_self_clearance_pass"]
        result["complete"] = True
        result["status"] = "sufficient_condition_pass" if result["qualified"] else "unproven_clearance"
    except Unproven as exc:
        result["status"] = "unproven_"+exc.code
        result["failures"].append({"code": exc.code, **exc.details})
    except (MemoryError, FloatingPointError, OverflowError, ValueError) as exc:
        result["status"] = "unproven_execution_error"
        result["failures"].append({"code": type(exc).__name__, "message": str(exc)})
    result["elapsed_seconds"] = time.monotonic()-budget.start
    result["process_peak_rss_mib"] = portable_resources.getrusage(portable_resources.RUSAGE_SELF).ru_maxrss/1024
    result["rss_accounting"] = portable_resources.accounting_receipt()
    result["source_hashes"] = {p.name: sha256_file(p) for p in (HERE/"native_curve_gate.py", HERE/"certified_segments.py")}
    result["audit_hashes"] = {name: sha256_file(AUDIT/name) for name in
                             ("SELF_CLEARANCE_AUDIT.md", "self_clearance_reference.py",
                              "cycles_representation_receipt.json", "cycles_basis_algebra_receipt.json") if (AUDIT/name).is_file()}
    return result


def qualify_npz(input_path, *, identity_world=False, config=None):
    config = GateConfig() if config is None else config
    path = Path(input_path)
    try:
        # Check total uncompressed size before np.load to bound ZIP expansion.
        with zipfile.ZipFile(path) as archive:
            if sum(info.file_size for info in archive.infolist()) > config.max_input_bytes:
                raise Unproven("npz_uncompressed_byte_budget_exceeded")
        hash_before = sha256_file(path)
        with np.load(path, allow_pickle=False) as values:
            required = ("body", "wraps", "body_radius", "wrap_radius")
            missing = [name for name in required if name not in values]
            if missing:
                raise Unproven("missing_input_keys", keys=missing)
            result = qualify_native(*(values[name] for name in required), identity_world=identity_world, config=config)
        result["input_npz_sha256"] = sha256_file(path)
        if result["input_npz_sha256"] != hash_before:
            result["qualified"] = False
            result["complete"] = False
            result["status"] = "unproven_input_file_changed_during_qualification"
            result["failures"].append({"code": "input_file_changed_during_qualification"})
        result["input_filename"] = path.name
        return result
    except (Unproven, OSError, ValueError, zipfile.BadZipFile) as exc:
        return {"qualified": False, "complete": False, "status": "unproven_input_file",
                "failures": [{"code": getattr(exc, "code", type(exc).__name__),
                              "message": str(exc), **getattr(exc, "details", {})}]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_npz", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--identity-world", action="store_true", help="Assert documented static identity-world native input contract")
    parser.add_argument("--wall-seconds", type=float, default=180.)
    parser.add_argument("--chord-tolerance-nm", type=float, default=2.5)
    parser.add_argument("--max-leaves", type=int, default=1_000_000)
    parser.add_argument("--_worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    config = GateConfig(wall_seconds=args.wall_seconds, chord_tolerance_m=args.chord_tolerance_nm*1e-9,
                        max_leaves=args.max_leaves)
    def wall_alarm(signum, frame):
        raise Unproven("cli_hard_wall_budget_exceeded")
    if os.name != 'nt':
        signal.signal(signal.SIGALRM, wall_alarm)
        signal.setitimer(signal.ITIMER_REAL, max(.1, config.wall_seconds))
    try:
        result = qualify_npz(args.input_npz, identity_world=args.identity_world, config=config)
    except Unproven as exc:
        result = {"qualified": False, "complete": False, "status": "unproven_"+exc.code,
                  "failures": [{"code": exc.code}]}
    finally:
        if os.name != 'nt':
            signal.setitimer(signal.ITIMER_REAL, 0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"qualified": result["qualified"], "status": result["status"], "receipt": str(args.output)}))
    return 0 if result["qualified"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
