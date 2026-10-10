#!/usr/bin/env python3
"""Small analytic/regression tests, never a carpet material geometry fixture."""
from dataclasses import replace
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from native_curve_gate import (Budget, GateConfig, Unproven, chord_arc_bounds,
    continuous_bounds, native_beziers, owner_curvature, owner_leaves,
    qualify_native, sha256_file, split_beziers)
from certified_segments import projection_distance_lower

HERE = Path(__file__).resolve().parent


def line_bundle(n=5, m=3):
    """Temporary source-test inputs only. No fixture NPZ is delivered."""
    body, wraps = np.empty((187, n, 3), "f4"), np.empty((6, m, 3), "f4")
    for i, keys in enumerate(body):
        keys[:] = np.c_[np.full(n, (i%14)*.0001), np.full(n, (i//14)*.0001), np.linspace(0, .004, n)]
    for i, keys in enumerate(wraps):
        keys[:] = np.c_[np.full(m, .002+i*.0001), np.zeros(m), np.linspace(0, .004, m)]
    return body, wraps, np.float32(15e-6), np.float32(16e-6)


def summary(result):
    return {key: result[key] for key in ("qualified", "complete", "status", "elapsed_seconds",
            "process_peak_rss_mib", "physical_nonpenetration_pass", "construction_gap_preserved_pass")
            if key in result} | ({"contacts": result["contacts"]} if "contacts" in result else {"failures": result["failures"]})


def run():
    results = {}
    config = GateConfig(wall_seconds=60.)
    print(json.dumps({'stage':'analytic_tests_started'}),flush=True)
    bundle = line_bundle()
    before = [np.asarray(x).copy() for x in bundle]
    result = qualify_native(*bundle, identity_world=True, config=config)
    assert result["qualified"] and result["construction_gap_preserved_pass"]
    assert result["contacts"]["same_owner_local_pairs"]>0
    assert result["contacts"]["same_owner_pairs_dropped_by_owner_only"] == 0
    assert all(np.array_equal(x, y) for x, y in zip(bundle, before))
    results["193_open_straight_fibres_input_immutability"] = summary(result)

    # A native span's endpoint tangents are exactly the audited half differences.
    keys = np.array([[0., 0, 0], [1., .3, 0], [2., -.1, .4], [3., 1, 0]], "f4")
    bez = native_beziers(keys)
    assert np.allclose(3*(bez[0, 1]-bez[0, 0]), (keys[1].astype(float)-keys[0])*.5, atol=1e-15)
    assert np.allclose(3*(bez[-1, 3]-bez[-1, 2]), (keys[-1].astype(float)-keys[-2])*.5, atol=1e-15)
    assert np.allclose(3*(bez[:-1, 3]-bez[:-1, 2]), 3*(bez[1:, 1]-bez[1:, 0]), atol=2e-15)
    results["repeated_endpoints_and_C1_native_joins"] = {"passed": True}

    # gamma(t)=(t,t^2,0): speed=sqrt(1+4t^2), max curvature=2,
    # length=sqrt(5)/2+asinh(2)/4; known exact continuous quantities.
    parabola = np.array([[[0., 0, 0], [1/3, 0, 0], [2/3, 1/3, 0], [1., 1., 0]]])
    speed, curvature = continuous_bounds(parabola)
    error, arc = chord_arc_bounds(parabola)
    exact_length = np.sqrt(5)/2+np.arcsinh(2)/4
    assert 0 < speed[0] <= 1 and curvature[0] >= 2
    assert arc[0] >= exact_length and error[0] >= 1/(4*np.sqrt(2))
    restricted = parabola.copy()
    for _ in range(4):
        restricted = split_beziers(restricted)
    fine_error, fine_arc = chord_arc_bounds(restricted)
    assert exact_length <= fine_arc.sum() < arc[0]
    assert fine_error.max()<error[0]
    short_config = replace(config, chord_tolerance_m=1e-5)
    leaves, _ = owner_leaves(np.concatenate((parabola, parabola+np.array([1., 1., 0.]))),
                              0, .4, short_config, Budget(short_config), 0., 10000)
    assert np.all(leaves["arc_upper"] <= .2)
    assert np.array_equal(np.lexsort((leaves["t0"], leaves["native_span"])), np.arange(len(leaves["a"])))
    results["continuous_parabola_bounds_and_ordered_arc_refinement"] = {
        "passed": True, "speed_lower": float(speed[0]), "curvature_upper": float(curvature[0]),
        "analytic_maximum_curvature": 2., "control_polygon_arc_upper": float(arc[0]),
        "analytic_arc_length": float(exact_length), "refined_arc_upper": float(fine_arc.sum()),
        "leaf_count": len(leaves["a"])}

    # This cubic has gamma'(1/2)=0. Refinement cannot legalize zero speed.
    cusp = np.array([[[0., 0, 0], [1/3, 0, 0], [0., 0, 0], [1/3, 0, 0]]])
    try:
        owner_curvature(cusp, 1e-5, replace(config, max_depth=8), Budget(config), 0., 0)
        raise AssertionError("zero speed was admitted")
    except Unproven as exc:
        assert exc.code == "curvature_or_speed_unresolved"
    results["zero_speed_fails_closed"] = {"passed": True}

    # Parabola curvature attains 2 at its endpoint: r=1/2 is focal equality,
    # which must fail the strict K*r<1 condition, not become a tolerance pass.
    try:
        owner_curvature(parabola, .5, replace(config, max_depth=8), Budget(config), 0., 0)
        raise AssertionError("focal equality was admitted")
    except Unproven as exc:
        assert exc.code == "curvature_or_speed_unresolved"
    results["strict_focal_equality_fails_closed"] = {"passed": True, "analytic_maximum_kappa_radius": 1.}

    margin_bundle = line_bundle()
    body, wraps, br, wr = margin_bundle
    body[1] = body[0]
    body[1, :, 0] = np.float32(2*float(br)+.5e-6)
    result = qualify_native(*margin_bundle, identity_world=True, config=config)
    assert result["qualified"] and result["physical_nonpenetration_pass"]
    assert not result["construction_gap_preserved_pass"]
    assert 0 < result["contacts"]["interfibre_surface_clearance_lower_m"] < 1.95e-6
    results["physical_clearance_pass_construction_margin_fail"] = summary(result)

    # Exactly tangent stored float32 straight curves: radius*2 is representable.
    body[1, :, 0] = np.float32(2*br)
    result = qualify_native(*margin_bundle, identity_world=True, config=config)
    assert not result["qualified"] and not result["physical_nonpenetration_pass"]
    assert result["contacts"]["interfibre_physical_unresolved_pairs"]>0
    results["strict_physical_tangency_rejected"] = summary(result)

    # Different owner families and key counts must use both actual radii.
    body, wraps, br, wr = line_bundle()
    wraps[0, :, 0] = np.float32(float(br)+float(wr)+.5e-6)
    result = qualify_native(body, wraps, br, wr, identity_world=True, config=config)
    assert result["qualified"] and not result["construction_gap_preserved_pass"]
    examples = result["contacts"]["construction_gap_failure_examples"]
    assert any(x["first"]["kind"] == "body" and x["second"]["kind"] == "wrap" for x in examples)
    results["unequal_radii_and_key_counts_body_wrap_clearance"] = summary(result)

    # Almost complete circle: its two ends are nonlocal and closer than 2r.
    # This tests the actual native cubic gate and the forbidden owner-only drop.
    body, wraps, br, wr = line_bundle(65, 5)
    body[:, :, 0] += .01
    wraps[:, :, 0] += .02
    angle = np.linspace(.1, 2*np.pi-.1, 65)
    body[0] = np.c_[1e-4*np.cos(angle), 1e-4*np.sin(angle), np.zeros(65)].astype("f4")
    result = qualify_native(body, wraps, br, wr, identity_world=True, config=config)
    assert not result["qualified"] and result["complete"]
    assert result["maximum_kappa_radius_upper"] < 1
    assert result["contacts"]["same_owner_nonlocal_pairs_tested"]>0
    assert result["contacts"]["same_fibre_self_unresolved_pairs"]>0
    assert result["contacts"]["interfibre_physical_unresolved_pairs"] == 0
    results["native_nonlocal_same_owner_near_return_rejected"] = summary(result)
    print(json.dumps({'stage':'geometric_cases_passed','groups':len(results)}),flush=True)

    crossings = []
    for angle in (2e-7, 1e-6, 1e-5, np.pi/2):
        for scale in (1e-9, 1e-6, 1., 1e6):
            a = np.zeros(3)
            b = np.array([1., 0., 0.])
            v = np.array([np.cos(angle), np.sin(angle), 0.])
            c = .37*b-.61*v
            d = c+v
            lower = float(projection_distance_lower(a*scale, b*scale, c*scale, d*scale))
            assert lower <= 1e-14*scale
            crossings.append({"angle_rad": angle, "scale": scale, "lower_bound": lower})
    results["r4_projection_near_parallel_and_scale_crossings"] = {"passed": True, "cases": crossings}

    # Guards use actual production API: no implicit cast or missing contract.
    for label, changed, options, expected in (
        ("float64_rejected", (bundle[0].astype("f8"), *bundle[1:]), {}, "stored_float32_required"),
        ("radius_float64_rejected", (*bundle[:2], float(bundle[2]), bundle[3]), {}, "constant_float32_radius_required"),
        ("missing_identity_rejected", bundle, {"identity_world": False}, "identity_world_contract_not_asserted"),
        ("leaf_budget_rejected", bundle, {"config": replace(config, max_leaves=100)}, "leaf_budget_exceeded"),
        ("wall_budget_rejected", bundle, {"config": replace(config, wall_seconds=1e-12)}, "wall_budget_exceeded"),
        ("pair_budget_rejected", bundle, {"config": replace(config, max_candidate_pairs=1)}, "candidate_pair_budget_exceeded"),
        ("neighbor_budget_rejected", bundle, {"config": replace(config, max_neighbors_per_query=1)}, "neighbor_budget_exceeded"),
        ("memory_budget_rejected", bundle, {"config": replace(config, memory_limit_mib=1)}, "memory_budget_exceeded"),
    ):
        args = {"identity_world": True, "config": config} | options
        result = qualify_native(*changed, **args)
        assert not result["qualified"] and not result["complete"]
        assert result["failures"][0]["code"] == expected, (label, result)
        results[label] = {"passed": True, "status": result["status"]}
    nonfinite = line_bundle()
    nonfinite[0][0, 0, 0] = np.nan
    result = qualify_native(*nonfinite, identity_world=True, config=config)
    assert result["failures"][0]["code"] == "nonfinite_input"
    results["nonfinite_rejected"] = {"passed": True}

    # Verify CLI, immutable NPZ hash, JSON finiteness, and fail exit status.
    print(json.dumps({'stage':'import_admission_and_budget_cases_passed','groups':len(results)}),flush=True)
    with tempfile.TemporaryDirectory(prefix="native-gate-source-test-", dir=HERE) as directory:
        path = Path(directory)/"analytic_only.npz"
        np.savez(path, body=bundle[0], wraps=bundle[1], body_radius=bundle[2], wrap_radius=bundle[3])
        original_hash = sha256_file(path)
        receipt = Path(directory)/"receipt.json"
        command = [sys.executable, str(HERE/"native_curve_gate.py"), str(path), "--output", str(receipt)]
        yes = subprocess.run(command+["--identity-world"], capture_output=True, text=True, timeout=20)
        print(json.dumps({'stage':'positive_cli_finished','returncode':yes.returncode}),flush=True)
        assert yes.returncode == 0, yes.stderr
        output = json.loads(receipt.read_text())
        assert output["qualified"] and output["input_npz_sha256"] == original_hash
        assert sha256_file(path) == original_hash
        no = subprocess.run(command, capture_output=True, text=True, timeout=20)
        assert no.returncode == 2
        assert not json.loads(receipt.read_text())["qualified"]
    results["cli_npz_immutability_hash_and_fail_exit"] = {"passed": True}

    report = {"all_assertions_passed": True,
              "scope": "Analytic native-gate source tests only. Temporary synthetic arrays were not delivered as carpet, and no production field was qualified.",
              "formal_interval_arithmetic": False,
              "failure_interpretation": "An unproven conservative bound is not an observed intersection.",
              "test_groups": len(results), "results": results,
              "source_sha256": {name: sha256_file(HERE/name) for name in
                                ("native_curve_gate.py", "certified_segments.py", "run_analytic_tests.py")}}
    (HERE/"analytic_tests_receipt.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"all_assertions_passed": True, "test_groups": len(results),
                      "receipt": str(HERE/"analytic_tests_receipt.json")}))


if __name__ == "__main__":
    run()
