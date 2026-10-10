"""Explicit native curvature and geometry requirements before any export."""
LOCAL_REQUIRED=('strict_nonfold_pass','positive_physical_clearance_pass',
                'same_fibre_self_clearance_pass','continuous_containment_pass',
                'original30um_transverse_disks_pass','roots_preserved',
                'unchanged_history_preserved','candidate_hash_matches_gate')

def decision(gates):
    missing=[name for name in LOCAL_REQUIRED if gates.get(name) is not True]
    local=not missing
    export=local and gates.get('construction_reserve_pass') is True and gates.get('production_loop_complete') is True and gates.get('actual_scene_contract_verified') is True
    promotion=export and gates.get('three_loop_qualified') is True and gates.get('matched_visual_review_pass') is True
    return {'local_repair_accepted':local,'production_export_allowed':export,'promotion_allowed':promotion,
            'failed_or_unproven_local_requirements':missing,
            'construction_reserve_pass':gates.get('construction_reserve_pass') is True,
            'positive_physical_clearance_is_separate_from_construction_reserve':True}
