"""Read-only preferred presentation references, separate from numerical experiments.

A source integrity pass never substitutes for a matched visual comparison.
This module does not mutate the registry or select the newest file in a directory.
"""
from .core import ContractError

def reference_identity(reference):
    for key in ('sha256', 'canonical_render_sha256', 'labeled_image_sha256'):
        if reference.get(key):
            return 'sha256', reference[key]
    if reference.get('library_file_id'):
        return 'library_file_id', reference['library_file_id']
    raise ContractError('Showcase reference needs a pinned content hash or private legacy identity')

def preferred_reference(material):
    reference=material.get('preferred_visual_reference')
    if reference is not None:
        if reference.get('kind')!='complete_panel':
            raise ContractError('Preferred showcase reference must be a confirmed complete panel')
        if reference.get('status') not in ('restored_visual_baseline','verified_visible_improvement','verified_narrow_visual_improvement'):
            raise ContractError('Diagnostic or numeric-only evidence cannot be a preferred showcase reference')
        key, identity = reference_identity(reference)
        result = dict(reference)
        if key == 'sha256': result.setdefault('sha256', identity)
        return result
    selected=material.get('selected') or {}
    return dict(selected.get('appearance_anchor') or {})

def validate_visual_promotion(material,candidate,comparison):
    """Validate the explicit review contract; this does not automate aesthetic judgment."""
    baseline=preferred_reference(material)
    if candidate.get('kind')!='complete_panel':
        raise ContractError('Detail studies cannot replace the main panel')
    baseline_key, baseline_id = reference_identity(baseline)
    candidate_key, candidate_id = reference_identity(candidate)
    if comparison.get('baseline_' + baseline_key)!=baseline_id:
        raise ContractError('Comparison must use the pinned preferred baseline')
    if comparison.get('candidate_' + candidate_key)!=candidate_id:
        raise ContractError('Comparison must identify the actual persisted candidate')
    if comparison.get('result')!='clear_visual_improvement':
        raise ContractError('Numerical correctness or marginal change does not promote a showcase image')
    if not comparison.get('native_pixels_reviewed') or not comparison.get('complete_view_compared'):
        raise ContractError('Native complete-view comparison is required')
    if not comparison.get('review_receipt') or not comparison.get('review_receipt_sha256'):
        raise ContractError('Explicit visual review evidence is required')
    return True
