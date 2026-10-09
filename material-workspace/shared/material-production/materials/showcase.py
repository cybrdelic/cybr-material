"""Read-only preferred presentation references, separate from numerical experiments.

A source integrity pass never substitutes for a matched visual comparison.
This module does not mutate the registry or select the newest file in a directory.
"""
from .core import ContractError

def preferred_reference(material):
    reference=material.get('preferred_visual_reference')
    if reference is not None:
        if reference.get('kind')!='complete_panel' or not reference.get('library_file_id'):
            raise ContractError('Preferred showcase reference must be a confirmed complete panel')
        if reference.get('status') not in ('restored_visual_baseline','verified_visible_improvement'):
            raise ContractError('Diagnostic or numeric-only evidence cannot be a preferred showcase reference')
        return dict(reference)
    selected=material.get('selected') or {}
    return dict(selected.get('appearance_anchor') or {})

def validate_visual_promotion(material,candidate,comparison):
    """Validate the explicit review contract; this does not automate aesthetic judgment."""
    baseline=preferred_reference(material)
    if candidate.get('kind')!='complete_panel':
        raise ContractError('Detail studies cannot replace the main panel')
    if comparison.get('baseline_library_file_id')!=baseline.get('library_file_id') or not baseline.get('library_file_id'):
        raise ContractError('Comparison must use the pinned preferred baseline')
    if comparison.get('candidate_library_file_id')!=candidate.get('library_file_id') or not candidate.get('library_file_id'):
        raise ContractError('Comparison must identify the actual persisted candidate')
    if comparison.get('result')!='clear_visual_improvement':
        raise ContractError('Numerical correctness or marginal change does not promote a showcase image')
    if not comparison.get('native_pixels_reviewed') or not comparison.get('complete_view_compared'):
        raise ContractError('Native complete-view comparison is required')
    if not comparison.get('review_receipt') or not comparison.get('review_receipt_sha256'):
        raise ContractError('Explicit visual review evidence is required')
    return True
