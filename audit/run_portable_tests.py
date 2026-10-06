"""Run only source-only contracts; native-scene integration remains a separate gate."""
from pathlib import Path
import os
import sys
import unittest
ROOT = Path(__file__).resolve().parents[1]
P = ROOT / 'material-workspace/shared/material-production'
os.chdir(P)
sys.path[:0] = [str(P), str(P/'tests')]
from test_integration_contracts import Contracts
from test_path_config import Paths
# This named historical case intentionally validates a real selected .blend.
# It remains in the source suite and is NOT reported as passing without that asset.
asset_case = 'test_registry_35_unique_and_leather_explicitly_unqualified'
names = [n for n in unittest.defaultTestLoader.getTestCaseNames(Contracts) if n != asset_case]
suite = unittest.TestSuite([Contracts(n) for n in names])
suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(Paths))
from test_native_png_and_experiments import PNGContracts, ExperimentContracts
from test_showcase_guard import ShowcaseGuard
for case in (PNGContracts, ExperimentContracts, ShowcaseGuard):
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(case))
from materials.core import registry
assert len(registry()['materials']) == 35
result = unittest.TextTestRunner(verbosity=2).run(suite)
print('Separate native-scene gate requires restored selected scene/map artifacts:', asset_case)
if not result.wasSuccessful():
    raise SystemExit(1)
import subprocess
C = ROOT / 'material-workspace/shared/cloud-eevee'
for relative in ('audit-20261006/cache_guard/test_cache_validation.py',
                 'audit-20261006/cache_guard/test_queue_paths.py',
                 'audit-20261006/resource_peak/test_resource_contract.py'):
    subprocess.run([sys.executable, str(C / relative)], check=True)
print('All portable source and renderer contracts passed.')
