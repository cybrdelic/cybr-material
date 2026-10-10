"""Bounded source contracts only. Do not execute held research solvers."""
from pathlib import Path
import os,subprocess,sys,tempfile,shutil
ROOT=Path(__file__).resolve().parents[1]
research=ROOT/'research/organic/research-20261006'
subprocess.run([sys.executable,str(research/'test_research_inputs.py')],check=True)
# The visual-model smoke writes tiny diagnostic images; isolate these from the
# curated source inventory rather than weakening the integrity checker.
source=ROOT/'research/asphalt/binder-refinement-20261009'
with tempfile.TemporaryDirectory(prefix='cybr-binder-source-') as tmp:
 target=Path(tmp)/'recipe';shutil.copytree(source,target,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 subprocess.run([sys.executable,str(target/'source/test_refinement.py')],check=True)
preset=ROOT/'material-workspace/shared/material-production/preview_presets'
subprocess.run([sys.executable,'-m','unittest','discover','-s',str(preset),'-v'],check=True)
print('Held research input contracts, isolated binder source and default preview preset checks passed.')
