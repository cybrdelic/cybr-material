"""Explicit non-selected experiment adapter; never participates in selection."""
import re
from pathlib import Path
from .core import ROOT,ContractError,read,resolve,check_hash,sha,material

def experiment(name):
 if not re.fullmatch(r'[a-z][a-z0-9_]{0,63}',name):raise ContractError('Invalid experiment ID')
 path=ROOT/'experiments'/name/'experiment.json';r=read(path)
 if r.get('schema_version')!=1 or r.get('experiment_id')!=name:raise ContractError('Experiment schema/ID mismatch')
 if r.get('selected') is not False or r.get('quality_acceptance') is not False:raise ContractError('Experiments cannot promote themselves')
 material(r['material_id'])
 if not r.get('purpose') or not r.get('pass_gates'):raise ContractError('Experiment purpose and pass gates are required')
 check_hash(resolve(r['source_scene']),r['source_sha256'])
 if not isinstance(r.get('dependencies'),list):raise ContractError('Declared dependency list is required')
 for d in r['dependencies']:check_hash(resolve(d['path']),d['sha256'])
 if any(k in r['capture'] for k in ['source','source_sha256','output']):raise ContractError('Capture cannot override protected input/output')
 return r,path

def verify_experiment_build(b):
 r,path=experiment(b['experiment_id'])
 if sha(path)!=b['experiment_manifest_sha256']:raise ContractError('Experiment changed since build; rebuild explicitly')
 if b['source']['scene_sha256']!=r['source_sha256']:raise ContractError('Experiment source changed')
 for im in b['inspection']['images']:
  if not im['packed']:check_hash(im['path'],im['sha256'])
 return r
