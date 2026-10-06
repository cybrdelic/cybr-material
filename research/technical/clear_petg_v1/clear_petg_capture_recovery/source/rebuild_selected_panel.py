"""Path-only replay of the preserved slab builder; reconstructed binary is explicit."""
from pathlib import Path
import sys,json,hashlib,resource,time
R=Path(__file__).resolve().parents[1];P=R.parents[1];BASE=Path('/workspace/scratch/b4387906eb93')
source=R/'source/original_build_petg_native4k.py';text=source.read_text();start=time.monotonic()
old="R=Path('/workspace/shared/material-slabs');NATIVE=Path('/workspace/shared/material-native4k');sys.path.insert(0,str(R/'source'));from studio import configure"
new="R=Path("+repr(str(R))+");sys.path.insert(0,"+repr(str(P.parent/'material-slabs/source'))+");from studio import configure"
assert text.count(old)==1;text=text.replace(old,new)
old="SRC=Path('/workspace/shared/material-technical-rebuild/scenes')/(cid+'_r5_hero.blend');P=NATIVE/'sets'/cid"
new="SRC=Path("+repr(str(BASE/'technical-recovery-20261006/restored/CYBR_technical_cloud_review/scenes'))+")/(cid+'_r5_hero.blend');P=R/'maps'"
assert text.count(old)==1;text=text.replace(old,new)
sys.argv=['rebuild_selected_panel.py','--','23_petg_transparent'];exec(compile(text,str(source),'exec'),{'__name__':'__main__','__file__':str(source)})
report=json.loads((R/'receipts/23_petg_transparent_native4k_slab.json').read_text());report.update(recovery='Exact source recipe and byte-exact native bound maps; newly saved reconstructed scene, not historical binary identity',builder_snapshot_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),wrapper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),wall_seconds=time.monotonic()-start,peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
(R/'receipts/reconstruction.json').write_text(json.dumps(report,indent=2)+'\n');print('CLEAR_SLAB_RECONSTRUCTED',report['scene_sha256'],report['peak_rss_mib'],flush=True)
