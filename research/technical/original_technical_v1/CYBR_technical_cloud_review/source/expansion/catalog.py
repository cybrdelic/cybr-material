"""Authored physical estimates in SI units, never measured claims."""
from copy import deepcopy
import json
from pathlib import Path

def spec(id, family, variant, tile=.12, relief=.001, rough=.5, color=(.4,.4,.4), **kw):
    return dict(id=id,name=f'{family} / {variant}',family=family,variant=variant,tile_m=tile,
                height_scale_m=relief,roughness=rough,color=color,seed=601+len(id)*71,
                ior=1.5,metallic=0,representation='solid displaced specimen',**kw)

SPECS=[
 spec('11_bark','Bark','Longitudinal plates',.18,.008,.88,(.28,.16,.075)),
 spec('12_lattice_wire','Lattice','Woven stainless wire',.06,.00008,.28,(.56,.6,.64),pitch_m=.006,wire_radius_m=.00042),
 spec('13_lattice_expanded','Lattice','Expanded diamond metal',.10,.00012,.32,(.4,.46,.5),pitch_m=.01,wire_radius_m=.00065),
 spec('14_chainmail','Chainmail','Japanese four in one',.08,.00004,.26,(.5,.56,.6),pitch_m=.010,ring_radius_m=.004,wire_radius_m=.00045),
 spec('15_future_ceramic','Engineered','Ceramic hex armor',.18,.0002,.25,(.19,.25,.29),panel_pitch_m=.016,coat=.25),
 spec('16_future_ribbed','Engineered','Anodized cooling fins',.12,.00008,.32,(.08,.15,.2),fin_pitch_m=.003,fin_height_m=.007),
 spec('17_blanket','Blanket','Soft brushed fleece',.12,.00008,.86,(.65,.47,.3),fiber_length_m=.005,fiber_radius_m=.000045,fiber_count=18000,sheen=.5),
 spec('18_lcd','LCD','RGB stripe display study',.035,.000003,.18,(.012,.012,.012),pixel_pitch_m=.00027,glass_thickness_m=.0007),
 spec('19_crt','CRT','RGB phosphor grille study',.035,.000004,.24,(.018,.018,.018),pixel_pitch_m=.00065,glass_thickness_m=.002),
 spec('20_plastic_smooth','Plastic','Molded satin',.12,.000018,.28,(.035,.18,.21),coat=.15),
 spec('21_plastic_texture','Plastic','Molded stipple',.12,.00012,.51,(.11,.14,.16)),
 spec('22_petg','PETG','Opaque printed wall',.02,.00016,.31,(.035,.32,.4),layer_height_m=.0002,bead_width_m=.00045,wall_thickness_m=.0012,anisotropy=.38),
 spec('23_petg_transparent','PETG','Clear printed wall',.02,.00016,.17,(.91,.97,.98),layer_height_m=.0002,bead_width_m=.00045,wall_thickness_m=.0012,transmission=.96,anisotropy=.38),
 spec('24_carpet','Carpet','Loop wool pile',.12,.0001,.89,(.28,.32,.26),fiber_length_m=.003,fiber_radius_m=.000065,fiber_count=13000,sheen=.25),
 spec('25_tile_ceramic','Tile','Glazed ceramic and grout',.30,.0005,.19,(.76,.78,.69),tile_pitch_m=.10,grout_width_m=.003,coat=.35),
 spec('26_tile_terrazzo','Tile','Honed terrazzo and grout',.30,.0004,.35,(.59,.56,.5),tile_pitch_m=.10,grout_width_m=.002),
 spec('27_grass','Grass','Blade and thatch tuft',.14,.0001,.73,(.12,.27,.045),blade_length_m=.04,blade_width_m=.0018,blade_count=1600),
 spec('28_asphalt','Asphalt','Rolled mineral aggregate',.25,.004,.88,(.048,.049,.045)),
 *[spec(f'{29+i}_concrete_{v}','Concrete',v.replace('_',' ').title(),.32,
        .004 if v in ('board_formed','exposed_aggregate','weathered') else .0012,
        .28 if v=='polished' else .77,(.49,.48,.44),finish=v)
   for i,v in enumerate(('cast','board_formed','brushed','polished','exposed_aggregate','weathered'))],
 spec('35_cement_paste','Cement paste','Troweled neat binder',.32,.0007,.72,(.52,.51,.47),finish='paste'),
]
GEOMETRY_FAMILIES={'Lattice','Chainmail','Blanket','Carpet','Grass','PETG','LCD','CRT','Engineered'}
for index,s in enumerate(SPECS):
    s['seed']=601+index*71
    if s['family'] in GEOMETRY_FAMILIES:s['representation']='explicit procedural geometry plus surface maps'
    s['physical_values']='authored approximations, not measured material data'
    s['tiling']='finite specimen' if s['family'] in GEOMETRY_FAMILIES else 'periodic field'
    if s['family']=='Bark':s['tiling']='finite specimen'
    s['ior']=1.57 if s['family']=='PETG' else 1.5
    if s['family'] in {'Lattice','Chainmail'}:s['metallic']=1
    if s['id']=='16_future_ribbed':s['metallic']=1

def get(id,parameters=None):
    s=deepcopy(next(s for s in SPECS if s['id']==id))
    if parameters:
        overrides=json.loads(Path(parameters).read_text()) if isinstance(parameters,(str,Path)) else parameters
        s.update(overrides.get(id,{}))
    assert .001<=s['tile_m']<=5 and s['height_scale_m']>0,'Invalid physical extent/relief'
    if s['family']=='PETG':
        assert .00008<=s['layer_height_m']<=.0006,'Layer pitch outside supported authoring range'
        assert .0004<=s['wall_thickness_m']<=.004,'Unsupported wall thickness'
        assert .0002<=s['bead_width_m']<=.001,'Unsupported bead width'
    if s['family']=='Chainmail':
        p=s['pitch_m'];r=s['ring_radius_m'];w=s['wire_radius_m']
        assert p>2*(r+w),'Main rings intersect'
        assert min(r-.15*p,.85*p-r)>2*w,'Bridge rings do not clear main links'
    return s
