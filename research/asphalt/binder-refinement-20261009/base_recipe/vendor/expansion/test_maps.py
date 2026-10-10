"""CPU numerical checks, authored parameter extremes, deterministic export validation."""
import hashlib,json,tempfile,unittest
from pathlib import Path
import sys
import numpy as np
from PIL import Image
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from expansion.catalog import SPECS,get
from expansion.maps import generate,normals,write

class Recipes(unittest.TestCase):
    def test_all_repeatable_finite_range_and_normal_unit(self):
        for s in SPECS:
            with self.subTest(recipe=s['id']):
                a=generate(s,128);b=generate(s,128)
                for key in a:
                    self.assertTrue(np.array_equal(a[key],b[key]),key)
                    self.assertTrue(np.isfinite(a[key]).all(),key)
                    self.assertGreaterEqual(float(a[key].min()),0);self.assertLessEqual(float(a[key].max()),1)
                v=normals(a['Height'],s)*2-1
                self.assertLess(float(abs(np.linalg.norm(v,axis=-1)-1).max()),2e-6)
    def test_petg_boundary_changes_roughness_separate_transmission(self):
        for id in ('22_petg','23_petg_transparent'):
            s=get(id)
            for pitch in (.0001,.0004):
                s['layer_height_m']=pitch;a=generate(s,256)
                self.assertGreater(np.corrcoef(a['PrintBoundary'].ravel(),a['Roughness'].ravel())[0,1],.98)
                self.assertEqual(float(a['Transmission'][0,0]),float(np.float32(s.get('transmission',0))))
                self.assertTrue(np.all(a['Opacity']==1),'glass is not alpha transparency')
    def test_display_color_and_emission_are_separate(self):
        for id in ('18_lcd','19_crt'):
            a=generate(get(id),256)
            self.assertLess(float(a['BaseColor'].max()),.04)
            self.assertGreater(float(a['Emission'].max()),.8)
            self.assertTrue(np.all(np.sum(a['Emission']>0,axis=-1)<=1))
    def test_concrete_paste_vs_aggregate(self):
        paste=generate(get('35_cement_paste'),128);stone=generate(get('33_concrete_exposed_aggregate'),128)
        self.assertEqual(float(paste['AggregateMask'].max()),0)
        self.assertGreater(float(stone['AggregateMask'].max()),.9)
    def test_export_bit_depth_metadata_and_dx(self):
        root=Path(__file__).resolve().parents[2]/'expansion/test_exports'
        s=get('22_petg');m=write(s,128,root);out=root/s['id']
        with Image.open(out/'Height.png') as im:self.assertEqual(im.mode,'I;16')
        with Image.open(out/'Normal_OpenGL.png') as im:gl=np.asarray(im)
        with Image.open(out/'Normal_DirectX.png') as im:dx=np.asarray(im)
        self.assertTrue(np.array_equal(gl[...,0],dx[...,0]));self.assertTrue(np.all(abs(gl[...,1].astype(int)+dx[...,1].astype(int)-255)<=1))
        self.assertEqual(m['maps']['BaseColor.png']['color_space'],'sRGB')
        self.assertEqual(m['maps']['Transmission.png']['color_space'],'linear / Non-Color')
        self.assertEqual(m['thickness_scale_m'],.0012)
        self.assertEqual(len(m['checks']),len(list(out.glob('*.png'))))

    def test_json_parameters_and_invalid_intersecting_chain(self):
        s=get('22_petg',{'22_petg':{'layer_height_m':.0004}})
        self.assertEqual(s['layer_height_m'],.0004)
        with self.assertRaises(AssertionError):get('14_chainmail',{'14_chainmail':{'pitch_m':.005}})

if __name__=='__main__':unittest.main()
