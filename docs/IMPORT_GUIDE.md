# Import and physical setup

## Blender

Open `blender/CYBR_Cycles_Details.blend` for ten displaced close-ups, or
`blender/CYBR_Cycles_Still_Lifes.blend` for the three scenes. Select the desired
scene and press F12. The materials and texture files are editable; keep the
`materials/` directory beside `blender/` so the relative references resolve.

For reusable material assets, add the `blender/` directory in Preferences →
File Paths → Asset Libraries, then use `CYBR_Material_Atelier.blend` in the Asset
Browser. The asset library has family catalogs and rendered preview thumbnails.

The group exposes tint, roughness offset, normal strength, wear tint and UV
scale. Default normal strength is 1. The two woods are finite procedural board cuts; their
default extents are 550 mm (walnut) and 450 mm (oak). Their shader derives normals
from physical height to handle fitted UV cuts. They are not certified seamless. Map at the physical
tile width in
`material.json` before judging the appearance.

For actual displacement, subdivide the mesh adequately, use `Height_Macro.png`
at scale `height_scale_m` and midlevel 0.5, and enable **Displacement Mode** in
the group. That selects the residual `Normal_Micro_OpenGL.png`. A regular
undisplaced mesh uses the full `Normal_OpenGL.png`.

`Opacity.png` is connected to Principled Alpha. Linen uses it for the openings
between yarns. Other materials have an all-white opacity map. A cut fabric edge
or close-up free-standing fiber still needs actual geometry.

## Unreal / glTF / other metallic-roughness engines

Set Base Color to sRGB. Disable sRGB for roughness, metalness, AO, height,
opacity and packed data. Use the normal format expected by the engine:
DirectX (-Y) for Unreal; OpenGL (+Y) for glTF and many other tools. Do not flip
the green channel twice.

The ORM texture packs AO in red, roughness in green and metallic in blue.
Normals and roughness should not share an sRGB decode with base color. Connect
linen opacity using the engine's suitable mask or transparency mode; an opaque
shader will close all weave holes.

Displacement support varies. The full normal is the default for a mesh without
real displacement. Parallax cannot reproduce a displaced silhouette or undercut
stone cavities.

## Unity HDRP

`exports/unity_hdrp/<material>/HDRP_MaskMap.png` uses:

| Channel | Value |
|---|---|
| R | Metallic |
| G | Ambient occlusion |
| B | Detail mask = 1 |
| A | Smoothness = 1 - roughness |

Treat the mask map as linear data. Base color remains sRGB. Choose normal map
import settings consistent with the supplied normal convention. Connect linen
opacity separately; it is not packed in the mask map.

## Unity URP

`URP_MetallicSmoothness.png` stores metallic in red and smoothness in alpha;
green and blue are 1. `URP_Occlusion.png` stores AO in green, with red and blue
equal to 1. Both are linear data. Linen opacity is a separate map.

## Color and AO

The saved scene previews use AgX. They show lit appearance, not an sRGB swatch
matched to a flat RGB value. Broad shadows and highlights in the renders must
not be copied back into Base Color. AO is a local cavity approximation and is
not multiplied into the Blender base color.

## Wear and geometry

Tiles describe surface history in the plane of the material. They cannot infer
which edge of a new object is handled or which part of a board is cut across
the grain. Add suitable end grain, edge wear, cut-hide fibers and broken-thread
geometry when the camera requires them. The supplied steel plate demonstrates
an exposed-metal bevel rather than pretending a generic tile knows its edges.
