# Material gallery

The gallery is an additive presentation of existing evidence. It never calls the
material generator, replaces native maps/renders, edits the frozen engineering
registry, or republishes the replaceable suite ZIP. Original SHA256SUMS remains
the suite receipt; `gallery/catalog.json` is the separate gallery image receipt.

## Browse and host

```sh
python -m http.server 8000
# Open http://localhost:8000/gallery/
```

All pages are static HTML. JavaScript only adds the material-family filter;
navigation and individual native PNG links work without it, including from a
local checkout. Serve the repository root, not the gallery directory alone:
native assets intentionally stay in their original sibling directories.

GitHub's repository viewer displays HTML as source. The README's
`docs/MATERIAL_GALLERY.md` index therefore also links every native base-color
image and Cycles specimen directly for people browsing on GitHub.

The separate gallery workflow validates pull requests and deploys a static
directory to GitHub Pages on master changes affecting the gallery. Pages must
use **GitHub Actions** as its publishing source in repository Settings → Pages.
The deployment URL is reported by the workflow, not assumed from a hostname.
Gallery-workflow changes and changes to the suite workflow alone are excluded
from automatic suite regeneration, including this gallery's scheduling edit.
Use the existing manual suite rebuild action when intentionally changing only
its workflow configuration. Source and Blender-installer changes keep their
normal automatic full-build behavior.
The deployment packages original native files without rerendering or changing
the downloadable suite. `python tools/gallery/package.py` builds the same
directory under `build/material-gallery` for another static host.

## Rebuild display derivatives

```sh
python -m pip install Pillow
python tools/gallery/build.py
python tools/gallery/verify.py
```

The builder validates original specimen/proof hashes against their existing
manifests. The reviewed supplemental plastic image is pinned to the exact
recorded reference hashes. Missing or changed accepted images fail the build.

Native maps stay 4096 × 4096. The generator makes 360/800-pixel WebP display
derivatives and never overwrites source images. HTML declares dimensions,
responsive srcsets, lazy loading and asynchronous decode. Native PNGs are
fetched only when their links are opened. Sixteen-bit height previews are
display approximations; the original Gray16 data remains in the linked PNG.

## Rendering scope

The ten rebuilt specimens and ten matched proofs are 1024 × 1024. Existing
still lifes and architecture retain their exact recorded dimensions and
render settings. Plastic's original render area is 960 × 960; its labeled
presentation canvas is larger. Its labels are preserved rather than
misrepresented as additional rendered pixels.

The original editable detail scene, native map dependencies and an additive
`tools/gallery/render_variant.py` recipe are available for a new native render.
The recipe rotates the existing lights by 90 degrees, preserves geometry and
shaders, renders to a new path and records scene/dependency hashes, camera,
lights, actual resolution, sample budget and denoising. It refuses to overwrite
a previous supplementary image. Its default is 4096 × 4096 / 128 maximum / 32
minimum samples, with guided OIDN. This is a supplementary capture setting,
not the accepted baseline's raw 768/128 setting.

```sh
blender -b blender/CYBR_Cycles_Details.blend \
  --python tools/gallery/render_variant.py -- \
  --material 05_champagne_brass --size 4096 --samples 128 --minimum 32
```

No new render from that recipe is included at this checkpoint: the official
Blender 4.5.3 archive passed its published SHA-256 verification, but its executable
failed at startup with SIGBUS in the implementation environment. The existing
public repository renders are reused. A future successful
supplementary capture needs visual review and an explicit gallery addition;
neither a recipe nor a numerical check silently promotes an appearance.

## Selection boundaries

The rebuilt ten-material suite and the separate engineering registry are
different evidence lines. The gallery does not overwrite either selection.
The registry started with 35 identities and now has a separate Victorville
entry. The plastic panel is already public in this repository and is directly displayed.
LCD/CRT and clear PETG review provenance is recorded, but their exact images
and scenes require archived inputs. This presentation uses files already published
in this repository; additional archived outputs need a separate publication step.
Concrete, new carpet constructions and organic
formation experiments stay unqualified. Archive IDs remain recovery references,
never fake public image URLs.

The original distribution receipt, baselines, maps, editable scenes, renderer
settings and registry bytes remain unchanged by the presentation build.

## Browser checks

`tools/gallery/browser-check.cjs` exercises 390 px mobile and 1440 px desktop
viewports, family filtering, original 4096 PNG response headers, zero initial
native-PNG requests, a 1 MiB initial-transfer budget, horizontal overflow and
navigation with JavaScript disabled. The gallery workflow installs an isolated
Playwright browser and runs these checks with screenshots.

The implementation host blocked local browser startup at the Unix socket
operation and also blocked a local HTTP listener. Local browser checks therefore
did not pass or produce screenshots; CI must establish those results on an
ordinary runner. The static verifier checks every HTML asset/link and original
image hash without requiring a browser.
