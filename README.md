# BiofilmAlgorithm

MATLAB algorithm for measuring biofilm-induced PDMS deformation using Newton's rings interference pattern analysis from time-lapse microscopy images.

## Overview

When a biofilm grows on a PDMS (polydimethylsiloxane) surface, it mechanically deforms the substrate. This deformation produces Newton's rings — concentric circular interference fringes visible under reflected light microscopy. By detecting these rings and measuring their diameters, the algorithm estimates the **radius of curvature** of the PDMS deformation over time, providing a quantitative measure of the mechanical forces exerted by the growing biofilm.

## How It Works

1. **Image loading** — Reads sequential time-lapse `.tif` frames (both raw microscopy and binarized versions).
2. **Polar transformation** — Converts each frame from Cartesian to polar coordinates using `ImToPolar`, centering on the biofilm. This turns concentric rings into horizontal lines, making them easier to detect.
3. **Pillar removal** — Masks out the central pillar region (a known artifact) so it isn't falsely detected as a ring.
4. **Biofilm radius estimation** — Models the biofilm as a circle from the binary images and computes its effective radius via area measurement (`bwarea`). Supports both full-circle and half-circle geometries.
5. **Dynamic masking** — Uses scaling parameters (`radiusparam`, `radiusparam0`) to crop the polar image to the biofilm boundary, isolating the Newton's rings from background noise.
6. **Ring detection** — Applies adaptive histogram equalization and Gaussian filtering, then thresholds the polar image to identify bright fringes. Extracts the radial positions (diameters) of detected rings.
7. **Interference condition matching** — Compares all pairs of detected ring diameters against the constructive interference condition:

$$\frac{D_m^2}{D_n^2} = \sqrt{\frac{m - 0.5}{n + 0.5}}$$

where *m* and *n* are ring order numbers. Matching pairs confirm valid Newton's rings.

8. **Radius of curvature calculation** — For each valid ring pair, computes the PDMS radius of curvature:

$$R = \frac{r^2}{\lambda (n - 0.5)}$$

where *r* is the ring radius, *&lambda;* is the illumination wavelength (default 500 nm), and *n* is the ring order. Averages across all valid pairs per frame.

9. **Time series output** — Tracks the biofilm radius and average PDMS radius of curvature across all frames, producing a scatter plot of *R* vs. frame number.

## Usage

```matlab
[biofilmradius, Rframelist, Raverage] = biofilm_PDMS( ...
    pillarlength, biofilmlocation, biofilmlocationBW, ...
    numframes, circletype, radiusparam, radiusparam0)
```

### Parameters

| Parameter | Type | Description |
|---|---|---|
| `pillarlength` | integer | Length of the central pillar in pixels (masked out to avoid false ring detection) |
| `biofilmlocation` | string | Path prefix to raw time-lapse images (e.g., `'stack19/overnight6_TL_'`) |
| `biofilmlocationBW` | string | Path prefix to binarized images (e.g., `'190319_radius/radiusfilm'`) |
| `numframes` | integer | Number of frames to process |
| `circletype` | `'full'` or `'half'` | Whether the biofilm image shows a full circle or half circle |
| `radiusparam` | float | Dynamic radius scaling factor — pixels beyond `radius * radiusparam` are masked |
| `radiusparam0` | float | Initial radius scaling factor — used until the biofilm grows past this boundary. Set equal to `radiusparam` if Newton's rings are not stationary |

### Outputs

| Output | Type | Description |
|---|---|---|
| `biofilmradius` | array | Estimated biofilm radius (pixels) for each frame |
| `Rframelist` | array | Frame indices corresponding to `Raverage` |
| `Raverage` | array | Average PDMS radius of curvature per frame |

### Example

```matlab
[biofilmradius, Rframelist, Raverage] = biofilm_PDMS( ...
    80, 'stack19/overnight6_oldChannels_37C_30C_3_190319_TL_', ...
    '190319_radius/radiusfilm', 294, 'half', 1.4, 2.0);
```

## Image Naming Convention

The algorithm expects sequentially numbered `.tif` files with zero-padded frame numbers appended to the path prefix:

```
<biofilmlocation>0001.tif
<biofilmlocation>0002.tif
...
<biofilmlocation>0294.tif
```

The same convention applies to the binarized (`BW`) images.

## Dependencies

- **MATLAB** with the **Image Processing Toolbox** (`adapthisteq`, `im2gray`, `imread`, `bwarea`, `bwareafilt`, `imresize`, `fspecial`, `montage`)
- [**ImToPolar**](https://www.mathworks.com/matlabcentral/fileexchange/17933-polar-to-from-rectangular-transform-of-images) — Cartesian-to-polar image transform (available on MATLAB File Exchange)

## Notes

- The algorithm displays a montage of the processed polar image alongside the original for each frame, allowing visual verification that ring detection is working correctly. Use this to tune `radiusparam` and `radiusparam0`.
- The tolerance for ring-pair matching (`epsilon`) is set to `0.01` by default. Adjust if too few or too many matches are found.
- The wavelength `lambda` is set to `500 nm`. Change this to match your actual illumination source.
