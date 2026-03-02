# BiofilmAlgorithm

Companion code for the paper:

> **[Pushing Matters: Growth Dynamics of Confined Biofilms](Pushing_Matters_Growth_Dynamics_of_Confined_Biofilms.pdf)**
> Twana Cheragwandi, George Fortune, and Nuno M. Oliveira
> *KTH Royal Institute of Technology & University of Cambridge*

## Background

Bacterial biofilms — surface-bound colonies enclosed in a self-produced extracellular matrix — have been observed to exhibit oscillatory growth when vertically confined ([Liu et al., *Nature*, 2015](https://doi.org/10.1038/nature14660)). Previous work attributed this to metabolic co-dependence (nutrient–ammonium exchange between peripheral and interior cells), but did not account for the physical effects of the confinement itself.

In the experimental setup studied here, a *Bacillus subtilis* biofilm grows sandwiched between a glass plate and a polydimethylsiloxane (PDMS) sheet. As the biofilm expands vertically, it pushes the elastic PDMS upward, creating a curved deformation. This curvature produces **Newton's rings** — concentric interference fringes visible in top-view microscopy footage. By tracking these rings over time, we can indirectly measure the PDMS radius of curvature and quantify how vertical expansion couples to horizontal growth.

The theoretical framework by [Fortune et al., *Phys. Rev. Lett.*, 2022](https://doi.org/10.1103/PhysRevLett.128.178102) predicts that the PDMS curvature relates to the biofilm radius as:

$$R_P = \frac{1}{2\lambda} \left[ R_B^2(1 - \Xi) + R_B \Xi \right]$$

where $\lambda$ is a flatness constant and $\Xi$ depends on the mechanical properties of the system. Our image analysis confirms this relationship across multiple independent experiments, and further derives a **growth rate equation**:

$$\log R_P R_B^{-4} = \log\frac{1}{2\lambda} - gt$$

enabling extraction of the biofilm growth rate $g$ from the slope of a linear fit. Across three experimental setups, we measured an average growth rate of $g = 0.0046\ \text{min}^{-1}$.

## What This Code Does

This MATLAB algorithm processes time-lapse microscopy images of confined biofilms to:

1. **Extract the biofilm radius** $R_B$ by modelling the biofilm as a circle and computing $r = \sqrt{A/\pi}$ from binarized images.
2. **Detect Newton's rings** in the interference pattern surrounding the biofilm.
3. **Compute the PDMS radius of curvature** $R_P$ from the ring radii using $r_N = [\lambda R(N - 1/2)]^{1/2}$, validated against the constructive interference condition.
4. **Track both quantities over time**, producing the datasets needed to verify the theoretical model and extract growth rates.

## How It Works

1. **Image loading** — Reads sequential time-lapse `.tif` frames (both raw microscopy and binarized versions).
2. **Polar transformation** — Converts each frame from Cartesian to polar coordinates using `ImToPolar`, centering on the biofilm. This turns concentric rings into horizontal lines, making them easier to detect.
3. **Pillar removal** — Masks out the central pillar region (a known artifact) so it isn't falsely detected as a ring.
4. **Biofilm radius estimation** — Models the biofilm as a circle from the binary images and computes its effective radius via area measurement (`bwarea`). Supports both full-circle and half-circle geometries.
5. **Dynamic masking** — Uses scaling parameters (`radiusparam`, `radiusparam0`) to crop the polar image to the biofilm boundary, isolating the Newton's rings from background noise.
6. **Ring detection** — Applies adaptive histogram equalization and Gaussian filtering, then thresholds the polar image to identify bright fringes. Extracts the radial positions (diameters) of detected rings.
7. **Interference condition matching** — For each pair of detected ring radii, verifies the constructive interference condition:

$$\frac{r_{N+1}}{r_N} = \sqrt{\frac{N + 1/2}{N - 1/2}}$$

Only ring pairs satisfying this within tolerance $\epsilon = 0.01$ are accepted.

8. **Radius of curvature calculation** — For each valid ring, computes the PDMS radius of curvature:

$$R = \frac{r^2}{\lambda (N - 1/2)}$$

where $r$ is the ring radius, $\lambda$ is the illumination wavelength (default 500 nm), and $N$ is the ring order. Averages across all valid pairs per frame.

9. **Time series output** — Tracks $R_B$ and $R_P$ across all frames, producing the data for the $R_P$ vs $R_B$ and growth rate plots shown in the paper (Figures 2 and 3).

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
| `biofilmradius` | array | Estimated biofilm radius $R_B$ (pixels) for each frame |
| `Rframelist` | array | Frame indices corresponding to `Raverage` |
| `Raverage` | array | Average PDMS radius of curvature $R_P$ per frame |

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

## References

1. Fortune et al., "Biofilm growth under elastic confinement," *Physical Review Letters*, 128, 178102, 2022.
2. Liu et al., "Metabolic co-dependence gives rise to collective oscillations within biofilms," *Nature*, 523, 550–554, 2015.
3. Gunka et al., "Metabolic co-dependence gives rise to collective oscillations within biofilms," *MicroReview*, 85, 213–224, 2012.
