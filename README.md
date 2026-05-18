# frackit2porepy

Generate 3D fracture networks with [**Frackit**](https://git.iws.uni-stuttgart.de/tools/frackit) and export them to [**PorePy-compatible CSV**](https://github.com/pmgbergen/porepy) files.

This repository provides a reproducible **Docker-based** workflow: the custom Frackit application `frackit2porepy` is compiled into Frackit’s `appl/` tree at image build time, and the container writes all runtime outputs to a mounted `shared/` directory.

## Credits

- **Frackit repository:** https://git.iws.uni-stuttgart.de/tools/frackit  
  Frackit is developed at the University of Stuttgart (Institute for Modelling Hydraulic and Environmental Systems, IWS).
- **Frackit paper:** Gläser, Dennis, Bernd Flemisch, Holger Class, and Rainer Helmig. "Frackit: a framework for stochastic fracture network generation and analysis." Journal of Open Source Software 5, no. 56 (2020): 2291. DOI: [10.21105/joss.02291](https://doi.org/10.21105/joss.02291)

If you use this project in academic work, please cite Frackit appropriately.

## Quick start

From the repository root:

```bash
docker build -t frackit2porepy:latest .
docker run --rm -v ".\shared:/frackit/shared" frackit2porepy:latest
```

### Line endings note

If you see (in Windows):

```
/usr/bin/env: ‘bash\r’: No such file or directory
```

then some files may have CRLF line endings. Convert it to **LF** and rebuild.

## The `shared/` folder (runtime I/O contract)

The folder `shared/` is the only directory you normally need to interact with.

- It is **mounted into the container** at `/frackit/shared`.
- It contains the **runtime configuration** (`shared/config.toml`).
- It receives the **generated outputs** (CSV + visualization artifacts).

Workflow: **edit inputs in `shared/` → run the container → read outputs from `shared/`.**

### Expected files in `shared/`

#### Inputs

- `config.toml` — controls domain, constraints, and all fracture-family sampler parameters.

#### Outputs

- `disks.csv` — **PorePy elliptic fracture format** (domain line + one line per ellipse).  
  Intended to be loaded with:
  `porepy.fracs.fracture_importer.network_from_csv(..., has_domain=True)`.

  Format (3D):
  - First line (domain): `X_MIN, Y_MIN, Z_MIN, X_MAX, Y_MAX, Z_MAX`
  - Each ellipse line (8 floats):
    `CENTER_X, CENTER_Y, CENTER_Z, MAJOR_AXIS, MINOR_AXIS, MAJOR_AXIS_ANGLE, STRIKE_ANGLE, DIP_ANGLE`

  Angle conventions follow PorePy’s documentation:
  - `major_axis_angle`: rotation of the ellipse major axis from the +x axis in radians, measured **before** strike–dip
  - `strike_angle`: clockwise from +y (north) in the horizontal plane (radians)
  - `dip_angle`: rotation around the strike direction (radians)

- `quads.csv` — **polygonal (point-based) fracture format** (domain line + one line per quad).
  Intended to be loaded with:
  `porepy.fracs.fracture_importer.network_from_csv(..., has_domain=True)`.

  Format (3D):
  - First line (domain): `X_MIN, Y_MIN, Z_MIN, X_MAX, Y_MAX, Z_MAX`
  - Each quad line (12 floats):
    `P0_X, P0_Y, P0_Z, P1_X, P1_Y, P1_Z, P2_X, P2_Y, P2_Z, P3_X, P3_Y, P3_Z`

  Notes:
  - This file is separate from `disks.csv`.
  - How you import polygonal fractures into PorePy depends on your workflow; see placeholder section below.

- `families.csv` — per-family metadata written by the generator:
  `family_id,type,target_num,sampled_num`

- `network.geo` — optional Gmsh geometry file for visualization.
- `network.brep` — optional BREP file (if produced by the underlying workflow).

## Configuration (`shared/config.toml`)

The TOML file configures:

1. **Domain** (`[domain]`) and **Subdomain** (`[subdomain]`)
2. **Constraints** (`[constraints]`)
3. **A dynamic number of fracture families** (`[sampler] num = N` and `[sampler.i]` blocks)
4. **Output paths** (`[output]`)

### 1) Domain and Subdomain

The domain (3D box) has no influence on the workflow and is merely used as domain in the PorePy CSV file.

```toml
[domain]
xmin = 0.0
ymin = 0.0
zmin = 0.0
xmax = 100.0
ymax = 100.0
zmax = 100.0
```

The subdomain (3D box) defines the sampling box for points during fracture generation. PorePy will throw errors if fractures lie outside of the domain, such that the subdomain offers some control on fracture placement wrt. domains.

```toml
[subdomain]
xmin = 25.0
ymin = 25.0
zmin = 25.0
xmax = 50.0
ymax = 50.0
zmax = 50.0
```

### 2) Constraints

Constraints are applied during sampling via Frackit’s entity-network constraint system.

```toml
[constraints]
min_distance = 0.05
min_intersecting_angle_deg_self = 30.0
min_intersecting_angle_deg_other = 40.0
min_intersection_magnitude = 0.05
min_intersection_distance = 0.05
```

Notes:
- `*_self` applies within the same family.
- `*_other` applies between different families.

### 3) Fracture families (dynamic)

The number of families is controlled by:

```toml
[sampler]
num = 2
```

Each family is defined in its own table `[sampler.1]`, `[sampler.2]`, …, `[sampler.N]`.

Common keys:

- `type` — geometry type (`"disks"` or `"quads"`)
- `target_num` — number of accepted fractures to generate for this family

#### 3a) Disk family (`type = "disks"`)

```toml
[sampler.1]
type = "disks"
target_num = 50

# Axis lengths are sampled from Normal(mean, stddev)
[sampler.1.major_axis_length]
mean = 30.0
stddev = 6.5

[sampler.1.minor_axis_length]
mean = 24.0
stddev = 4.5

# Orientation is sampled as 3 independent Normal rotations (degrees)
[sampler.1.rotation_deg.x]
mean = 0.0
stddev = 7.5

[sampler.1.rotation_deg.y]
mean = 0.0
stddev = 7.5

[sampler.1.rotation_deg.z]
mean = 0.0
stddev = 7.5
```

#### 3b) Quadrilateral family (`type = "quads"`)

```toml
[sampler.2]
type = "quads"
target_num = 25

# Orientation (degrees) sampled from Normal(mean, stddev)
[sampler.2.strike_deg]
mean = 45.0
stddev = 5.0

[sampler.2.dip_deg]
mean = 90.0
stddev = 5.0

# Edge lengths sampled uniformly from [min, max]
[sampler.2.strike_length]
min = 30.0
max = 60.0

[sampler.2.dip_length]
min = 30.0
max = 60.0
```

### 4) Output paths

```toml
[output]
disks_csv = "disks.csv"
quads_csv = "quads.csv"
families_csv = "families.csv"
```
Paths are interpreted relative to the working directory inside the container; the entrypoint copies known outputs back to `/frackit/shared`.

## PorePy integration

A minimal example script (`porepy_example.py`) is provided that:

1. **Loads the fracture network** from `shared/disks.csv` (or `shared/quads.csv`) using PorePy's fracture importer.
2. **Creates a mixed-dimensional grid** (fractures + matrix).
3. **Runs a single-phase flow simulation** on the network.
4. **Exports results** to Paraview format (`.vtu` files).

### Usage

PorePy needs to be installed in order to run the example. This repo and associated docker image do not come with a working PorePy version.

### How it works

The `ImportedGeometry` class reads the CSV:

```python
self.fracture_network = pp.fracture_importer.network_from_csv(
    Path("shared/disks.csv"),
    has_domain=True,
    tol=expected_domain_size * 1e-6,
)
```

The first line of `disks.csv` defines the domain; subsequent lines define elliptical fractures. For **polygonal fractures** from `quads.csv`, adapt the script similarly (PorePy also supports polygon import via `network_from_csv`).

The `SinglePhaseFlowGeometry` class combines the imported geometry with PorePy's `SinglePhaseFlow` model. Running `pp.ModelRunner(model).run()` meshes, assembles, solves, and exports results.

### Customization

- **Mesh size:** Adjust `cell_size` and `cell_size_fracture` in the `meshing_arguments()` method.
- **Physics model:** Replace `SinglePhaseFlow` with a different PorePy model (e.g., `BiotContactMechanics`) for coupled simulations.
- **Boundary conditions & parameters:** Override `set_*` methods in the geometry/model class.

### Known issues

- If fractures intersect the domain boundary, the importer may fail. Ensure fractures are strictly interior.
- Requires PorePy with Gmsh ≥ 4.8 (OpenCascade ≥ 7.5 recommended for robust geometry handling).
