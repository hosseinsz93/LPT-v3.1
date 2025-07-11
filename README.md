# LPT-v3.1 - Lagrangian Particle Tracking

A computational fluid dynamics (CFD) code for Lagrangian particle tracking, designed to work with VFS Geophysics flow simulations. This code tracks the movement of particles through fluid flows using Lagrangian methods, particularly useful for studying particle transport phenomena in biological and engineering applications.

## Overview

LPT-v3.1 is a C++ based particle tracking simulation that includes:
- **Breath particle generation** for respiratory flow studies
- **Particle acceleration calculations** with drag forces
- **3D geometric operations** for complex geometries
- **Spatial indexing** for efficient collision detection
- **Data extraction tools** for post-processing

## Features

- **Multi-platform support**: Configured for various HPC systems (Zagros, SeaWulf, etc.)
- **Breath particle simulation**: Specialized modules for respiratory particle tracking
- **Advanced geometry handling**: Comprehensive 3D geometric operations and spatial queries
- **Efficient algorithms**: R-tree spatial indexing for performance optimization
- **Flexible configuration**: Parameter control through `.dat` files
- **Python integration**: Particle generation scripts in Python

## Key Components

### Core Modules
- `BreathParticles.cpp/h` - Breath particle dynamics and generation
- `calcParticleAcceleration.cpp/h` - Particle physics calculations
- `DragCalcDebug.cpp/h` - Drag force computations with debugging

### Geometry Engine
- `GeomRtree.cpp/h` - R-tree spatial indexing
- `GeomBoundingBox.cpp/h` - Bounding box operations
- `GeomSpatialIndex.cpp/h` - Spatial query optimization
- Various geometric primitives (Point3, Line3, Plane3, Sphere3, etc.)

### Utilities
- `BreathParticleGeneration.py` - Python script for particle initialization
- `dataPartExtr.cpp` - Data extraction and post-processing
- `control.dat` - Simulation parameters and configuration

## Building the Project

### Prerequisites
- C++ compiler with C++11 support
- PETSc library (version 3.1 or higher)
- MPI implementation (OpenMPI recommended)
- ACML or equivalent math library
- Python 3.x for particle generation scripts

### Compilation

The project includes build scripts for different systems:

```bash
# Standard build
./compile

# Debug build
./compile.debug
```

The build system automatically detects the host system and configures appropriate compiler flags and library paths.

## Configuration

Simulation parameters are controlled through `control.dat`:

```
-ren 25000          # Reynolds number
-cfl 0.2            # CFL number for stability
-dt 0.01            # Non-dimensional time step
-totalsteps 10      # Total simulation steps
-chact_leng 25.4    # Characteristic length
```

## Usage

### Basic Workflow

1. **Configure simulation parameters** in `control.dat`
2. **Generate particles** using `BreathParticleGeneration.py`
3. **Compile the code** using the provided build scripts
4. **Run the simulation** with the compiled executable
5. **Extract data** using the provided extraction tools

### Particle Generation

```bash
python3 BreathParticleGeneration.py
```

This generates initial particle distributions saved to `BreathParticleList.dat`.

### Running Simulations

After compilation, run the tracking simulation:

```bash
./tracking
```

## File Structure

- **Source files**: `*.cpp`, `*.h` - Core C++ implementation
- **Python scripts**: `*.py` - Particle generation and utilities
- **Data files**: `*.dat` - Configuration and particle data
- **Build scripts**: `compile*` - Compilation automation
- **Object files**: `*.o` - Compiled objects (generated during build)

## Contributing

This is a research code developed for computational fluid dynamics applications. When making modifications:

1. Maintain compatibility with the existing parameter system
2. Follow the established naming conventions
3. Update version numbers appropriately (e.g., `v04.py` for new versions)
4. Test on target HPC systems before committing changes


## Contact

For questions about the code or collaboration opportunities, please contact the development team through the VFS Geophysics organization.
