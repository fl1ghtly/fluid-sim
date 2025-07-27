# Fluid Simulation

An implementation of the Smoothed Particle Hydronamics (SPH) Algorithm as described by Ihmsen et al and Müller et al. [\[1](#1), [3\]](#3).

Static rigidbody-fluid interactions are modeled from the methodology of Akinci et al. [[2]](#2)

### Running
A sandbox and visual demo executable is available in the releases tab.

## Getting Started

### Prerequisites
Before building this repository, make sure to have the following installed.
* [CMake](https://cmake.org/)
* [oneTBB](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html)

### Building
To build, download the repository and run either 
```bash
cmake --build build --target all
```

```bash
cmake --build build --target sandbox
```

```bash
cmake --build build --target demo
```

to get the library, sandbox, or demo binaries respectively.

## Usage
To import and use this project for yourself add the following to your CMakeList.txt file.
```CMake
FetchContent_Declare(FluidSimulation
    GIT_REPOSITORY https://github.com/fl1ghtly/fluid-sim.git
    GIT_TAG 1.0.0
    GIT_SHALLOW ON
    EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(FluidSimulation)

target_link_libraries(YOUR_TARGET PRIVATE
    FluidSimulation
)
```

## References
<a id="1">[1]</a>
Ihmsen, M., Orthmann, J., Solenthaler, B., Kolb, A., & Teschner, M. (2014). SPH fluids in computer graphics. https://doi.org/10.2312/egst.20141034
<a id="2">[2]</a>
Akinci, N., Ihmsen, M., Akinci, G., Solenthaler, B., & Teschner, M. (2012). Versatile rigid-fluid coupling for incompressible SPH. ACM Transactions on Graphics (TOG), 31(4), 1-8. https://dl.acm.org/doi/10.1145/2185520.2185558
<a id="3">[3]</a>
Müller, M., Charypar, D., & Gross, M. (2003). Particle-based fluid simulation for interactive applications. In Proceedings of the 2003 ACM SIGGRAPH/Eurographics symposium on Computer animation (pp. 154-159). https://dl.acm.org/doi/10.5555/846276.846298