cmake_minimum_required(VERSION 3.28)
project(FluidSimulation LANGUAGES CXX VERSION 1.0.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

# Enable tracy by default on debug builds
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(ENABLE_TRACY_DEFAULT ON)
else()
    set(ENABLE_TRACY_DEFAULT OFF)
endif()

option(
    ENABLE_TRACY 
    "Enable Tracy Profiling"
    ${ENABLE_TRACY_DEFAULT}
)

message(STATUS "ENABLE_TRACY: ${ENABLE_TRACY}")

include(FetchContent)

find_package(OpenMP)
find_package(TBB REQUIRED)

FetchContent_Declare(SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 3.0.1
    GIT_SHALLOW ON
    EXCLUDE_FROM_ALL
    SYSTEM)
FetchContent_MakeAvailable(SFML)

if(ENABLE_TRACY)
    FetchContent_Declare (tracy
        GIT_REPOSITORY https://github.com/wolfpld/tracy.git
        GIT_TAG v0.12.2
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    FetchContent_MakeAvailable(tracy)
endif()

# Create library
add_library(FluidSimulation src/Simulation.cpp)

# Include directory for headers
target_include_directories(FluidSimulation PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(FluidSimulation PRIVATE TBB::tbb)

if(OpenMP_CXX_FOUND)
target_link_libraries(FluidSimulation PRIVATE OpenMP::OpenMP_CXX)
endif()

if(ENABLE_TRACY)
    target_compile_definitions(FluidSimulation PRIVATE TRACY_ENABLED)
    target_link_libraries(FluidSimulation PRIVATE TracyClient)
endif()

# Add executable with all source files
add_executable(fluid_sim
    src/main.cpp
    src/Simulation.cpp
    src/Vector2f.cpp
    src/GridCell.cpp
    src/ParticleSystem.cpp
    src/FluidParameters.cpp
    src/Boundary.cpp
    src/PressureSystem.cpp
)

# Link dependencies only for executable
target_link_libraries(fluid_sim PRIVATE
    FluidSimulation
    SFML::Graphics
)

# Create targets to export
install(TARGETS FluidSimulation
    EXPORT FluidSimulationTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install all header files in include directory
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Export
install(EXPORT FluidSimulationTargets
    FILE FluidSimulationTargets.cmake
    NAMESPACE FluidSimulation::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FluidSimulation
)

# Create package config file
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FluidSimulationConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/FluidSimulationConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FluidSimulation
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/FluidSimulationConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/FluidSimulationConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FluidSimulation
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    FluidSimulationConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)