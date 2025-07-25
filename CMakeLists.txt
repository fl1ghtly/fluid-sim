cmake_minimum_required(VERSION 3.28)
project(FluidSimulation LANGUAGES CXX VERSION 1.0.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

find_package(OpenMP)
find_package(TBB REQUIRED)

# Create library
add_library(FluidSimulation 
    src/Simulation.cpp
    src/Vector2f.cpp
    src/GridCell.cpp
    src/FluidParameters.cpp
    src/Boundary.cpp
)

add_subdirectory(apps EXCLUDE_FROM_ALL)

# Include directory for headers
target_include_directories(FluidSimulation PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(FluidSimulation PRIVATE TBB::tbb)

if(OpenMP_CXX_FOUND)
target_link_libraries(FluidSimulation PRIVATE OpenMP::OpenMP_CXX)
endif()

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

include(CMakePackageConfigHelpers)

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

write_basic_package_version_file(
    FluidSimulationConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)