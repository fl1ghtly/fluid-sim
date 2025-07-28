cmake_minimum_required(VERSION 3.28)
project(FluidSimulation LANGUAGES CXX VERSION 1.0.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

find_package(OpenMP)
find_package(TBB)

set(ENABLE_MULTIPROCESSING "Enable multiprocessing support" ON)

# Create library
add_library(FluidSimulation 
    src/Simulation.cpp
    src/Vector2f.cpp
    src/GridCell.cpp
    src/FluidParameters.cpp
    src/Boundary.cpp
)

set(GCC_CLANG_WARNINGS
    -Wall
    -Wextra
    -Wpedantic
    -Wconversion
    -Wunused
    -Wsign-conversion
)

target_compile_options(FluidSimulation PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /WX>
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
        ${GCC_CLANG_WARNINGS}
    >
)

if(PROJECT_IS_TOP_LEVEL)
    add_subdirectory(apps EXCLUDE_FROM_ALL)
endif()

# Include directory for headers
target_include_directories(FluidSimulation PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Add multiprocessing support if packages are found
if(OpenMP_CXX_FOUND AND TBB_FOUND AND ENABLE_MULTIPROCESSING)
    target_compile_definitions(FluidSimulation PRIVATE MULTIPROCESSING_ENABLED)
    target_link_libraries(FluidSimulation PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(FluidSimulation PRIVATE TBB::tbb)
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