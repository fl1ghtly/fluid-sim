cmake_minimum_required(VERSION 3.28)
project(FluidSimulation LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 3.0.1
    GIT_SHALLOW ON
    EXCLUDE_FROM_ALL
    SYSTEM)
FetchContent_MakeAvailable(SFML)

FetchContent_Declare (tracy
    GIT_REPOSITORY https://github.com/wolfpld/tracy.git
    GIT_TAG v0.12.2
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(tracy)

# Add executable with all source files
add_executable(fluid_sim
    src/main.cpp
    src/Simulation.cpp
    src/Vector2f.cpp
    src/GridCell.cpp
    src/ParticleSystem.cpp
    src/FluidParameters.cpp
    src/PressureSystem.cpp
)

# Include directory for headers
target_include_directories(fluid_sim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(fluid_sim PRIVATE SFML::Graphics)
target_link_libraries(fluid_sim PRIVATE TracyClient)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
target_link_libraries(fluid_sim PRIVATE OpenMP::OpenMP_CXX)
endif()
find_package(TBB REQUIRED)
target_link_libraries(fluid_sim PRIVATE TBB::tbb)