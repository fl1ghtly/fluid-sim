cmake_minimum_required(VERSION 3.10)
project(FluidSimulation LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 3.0.1
    GIT_SHALLOW ON
    EXCLUDE_FROM_ALL
    SYSTEM)
FetchContent_MakeAvailable(SFML)

# Add executable with all source files
add_executable(fluid_sim
    src/main.cpp
    src/fluid.cpp
    src/Vector2f.cpp
    src/GridCell.cpp
)

# Include directory for headers
target_include_directories(fluid_sim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(fluid_sim PRIVATE SFML::Graphics)