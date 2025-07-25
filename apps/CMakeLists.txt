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

include(FetchContent)

FetchContent_Declare(SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 3.0.1
    GIT_SHALLOW ON
    EXCLUDE_FROM_ALL
    SYSTEM
)
FetchContent_MakeAvailable(SFML)

FetchContent_Declare (tracy
	GIT_REPOSITORY https://github.com/wolfpld/tracy.git
	GIT_TAG v0.12.2
	GIT_SHALLOW TRUE
	GIT_PROGRESS TRUE
	EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(tracy)

# Add executable with all source files
add_executable(fluid_sim
    EXCLUDE_FROM_ALL
    app.cpp
    ../src/Simulation.cpp
    ../src/Vector2f.cpp
    ../src/GridCell.cpp
    ../src/ParticleSystem.cpp
    ../src/FluidParameters.cpp
    ../src/Boundary.cpp
    ../src/PressureSystem.cpp
)

# Link dependencies only for executable
target_link_libraries(fluid_sim PRIVATE
    FluidSimulation
    SFML::Graphics
)

if(ENABLE_TRACY)
    target_compile_definitions(fluid_sim PRIVATE TRACY_ENABLED)
    target_link_libraries(fluid_sim PRIVATE TracyClient)
endif()