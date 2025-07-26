# Add executable target
add_executable(sandbox
    EXCLUDE_FROM_ALL
    src/Sandbox.cpp
    src/Visualization/ParticleSystem.cpp
)

add_executable(demo
    EXCLUDE_FROM_ALL
    src/Demo.cpp
    src/Visualization/ParticleSystem.cpp
)

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

target_include_directories(sandbox PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include # Apps/Include
)

# Link for sandbox target
# Link dependencies only for executable
target_link_libraries(sandbox PRIVATE
    FluidSimulation
    SFML::Graphics
)

if(ENABLE_TRACY)
    target_compile_definitions(sandbox PRIVATE TRACY_ENABLED)
    target_link_libraries(sandbox PRIVATE TracyClient)
endif()

# Link for demo target
target_include_directories(demo PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include # Apps/Include
)

# Link dependencies only for executable
target_link_libraries(demo PRIVATE
    FluidSimulation
    SFML::Graphics
)

if(ENABLE_TRACY)
    target_compile_definitions(demo PRIVATE TRACY_ENABLED)
    target_link_libraries(demo PRIVATE TracyClient)
endif()