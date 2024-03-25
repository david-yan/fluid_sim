cmake_minimum_required(VERSION 3.16)
project(FluidSim LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

include(FetchContent)
FetchContent_Declare(SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 2.6.x)
FetchContent_MakeAvailable(SFML)

include_directories(header)

add_library(Vector2f src/Vector2f.cpp)
add_executable(FluidSim src/main.cpp)
target_link_libraries(Vector2f PRIVATE sfml-graphics)
target_link_libraries(FluidSim PRIVATE sfml-graphics Vector2f)
target_compile_features(FluidSim PRIVATE cxx_std_17)

if(WIN32)
    add_custom_command(
        TARGET FluidSim
        COMMENT "Copy OpenAL DLL"
        PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${SFML_SOURCE_DIR}/extlibs/bin/$<IF:$<EQUAL:${CMAKE_SIZEOF_VOID_P},8>,x64,x86>/openal32.dll $<TARGET_FILE_DIR:FluidSim>
        VERBATIM)
endif()

install(TARGETS FluidSim)
