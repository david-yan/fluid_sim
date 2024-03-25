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

add_library(vector2f src/vector2f.cpp)
add_library(circle src/circle.cpp)
add_executable(fluid_sim src/main.cpp)
target_link_libraries(vector2f PRIVATE sfml-graphics)
target_link_libraries(circle PRIVATE sfml-graphics)
target_link_libraries(fluid_sim PRIVATE sfml-graphics vector2f circle)
target_compile_features(fluid_sim PRIVATE cxx_std_17)

if(WIN32)
    add_custom_command(
        TARGET fluid_sim
        COMMENT "Copy OpenAL DLL"
        PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${SFML_SOURCE_DIR}/extlibs/bin/$<IF:$<EQUAL:${CMAKE_SIZEOF_VOID_P},8>,x64,x86>/openal32.dll $<TARGET_FILE_DIR:fluid_sim>
        VERBATIM)
endif()

install(TARGETS fluid_sim)
