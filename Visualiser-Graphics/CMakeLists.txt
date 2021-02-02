add_library(
    gfx STATIC
    src/MainGFX/Graphics.hpp
    src/MainGFX/Graphics.cpp
)

target_include_directories(gfx PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")  