add_library(
    gfx STATIC
    src/MainGFX/Graphics.hpp
    src/MainGFX/Graphics.cpp
    src/MainGFX/shader.hpp
    src/MainGFX/shader.cpp
    src/Shaders/SimpleVertexShader.vertexshader
    src/Shaders/SimpleFragmentShader.fragmentshader
)

target_include_directories(gfx PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")  