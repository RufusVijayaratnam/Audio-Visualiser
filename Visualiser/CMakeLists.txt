add_executable(Visualise main.cpp)


find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)

INCLUDE(FindPkgConfig)

PKG_SEARCH_MODULE(SDL2 REQUIRED sdl2)
PKG_SEARCH_MODULE(SDL2MIXER REQUIRED SDL2_mixer>=2.0.0)

INCLUDE_DIRECTORIES(${SDL2_INCLUDE_DIRS} ${SDL2IMAGE_INCLUDE_DIRS})


target_link_libraries(Visualise PRIVATE Audio)
target_link_libraries(Visualise PRIVATE gfx)
target_link_libraries(Visualise PRIVATE ${SDL2_LIBRARIES} ${SDL2MIXER_LIBRARIES})
target_link_libraries(Visualise PRIVATE GLEW glfw OpenGL::GL)