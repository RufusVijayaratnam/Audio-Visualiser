add_executable(Visualise main.cpp)

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)
target_link_libraries(Visualise PRIVATE Audio)
target_link_libraries(Visualise PRIVATE gfx)
target_link_libraries(Visualise PRIVATE GLEW glfw OpenGL::GL)