//
//  Graphics.hpp
//  Visualiser
//
//  Created by Rufus Vijayaratnam on 21/01/2021.
//  Copyright Â© 2021 Rufus Vijayaratnam. All rights reserved.
//

#ifndef Graphics_hpp
#define Graphics_hpp
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <stdio.h>
#include <iostream>
#include "shader.hpp"
#include <string>
#include <vector>
#include <cmath>
#include <SDL2/SDL.h>
#include <SDL2/SDL_mixer.h>




namespace gfx {
    bool InitialiseGLFW();
    GLFWwindow* OpenWindow(const char * windowName, bool &windowOpened);
    void Main(GLFWwindow* window, std::vector<std::vector<double>> &spectrumMagnitudes, double &frameTime, Mix_Chunk* sound);
    std::vector<glm::vec3> GetFrameVertices(std::vector<double> &magnitudes);
}
#endif