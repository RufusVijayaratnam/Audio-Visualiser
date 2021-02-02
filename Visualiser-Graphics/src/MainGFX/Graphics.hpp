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



namespace gfx {
    bool InitialiseGLFW();
    GLFWwindow* OpenWindow(const char * windowName, bool &windowOpened);
    void MainLoop(GLFWwindow* window);
}
#endif