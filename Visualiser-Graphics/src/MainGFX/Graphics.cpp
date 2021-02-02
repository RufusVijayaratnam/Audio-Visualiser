//
//  Graphics.cpp
//  Visualiser
//
//  Created by Rufus Vijayaratnam on 21/01/2021.
//  Copyright Â© 2021 Rufus Vijayaratnam. All rights reserved.
//

#include "Graphics.hpp"
//extern GLFWwindow * window;

bool gfx::InitialiseGLFW() {
    // Initialise GLFW
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        getchar();
        return false;
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    return true;
}

GLFWwindow* gfx::OpenWindow(const char * windowName, bool &windowOpened) {
    GLFWwindow* window;
    window = glfwCreateWindow(1024, 720, windowName, NULL, NULL);
    std::cout << "here window is: " << window << std::endl;

    if (window == NULL) {
        printf("Could not open a window, terminating... \n");
        //getchar();
        glfwTerminate();
        windowOpened = false;
        return window;
    }
    glfwMakeContextCurrent(window);

    glewExperimental = true;
    if (glewInit() != GLEW_OK) {
        printf("It is not ok... \n");
        getchar();
        glfwTerminate();
        windowOpened = false;
        return window;
    }

    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwPollEvents();
    glfwSetCursorPos(window, 1024 / 2, 720 / 2);

    //Set coloured background, Disable if some fancier method desired.
    glClearColor(0.0f, 0.0f, 0.4f, 0.0f);
    windowOpened = true;
    return window;
}

void gfx::MainLoop(GLFWwindow* window) {
    do{
    // Clear the screen. It's not mentioned before Tutorial 02, but it can cause flickering, so it's there nonetheless.
    glClear( GL_COLOR_BUFFER_BIT );

    // Draw nothing, see you in tutorial 2 !

    // Swap buffers
    glfwSwapBuffers(window);
    glfwPollEvents();

    } while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
       glfwWindowShouldClose(window) == 0 );
}