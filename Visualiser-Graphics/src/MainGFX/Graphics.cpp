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

void gfx::Main(GLFWwindow* window) {

    GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

    const char* vertexPath = "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser-Graphics/src/Shaders/SimpleVertexShader.vertexshader";
    const char* fragmentPath = "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser-Graphics/src/Shaders/SimpleFragmentShader.fragmentshader";

    GLuint programID = LoadShaders("/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser-Graphics/src/Shaders/SimpleVertexShader.vertexshader", "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser-Graphics/src/Shaders/SimpleFragmentShader.fragmentshader");
    //GLuint programID = LoadShaders(vertexPath, fragmentPath);


	static const GLfloat g_vertex_buffer_data[] = { 
		-1.0f, -1.0f, 0.0f,
		 1.0f, -1.0f, 0.0f,
		 0.0f,  1.0f, 0.0f,
	};

	GLuint vertexbuffer;
	glGenBuffers(1, &vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

    do{
    // Clear the screen. It's not mentioned before Tutorial 02, but it can cause flickering, so it's there nonetheless.
    glClear( GL_COLOR_BUFFER_BIT );

    glUseProgram(programID);

    // Draw nothing, see you in tutorial 2 !
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(
        0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
		3,                  // size
		GL_FLOAT,           // type
		GL_FALSE,           // normalized?
		0,                  // stride
		(void*)0            // array buffer offset
    );

    glDrawArrays(GL_TRIANGLES, 0, 3);
    glDisableVertexAttribArray(0);

    // Swap buffers
    glfwSwapBuffers(window);
    glfwPollEvents();

    } while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
       glfwWindowShouldClose(window) == 0 );

    glDeleteBuffers(1, &vertexbuffer);
	glDeleteVertexArrays(1, &VertexArrayID);
	glDeleteProgram(programID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

}