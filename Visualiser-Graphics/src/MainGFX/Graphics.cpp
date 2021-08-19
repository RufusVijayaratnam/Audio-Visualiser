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
    window = glfwCreateWindow(1080, 720, windowName, NULL, NULL);
    std::cout << "here window is: " << window << std::endl;

    if (window == NULL) {
        printf("Could not open a window, terminating... \n");
        //getchar();
        glfwTerminate();
        windowOpened = false;
        return window;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(0);

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
    //glfwSetCursorPos(window, 1024 / 2, 720 / 2);

    //Set coloured background, Disable if some fancier method desired.
    glClearColor(0.0f, 0.0f, 0.4f, 0.0f);
    windowOpened = true;
    return window;
}

void gfx::Main(GLFWwindow* window, std::vector<std::vector<double>> &spectrumMagnitudes, double &frameTime, Mix_Chunk* sound) {

    int channel;
    channel = Mix_PlayChannel(-1, sound, 0); 
    if(channel == -1) {
        fprintf(stderr, "Unable to play WAV file: %s\n", Mix_GetError()); 
    }
    
    GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

    const char* vertexPath = "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser-Graphics/src/Shaders/SimpleVertexShader.vertexshader";
    const char* fragmentPath = "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser-Graphics/src/Shaders/SimpleFragmentShader.fragmentshader";

    GLuint programID = LoadShaders("/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser-Graphics/src/Shaders/SimpleVertexShader.vertexshader", "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser-Graphics/src/Shaders/SimpleFragmentShader.fragmentshader");
    //GLuint programID = LoadShaders(vertexPath, fragmentPath);

	GLuint vertexbuffer;
	glGenBuffers(1, &vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    
    glfwSetTime(0);
    int currentFrame = 0;
    bool began = false;
    double prevTime = 0;
    double currentTime = 0;

    const double waitTime = 0.45;
    //int currentFrame;
    do  {

        if (Mix_Playing(channel) == 1) {
        //printf("current frame is: %i", currentFrame);

        if (!began) {
            began = true;
            glfwSetTime(0);
            while (waitTime > glfwGetTime()) {
                SDL_Delay(1);
            }
            glfwSetTime(0);
        }

        
        
        std::vector<glm::vec3> g_vertex_buffer_data = gfx::GetFrameVertices(spectrumMagnitudes[currentFrame]);
        //glBufferData(GL_ARRAY_BUFFER, g_vertex_buffer_data.size() * sizeof(glm::vec3), float(g_vertex_buffer_data), GL_STATIC_DRAW);
        glBufferData(GL_ARRAY_BUFFER, g_vertex_buffer_data.size() * sizeof(glm::vec3), &g_vertex_buffer_data[0], GL_STATIC_DRAW);

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

        glDrawArrays(GL_TRIANGLES, 0, g_vertex_buffer_data.size() * 3);
        glDisableVertexAttribArray(0);

        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();

        double timePassed = glfwGetTime();
        while (timePassed < frameTime) {
            SDL_Delay(1);
            timePassed = glfwGetTime();
        }
        currentFrame++;
        glfwSetTime(0);
        
        //currentFrame++; 
        } else {
            //SDL_Delay(1);
            glfwSetTime(0);
        }
    

    } while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS);

    glDeleteBuffers(1, &vertexbuffer);
	glDeleteVertexArrays(1, &VertexArrayID);
	glDeleteProgram(programID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

}

// XYZ convention
std::vector<glm::vec3> gfx::GetFrameVertices(std::vector<double> &magnitudes) {
    std::vector<glm::vec3> frameVertices;
    int nBars = magnitudes.size();
    //Every bar will have six vertices, each of those with 3 coordinates, z is always 0
    //Length of frameVertices will be nBars * 6 * 3
    const float maxBarWidth = 2.0 / nBars;
    const float spacing = 0.5 * maxBarWidth;
    const float barWidth = spacing;
    for (int i = 0; i < nBars; i++) { //For every bar
        //Triangle points top to bottom, left to right
        const float y = float(magnitudes[i]);
        const float x1 = spacing + spacing * i + i * barWidth - 1.0;
        const float x2 = x1 + barWidth;
        const float z = 0.0f;
        const float bottom = -1.0f;

        const glm::vec3 p1 = {x1, y, z};
        const glm::vec3 p2 = {x1, bottom, z};
        const glm::vec3 p3 = {x2, bottom, z};
        const glm::vec3 p4 = p1;
        const glm::vec3 p5 = {x2, y, z};
        const glm::vec3 p6 = p3;

        frameVertices.push_back(p1);
        frameVertices.push_back(p2);
        frameVertices.push_back(p3);
        frameVertices.push_back(p4);
        frameVertices.push_back(p5);
        frameVertices.push_back(p6);
    }
    return frameVertices;
}

