// This example is heavily based on the tutorial at https://open.gl

// OpenGL Helpers to reduce the clutter
#include "Helpers.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#else
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#endif


// Linear Algebra Library
#include <Eigen/Core>
using namespace Eigen;

// Timer
#include <chrono>

#include <iostream>

// VertexBufferObject wrapper
VertexBufferObject line_VBO;
VertexBufferObject tri_VBO;

// Contains the vertex positions
Eigen::MatrixXf line_V(2,1);
Eigen::MatrixXf tri_V(2,1);

// bool tri_first = true;
// bool tri_complete = false;
bool tri_insert_mode = false;
int click_count = 0;
int num_triangles = 0;
bool mouse_move = false;

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    //The callback functions receives the cursor position
    //measured in screen coordinates but relative to the top-left corner of the window content area.
    if(tri_insert_mode == true)
    {
        //after the first click, draw a segment from first click to wherever the mouse cursor is
        //(continuously set the last column of the line matrix to mouse cursor position)
        if(click_count>0)
        {
            // Get the size of the window
            int width, height;
            glfwGetWindowSize(window, &width, &height);

            // Convert screen position to world coordinates
            double xworld = ((xpos/double(width))*2)-1;
            double yworld = (((height-1-ypos)/double(height))*2)-1; // NOTE: y axis is flipped in glfw

            //set line matrix column to mouse cursor position
            line_V.col(click_count) << xworld, yworld;

            //update VBO
            line_VBO.update(line_V);

            //mouse is now moving, and the line can now be rendered
            mouse_move = true;
        }

    }
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if(tri_insert_mode == true)
    {
        // Get the position of the mouse in the window
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        // Get the size of the window
        int width, height;
        glfwGetWindowSize(window, &width, &height);

        // Convert screen position to world coordinates
        double xworld = ((xpos/double(width))*2)-1;
        double yworld = (((height-1-ypos)/double(height))*2)-1; // NOTE: y axis is flipped in glfw

        // Add mouse click coordinates to V if the left button is pressed
        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
        {
            //at every click expand the line matrix by one column
            line_V.conservativeResize(NoChange ,line_V.cols()+1);

            switch (click_count)
            {
                case  0:
                    //if it's the first click set first column of line matrix to click position
                    line_V.col(click_count) << xworld, yworld;
                    click_count = click_count+1;
                    mouse_move=false;
                    break;
                case 1:
                    //on the second click, add click position to line matrix
                    //then draw a line strip using the following line matrix
                    //first click, second click, wherever the mouse cursor is, and first click
                    line_V.col(click_count) << xworld, yworld;
                    
                    //third vertex in matrix starts as the same as second vertex as a placeholder
                    line_V.col(line_V.cols()-1) << xworld, yworld;

                    line_V.conservativeResize(NoChange ,line_V.cols()+1);

                    //last vertex of line strip is same as first vertex of line strip (to just preview the triangle not draw one)
                    line_V.col(line_V.cols()-1) << line_V.col(0);
                    click_count++;
                    break;
                case 2:
                    //add all three click positions to triangle matrix
                    int insert_start = tri_V.cols()-1;
                    tri_V.conservativeResize(NoChange ,tri_V.cols()+3);
                    for(unsigned i=0; i<3; i++)
                    {
                        tri_V.col(insert_start+i) << line_V.col(i);
                    }

                    //update triangle VBO
                    tri_VBO.update(tri_V);

                    //increment triangle count
                    num_triangles++;

                    //clear line matrix
                    line_V = Eigen::MatrixXf::Zero(2,1),
                    //reset click count
                    click_count=0;
                    break;
            }   
        }
        //update line VBO
        // Upload the change to the GPU
        line_VBO.update(line_V);
    }
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // Update the position of the first vertex if the keys 1,2, or 3 are pressed
    switch (key)
    {
        case  GLFW_KEY_I:
            tri_insert_mode = true;
            break;
        default:
            break;
    }

    // Upload the change to the GPU
    line_VBO.update(line_V);
}

int main(void)
{
    GLFWwindow* window;

    // Initialize the library
    if (!glfwInit())
        return -1;

    // Activate supersampling
    glfwWindowHint(GLFW_SAMPLES, 8);

    // Ensure that we get at least a 3.2 context
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);

    // On apple we have to load a core profile with forward compatibility
    #ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif

    // Create a windowed mode window and its OpenGL context
    window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    #ifndef __APPLE__
      glewExperimental = true;
      GLenum err = glewInit();
      if(GLEW_OK != err)
      {
        /* Problem: glewInit failed, something is seriously wrong. */
       fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
      }
      glGetError(); // pull and savely ignonre unhandled errors like GL_INVALID_ENUM
      fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
    #endif

    int major, minor, rev;
    major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
    minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
    rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
    printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
    printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
    printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));

    // Initialize the VAO
    // A Vertex Array Object (or VAO) is an object that describes how the vertex
    // attributes are stored in a Vertex Buffer Object (or VBO). This means that
    // the VAO is not the actual object storing the vertex data,
    // but the descriptor of the vertex data.
    VertexArrayObject line_VAO;
    line_VAO.init();
    line_VAO.bind();

    // VertexArrayObject tri_VAO;
    // tri_VAO.init();
    // tri_VAO.bind();

    // Initialize the VBO with the vertices data
    // A VBO is a data container that lives in the GPU memory
    line_VBO.init();
    line_V.resize(2,1);
    line_VBO.update(line_V);

    tri_VBO.init();
    tri_V.resize(2,1);
    tri_VBO.update(tri_V);

    // Initialize the OpenGL Program
    // A program controls the OpenGL pipeline and it must contains
    // at least a vertex shader and a fragment shader to be valid
    Program program;
    const GLchar* vertex_shader =
            "#version 150 core\n"
                    "in vec2 position;"
                    "void main()"
                    "{"
                    "    gl_Position = vec4(position, 0.0, 1.0);"
                    "}";
    const GLchar* fragment_shader =
            "#version 150 core\n"
                    "out vec4 outColor;"
                    "uniform vec3 triangleColor;"
                    "void main()"
                    "{"
                    "    outColor = vec4(triangleColor, 1.0);"
                    "}";

    // Compile the two shaders and upload the binary to the GPU
    // Note that we have to explicitly specify that the output "slot" called outColor
    // is the one that we want in the fragment buffer (and thus on screen)
    program.init(vertex_shader,fragment_shader,"outColor");
    program.bind();

    // The vertex shader wants the position of the vertices as an input.
    // The following line connects the VBO we defined above with the position "slot"
    // in the vertex shader
    //program.bindVertexAttribArray("position",line_VBO);

    // Save the current time --- it will be used to dynamically change the triangle color
    auto t_start = std::chrono::high_resolution_clock::now();

    // Register the keyboard callback
    glfwSetKeyCallback(window, key_callback);

    //Register the mouse cursor position callback
    glfwSetCursorPosCallback(window, cursor_position_callback);

    // Register the mouse callback
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Update viewport
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window))
    {
        line_VAO.bind();

        // Bind your program
        program.bind();

        // Set the uniform value depending on the time difference
        auto t_now = std::chrono::high_resolution_clock::now();
        float time = std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
        // glUniform3f(program.uniform("triangleColor"), (float)(sin(time * 4.0f) + 1.0f) / 2.0f, 0.0f, 0.0f);

        // Clear the framebuffer
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);


        // Draw all complete triangles
        if(num_triangles>0)
        {
            //bind triangle VAO
            //tri_VAO.bind();
            program.bindVertexAttribArray("position",tri_VBO);
            //draw triangles
            glDrawArrays(GL_TRIANGLES, 0, 3*num_triangles);
            //tri_VAO.free();
        }

        
        //if a line is being drawn
        if(tri_insert_mode)
        {
            program.bindVertexAttribArray("position",line_VBO);
            if(click_count>0)
            {
                if(mouse_move)
                {
                    glDrawArrays(GL_LINE_STRIP,0,line_V.cols());
                }
                else
                {
                    glDrawArrays(GL_LINE_STRIP,0,line_V.cols()-1);
                }   
            }
            
        }

        // Swap front and back buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();
    }

    // Deallocate opengl memory
    program.free();

    line_VBO.free();
    tri_VBO.free();

    // Deallocate glfw internals
    glfwTerminate();
    return 0;
}
