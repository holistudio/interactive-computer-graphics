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
using namespace std;

// Timer
#include <chrono>

#include <iostream>
#include <limits>

// VertexBufferObject wrapper
VertexBufferObject line_VBO;
VertexBufferObject tri_VBO;

// Contains the vertex positions
Eigen::MatrixXf line_V(5,1);
Eigen::MatrixXf tri_V(5,1);

Matrix2f view_scale;
Vector2f view_pos;

// bool tri_first = true;
// bool tri_complete = false;
char mode = ' ';

int click_count = 0;
int num_triangles = 0;
bool mouse_move = false;

Vector2d start_click;

class color
{
    public:

        float r;
        float g;
        float b;
        color()
        {
            r = 0.0;
            g = 0.0;
            b = 0.0;
        }
        color(float r1, float g1, float b1)
        {
            r = r1;
            g = g1;
            b = b1;
        }

};

class point
{
    public:
        float x;
        float y;
        color rgb;
        point()
        {
            x=0.0;
            y=0.0;
            rgb= color();
        }
        point(float x1, float y1, color c1)
        {
            x = x1;
            y = y1;
            rgb = c1;
        }
};
class triangle
{
    public:
        vector<point> v; //vector of vertices
        point center;

        bool clicked = false;
        vector<point> clicked_v;
};

bool tri_clicked = false;
int clicked_index = 0;
int v_clicked = 0;
vector<triangle> triangles; //TODO: maybe not necessary. only one triangle needs to be stored, the one that selected. 

vector<color> colors;


void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

point screen_to_world(GLFWwindow* window)
{
    // Get the position of the mouse in the window
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    Matrix2f inv_scale;
    inv_scale << 1/view_scale.coeff(0,0),0,0,1/view_scale.coeff(1,1);

    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    // Convert screen position to world coordinates
    double xworld = ((xpos/double(width))*2)-1;
    double yworld = (((height-1-ypos)/double(height))*2)-1; 

    Vector2f screen_pos;
    screen_pos << float(xworld), float(yworld);
    screen_pos = inv_scale * (screen_pos - view_pos);

    xworld = double(screen_pos.coeff(0));
    yworld = double(screen_pos.coeff(1));
    return point(xworld, yworld, color());
}

bool click_triangle(point click_point, triangle test_triangle)
{
    point a = test_triangle.v[0];
    point b = test_triangle.v[1];
    point c = test_triangle.v[2];
    float ya_minus_yb = a.y - b.y;
    float xb_minus_xa = b.x - a.x;
    float xayb_minus_xbya = a.x*b.y - b.x*a.y;
    float ya_minus_yc = a.y - c.y;
    float xc_minus_xa = c.x - a.x;
    float xayc_minus_xcya = a.x*c.y - c.x*a.y;

    float gamma = (ya_minus_yb*click_point.x + xb_minus_xa*click_point.y + xayb_minus_xbya) / (ya_minus_yb*c.x + xb_minus_xa*c.y + xayb_minus_xbya);
    
    if(gamma < 1 && gamma >=0)
    {
        float beta = (ya_minus_yc*click_point.x + xc_minus_xa*click_point.y + xayc_minus_xcya) / (ya_minus_yc*b.x + xc_minus_xa*b.y + xayc_minus_xcya);
        if(beta < 1 && beta >=0)
        {
            float alpha = 1-beta-gamma;
            if(alpha < 1 && alpha >0)
            {
                return true;
                //TODO: turn all vertices colors blue
            }
            
        }
    }
    return false;
}

void transform_triangle(GLFWwindow* window, triangle sel_triangle, Matrix2f transform)
{
    vector<Vector2f> new_vertices;
    //find barycentric center
    Vector2f v1 = Vector2f(sel_triangle.v[0].x,sel_triangle.v[0].y);
    Vector2f v2 = Vector2f(sel_triangle.v[1].x,sel_triangle.v[1].y);
    Vector2f v3 = Vector2f(sel_triangle.v[2].x,sel_triangle.v[2].y);
    Vector2f bary_center = (v1 + v2 + v3)/3;

    // subtract the coordinates of the centre from each vertex
    v1 = v1 - bary_center;
    v2 = v2 - bary_center;
    v3 = v3 - bary_center;

    // Then multiply each new vertex by the rotation matrix
    v1 = transform * v1;
    v2 = transform * v2;
    v3 = transform * v3;

    // add back the centre coordinates to each rotated vertex
    v1 = v1 + bary_center;
    v2 = v2 + bary_center;
    v3 = v3 + bary_center;
    new_vertices.push_back(v1);
    new_vertices.push_back(v2);
    new_vertices.push_back(v3);

    for(unsigned i = 0; i<new_vertices.size(); i++)
    {
        triangles[clicked_index].v[i].x = new_vertices[i].coeff(0);
        triangles[clicked_index].v[i].y = new_vertices[i].coeff(1);
        tri_V.col(clicked_index*3+i).coeffRef(0) = new_vertices[i].coeff(0);
        tri_V.col(clicked_index*3+i).coeffRef(1) = new_vertices[i].coeff(1);
    }

    triangles[clicked_index].clicked_v = triangles[clicked_index].v;
    tri_VBO.update(tri_V);

    point world_click = screen_to_world(window);
    start_click << world_click.x, world_click.y;
}

void removeColumn(MatrixXf& matrix, unsigned int colToRemove)
{
    //Eigen matrix column remover function courtesy of https://stackoverflow.com/questions/13290395/how-to-remove-a-certain-row-or-column-while-using-eigen-library-c
    //Thank the gods for StackOverflow
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    point world_click = screen_to_world(window);

    switch (mode)
    {
        case 'i':
            //after the first click, draw a segment from first click to wherever the mouse cursor is
            //(continuously set the last column of the line matrix to mouse cursor position)
            if(click_count>0)
            {

                //set line matrix column to mouse cursor position
                line_V.col(click_count) << world_click.x, world_click.y , 0.0f, 0.0f, 0.0f;

                //update VBO
                line_VBO.update(line_V);

                //mouse is now moving, and the line can now be rendered
                mouse_move = true;
            }
            break;
        case 'm':
            if(tri_clicked)
            {
                //set current mouse position to vector mouse_pos
                Vector2d mouse_pos;
                mouse_pos << world_click.x, world_click.y;
                
                //calculated difference btw start_click and mouse_pos
                Vector2d tr = mouse_pos - start_click;
                //translate all vertices of the clicked triangle
                for(unsigned i = 0; i < triangles[clicked_index].v.size(); i++)
                {
                    triangles[clicked_index].v[i].x = triangles[clicked_index].clicked_v[i].x + tr.coeffRef(0);
                    triangles[clicked_index].v[i].y = triangles[clicked_index].clicked_v[i].y + tr.coeffRef(1);

                    tri_V.col(clicked_index*3+i).coeffRef(0)=triangles[clicked_index].v[i].x;
                    tri_V.col(clicked_index*3+i).coeffRef(1)=triangles[clicked_index].v[i].y;
                }

                //update triangle VBO
                tri_VBO.update(tri_V);
            }
            break;
        default:
            break;
    }
    
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    point world_click = screen_to_world(window);

    // Add mouse click coordinates to V if the left button is pressed
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        //TODO: switch on a mode string variable
        switch (mode)
        {
            case 'i':
            {
                //at every click expand the line matrix by one column
                line_V.conservativeResize(NoChange ,line_V.cols()+1);

                switch (click_count)
                {
                    case  0:
                    {
                        //if it's the first click set first column of line matrix to click position
                        line_V.col(click_count) << world_click.x, world_click.y, 0.0f, 0.0f, 0.0f;
                        click_count = click_count+1;
                        mouse_move=false;
                        break;
                    }
                    case  1:
                    {
                        //on the second click, add click position to line matrix
                        //then draw a line strip using the following line matrix
                        //first click, second click, wherever the mouse cursor is, and first click
                        line_V.col(click_count) << world_click.x, world_click.y, 0.0f, 0.0f, 0.0f;
                        
                        //third vertex in matrix starts as the same as second vertex as a placeholder
                        line_V.col(line_V.cols()-1) << world_click.x, world_click.y, 0.0f, 0.0f, 0.0f;;

                        line_V.conservativeResize(NoChange ,line_V.cols()+1);

                        //last vertex of line strip is same as first vertex of line strip (to just preview the triangle not draw one)
                        line_V.col(line_V.cols()-1) << line_V.col(0);
                        click_count++;
                        break;
                    }
                    case  2:
                    {
                        //add all three click positions to triangle matrix
                        int insert_start = tri_V.cols()-1;
                        tri_V.conservativeResize(NoChange ,tri_V.cols()+3);
                        triangle new_triangle;
                        point new_vertex;

                        for(unsigned i=0; i<3; i++)
                        {
                            line_V.col(i).coeffRef(2) = 1.0f;
                            tri_V.col(insert_start+i) << line_V.col(i);

                            new_vertex.x = line_V.col(i).coeffRef(0);
                            new_vertex.y = line_V.col(i).coeffRef(1);
                            new_triangle.v.push_back(new_vertex);
                        }

                        triangles.push_back(new_triangle);

                        //update triangle VBO
                        tri_VBO.update(tri_V);

                        //increment triangle count
                        num_triangles++;

                        //clear line matrix
                        line_V = MatrixXf::Zero(5,1);
                        //reset click count
                        click_count=0;
                        break;
                    }
                    default:
                        break;
                }
                //update line VBO
                line_VBO.update(line_V);
                break;
            }
            case 'm':
            {
                if(click_count>0)
                {
                    //reset color to red
                    //TODO: revert to original colors before click
                    //Set vertex colors to blue
                    for(unsigned i = 0; i<3; i++)
                    {
                        triangles[clicked_index].v[i].rgb.r = 1.0f;
                        triangles[clicked_index].v[i].rgb.g = 0.0f;
                        triangles[clicked_index].v[i].rgb.b = 0.0f;

                        tri_V.col(clicked_index*3+i).coeffRef(2) = 1.0f;
                        tri_V.col(clicked_index*3+i).coeffRef(3) = 0.0f;
                        tri_V.col(clicked_index*3+i).coeffRef(4) = 0.0f;
                    }
                    tri_VBO.update(tri_V);

                    tri_clicked = false;
                    clicked_index = 0;
                    click_count = 0;
                }
                else
                {
                    //check if click coordinates is in any triangles

                    for(unsigned i=0; i<triangles.size(); i++)
                    {
                        if(click_triangle(world_click,triangles[i]))
                        {
                            tri_clicked = true;
                            clicked_index = i;
                            start_click << world_click.x, world_click.y;
                            triangles[i].clicked_v = triangles[i].v;
                            
                            click_count++;
                        }
                    }
                    //Set vertex colors to blue
                    for(unsigned i = 0; i<3; i++)
                    {
                        triangles[clicked_index].v[i].rgb.r = 0.0f;
                        triangles[clicked_index].v[i].rgb.g = 0.0f;
                        triangles[clicked_index].v[i].rgb.b = 1.0f;

                        tri_V.col(clicked_index*3+i).coeffRef(2) = 0.0f;
                        tri_V.col(clicked_index*3+i).coeffRef(3) = 0.0f;
                        tri_V.col(clicked_index*3+i).coeffRef(4) = 1.0f;
                    }
                    
                    tri_VBO.update(tri_V);
                }
                break;
            }
            case 'd':
            {
                //check if click coordinates is in any triangles

                for(unsigned i=0; i<triangles.size(); i++)
                {
                    if(click_triangle(world_click,triangles[i]))
                    {
                        triangles.erase(triangles.begin()+i);
                        removeColumn(tri_V,i*3);
                        removeColumn(tri_V,i*3);
                        removeColumn(tri_V,i*3);
                        i--;
                    }
                }
                tri_VBO.update(tri_V);
                break;
            }
            case 'c':
            {
                float v_dist; 
                float min_dist = numeric_limits<float>::infinity(); //set to infinity
                
                //for all vertices in tri_V
                for(unsigned i=0; i<tri_V.cols()-1; i++)
                {
                    //calculate distance between mouse click and vertex, v_dist
                    v_dist = (tri_V.col(i).coeff(0) - world_click.x)*(tri_V.col(i).coeff(0) - world_click.x) + (tri_V.col(i).coeff(1) - world_click.y)*(tri_V.col(i).coeff(1) - world_click.y);

                    //if distance is less than min_dist
                    if(v_dist < min_dist)
                    {
                        //set to v_dist as new min_dist
                        min_dist = v_dist;

                        //store index of vertex
                        v_clicked = i;
                    }
                }

                //change closest vertex color to blue
                tri_V.col(v_clicked).coeffRef(2) = 0.0f;
                tri_V.col(v_clicked).coeffRef(3) = 0.0f;
                tri_V.col(v_clicked).coeffRef(4) = 1.0f;
                tri_VBO.update(tri_V);
                click_count++;
                break;
            }
            default:
                break;
        }
    }
    // The follow commented out section of code is meant only for Task 1.1
    // where releasing the left mouse button "releases" the clicked triangle.
    // For Task 1.2 onwards, the triangle remains "selected" after mouse button release
    // to allow for rotation and scaling
    // else
    // {
    //     if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    //     {
    //         if(tri_move_mode)
    //         {
    //             tri_clicked = false;
    //             clicked_index = 0;
    //         } 
    //     }
    // }
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
        // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    //only perform action on key press, not key release
    if(action == GLFW_PRESS)
    {
        //TODO: re-organize switch and boolean if statements
        switch (key)
        {
            case  GLFW_KEY_I:
                mode='i';
                click_count = 0;
                std::cout << "Triangle Insert Mode" << std::endl;
                break;
            case  GLFW_KEY_O:
                mode='m';
                click_count = 0;
                std::cout << "Triangle Move Mode" << std::endl;
                break;
            case  GLFW_KEY_P:
                mode='d';
                click_count = 0;
                std::cout << "Triangle Delete Mode" << std::endl;
                break;
            case  GLFW_KEY_C:
                mode = 'c';
                click_count = 0; 
                std::cout << "Vertex Color Mode" << std::endl;
                break;
            case  GLFW_KEY_EQUAL: //TODO: see if you can combine this with SHIFT key so it's actually '+' sign
                view_scale = view_scale + 0.2*MatrixXf::Identity(2,2);
                break;
            case  GLFW_KEY_MINUS:
                view_scale = view_scale - 0.2*MatrixXf::Identity(2,2);
                break;
            case  GLFW_KEY_W:
                view_pos = view_pos - (0.2*2/view_scale.coeff(1,1))*Vector2f::UnitY();
                break;
            case  GLFW_KEY_S:
                view_pos = view_pos + (0.2*2/view_scale.coeff(1,1))*Vector2f::UnitY();
                break;
            default:
                break;
        }
        if(tri_clicked)
        {
            switch (key)
            {
                case  GLFW_KEY_H:
                {
                    //rotation matrix
                    float radians = -10 * 3.141592f / 180;
                    Matrix2f rotation;
                    rotation << cos(radians), sin(radians), -sin(radians), cos(radians);
                    transform_triangle(window, triangles[clicked_index], rotation);
                    break;
                }
                case  GLFW_KEY_J:
                {
                    //rotation matrix
                    float radians = 10 * 3.141592f / 180;
                    Matrix2f rotation;
                    rotation << cos(radians), sin(radians), -sin(radians), cos(radians);
                    transform_triangle(window, triangles[clicked_index], rotation);
                    break;
                }
                case  GLFW_KEY_K:
                {
                    Matrix2f scale;
                    scale << 1.25, 0, 0, 1.25;
                    transform_triangle(window, triangles[clicked_index], scale);
                    break;
                }
                case  GLFW_KEY_L:
                {
                    Matrix2f scale;
                    scale << 0.75, 0, 0, 0.75;
                    transform_triangle(window, triangles[clicked_index], scale);
                    break;
                }
                default:
                    break;
            }
        }
        if(mode=='c')
        {
            if(click_count>0)
            {
                int color_index = 0;
                switch (key)
                {
                    case  GLFW_KEY_1:
                        color_index = 0;
                        break;
                    case  GLFW_KEY_2:
                        color_index = 1;
                        break;
                    case  GLFW_KEY_3:
                        color_index = 2;
                        break;
                    case  GLFW_KEY_4:
                        color_index = 3;
                        break;
                    case  GLFW_KEY_5:
                        color_index = 4;
                        break;
                    case  GLFW_KEY_6:
                        color_index = 5;
                        break;
                    case  GLFW_KEY_7:
                        color_index = 6;
                        break;
                    case  GLFW_KEY_8:
                        color_index = 7;
                        break;
                    case  GLFW_KEY_9:
                        color_index = 8;
                        break;
                    default:
                        break;
                }

                tri_V.col(v_clicked).coeffRef(2) = colors[color_index].r;
                tri_V.col(v_clicked).coeffRef(3) = colors[color_index].g;
                tri_V.col(v_clicked).coeffRef(4) = colors[color_index].b;

                tri_VBO.update(tri_V);
            }
        }
    }
}

int main(void)
{
    colors.push_back(color(0.0f,0.0f,1.0f));
    colors.push_back(color(0.0f,float(195./255.),1.0f));
    colors.push_back(color(0.0f,1.0f,float(187./255.)));
    colors.push_back(color(float(149./255.),1.0f,0.));
    colors.push_back(color(1.0f,float(242./255.),0.));
    colors.push_back(color(1.0f,float(106./255.),0.));
    colors.push_back(color(1.0f,0.,float(81./255.)));
    colors.push_back(color(1.0f,0.,float(234./255.)));
    colors.push_back(color(float(140./255.),0.,float(227./255.)));

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
    line_V.resize(5,1);
    line_VBO.update(line_V);

    tri_VBO.init();
    tri_V.resize(5,1);
    tri_VBO.update(tri_V);

    view_scale << 1, 0, 0, 1;
    view_pos << 0, 0;

    // Initialize the OpenGL Program
    // A program controls the OpenGL pipeline and it must contains
    // at least a vertex shader and a fragment shader to be valid
    Program program;
    const GLchar* vertex_shader =
            "#version 150 core\n"
                    "in vec2 position;"
                    "in vec3 inColor;"
                    "uniform mat2 scale;"
                    "uniform vec2 translation;"
                    "out vec3 vertexColor;"
                    "void main()"
                    "{"
                    "    vertexColor = inColor;"
                    "    gl_Position = vec4(scale * position + translation, 0.0, 1.0);"
                    "}";
    const GLchar* fragment_shader =
            "#version 150 core\n"
                    "in vec3 vertexColor;"
                    "out vec4 fragmentColor;"
                    "void main()"
                    "{"
                    "    fragmentColor = vec4(vertexColor, 1.0);"
                    "}";

    // Compile the two shaders and upload the binary to the GPU
    program.init(vertex_shader,fragment_shader,"fragmentColor");
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

        // Set the uniform view matrix and translation vectors
        glUniformMatrix2fv(program.uniform("scale"),1, GL_FALSE, view_scale.data());
        glUniform2fv(program.uniform("translation"),1, view_pos.data());
        
        // Clear the framebuffer
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);


        // Draw all complete triangles
        //TODO: remove num_triangles from code
        if(num_triangles>0)
        {
            tri_VBO.bind();
            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), 0);

            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), (const GLvoid *)8);
            for(unsigned i=0; i<triangles.size(); i++)
            {
                //draw triangles
                glDrawArrays(GL_TRIANGLES, i*3, 3);
            }

        }

        
        //if a line is being drawn
        if(mode == 'i')
        {
            line_VBO.bind();

            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), (const GLvoid *)8);

            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), 0);
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
