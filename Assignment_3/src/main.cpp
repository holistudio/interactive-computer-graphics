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
#include <Eigen/Geometry>
using namespace Eigen;
using namespace std;
#include <chrono>
#include <iostream>
#include <cmath>
#include <fstream>
#include <limits>

// VertexBufferObject wrapper
VertexBufferObject line_VBO;
VertexBufferObject tri_VBO;
VertexBufferObject mesh_VBO;

// Viewing Transformation Matrices
float near = -0.5;
float far = -5;
float l = -1.0;
float r = 1.0;
float t = 1.0;
float b = -1.0;
Matrix4f M_vp;

Matrix4f M_orth;
Matrix4f P;
Matrix4f M_proj;

Matrix4f M_cam;
Matrix4f M_model;

Vector3f eye_pos(0,0,1.5);

// Contains the vertex positions of the lines and triangles
MatrixXf line_V(5,1);
MatrixXf tri_V(5,1);
MatrixXf mesh_V(6,1);

// variables for tracking which mode the drawing application is in
char mode = ' ';
bool animation_mode = false;

// variables for tracking mouse states
Vector2d start_click;
int click_count = 0;
bool mouse_move = false;

//View control variables that link to shader uniforms
Matrix2f view_scale;
Vector2f view_pos;

//Light Class
class light
{
    public:
        Vector3f position;
        Vector3f color;
        Vector3f intensity;
};


// Container for tracking the triangle vertex positions and colors at each key frame of an animation
Eigen::MatrixXf anim_tri_V(5,1);
// Number of triangle vertices throughout the animation
int num_anim_V=0;
// Number of animation keyframes
int num_keyframes = 0;
// Time step between animation key frames
float time_step = 1.0; 
// Variable for recording when the animation starts
chrono::time_point<chrono::high_resolution_clock> t_start;

// Color, Point, and Triangle Classes for easier-to-read code
class color
{
    public:

        float r;
        float g;
        float b;
        color()
        {
            r = 0.0; //values between 0 and 1 for convenient use with OpenGL shaders
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
        point()
        {
            x=0.0;
            y=0.0;
        }
        point(float x1, float y1)
        {
            x = x1;
            y = y1;
        }
};

class triangle
{
    public:
        vector<point> v; //vector of vertices
        bool clicked = false; // tracks whether there is a triangle selected on the screen or not
        int clicked_index; //start index on tri_V matrix
};

struct Vertex 
{
  Vector3f position;
  Vector3f normal;
  Vector3f color;
};

class tri_mesh
{
    //basic triangular mesh
    public:
        //V is a matrix (dimension 6 x #V where #V is the number of vertices)
        // contains the positions and normals of the vertices of the mesh 
        MatrixXf V;
        //F is an matrix (dimension 6 x #F where #F is the number of faces) which contains the descriptions of the triangles in the mesh. 
        MatrixXf F;
        Vector3f diff_color;
        Vector3f ka;
        Vector3f ks;
        char shader_type;
        float phong_exp;
        bool mirror=false;
};

// global variable for storing the clicked triangle's properties
triangle clicked_triangle;
tri_mesh clicked_mesh;

// global variable for storing the clicked vertex's index position in tri_V
int v_clicked = 0;

// vector of colors for shading vertices in Task 1.3
vector<color> colors;

vector<tri_mesh> meshes;

vector<string> space_sep_string(string line)
{
    //separates string into a vector of substring as space dividers
    vector<int> split_index;
    vector<string> substrings;
    for (unsigned i=0;i<line.size();i++)
    {
        if(line[i]==' ')
        {
            split_index.push_back(i);
        }
    }

    int split_start = 0;
    int substr_len;
    for (unsigned i=0;i<split_index.size()+1;i++)
    {
        if(i<split_index.size())
        {
            substr_len= split_index[i] - split_start;
        }
        else
        {
            substr_len= line.size() - split_start;
        }
        substrings.push_back(line.substr(split_start,substr_len));
        if(i<split_index.size())
        {
            split_start = split_index[i]+1;
        }
        
    }
    return substrings;
}

tri_mesh load_mesh(string off_filepath, Vector3d position, double scale)
{
    // returns a triangule mesh object at the file path for the OFF file
    tri_mesh mesh_structure;
    mesh_structure.diff_color = Vector3f(255,255,255);
    mesh_structure.ka = mesh_structure.diff_color * .15;
    mesh_structure.ks =  Vector3f(0.0,0.0,0.0);
    mesh_structure.phong_exp = 1.0;
    mesh_structure.shader_type = 'f';

    // mesh vertices
    int num_vertices;
    // mesh faces
    int num_faces;
    ifstream in(off_filepath);

    char str[255];
    bool summary_line_read = false;
    vector<string> substrings;

    //read the first line
    in.getline(str,255);

    //if the first line says 'OFF' continue
    if(string(str).compare("OFF")==0)
    {  
        std::cout << "Valid OFF File" << std::endl;
        while (summary_line_read == false && in.getline(str,255)) 
        {
            if(str[0]!='#')
            {
                //then read the next line that doesn't start with '#'
                
                substrings = space_sep_string(string(str));
                // std::cout << stoi(substrings[1]) << std::endl;
                num_vertices = stoi(substrings[0]);
                num_faces = stoi(substrings[1]);
                summary_line_read = true;
            }
        }

        //vertices.resize(num_vertices,3);
        mesh_structure.V.resize(6,num_vertices);
        mesh_structure.F.resize(6,num_faces);



        //load matrix V with vertices
        for (unsigned i=0;i<num_vertices;i++)
        {
            in.getline(str,255);
            substrings = space_sep_string(string(str));
            mesh_structure.V.col(i) << stof(substrings[0])*scale+position.coeff(0),
            stof(substrings[1])*scale+position.coeff(1),
            stof(substrings[2])*scale+position.coeff(2),
            0.0f, 0.0f, 0.0f;
        }

        //load matrix F with vertices of each face
        for (unsigned i=0;i<num_faces;i++)
        {
            in.getline(str,255);
            substrings = space_sep_string(string(str));
            
            //calculate face normals un-normalized
            //which will weight the vertex normal average based on each face's area
            Vector3f v1 = (mesh_structure.V.block(0,stoi(substrings[1]),3,1));
            Vector3f v2 = (mesh_structure.V.block(0,stoi(substrings[2]),3,1));
            Vector3f v3 = (mesh_structure.V.block(0,stoi(substrings[3]),3,1));

            Vector3f face_normal = (v2-v1).cross(v3-v1);

            mesh_structure.F.col(i) << stoi(substrings[1]),
            stoi(substrings[2]),
            stoi(substrings[3]),
            0.0f, 0.0f, 0.0f;

            for (unsigned j=0; j<face_normal.size(); j++)
            {
                mesh_structure.F.coeffRef(j+3,i) = face_normal.coeffRef(j);
            }

        }
        
        //iterate through matrix V
        for(unsigned i=0; i<mesh_structure.V.cols(); i++)
        {
            //for each face 
            for(unsigned j=0; j<mesh_structure.F.cols(); j++)
            {
                //for each face vertex
                for(unsigned k=0; k<3; k++)
                {
                    int face_vertex_index = int(mesh_structure.F.coeffRef(k,j));
                    //if face vertex index == i
                    if(face_vertex_index == i)
                    {
                        //add face normal to vertex normal
                        mesh_structure.V.block(3,i,3,1) = mesh_structure.V.block(3,i,3,1) + mesh_structure.F.block(3,j,3,1);
                    }

                }
            }
            //normalize vertex (3,4,5);
            Vector3f vertex_normal = mesh_structure.V.block(3,i,3,1);
            vertex_normal = vertex_normal.normalized();
            mesh_structure.V.block(3,i,3,1) << vertex_normal;
        }
        return mesh_structure;
    }
    else
    {
        //otherwise exit with output "Invalid file (not OFF file)"
        std::cout << "Invalid file (not OFF file)" << std::endl;
        return mesh_structure;
    }
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);

    M_vp << width/2, 0, 0, (width - 1)/2,
    0, height/2, 0, (height-1)/2,
    0, 0, 1, 0,
    0, 0, 0, 1;

    // M_vp << 2/width, 0, 0, 1,
    // 0, 2/height, 0, 1,
    // 0, 0, 1, 0,
    // 0, 0, 0, 1;
}

point screen_to_world(GLFWwindow* window)
{
    // Given the mouse's screen coordinates
    // the corresponding point in world coordinates is returned
    // accounting for screen size and view controls

    // Get the position of the mouse in the window
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    // Convert screen position to world coordinates
    double xworld = ((xpos/double(width))*2)-1;
    double yworld = (((height-1-ypos)/double(height))*2)-1; 

    // adjust position accounting for view controls
    Matrix2f inv_scale;
    inv_scale << 1/view_scale.coeff(0,0),0,0,1/view_scale.coeff(1,1);
    
    Vector2f adj_pos;
    adj_pos << float(xworld), float(yworld);
    adj_pos = inv_scale * adj_pos - view_pos;

    xworld = double(adj_pos.coeff(0));
    yworld = double(adj_pos.coeff(1));
    return point(xworld, yworld);
}

bool click_triangle(point click_point, triangle test_triangle)
{
    // Given a click point and a triangle
    // return whether click point is inside triangle
    // using barycentric coordinates

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
            }
            
        }
    }
    return false;
}

void transform_triangle(GLFWwindow* window, triangle sel_triangle, Matrix2f transform)
{
    // Given a triangle and a transformation matrix,
    // compute transformed coordinates of its vertices

    vector<Vector2f> new_vertices;

    //find barycentric center, referencing the current rendered positions of the triangles
    Vector2f v1 = Vector2f(tri_V.col(sel_triangle.clicked_index).coeff(0),tri_V.col(sel_triangle.clicked_index).coeff(1));
    Vector2f v2 = Vector2f(tri_V.col(sel_triangle.clicked_index+1).coeff(0),tri_V.col(sel_triangle.clicked_index+1).coeff(1));
    Vector2f v3 = Vector2f(tri_V.col(sel_triangle.clicked_index+2).coeff(0),tri_V.col(sel_triangle.clicked_index+2).coeff(1));
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

    // update matrices and VBOs
    // line vertices and VBO are also updated
    // given the transformed triangle is selected with a white border
    for(unsigned i = 0; i<new_vertices.size(); i++)
    {
        clicked_triangle.v[i].x = new_vertices[i].coeff(0);
        clicked_triangle.v[i].y = new_vertices[i].coeff(1);
        tri_V.col(clicked_triangle.clicked_index+i).coeffRef(0) = new_vertices[i].coeff(0);
        tri_V.col(clicked_triangle.clicked_index+i).coeffRef(1) = new_vertices[i].coeff(1);
        line_V.col(i) << tri_V.col(clicked_triangle.clicked_index+i).coeff(0),
                                        tri_V.col(clicked_triangle.clicked_index+i).coeff(1), 
                                        1.0,1.0,1.0;
    }

    tri_VBO.update(tri_V);
    line_VBO.update(line_V);

    // update starting position of the triangle to where the current mouse cursor is
    point world_click = screen_to_world(window);
    start_click << world_click.x, world_click.y;
}

Matrix4f camera_matrix(Vector3f eye_pos)
{
    Vector3f gaze=(Vector3f(0,0,0) - eye_pos);
    Vector3f view_up(0,1,0);

    Vector3f cam_u;
    Vector3f cam_v;
    Vector3f cam_w;

    cam_w = -(gaze.normalized());
    cam_u = view_up.cross(cam_w).normalized();
    cam_v = cam_w.cross(cam_u);

    Matrix4f R;
    R << cam_u.coeff(0),cam_u.coeff(1),cam_u.coeff(2), 0,
    cam_v.coeff(0),cam_v.coeff(1),cam_v.coeff(2), 0,
    cam_w.coeff(0),cam_w.coeff(1),cam_w.coeff(2), 0,
    0,0,0,1;

    R.transposeInPlace();

    Matrix4f T;
    T << 1, 0, 0, -eye_pos.coeff(0),
    0, 1, 0, -eye_pos.coeff(1),
    0, 0, 1, -eye_pos.coeff(2),
    0, 0, 0, 1;

    Matrix4f camera = R * T;

    return camera;
}

void removeColumn(MatrixXf& matrix, unsigned int colToRemove)
{
    // Eigen matrix column remover function courtesy of https://stackoverflow.com/questions/13290395/how-to-remove-a-certain-row-or-column-while-using-eigen-library-c
    // Thank the gods for StackOverflow
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    // Get mouse cursor position
    point world_click = screen_to_world(window);

    switch (mode)
    {
        default:
            break;
    }
    
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    // Get mouse cursor position
    point world_click = screen_to_world(window);

    // When left mouse button is pressed
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        switch (mode)
        {
            default:
                break;
        }
    }
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    //only perform action on key press, not key release
    if(action == GLFW_PRESS)
    {
        switch (key)
        {
            case GLFW_KEY_P:
            {
                //perspective projection matrix
                M_proj << 2*abs(near)/(r-l), 0, (r+l)/(r-l), 0,
                0, abs(near)*2/(t-b), (t+b)/(t-b), 0,
                0, 0, (abs(near)+abs(far))/(abs(near)-abs(far)), 2*abs(far)*abs(near)/(abs(near)-abs(far)),
                0, 0, -1, 0;
                break; 
            }
            case GLFW_KEY_O:
            {
                M_proj= M_orth;
                break; 
            }
            case GLFW_KEY_1:
            {                
                tri_mesh cube = load_mesh("../data/cube.off",Vector3d(0,0,0),1);
                meshes.push_back(cube);
                int insert_start;
                if(mesh_V.cols()==1)
                {
                    insert_start = 0;
                    mesh_V.conservativeResize(NoChange, cube.F.cols()*3);
                }
                else
                {
                    insert_start = mesh_V.cols(); 
                    mesh_V.conservativeResize(NoChange, mesh_V.cols()+cube.F.cols()*3); 
                }

                for(unsigned i=0; i<cube.F.cols();i++)
                {
                    for(unsigned j=0; j<3;j++)
                    {
                        mesh_V.col(insert_start) << cube.V.col((int)cube.F.coeffRef(j,i));
                        insert_start++;
                    }
                }

                insert_start = 0;
                for(unsigned i=0; i < cube.F.cols();i++)
                {
                    for(unsigned j=0; j<3;j++)
                    {
                        mesh_V.block(3,insert_start,3,1) = cube.F.block(3,i,3,1);
                        insert_start++;
                    }
                }

                //update VBO
                mesh_VBO.update(mesh_V);
                break;
            }
            case GLFW_KEY_2:
            {
                tri_mesh cube = load_mesh("../data/bumpy_cube.off",Vector3d(0,0,0),0.1);
                meshes.push_back(cube);
                int insert_start;
                if(mesh_V.cols()==1)
                {
                    insert_start = 0;
                    mesh_V.conservativeResize(NoChange, cube.F.cols()*3);
                }
                else
                {
                    insert_start = mesh_V.cols(); 
                    mesh_V.conservativeResize(NoChange, mesh_V.cols()+cube.F.cols()*3); 
                }

                for(unsigned i=0; i<cube.F.cols();i++)
                {
                    for(unsigned j=0; j<3;j++)
                    {
                        mesh_V.col(insert_start) << cube.V.col((int)cube.F.coeffRef(j,i));
                        insert_start++;
                    }
                }

                //update VBO
                mesh_VBO.update(mesh_V);
                break;
            } 
            case GLFW_KEY_3:
            {
                load_mesh("../data/bunny.off",Vector3d(0,0,0),1);
                break;
            }
            case  GLFW_KEY_W:
                eye_pos = eye_pos - 0.1*Vector3f::UnitZ();
                M_cam = camera_matrix(eye_pos);
            case  GLFW_KEY_S:
                eye_pos = eye_pos + 0.1*Vector3f::UnitZ();
                M_cam = camera_matrix(eye_pos);
                break;
            case  GLFW_KEY_A:
                eye_pos = eye_pos - 0.1*Vector3f::UnitX();
                M_cam = camera_matrix(eye_pos);
                break;
            case  GLFW_KEY_D:
                eye_pos = eye_pos + 0.1*Vector3f::UnitX();
                M_cam = camera_matrix(eye_pos);
                break;
            case  GLFW_KEY_Q:
                eye_pos = eye_pos + 0.1*Vector3f::UnitY();
                M_cam = camera_matrix(eye_pos);
                break;
            case  GLFW_KEY_Z:
                eye_pos = eye_pos - 0.1*Vector3f::UnitY();
                M_cam = camera_matrix(eye_pos);
                break;
            default:
                break;
        }

    }
}

int main(void)
{
    // add color palette to colors vector


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
    window = glfwCreateWindow(640, 480, "Triangles!", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);
    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);

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
    VertexArrayObject VAO;
    VAO.init();
    VAO.bind();

    // Initialize the VBOs
    line_VBO.init();
    line_V.resize(5,1);
    line_VBO.update(line_V);

    tri_VBO.init();
    tri_V.resize(5,1);
    tri_VBO.update(tri_V);

    mesh_VBO.init();
    mesh_V.resize(6,1);
    mesh_VBO.update(mesh_V);

    // Initialize view scale and position
    view_scale << 1, 0, 0, 1;
    view_pos << 0, 0;

    light spotlight;
    spotlight.position << 0, 5.0, 5.0;
    spotlight.color << 1.0, 1.0, 1.0;
    spotlight.intensity << 1.0, 1.0, 1.0;

    //viewport transformation matrix
    // M_vp << width/2, 0, 0, (width-1)/2,
    // 0, height/2, 0, (height-1)/2,
    // 0, 0, 1, 0,
    // 0, 0, 0, 1;

    // M_vp = M_vp.inverse().eval();

    M_vp << 2/width, 0, 0, 1,
    0, 2/height, 0, 1,
    0, 0, 1, 0,
    0, 0, 0, 1;

    //orthographic projection matrix
    M_orth << 2/(r-l), 0, 0, -(r+l)/(r-l),
    0, 2/(t-b), 0, -(t+b)/(t-b),
    0, 0, 2/(near-far), -(near+far)/(near-far),
    0, 0, 0, 1;

    //perspective project matrix
    // P << near, 0, 0, 0,
    // 0, near, 0, 0,
    // 0, 0, near+far, -far*near,
    // 0, 0, -1, 0;

    P << 2*abs(near)/(r-l), 0, (r+l)/(r-l), 0,
                0, abs(near)*2/(t-b), (t+b)/(t-b), 0,
                0, 0, (abs(near)+abs(far))/(abs(near)-abs(far)), 2*abs(far)*abs(near)/(abs(near)-abs(far)),
                0, 0, -1, 0;

    M_proj = M_orth;
    //camera view matrix
    M_cam = camera_matrix(eye_pos);

    //model to world frame transformation matrix
    Matrix4f M_model = Matrix<float, 4, 4>::Identity();

    Vector4f test(0.5,0.0,0.0,1.0);
    cout << M_cam * M_model * test << endl;
    cout << "---" << endl;
    cout << M_orth * M_cam * M_model * test << endl;
    cout << "---" << endl;
    Vector4f persp_test =  P * M_cam * M_model * test;
    persp_test = persp_test / persp_test.coeff(3);
    cout << persp_test << endl;

    Matrix4f M_comb = M_cam * M_model;
    Matrix4f M_normal = M_comb.inverse().transpose();

    // Initialize the OpenGL Program
    Program program;
    const GLchar* vertex_shader =
            "#version 330 core\n"
                    "layout(location=0) in vec3 position;"
                    "layout(location=1) in vec3 inNormal;"
                    "out vec4 normal;"
                    "out vec3 halfVec;"
                    "out vec3 lightDir;"
                    "uniform mat4 viewportMatrix;"
                    "uniform mat4 projMatrix;"
                    "uniform mat4 camMatrix;"
                    "uniform mat4 modelMatrix;"
                    "uniform mat4 normalMatrix;"
                    "uniform int shaderMode;"
                    "uniform mat2 scale;"
                    "uniform vec2 translation;"
                    "uniform vec3 lightPosition;"
                    "void main()"
                    "{"
                    "    vec4 pos = camMatrix * modelMatrix * vec4(position, 1.0);"
                    "    if(shaderMode != 0){"
                    "       vec4 lightPos = camMatrix * vec4(lightPosition, 1.0);"
                    "       normal = normalMatrix * vec4(inNormal, 0.0);"
                    "       vec3 v = normalize(-pos.xyz);"
                    "       lightDir = normalize(lightPos.xyz - pos.xyz);"
                    "       halfVec = normalize(v + lightDir);"
                    "    }"
                    "    gl_Position = projMatrix * pos;"
                    "}";
    const GLchar* fragment_shader =
            "#version 330 core\n"
                    "in vec4 normal;"
                    "in vec3 halfVec;"
                    "in vec3 lightDir;"
                    "uniform vec3 lightIntensity;"
                    "uniform int shaderMode;"
                    "uniform vec3 Ia;"
                    "uniform vec3 ka, kd, ks;"
                    "uniform float phongExp;"
                    "out vec4 fragmentColor;"
                    "void main()"
                    "{"
                    "    if (shaderMode == 0) {"
                    "       fragmentColor = vec4(kd, 1.0);"
                    "    }"
                    "    else {"
                    "       vec3 n = normalize(normal.xyz);"
                    "       vec3 h = normalize(halfVec);"
                    "       vec3 l = normalize(lightDir);"
                    "       vec3 intensity = ka * Ia + kd * lightIntensity * max(0.0, dot(n,l)) + ks * lightIntensity * pow(max(0.0,dot(n,h)),phongExp);"
                    "       fragmentColor = vec4(intensity, 1.0);"
                    "    }"
                    "}";

    // Compile the two shaders and upload the binary to the GPU
    program.init(vertex_shader,fragment_shader,"fragmentColor");
    program.bind();

    // Register the keyboard callback
    glfwSetKeyCallback(window, key_callback);

    //Register the mouse cursor position callback
    glfwSetCursorPosCallback(window, cursor_position_callback);

    // Register the mouse callback
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Update viewport
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    
    // glEnable(GL_DEPTH_TEST);
    // glDepthFunc(GL_LESS);
    
    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window))
    {
        VAO.bind();

        // Bind program
        program.bind();

        // Set the uniform view matrix and translation vectors
        glUniformMatrix2fv(program.uniform("scale"),1, GL_FALSE, view_scale.data());
        glUniform2fv(program.uniform("translation"),1, view_pos.data());

        //uniforms for transformation matrices
        glUniformMatrix4fv(program.uniform("viewportMatrix"),1, GL_FALSE, M_vp.data());
        glUniformMatrix4fv(program.uniform("camMatrix"),1, GL_FALSE, M_cam.data());
        glUniformMatrix4fv(program.uniform("modelMatrix"),1, GL_FALSE, M_model.data());
        glUniformMatrix4fv(program.uniform("projMatrix"),1, GL_FALSE, M_proj.data());
        glUniformMatrix4fv(program.uniform("normalMatrix"),1, GL_FALSE, M_normal.data());
        
        //light uniforms
        glUniform3fv(program.uniform("lightPosition"),1, spotlight.position.data());
        glUniform3fv(program.uniform("lightIntensity"),1, spotlight.intensity.data());
        Vector3f ambient_intensity(1.0,1.0,1.0);
        glUniform3fv(program.uniform("Ia"),1,ambient_intensity.data());

        glUniform1i(program.uniform("shaderMode"),0);

        // Clear the framebuffer
        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

        glClear(GL_COLOR_BUFFER_BIT);
        //TODO: glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_TEST);

        if(mesh_V.cols()>1)
        {
            mesh_VBO.bind();
            //for each mesh
            for(unsigned i = 0; i<meshes.size(); i++)
            {
                //draw mesh elements
                Vector3f color_gl = meshes[i].diff_color/255;
                glUniform3fv(program.uniform("kd"),1, color_gl.data());
                
                glEnableVertexAttribArray(0);
                glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), 0);

                if(meshes[i].shader_type == 'p' || meshes[i].shader_type == 'f')
                {
                    glUniform1i(program.uniform("shaderMode"),2);
                    Vector3f ka_gl = meshes[i].ka/255;
                    glUniform3fv(program.uniform("ka"),1, ka_gl.data());
                    glUniform3fv(program.uniform("ks"),1, meshes[i].ks.data());
                    glUniform1f(program.uniform("phongExp"), meshes[i].phong_exp);

                    glEnableVertexAttribArray(1);
                    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (const GLvoid *)12);
                }
                else if (meshes[i].shader_type == 'w')
                {
                    glUniform1i(program.uniform("shaderMode"),0);
                }

                for(unsigned j=0; j<mesh_V.cols()/3; j++)
                {
                    //draw triangles
                    if(meshes[i].shader_type == 'w')
                    {
                        glDrawArrays(GL_LINE_LOOP,j*3,3);
                    }
                    else
                    {
                        glDrawArrays(GL_TRIANGLES, j*3, 3);
                    }
                    
                }
                if(meshes[i].shader_type == 'f')
                {
                    glUniform3fv(program.uniform("kd"),1, Vector3f(0,0,0).data());
                    for(unsigned j=0; j<mesh_V.cols()/3; j++)
                    {
                        glDrawArrays(GL_LINE_LOOP,j*3,3);
                    }
                }
            }
  
        }
            

        //if triangle is selected
        if(mode == 'm')
        {
            
        }
        
        //if a triangle is being inserted
        if(mode == 'i')
        {

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
    mesh_VBO.free();

    // Deallocate glfw internals
    glfwTerminate();
    return 0;
}
