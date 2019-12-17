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
#include <algorithm> 
#include <cmath>
#include <fstream>
#include <string>
#include <limits>

// VertexBufferObject wrappers
VertexBufferObject line_VBO;
VertexBufferObject tri_VBO;
VertexBufferObject mesh_VBO;

// Viewing Volume (camera coordinates)
float near = 1;
float far = -20;
float l = -2.0;
float r = 2.0;
float t = 2.0;
float b = -2.0;

// Viewing Transformation Matrices
Matrix4f M_vp;

Matrix4f M_orth;
Matrix4f P;
Matrix4f M_proj;

Matrix4f M_cam;
Matrix4f M_comb;
Matrix4f M_normal;

// Camera initial position
Vector3f eye_pos(0,0,3);

// Matrices contain the vertex positions of the lines and triangles
MatrixXf line_V(6,1);
MatrixXf tri_V(6,1);
MatrixXf mesh_V(6,1);

// variables for tracking which mode the drawing application is in
char mode = ' ';


//Light Class
class light
{
    public:
        Vector3f position;
        Vector3f color;
        Vector3f intensity;
};

// Human Pose Parameters
// Vector of pose coordinates for each animation frame
vector<MatrixXf> poses;

// The following code is for rendering the lines of the pose skeleton
// MatrixXf pose_V(6,2);
// VertexBufferObject pose_VBO;

//start and end point indices of Pose_3D keypoints for drawing pose skeleton
vector<int> point_i{0,1,2,0,6,7, 0,13,13,17,18,13,25,26};
vector<int> point_j{1,2,3,6,7,8,13,15,17,18,19,25,26,27};

class tri_mesh
{
    //basic triangular mesh
    public:
        //V is a matrix (dimension 6 x #V where #V is the number of vertices)
        // contains the positions and normals of the vertices of the mesh 
        MatrixXf V;
        //F is an matrix (dimension 6 x #F where #F is the number of faces) which contains the descriptions of the triangles in the mesh. 
        MatrixXf F;
        Matrix4f M_model;
        Matrix4f M_model_clicked;
        Vector3f diff_color;
        Vector3f ka;
        float ks;
        char shader_type;
        float phong_exp;
};

//Store all mesh objects in vector, 'meshes'
vector<tri_mesh> meshes;
// Number of animation keyframes
int num_keyframes = 0;
// Time step between animation key frames
float time_step = 1.0; 
// Variable for recording when the animation starts
chrono::time_point<chrono::high_resolution_clock> t_start;

// Color, Point, and Triangle Classes for easier-to-read code

class point
{
    public:
        float x;
        float y;
        float z;
        point()
        {
            x=0.0;
            y=0.0;
            z=0.0;
        }
        point(float x1, float y1, float z1)
        {
            x = x1;
            y = y1;
            z = z1;
        }
};

//Useful functions for splitting string based on space and comma delimiters
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
vector<string> comma_sep_string(string line)
{
    //separates string into a vector of substring as delimited by commas
    vector<int> split_index;
    vector<string> substrings;
    for (unsigned i=0;i<line.size();i++)
    {
        if(line[i]==',')
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

void load_pose(string csv_filepath, Vector3f position, float scale)
{
    //load pose given a csv file of keypoints
    string line;
    vector<string> substrings;
    ifstream in(csv_filepath);

    if(in.is_open())
    {
        while(getline(in,line))
        {
            MatrixXf vertices(3,1);
            unsigned i = 0;
            substrings = comma_sep_string(line);
            for (unsigned j=0;j<substrings.size();j++)
            {
                if(j>0)
                {
                    vertices.conservativeResize(3,vertices.cols()+1);
                }
                vertices.col(i) << stof(substrings[j])*scale+position.coeff(0), stof(substrings[j+2])*scale+position.coeff(1), stof(substrings[j+1])*-scale+position.coeff(2);
                
                j=j+2;
                i++;
            }
            poses.push_back(vertices);
        }
        in.close();
    }
}

Matrix4f cube_transform(point p1, point p2, float x_scale, float z_scale)
{
    //transform basic cube mesh to fit the body part between 2 points
    // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    Vector3f cube_vector(0,-1,0);
    Vector3f part_vector(p2.x-p1.x, p2.y-p1.y, p2.z-p1.z);
    float y_scale = part_vector.norm();

    part_vector = part_vector.normalized();

    Vector3f v = cube_vector.cross(part_vector);
    float s = v.norm();
    float c = cube_vector.dot(part_vector);

    Matrix3f v_x;
    v_x<< 0, -v.coeff(2), v.coeff(1), 
    v.coeff(2), 0, -v.coeff(0), 
    -v.coeff(1), v.coeff(0), 0;

    Matrix3f rotation;
    rotation = Matrix3f::Identity()+ v_x + v_x*v_x*(1-c)/(s*s);
    
    Matrix4f rotation_h;
    rotation_h << rotation.row(0), 0,
    rotation.row(1), 0,
    rotation.row(2), 0,
    0, 0, 0, 1;

    //scale/stretch cube
    Matrix4f transformation; 
    transformation << x_scale, 0, 0, 0,
    0, y_scale, 0, 0,
    0, 0, z_scale, 0,
    0, 0, 0, 1;

    //rotate stretched cube
    transformation = rotation_h * transformation;

    //translate stretched cube
    Matrix4f translate;
    translate << 1, 0, 0, p1.x,
    0, 1, 0, p1.y,
    0, 0, 1, p1.z,
    0, 0, 0, 1;

    transformation = translate * transformation;

    return transformation;
}

tri_mesh load_mesh(string off_filepath, Vector3f position, double scale)
{
    // returns a triangule mesh object at the file path for the OFF file
    tri_mesh mesh_structure;
    mesh_structure.diff_color = Vector3f(255,255,255);
    mesh_structure.ka = mesh_structure.diff_color * .15;
    mesh_structure.ks =  0.5;
    mesh_structure.phong_exp = 64;
    mesh_structure.shader_type = 'f';
    mesh_structure.M_model = Matrix<float, 4, 4>::Identity();

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
        // std::cout << "Valid OFF File" << std::endl;
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
    float f_width, f_height;
    f_width = float(width);
    f_height = float(height);
    float aspect_x = fmin(f_height/f_width,1.);
    float aspect_y = fmin(f_width/f_height,1.);
    M_vp << aspect_x, 0, 0, 0,
    0, aspect_y, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1;

    // M_vp << width/2, 0, 0, (width - 1)/2,
    // 0, height/2, 0, (height-1)/2,
    // 0, 0, 1, 0,
    // 0, 0, 0, 1;

    // M_vp << 2/width, 0, 0, 1,
    // 0, 2/height, 0, 1,
    // 0, 0, 1, 0,
    // 0, 0, 0, 1;
}

Matrix4f camera_matrix()
{
    Vector3f gaze=(Vector3f(0,0,0) - eye_pos);
    Vector3f view_up(0,1,0);

    Vector3f cam_u;
    Vector3f cam_v;
    Vector3f cam_w;

    cam_w = -(gaze.normalized());
    cam_u = view_up.cross(cam_w).normalized();
    cam_v = cam_w.cross(cam_u);

    Matrix4f camera;
    camera << cam_u.coeff(0),cam_v.coeff(0),cam_w.coeff(0), eye_pos.coeff(0),
    cam_u.coeff(1),cam_v.coeff(1),cam_w.coeff(1), eye_pos.coeff(1),
    cam_u.coeff(2),cam_v.coeff(2),cam_w.coeff(2), eye_pos.coeff(2),
    0,0,0,1;

    camera = camera.inverse().eval();

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


void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
 
}

void mesh_V_update(tri_mesh new_mesh)
{
    //given a new tri_mesh object, update the mesh_V containing vertices of all mesh triangles

    //insert mesh vertices at insert_start, which depends on existing mesh_V size
    int insert_start;

    if(mesh_V.cols()<2)
    {
        insert_start = 0;
        mesh_V.conservativeResize(NoChange, new_mesh.F.cols()*3);
    }
    else
    {
        insert_start = mesh_V.cols(); 
        mesh_V.conservativeResize(NoChange, mesh_V.cols()+new_mesh.F.cols()*3); 
    }

    //for each mesh face
    for(unsigned i=0; i < new_mesh.F.cols();i++)
    {
        //for each of mesh face's vertex
        for(unsigned j=0; j<3;j++)
        {
            //transform vertex based on model Matrix
            Vector4f vertex_trans;
            Vector4f normal_trans;
            vertex_trans << new_mesh.V.col((int)new_mesh.F.coeffRef(j,i)).head<3>(),1;
            normal_trans << new_mesh.F.block(3,i,3,1), 1;

            vertex_trans = new_mesh.M_model * vertex_trans;
            normal_trans = new_mesh.M_model * normal_trans;

            // new_mesh.V.col((int)new_mesh.F.coeffRef(j,i)) << vertex_trans.head<3>(), normal_trans.head<3>();

            // insert mesh vertex and vertex normals into the mesh_V matrix
            mesh_V.col(insert_start) << vertex_trans.head<3>(), 1, 1, 1;
            
            //vertex normals equal to face normals
            mesh_V.block(3,insert_start,3,1) = normal_trans.head<3>();
            insert_start++;
        }
    }
}

void tri_V_update(tri_mesh new_mesh)
{
    //given a new tri_mesh object, update the tri_V containing vertices of all mesh triangles

    //insert mesh vertices at insert_start, which depends on existing tri_V size
    int insert_start;

    if(tri_V.cols()<2)
    {
        insert_start = 0;
        tri_V.conservativeResize(NoChange, new_mesh.F.cols()*3);
    }
    else
    {
        insert_start = tri_V.cols(); 
        tri_V.conservativeResize(NoChange, tri_V.cols()+new_mesh.F.cols()*3); 
    }

    //for each mesh face
    for(unsigned i=0; i < new_mesh.F.cols();i++)
    {
        //for each of mesh face's vertex
        for(unsigned j=0; j<3;j++)
        {
            // insert mesh vertex and vertex normals into the tri_V matrix
            tri_V.col(insert_start) << new_mesh.V.col((int)new_mesh.F.coeffRef(j,i));
            
            //vertex normals equal to face normals
            tri_V.block(3,insert_start,3,1) = -new_mesh.F.block(3,i,3,1);
            insert_start++;
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
                load_pose("../data/vertices.csv",Vector3f(0,.917,0),0.00123);
                num_keyframes = poses.size();
                
                // The following code is for rendering the lines of the pose skeleton
                // for(unsigned i=0; i<point_i.size(); i++)
                // {
                //     if(i>0)
                //     {
                //         pose_V.conservativeResize(6,pose_V.cols()+2);
                //     }
                //     pose_V.col(i*2) << poses[0].coeffRef(0,point_i[i]), poses[0].coeffRef(1,point_i[i]), poses[0].coeffRef(2,point_i[i]), 1, 1, 1; 
                //     pose_V.col(2*i+1) << poses[0].coeffRef(0,point_j[i]), poses[0].coeffRef(1,point_j[i]), poses[0].coeffRef(2,point_j[i]), 1, 1, 1;
                // }
                cout << "Pose coordinates loaded" << endl;

                //for each pose keyframe
                for(unsigned j=0; j<poses.size(); j++) //j=0
                {
                    //for each body part
                    for(unsigned i=0; i<point_i.size(); i++) //i=0
                    {
                        tri_mesh cube = load_mesh("../data/cube.off",Vector3f(0,0,0),1); //cube V=8 F=12
                        point arm1 = point(poses[j].coeffRef(0,point_i[i]),poses[j].coeffRef(1,point_i[i]),poses[j].coeffRef(2,point_i[i]));
                        point arm2 = point(poses[j].coeffRef(0,point_j[i]),poses[j].coeffRef(1,point_j[i]),poses[j].coeffRef(2,point_j[i]));
                        
                        if(i==6)
                        {
                            cube.M_model = cube_transform(arm1,arm2,0.1,0.1)*cube.M_model;
                        }
                        else
                        {
                            cube.M_model = cube_transform(arm1,arm2,0.05,0.05)*cube.M_model;
                        }
                    
                        //update mesh_V
                        mesh_V_update(cube);
                        meshes.push_back(cube);
                    }
                }
                
                cout << "All " << poses.size() << " keyframe poses loaded" << endl;
                // cout << mesh_V.cols() <<endl;
                //update VBO
                mesh_VBO.update(mesh_V.block(0,0,6,504));
                break;
            }
            case  GLFW_KEY_SLASH:
            {
                cout << "*** Animation Now Playing! ***" << endl;
                // Store the current time for tracking animation time
                t_start = std::chrono::high_resolution_clock::now();

                // change application mode
                mode = 'a';
                break;
            }
            case  GLFW_KEY_W:
                eye_pos = eye_pos - 0.1*Vector3f::UnitZ();
                M_cam = camera_matrix();
                break;
            case  GLFW_KEY_S:
                eye_pos = eye_pos + 0.1*Vector3f::UnitZ();
                M_cam = camera_matrix();
                break;
            case  GLFW_KEY_A:
                eye_pos = eye_pos - 0.1*Vector3f::UnitX();
                M_cam = camera_matrix();
                break;
            case  GLFW_KEY_D:
                eye_pos = eye_pos + 0.1*Vector3f::UnitX();
                M_cam = camera_matrix();
                break;
            case  GLFW_KEY_Q:
                eye_pos = eye_pos + 0.1*Vector3f::UnitY();
                M_cam = camera_matrix();
                break;
            case  GLFW_KEY_Z:
                eye_pos = eye_pos - 0.1*Vector3f::UnitY();
                M_cam = camera_matrix();
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
    window = glfwCreateWindow(640, 480, "Wushu!", NULL, NULL);
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

    float f_width, f_height;
    f_width = float(width);
    f_height = float(height);

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
    cout << "====================" << endl;
    // Initialize the VAO
    VertexArrayObject VAO;
    VAO.init();
    VAO.bind();

    // Initialize the VBOs
    line_VBO.init();
    line_V.resize(6,6);
    line_V.col(0) << -100,0.,0.,1.,0.,0.;
    line_V.col(1) << 100,0.,0.,1.,0.,0.;
    line_V.col(2) << 0.,-100,0.,0.,0.,1.;
    line_V.col(3) << 0.,100,0.,0.,0.,1.;
    line_V.col(4) << 0,0.,-100.,0.,1.,0.;
    line_V.col(5) << 0,0.,100.,0.,1.,0.;

    line_VBO.update(line_V);

    //load floor
    tri_VBO.init();
    tri_V.resize(6,1);

    tri_mesh floor_mesh = load_mesh("../data/ground_plane.off",Vector3f(0,0,0),1.0);
    tri_V_update(floor_mesh);
    tri_VBO.update(tri_V);

    mesh_VBO.init();
    mesh_V.resize(6,1);
    mesh_VBO.update(mesh_V);

    // pose_VBO.init();
    // pose_V.resize(6,2);
    // pose_VBO.update(pose_V);

    light spotlight;
    spotlight.position << 0, 3, 3.0;
    spotlight.color << 1.0, 1.0, 1.0;
    spotlight.intensity << 1.0, 1.0, 1.0;

    //viewport transformation matrix
    float aspect_x = fmin(f_height/f_width,1.);
    float aspect_y = fmin(f_width/f_height,1.);
    M_vp << aspect_x, 0, 0, 0,
    0, aspect_y, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1;

    //orthographic projection matrix
    M_orth << 2/(r-l), 0, 0, -(r+l)/(r-l),
    0, 2/(t-b), 0, -(t+b)/(t-b),
    0, 0, 2/(far-near), -(near+far)/(far-near),
    0, 0, 0, 1;

    M_proj = M_orth;
    //camera view matrix
    M_cam = camera_matrix();

    // Vector4f test(1.,1.,0.,1.0);
    // cout << M_cam * test << endl;
    // cout << "---" << endl;
    // cout << M_orth * M_cam * test << endl;
    // cout << "---" << endl;
    // cout << M_vp * M_orth * M_cam * test << endl;
    // cout << "---" << endl;
    
    // Vector4f persp_test =   M_orth * P * M_cam * test;
    // persp_test = persp_test / persp_test.coeff(3);
    // cout << persp_test << endl;

    // Initialize simpler program for drawing axes
    Program axes_program;
    const GLchar* ax_vertex_shader =
            "#version 330 core\n"
                    "layout(location=0) in vec3 position;"
                    "layout(location=1) in vec3 inColor;"
                    "uniform mat4 viewportMatrix;"
                    "uniform mat4 projMatrix;"
                    "uniform mat4 camMatrix;"
                    "out vec3 vertexColor;"
                    "void main()"
                    "{"
                    "    vertexColor = inColor;"
                    "    vec4 pos = camMatrix * vec4(position, 1.0);"
                    "    gl_Position = viewportMatrix * projMatrix * pos;"
                    "}";
    const GLchar* ax_fragment_shader =
            "#version 330 core\n"
                    "in vec3 vertexColor;"
                    "out vec4 fragmentColor;"
                    "void main()"
                    "{"
                    "    fragmentColor = vec4(vertexColor, 1.0);"
                    "}";

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
                    "uniform mat4 normalMatrix;"
                    "uniform vec3 lightPosition;"
                    "void main()"
                    "{"
                    "   vec4 pos = camMatrix * vec4(position, 1.0);"
                    "   vec4 lightPos = camMatrix * vec4(lightPosition, 1.0);"
                    "   normal = normalMatrix * vec4(inNormal, 0.0);"
                    "   vec3 v = normalize(-pos.xyz);"
                    "   lightDir = normalize(lightPos.xyz - pos.xyz);"
                    "   halfVec = normalize(v + lightDir);"
                    "   gl_Position = viewportMatrix * projMatrix * pos;"
                    "}";
    const GLchar* fragment_shader =
            "#version 330 core\n"
                    "in vec4 normal;"
                    "in vec3 halfVec;"
                    "in vec3 lightDir;"
                    "uniform vec3 lightIntensity;"
                    "uniform vec3 Ia;"
                    "uniform vec3 ka;"
                    "uniform vec3 kd;"
                    "uniform float ks;"
                    "uniform float phongExp;"
                    "out vec4 fragmentColor;"
                    "void main()"
                    "{"
                    "   vec3 n = normalize(normal.xyz);"
                    "   vec3 l = normalize(lightDir);"
                    "   float diffuse =  max(0.0, dot(n,l));"
                    "   vec3 intensity = ka * Ia + kd * lightIntensity * (diffuse);"
                    "   fragmentColor = vec4(intensity, 1.0);"
                    "}";

    axes_program.init(ax_vertex_shader,ax_fragment_shader,"fragmentColor");
    axes_program.bind();
    // Compile the two shaders and upload the binary to the GPU
    program.init(vertex_shader,fragment_shader,"fragmentColor");
    program.bind();

    // Register the keyboard callback
    glfwSetKeyCallback(window, key_callback);

    // Register the mouse callback
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Update viewport
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    
    // Enable depth tests
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // Clear the framebuffer
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window))
    {
        //glClear(GL_COLOR_BUFFER_BIT);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        axes_program.bind();
        VAO.bind();

        glUniformMatrix4fv(axes_program.uniform("viewportMatrix"),1, GL_FALSE, M_vp.data());
        glUniformMatrix4fv(axes_program.uniform("camMatrix"),1, GL_FALSE, M_cam.data());
        glUniformMatrix4fv(axes_program.uniform("projMatrix"),1, GL_FALSE, M_proj.data());

        line_VBO.bind();

        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), 0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (const GLvoid *)(12));
        glDrawArrays(GL_LINES,0,line_V.cols());

        // The following code is for rendering the lines of the pose skeleton
        // if(pose_V.cols()>1)
        // {
        //     pose_VBO.bind();
        //     glEnableVertexAttribArray(0);
        //     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), 0);
        //     glEnableVertexAttribArray(1);
        //     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (const GLvoid *)(12));
        //     glDrawArrays(GL_LINES,0,pose_V.cols());
        // }
        // Bind program
        program.bind();

        //uniforms for transformation matrices
        glUniformMatrix4fv(program.uniform("viewportMatrix"),1, GL_FALSE, M_vp.data());
        glUniformMatrix4fv(program.uniform("camMatrix"),1, GL_FALSE, M_cam.data());

        glUniformMatrix4fv(program.uniform("projMatrix"),1, GL_FALSE, M_proj.data());
        //glUniformMatrix4fv(program.uniform("projMatrix"),1, GL_FALSE, &perspective_projection[0][0]);
        
        //light uniforms
        glUniform3fv(program.uniform("lightPosition"),1, spotlight.position.data());
        glUniform3fv(program.uniform("lightIntensity"),1, spotlight.intensity.data());
        Vector3f ambient_intensity(1.0,1.0,1.0);
        glUniform3fv(program.uniform("Ia"),1,ambient_intensity.data());

        glUniform1i(program.uniform("shaderMode"),0);


        //tri_V and tri_VBO models the floor
        tri_VBO.bind();

        M_comb = M_cam;
        M_normal = M_comb.inverse().transpose();
        
        glUniformMatrix4fv(program.uniform("normalMatrix"),1, GL_FALSE, M_normal.data());

        //draw mesh elements
        Vector3f color_gl(0.278,0.455,0.706);

        glUniform3fv(program.uniform("kd"),1, color_gl.data());
        
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), 0);

        Vector3f ka_gl(0.15,0.15,0.15);
        glUniform3fv(program.uniform("ka"),1, ka_gl.data());

        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (const GLvoid *)(12));

        glDrawArrays(GL_TRIANGLES, 0, 6);

        
        //mesh_V and mesh_VBO models the human pose
        if(mesh_V.cols()>1)
        {
            mesh_VBO.bind();

            M_comb = M_cam;
            M_normal = M_comb.inverse().transpose();
            
            glUniformMatrix4fv(program.uniform("normalMatrix"),1, GL_FALSE, M_normal.data());

            //draw mesh elements
            Vector3f color_gl(1,1,1);

            glUniform3fv(program.uniform("kd"),1, color_gl.data());
            
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), 0);

            Vector3f ka_gl(0.15,0.15,0.15);
            glUniform3fv(program.uniform("ka"),1, ka_gl.data());

            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (const GLvoid *)(12));

            glDrawArrays(GL_TRIANGLES, 0, 504);
        }

        if(mode == 'a')
        {
            MatrixXf mesh_interp;
            mesh_interp.resize(6,504);
            // track time elapsed since animation start
            auto t_now = std::chrono::high_resolution_clock::now();
            float time = std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();

            if(time > 0)
            {
                //for each frame
                // if it's not the end of the animation
                if(ceil(time/time_step)<num_keyframes)
                {
                    int ins = 0;
                    //for each vertex
                    for (unsigned i = 0; i<504; i++)
                    {
                        // keyframe pairs at time 
                        int k_0 = floor(time / time_step);
                        int k_1 = ceil(time / time_step);
                         
                        //same body part mesh in two different keyframe poses
                        // tri_mesh part_0 = meshes[k_0*point_i.size()];
                        // tri_mesh part_1 = meshes[k_1*point_i.size()];

                        //interpolate between vertices of body part in the two different keyframe poses
                        
                        for(unsigned j=0; j < mesh_interp.cols();j++)
                        {
                            VectorXf v_0_col(6,1);
                            v_0_col = mesh_V.col(k_0*504+j);
                            VectorXf v_1_col(6,1);
                            v_1_col = mesh_V.col(k_1*504+j);
                            mesh_interp.col(j) << v_0_col + (time-floor(time / time_step)) * (v_1_col - v_0_col);
                        }
                    }
                    mesh_VBO.update(mesh_interp);
                }
                else
                {
                    t_start = t_now;
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
    axes_program.free();
    line_VBO.free();
    tri_VBO.free();
    mesh_VBO.free();

    // Deallocate glfw internals
    glfwTerminate();
    return 0;
}
