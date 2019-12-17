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

// VertexBufferObject wrapper
VertexBufferObject line_VBO;
VertexBufferObject tri_VBO;
VertexBufferObject mesh_VBO;

// Viewing Transformation Matrices
float near = 20;
float far = -20;
float l = -3.0;
float r = 3.0;
float t = 3.0;
float b = -3.0;
Matrix4f M_vp;

Matrix4f M_orth;
Matrix4f P;
Matrix4f M_proj;

Matrix4f M_cam;
Matrix4f M_comb;
Matrix4f M_normal;

Vector3f eye_pos(0,0,3);

// Contains the vertex positions of the lines and triangles
MatrixXf line_V(6,1);
MatrixXf tri_V(5,1);
MatrixXf mesh_V(6,1);

// variables for tracking which mode the drawing application is in
char mode = ' ';
bool animation_mode = false;

// variables for tracking mouse states
Vector3f start_click;
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

// Human Pose Parameters
// Vector of pose coordinates for each animation frame
vector<MatrixXf> poses;

MatrixXf pose_V(6,2);
// MatrixXf pose_V(3,1);
VertexBufferObject pose_VBO;

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
        Matrix4f M_model;
        Matrix4f M_model_clicked;
        Vector3f diff_color;
        Vector3f ka;
        float ks;
        char shader_type;
        float phong_exp;
        float bound_radius;
        Vector3f bound_center;
        Vector3f bound_center_clicked;
        bool clicked = false; // tracks whether the mesh is selected
        int clicked_index; //start index on meshes vector
};

class ray
{
    public:
        Vector3f e;
        Vector3f d;
};
class sphere
{
    public:
        Vector3f center;
        double radius;
};

class hit_sphere
{
    public:
        bool hit = false;
        double t;
        sphere hit_obj;
};

// global variable for storing the clicked triangle's properties
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
vector<string> comma_sep_string(string line)
{
    //separates string into a vector of substring as space dividers
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
    cout << scale << endl;
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
        
        mesh_structure.bound_center << 0,0,0;
        

        //iterate through matrix V
        for(unsigned i=0; i<mesh_structure.V.cols(); i++)
        {
            for(unsigned j=0; j<3; j++)
            {
                mesh_structure.bound_center.coeffRef(j)+=mesh_structure.V.coeffRef(j,i);
            }
            
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

        mesh_structure.bound_center = mesh_structure.bound_center / num_vertices;
        mesh_structure.bound_radius = 0.0;

        for(unsigned i=0; i<mesh_structure.V.cols(); i++)
        {
            Vector3f vertex = mesh_structure.V.block(0,i,3,1);
            float cand_radius = (vertex - mesh_structure.bound_center).norm();
            if(cand_radius > mesh_structure.bound_radius)
            {
                mesh_structure.bound_radius = cand_radius;
            }
        }
        // cout << mesh_structure.bound_center << endl;
        // cout << mesh_structure.bound_radius << endl;
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

point screen_to_world(GLFWwindow* window)
{
    // Given the mouse's screen coordinates
    // the corresponding point in world coordinates is returned
    // accounting for screen size and view controls

    // Get the position of the mouse in the window
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    //printf("===\n%f, %f\n---\n", xpos, ypos);

    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    // Convert screen position to canonical
    double xworld = ((xpos/double(width))*2)-1;
    double yworld = (((height-1-ypos)/double(height))*2)-1; 
    double zworld = 0;
    //printf("%f, %f\n---\n", xworld, yworld);

    // adjust position accounting for view controls
    // Matrix2f inv_scale;
    // inv_scale << 1/view_scale.coeff(0,0),0,0,1/view_scale.coeff(1,1);
    
    Vector4f adj_pos;
    adj_pos << float(xworld), float(yworld), 0, 1;
    // adj_pos = inv_scale * adj_pos - view_pos;
    Matrix4f inv_proj = M_proj.inverse().eval();
    Matrix4f inv_cam = M_cam.inverse().eval();
    adj_pos = inv_cam * inv_proj * adj_pos;
    //printf("%f, %f\n===\n", adj_pos.coeff(0), adj_pos.coeff(1),adj_pos.coeff(2));

    xworld = adj_pos.coeff(0);
    yworld = adj_pos.coeff(1);
    zworld = adj_pos.coeff(2);
    return point(xworld, yworld, zworld);
}

hit_sphere ray_hit_sphere(ray r, double t_lower, sphere test_sphere)
{
    hit_sphere test;
    const double sphere_radius = test_sphere.radius;
    Vector3f sphere_center = test_sphere.center;

    //Calculate discriminant
    Vector3f e_c = r.e - sphere_center;
    //quadratic equation
    double a = r.d.dot(r.d);
    double b = r.d.dot(e_c);
    double c = e_c.dot(e_c) - sphere_radius * sphere_radius;

    //discriminant
    double discriminant = b*b - a*c;

    if(discriminant>=0)
    {
        test.hit = true;
        test.hit_obj = test_sphere;
        //if discriminant is greater than or equal to 0
        //Calculate sphere intersection parameter, test.t
        if(discriminant ==0)
        {
            test.t = -b/a;
        }
        else
        {
            test.t = (-b - sqrt(discriminant)) / a;
            if(test.t < 0)
            {
                test.t = (-b + sqrt(discriminant)) / a;
            }
            else
            {
                double t2 = (-b + sqrt(discriminant)) / a;
                if (t2>0 && t2<test.t)
                {
                    test.t = t2;
                }
            }
        }
    }
    return test;
}

void transform_triangle(GLFWwindow* window, tri_mesh sel_triangle, Matrix2f transform)
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

    tri_VBO.update(tri_V);
    line_VBO.update(line_V);

    // update starting position of the triangle to where the current mouse cursor is
    point world_click = screen_to_world(window);
    start_click << world_click.x, world_click.y, world_click.z;
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

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    // Get mouse cursor position
    point world_click = screen_to_world(window);

    if(clicked_mesh.clicked)
    {
        //set current mouse position to vector mouse_pos
        Vector3f mouse_pos;
        mouse_pos << world_click.x, world_click.y, world_click.z;

        //calculated difference btw start_click and mouse_pos
        Vector3f tr = mouse_pos - start_click;

        Matrix4f translation;
        translation << 0, 0, 0, tr.coeff(0),
        0, 0, 0, tr.coeff(1),
        0, 0, 0, tr.coeff(2),
        0, 0, 0, 0;

        // int insert_start = 0;
        // for(unsigned i = 0; i < clicked_mesh.clicked_index; i++)
        // {
        //     insert_start+= meshes[i].F.cols()*3;
        // }
        
        //translate all vertices of the clicked mesh

        clicked_mesh.M_model = clicked_mesh.M_model_clicked + translation;
        clicked_mesh.bound_center = clicked_mesh.bound_center_clicked + tr;
        //cout << clicked_mesh.M_model << endl;
        meshes[clicked_mesh.clicked_index] = clicked_mesh;

        // meshes[clicked_mesh.clicked_index].bound_center

    }
    
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    // Get mouse cursor position
    point world_click = screen_to_world(window);

    // When left mouse button is pressed
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        ray click_ray;
        click_ray.e << world_click.x,world_click.y,world_click.z;
        click_ray.d << 0,0,-1;

        Vector4f d_hom;
        d_hom << click_ray.d, 0;
        Matrix4f inv_proj = M_proj.inverse().eval();
        Matrix4f inv_cam = M_cam.inverse().eval();
        d_hom = inv_cam * inv_proj * d_hom;
        click_ray.d << (d_hom.head(3)).normalized();
        //cout<< click_ray.d <<endl;

        float min_dist = numeric_limits<float>::infinity(); 
        for(unsigned i=0; i<meshes.size(); i++)
        {
            sphere test_sphere;
            test_sphere.center = meshes[i].bound_center;
            test_sphere.radius = meshes[i].bound_radius;
            hit_sphere check = ray_hit_sphere(click_ray,0,test_sphere);
            if(check.hit)
            {
                if(check.t < min_dist)
                {
                    min_dist = check.t;
                    clicked_mesh = meshes[i];
                    clicked_mesh.clicked_index = i;
                    clicked_mesh.clicked = true;
                    cout<< i <<endl;
                }
                start_click << world_click.x, world_click.y, world_click.z;
                
            }
        }
        if(min_dist == numeric_limits<float>::infinity())
        {
            clicked_mesh.clicked_index = 0;
            clicked_mesh.clicked = false;
        }
        else
        {
            clicked_mesh.M_model_clicked = clicked_mesh.M_model;
            clicked_mesh.bound_center_clicked << world_click.x, world_click.y, world_click.z;
        }
        

 
    }
    else
    {
        {
            clicked_mesh.clicked = false;
            clicked_mesh.clicked_index = 0;
        }
    }
}

void mesh_V_update(tri_mesh new_mesh)
{
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

    
    for(unsigned i=0; i < new_mesh.F.cols();i++)
    {
        for(unsigned j=0; j<3;j++)
        {
            mesh_V.col(insert_start) << new_mesh.V.col((int)new_mesh.F.coeffRef(j,i));
            //vertex normals equal to face normals
            mesh_V.block(3,insert_start,3,1) = new_mesh.F.block(3,i,3,1);
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

                // M_proj << 2*near/(r-l), 0, (l+r)/(l-r), 0,
                // 0, near*2/(t-b), (b+t)/(b-t), 0,
                // 0, 0, (near+far)/(near-far), 2*far*near/(near-far),
                // 0, 0, 1, 0;

                // cout<< M_proj <<endl;

                // M_proj = M_orth * P ;
                // cout<< M_proj <<endl;
                break; 
            }
            case GLFW_KEY_O:
            {
                M_proj= M_orth;
                break; 
            }
            case GLFW_KEY_C:
            {
                load_pose("../data/vertices.csv",Vector3f(0,.917,0),0.00123);
                vector<int> point_i{0,1,2,0,6,7, 0,13,13,17,18,13,25,26};
                vector<int> point_j{1,2,3,6,7,8,13,15,17,18,19,25,26,27};
                
                for(unsigned i=0; i<point_i.size(); i++)
                {
                    if(i>0)
                    {
                        pose_V.conservativeResize(6,pose_V.cols()+2);
                    }
                    pose_V.col(i*2) << poses[0].coeffRef(0,point_i[i]), poses[0].coeffRef(1,point_i[i]), poses[0].coeffRef(2,point_i[i]), 1, 1, 1; 
                    pose_V.col(2*i+1) << poses[0].coeffRef(0,point_j[i]), poses[0].coeffRef(1,point_j[i]), poses[0].coeffRef(2,point_j[i]), 1, 1, 1;
                }
                // cout << pose_V << endl;
                pose_VBO.update(pose_V);
                break; 
            }
            case GLFW_KEY_1:
            {                
                vector<int> point_i{0,1,2,0,6,7, 0,13,13,17,18,13,25,26};
                vector<int> point_j{1,2,3,6,7,8,13,15,17,18,19,25,26,27};
                for(unsigned i=0; i<point_i.size(); i++)
                {
                    tri_mesh cube = load_mesh("../data/cube.off",Vector3f(0,0,0),1);
                    point arm1 = point(poses[0].coeffRef(0,point_i[i]),poses[0].coeffRef(1,point_i[i]),poses[0].coeffRef(2,point_i[i]));
                    point arm2 = point(poses[0].coeffRef(0,point_j[i]),poses[0].coeffRef(1,point_j[i]),poses[0].coeffRef(2,point_j[i]));
                    cube.M_model = cube_transform(arm1,arm2,0.05,0.05)*cube.M_model;
                    meshes.push_back(cube);
                    mesh_V_update(cube);
                }
                
                //update VBO
                mesh_VBO.update(mesh_V);
                break;
            }
            case GLFW_KEY_2:
            {
                tri_mesh cube2 = load_mesh("../data/bumpy_cube.off",Vector3f(0,0,0),0.11343);
                meshes.push_back(cube2);
                mesh_V_update(cube2);
                cout << meshes.size() << endl;
                //update VBO
                mesh_VBO.update(mesh_V);
                break;
            } 
            case GLFW_KEY_3:
            {
                tri_mesh bunny = load_mesh("../data/bunny.off",Vector3f(0.120779,-0.405439,-0.0351582),5*0.86);
                meshes.push_back(bunny);
                mesh_V_update(bunny);
                
                //update VBO
                mesh_VBO.update(mesh_V);
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
        if(clicked_mesh.clicked)
        {
            switch (key)
            {
                case  GLFW_KEY_J:
                {
                    clicked_mesh.shader_type = 'w';
                    
                    meshes[clicked_mesh.clicked_index].shader_type = 'w';
                    break;
                }
                case  GLFW_KEY_K:
                {
                    clicked_mesh.shader_type = 'f';
                    meshes[clicked_mesh.clicked_index].shader_type = 'f';

                    //vertex normals equal to face normals
                    int insert_start = 0;
                    for(unsigned i = 0; i < clicked_mesh.clicked_index; i++)
                    {
                        insert_start+= meshes[i].F.cols()*3;
                    }

                    for(unsigned i=0; i < clicked_mesh.F.cols();i++)
                    {
                        for(unsigned j=0; j<3;j++)
                        {
                            mesh_V.block(3,insert_start,3,1) = clicked_mesh.F.block(3,i,3,1);
                            insert_start++;
                        }
                    }

                    mesh_VBO.update(mesh_V);
                    break;
                }
                case  GLFW_KEY_L:
                {
                    clicked_mesh.shader_type = 'p';
                    meshes[clicked_mesh.clicked_index].shader_type = 'p';
                    int insert_start = 0;
                    for(unsigned i = 0; i < clicked_mesh.clicked_index; i++)
                    {
                        insert_start+= meshes[i].F.cols()*3;
                    }

                    //vertex normals equal to average normal of all faces
                    for(unsigned i=0; i<clicked_mesh.F.cols();i++)
                    {
                        for(unsigned j=0; j<3;j++)
                        {
                            mesh_V.col(insert_start) << clicked_mesh.V.col((int)clicked_mesh.F.coeffRef(j,i));
                            insert_start++;
                        }
                    }

                    mesh_VBO.update(mesh_V);
                    break;
                }
                case GLFW_KEY_SLASH:
                {
                    for(unsigned i = 0; i<3; i++)
                    {
                        meshes[clicked_mesh.clicked_index].M_model.coeffRef(i,i) = 
                        meshes[clicked_mesh.clicked_index].M_model.coeffRef(i,i)*1.25;
                    }
                    break;
                }
                case GLFW_KEY_PERIOD:
                {
                    for(unsigned i = 0; i<3; i++)
                    {
                        meshes[clicked_mesh.clicked_index].M_model.coeffRef(i,i) = 
                        meshes[clicked_mesh.clicked_index].M_model.coeffRef(i,i)*0.75;
                    }
                    break;
                }
                case GLFW_KEY_SEMICOLON:
                {
                    float radians = -10 * 3.141592f / 180;
                    Matrix4f rotation;
                    rotation << cos(radians), sin(radians), 0, 0,
                    -sin(radians), cos(radians), 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1;
                    meshes[clicked_mesh.clicked_index].M_model=rotation * meshes[clicked_mesh.clicked_index].M_model;

                    break;
                }
                case GLFW_KEY_APOSTROPHE:
                {
                    float radians = 10 * 3.141592f / 180;
                    Matrix4f rotation;
                    rotation << cos(radians), sin(radians), 0, 0,
                    -sin(radians), cos(radians), 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1;
                    meshes[clicked_mesh.clicked_index].M_model=rotation * meshes[clicked_mesh.clicked_index].M_model;

                    break;
                }
                case GLFW_KEY_DELETE:
                {
                    clicked_mesh.clicked = false;

                    int start = 0;
                    for(unsigned i = 0; i < clicked_mesh.clicked_index; i++)
                    {
                        start+= meshes[i].F.cols()*3;
                    }
                    for(unsigned i=0; i<clicked_mesh.F.cols()*3;i++)
                    {
                        removeColumn(mesh_V,start);
                    }

                    meshes.erase(meshes.begin()+clicked_mesh.clicked_index);
                    mesh_VBO.update(mesh_V);
                    break;
                }
                case GLFW_KEY_BACKSPACE:
                {
                    clicked_mesh.clicked = false;

                    int start = 0;
                    for(unsigned i = 0; i < clicked_mesh.clicked_index; i++)
                    {
                        start+= meshes[i].F.cols()*3;
                    }
                    for(unsigned i=0; i<clicked_mesh.F.cols()*3;i++)
                    {
                        removeColumn(mesh_V,start);
                    }

                    meshes.erase(meshes.begin()+clicked_mesh.clicked_index);
                    mesh_VBO.update(mesh_V);
                    break;
                }
                default:
                    break;
            }
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

    tri_VBO.init();
    tri_V.resize(5,1);
    tri_VBO.update(tri_V);

    mesh_VBO.init();
    mesh_V.resize(6,1);
    mesh_VBO.update(mesh_V);

    pose_VBO.init();
    pose_V.resize(6,2);
    pose_VBO.update(pose_V);
    // Initialize view scale and position
    view_scale << 1, 0, 0, 1;
    view_pos << 0, 0;

    light spotlight;
    spotlight.position << 0, 1, 1.0;
    spotlight.color << 1.0, 1.0, 1.0;
    spotlight.intensity << 1.0, 1.0, 1.0;

    //viewport transformation matrix
    // M_vp << width/2, 0, 0, (width-1)/2,
    // 0, height/2, 0, (height-1)/2,
    // 0, 0, 1, 0,
    // 0, 0, 0, 1;

    // M_vp = M_vp.inverse().eval();

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

    //perspective project matrix
    // P << near, 0, 0, 0,
    // 0, near, 0, 0,
    // 0, 0, near+far, -far*near,
    // 0, 0, 1, 0;

    // P << 2*abs(near)/(r-l), 0, (r+l)/(r-l), 0,
    //             0, abs(near)*2/(t-b), (t+b)/(t-b), 0,
    //             0, 0, (abs(near)+abs(far))/(abs(near)-abs(far)), 2*abs(far)*abs(near)/(abs(near)-abs(far)),
    //             0, 0, -1, 0;

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
                    "    gl_Position = viewportMatrix * projMatrix * pos;"
                    "}";
    const GLchar* fragment_shader =
            "#version 330 core\n"
                    "in vec4 normal;"
                    "in vec3 halfVec;"
                    "in vec3 lightDir;"
                    "uniform vec3 lightIntensity;"
                    "uniform int shaderMode;"
                    "uniform vec3 Ia;"
                    "uniform vec3 ka;"
                    "uniform vec3 kd;"
                    "uniform float ks;"
                    "uniform float phongExp;"
                    "out vec4 fragmentColor;"
                    "void main()"
                    "{"
                    "    if (shaderMode == 0) {"
                    "       fragmentColor = vec4(kd, 1.0);"
                    "    }"
                    "    else {"
                    "       if(shaderMode == 1){"
                    "       vec3 n = normalize(normal.xyz);"
                    "       vec3 l = normalize(lightDir);"
                    "       float diffuse =  max(0.0, dot(n,l));"
                    "       vec3 intensity = ka * Ia + kd * lightIntensity * (diffuse);"
                    "       fragmentColor = vec4(intensity, 1.0);"
                    "       }"
                    "       else{"
                    "       vec3 n = normalize(normal.xyz);"
                    "       vec3 h = normalize(halfVec);"
                    "       vec3 l = normalize(lightDir);"
                    "       float diffuse =  max(0.0, dot(n,l));"
                    "       float specular = ks * pow(max(0.0,dot(n,h)),phongExp);"
                    "       vec3 intensity = ka * Ia + kd * lightIntensity * (diffuse + specular);"
                    "       fragmentColor = vec4(intensity, 1.0);"
                    "       }"
                    "    }"
                    "}";

    axes_program.init(ax_vertex_shader,ax_fragment_shader,"fragmentColor");
    axes_program.bind();
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

        if(pose_V.cols()>1)
        {
            pose_VBO.bind();
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), 0);
            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (const GLvoid *)(12));
            glDrawArrays(GL_LINES,0,pose_V.cols());
        }
        // Bind program
        program.bind();

        // Set the uniform view matrix and translation vectors
        glUniformMatrix2fv(program.uniform("scale"),1, GL_FALSE, view_scale.data());
        glUniform2fv(program.uniform("translation"),1, view_pos.data());

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

        if(mesh_V.cols()>1)
        {
            mesh_VBO.bind();
            //for each mesh
            long start = 0;
                
            for(unsigned i = 0; i<meshes.size(); i++)
            {
                M_comb = M_cam * meshes[i].M_model;
                M_normal = M_comb.inverse().transpose();

                glUniformMatrix4fv(program.uniform("modelMatrix"),1, GL_FALSE, meshes[i].M_model.data());
                glUniformMatrix4fv(program.uniform("normalMatrix"),1, GL_FALSE, M_normal.data());

                //draw mesh elements
                Vector3f color_gl = meshes[i].diff_color/255;
                if(clicked_mesh.clicked && clicked_mesh.clicked_index == i)
                {
                    glUniform3fv(program.uniform("kd"),1, Vector3f(0,0,1).data());
                }
                else
                {
                    glUniform3fv(program.uniform("kd"),1, color_gl.data());
                }
                
                glEnableVertexAttribArray(0);
                glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (const GLvoid *)(24*start));

                if(meshes[i].shader_type == 'p')
                {
                    glUniform1i(program.uniform("shaderMode"),2);
                    Vector3f ka_gl = meshes[i].ka/255;
                    glUniform3fv(program.uniform("ka"),1, ka_gl.data());
                    glUniform1f(program.uniform("ks"), meshes[i].ks);
                    glUniform1f(program.uniform("phongExp"), meshes[i].phong_exp);

                    glEnableVertexAttribArray(1);
                    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (const GLvoid *)(24*start+12));
                }
                else if( meshes[i].shader_type == 'f')
                {
                    glUniform1i(program.uniform("shaderMode"),1);
                    Vector3f ka_gl = meshes[i].ka/255;
                    glUniform3fv(program.uniform("ka"),1, ka_gl.data());

                    glEnableVertexAttribArray(1);
                    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (const GLvoid *)(24*start+12));
                }
                else if (meshes[i].shader_type == 'w')
                {
                    glUniform1i(program.uniform("shaderMode"),0);
                }

                for(unsigned j=0; j<meshes[i].F.cols(); j++)
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
                    glUniform3fv(program.uniform("kd"),1, Vector3f(1,1,1).data());
                    glUniform1i(program.uniform("shaderMode"),0);
                    for(unsigned j=0; j<meshes[i].F.cols(); j++)
                    {
                        glDrawArrays(GL_LINE_LOOP,j*3,3);
                    }
                }
                start+= meshes[i].F.cols()*3;
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
