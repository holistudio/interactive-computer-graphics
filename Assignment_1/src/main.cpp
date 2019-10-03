// C++ include
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include<Eigen/Geometry>

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"

// Shortcut to avoid Eigen:: and std:: everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;

struct tri_mesh
{
    //basic triangular mesh
    //V is a matrix (dimension #V × 3 where #V is the number of vertices) that contains the positions of the vertices of the mesh
    MatrixXd V;
    //F is an matrix (dimension #faces × 3 where #F is the number of faces) which contains the descriptions of the triangles in the mesh. 
    MatrixXd F;
    vector<int> color;
};

std::ostream& operator<<(std::ostream &strm, const tri_mesh &a) {
  return strm << "Vertices\n" << a.V << "\nFaces\n" << a.F;
}

class ray
{
    public:
        Vector3d e;
        Vector3d d;
};

void part1()
{
    std::cout << "Part 1: Writing a grid png image" << std::endl;

    const std::string filename("part1.png");
    Eigen::MatrixXd M(800,800);
    //additional matrices for testing colors of checkerboards
    MatrixXd Z = MatrixXd::Zero(800,800);
    MatrixXd ones = MatrixXd::Ones(800,800);
    // Draw a grid, each square has a side of e pixels
    const int e = 50;
    const double black = 0;
    const double white = 1;

    for (unsigned wi = 0; wi<M.cols();++wi)
        for (unsigned hi = 0; hi < M.rows(); ++hi)
            M(hi,wi) = (lround(wi / e) % 2) == (lround(hi / e) % 2) ? black : white;

    // Write it in a png image. Note that the alpha channel is reversed to make the white (color = 1) pixels transparent (alhpa = 0)
    write_matrix_to_png(M,Z,Z,ones,filename);
}

void part2()
{
    std::cout << "Part 2: Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("part2.png");
    MatrixXd C = MatrixXd::Zero(800,800); // Store the color
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/C.cols(),0,0);
    Vector3d y_displacement(0,-2.0/C.rows(),0);

    // Single light source
    const Vector3d light_position(-1,1,1);

    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {
            // Prepare the ray
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            //Vector3d ray_direction = RowVector3d(0,0,-1);

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0),ray_origin(1));
            const double sphere_radius = 0.9;

            if (ray_on_xy.norm()<sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(ray_on_xy(0),ray_on_xy(1),sqrt(sphere_radius*sphere_radius - ray_on_xy.squaredNorm()));

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // Simple diffuse model
                C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i,j) = max(C(i,j),0.);

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C,C,C,A,filename);

}

void task1_1(vector<float> sphere_radii, MatrixXd sphere_centers)
{
    std::cout << "Task 1.1: Multiple Spheres" << std::endl;
    //camera coordinates
    double l = -1.0;
    double r = 1.0;
    double t = 1.0;
    double b = -1.0;

    const std::string filename("task1_1.png");
    MatrixXd C = MatrixXd::Zero(800,800); // Store the color
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d cam_origin(0,0,2);
    Vector3d x_displacement;
    Vector3d y_displacement;

    // Single light source
    const Vector3d light_position(-2.0,  2.0,  2.0);

    // Ray intersection parameter
    double t1 = 0;
    double t_new = 0;
    bool hit = false;

    //hit sphere (cx, cy, cz, R)
    double hit_sphere_radius;
    Vector3d hit_sphere_center;

    // For each ray with direction -w
    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {
            // Prepare the ray
            x_displacement = Vector3d(l+(r-l)*(i)/C.cols(),0,0);
            y_displacement = Vector3d(0,t-(t-b)*(j)/C.rows(),0);
            Vector3d ray_origin = cam_origin + x_displacement + y_displacement;
            Vector3d ray_direction = Vector3d(0,0,-1);
            hit = false;

            // For each sphere
            for (unsigned k=0;k<sphere_radii.size();k++)
            {

                const double sphere_radius = sphere_radii[k];
                Vector3d sphere_center(sphere_centers.coeffRef(k, 0),sphere_centers.coeffRef(k, 1),sphere_centers.coeffRef(k, 2));

                //Calculate discriminant
                Vector3d e_c = ray_origin - sphere_center;
                //quadratic equation
                double a = ray_direction.dot(ray_direction);
                double b = ray_direction.dot(e_c);
                double c = e_c.dot(e_c) - sphere_radius * sphere_radius;

                //discriminant
                double discriminant = b*b - a*c;

                if(discriminant>=0)
                {
                    //if discriminant is greater than or equal to 0
                    //Calculate sphere intersection parameter, t_new
                    if(discriminant ==0)
                    {
                        t_new = -b/a;
                    }
                    else
                    {
                        t_new = (-b - sqrt(discriminant)) / a;
                        if(t_new < 0)
                        {
                            t_new = (-b + sqrt(discriminant)) / a;
                        }
                        else
                        {
                            double t2 = (-b + sqrt(discriminant)) / a;
                            if (t2>0 && t2<t_new)
                            {
                                t_new = t2;
                            }
                        }
                    }

                    if(hit == false)
                    {
                        //this must be the first ray hit
                        t1 = t_new;
                        hit_sphere_center = sphere_center;
                        hit_sphere_radius = sphere_radius;
                        hit = true;
                    }
                    else
                    {
                        //this means there's already a hit from another sphere
                        //If that t is less than stored t, then replace stored t with the new t
                        if(t_new<t1)
                        {
                            t1 = t_new;
                            hit_sphere_center = sphere_center;
                            hit_sphere_radius = sphere_radius;
                        }
                    }

                }

            }
            //if hit is true, then do lighting calculations
            if (hit == true)
            {
                Vector3d ray_intersection = ray_origin + t1 * ray_direction;
                Vector3d ray_normal = (ray_intersection-hit_sphere_center)/hit_sphere_radius;

                // Simple diffuse model
                C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i,j) = max(C(i,j),0.);

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }

        }
    }
    write_matrix_to_png(C,C,C,A,filename);
}

void task1_2(vector<float> sphere_radii, MatrixXd sphere_centers, MatrixXd sphere_colors, vector<char> sphere_shading)
{
    //this function alternates between diffuse shading and specular shading
    //e.g sphere 1 gets diffuse shading, sphere 2 gets specular, sphere 3 gets diffuse...
    std::cout << "Task 1.2: Multiple Shades and Colors" << std::endl;
    //camera coordinates
    double l = -1.0;
    double r = 1.0;
    double t = 1.0;
    double b = -1.0;

    const std::string filename("task1_2.png");
    MatrixXd R = MatrixXd::Zero(800,800); // Store the red
    MatrixXd G = MatrixXd::Zero(800,800); // green
    MatrixXd B = MatrixXd::Zero(800,800); // blue
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d cam_origin(0,0,2);
    Vector3d x_displacement;
    Vector3d y_displacement;

    // Multiple Light Sources
    MatrixXd light_positions(2,3);
    light_positions <<  -0.8,  1.0,  1.2,
                         1.0, -1.0, -1.0;

    // Ray intersection parameter
    double t1 = 0;
    double t_new = 0;
    bool hit = false;

    //hit sphere properties
    double hit_sphere_radius;
    Vector3d hit_sphere_center;
    Vector3d hit_sphere_color;
    char hit_sphere_shading;

    // For each ray with direction -w
    for (unsigned i=0;i<R.cols();i++)
    {
        for (unsigned j=0;j<R.rows();j++)
        {
            // Prepare the ray
            x_displacement = Vector3d(l+(r-l)*(i)/R.cols(),0,0);
            y_displacement = Vector3d(0,t-(t-b)*(j)/R.rows(),0);
            Vector3d ray_origin = cam_origin + x_displacement + y_displacement;
            Vector3d ray_direction = Vector3d(0,0,-1);
            hit = false;

            // For each sphere
            for (unsigned k=0;k<sphere_radii.size();k++)
            {
                const double sphere_radius = sphere_radii[k];
                Vector3d sphere_center(sphere_centers.coeffRef(k, 0),sphere_centers.coeffRef(k, 1),sphere_centers.coeffRef(k, 2));

                //Calculate discriminant
                Vector3d e_c = ray_origin - sphere_center;
                //quadratic equation
                double a = ray_direction.dot(ray_direction);
                double b = ray_direction.dot(e_c);
                double c = e_c.dot(e_c) - sphere_radius * sphere_radius;

                //discriminant
                double discriminant = b*b - a*c;

                if(discriminant>=0)
                {
                    //if discriminant is greater than or equal to 0
                    //Calculate sphere intersection parameter, t_new
                    if(discriminant ==0)
                    {
                        t_new = -b/a;
                    }
                    else
                    {
                        t_new = (-b - sqrt(discriminant)) / a;
                        if(t_new < 0)
                        {
                            t_new = (-b + sqrt(discriminant)) / a;
                        }
                        else
                        {
                            double t2 = (-b + sqrt(discriminant)) / a;
                            if (t2>0 && t2<t_new)
                            {
                                t_new = t2;
                            }
                        }
                    }
                    if(hit == false)
                    {
                        //this must be the first ray hit
                        t1 = t_new;
                        hit_sphere_center = sphere_center;
                        hit_sphere_radius = sphere_radius;
                        hit_sphere_color = Vector3d(sphere_colors.row(k));
                        hit_sphere_shading = sphere_shading[k];
                        hit = true;
                    }
                    else
                    {
                        //this means there's already a hit from another sphere
                        //If that t is less than stored t, then replace stored t with the new t
                        if(t_new<t1)
                        {
                            t1 = t_new;
                            hit_sphere_center = sphere_center;
                            hit_sphere_radius = sphere_radius;
                            hit_sphere_color = Vector3d(sphere_colors.row(k));
                            hit_sphere_shading = sphere_shading[k];
                        }
                    }
                }

            }

            //if hit is true, then do lighting calculations
            if (hit == true)
            {
                Vector3d ray_intersection = ray_origin + t1 * ray_direction;
                Vector3d ray_normal = (ray_intersection-hit_sphere_center)/hit_sphere_radius;
                double color_intensity=0.0;

                //for each light source
                for(unsigned m=0;m<light_positions.rows();m++)
                {
                    // Simple diffuse model assuming light intensity I=1
                    Vector3d light_vector = (Vector3d(light_positions.row(m))-ray_intersection).normalized();
                    color_intensity = light_vector.transpose() * ray_normal;
                    // Clamp to zero
                    color_intensity = max(color_intensity,0.);

                    //Specular shading, assuming Phong exponent p=100
                    if(hit_sphere_shading=='s')
                    {
                        Vector3d half_vector=(-ray_direction+light_vector).normalized();
                        color_intensity = color_intensity + pow(max(ray_normal.dot(half_vector),0.),100);
                    }
                    R(i,j) = R(i,j)+color_intensity;
                    G(i,j) = G(i,j)+color_intensity;
                    B(i,j) = B(i,j)+color_intensity;
                }

                //adjust RGB based on the hit sphere's color
                //assume kd and ks are same, based on RGB values provided in arguments
                R(i,j)=R(i,j)*hit_sphere_color.coeff(0)/255;
                G(i,j)=G(i,j)*hit_sphere_color.coeff(1)/255;
                B(i,j)=B(i,j)*hit_sphere_color.coeff(2)/255;

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }

        }
    }
    write_matrix_to_png(R,G,B,A,filename);
}

void task1_3(vector<float> sphere_radii, MatrixXd sphere_centers, MatrixXd sphere_colors, vector<char> sphere_shading)
{
    //Duplicate of Task 1.2 but with perspective
    std::cout << "Task 1.3: Perspective" << std::endl;
    //camera coordinates
    double l = -1.0;
    double r = 1.0;
    double t = 1.0;
    double b = -1.0;

    const std::string filename("task1_3.png");
    MatrixXd R = MatrixXd::Zero(800,800); // Store the red
    MatrixXd G = MatrixXd::Zero(800,800); // green
    MatrixXd B = MatrixXd::Zero(800,800); // blue
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // Perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d cam_origin(0,0,2);
    Vector3d cam_d = Vector3d(0,0,-1.5);
    Vector3d x_displacement;
    Vector3d y_displacement;

    // Multiple Light Sources
    MatrixXd light_positions(2,3);
    light_positions <<  -0.8,  1.0,  1.2,
                         1.0, -1.0, -1.0;

    // Ray intersection parameter
    double t1 = 0;
    double t_new = 0;
    bool hit = false;

    //hit sphere properties
    double hit_sphere_radius;
    Vector3d hit_sphere_center;
    Vector3d hit_sphere_color;
    char hit_sphere_shading;

    // For each ray with direction -w
    for (unsigned i=0;i<R.cols();i++)
    {
        for (unsigned j=0;j<R.rows();j++)
        {
            // Prepare the ray, perspective camera
            x_displacement = Vector3d(l+(r-l)*(i)/R.cols(),0,0);
            y_displacement = Vector3d(0,t-(t-b)*(j)/R.rows(),0);
            Vector3d ray_origin = cam_origin;
            Vector3d ray_direction = cam_d+x_displacement+y_displacement;
            hit = false;

            // For each sphere
            for (unsigned k=0;k<sphere_radii.size();k++)
            {
                const double sphere_radius = sphere_radii[k];
                Vector3d sphere_center(sphere_centers.coeffRef(k, 0),sphere_centers.coeffRef(k, 1),sphere_centers.coeffRef(k, 2));

                //Calculate discriminant
                Vector3d e_c = ray_origin - sphere_center;
                //quadratic equation
                double a = ray_direction.dot(ray_direction);
                double b = ray_direction.dot(e_c);
                double c = e_c.dot(e_c) - sphere_radius * sphere_radius;

                //discriminant
                double discriminant = b*b - a*c;

                if(discriminant>=0)
                {
                    //if discriminant is greater than or equal to 0
                    //Calculate sphere intersection parameter, t_new
                    if(discriminant ==0)
                    {
                        t_new = -b/a;
                    }
                    else
                    {
                        t_new = (-b - sqrt(discriminant)) / a;
                        if(t_new < 0)
                        {
                            t_new = (-b + sqrt(discriminant)) / a;
                        }
                        else
                        {
                            double t2 = (-b + sqrt(discriminant)) / a;
                            if (t2>0 && t2<t_new)
                            {
                                t_new = t2;
                            }
                        }
                    }
                    if(hit == false)
                    {
                        //this must be the first ray hit
                        t1 = t_new;
                        hit_sphere_center = sphere_center;
                        hit_sphere_radius = sphere_radius;
                        hit_sphere_color = Vector3d(sphere_colors.row(k));
                        hit_sphere_shading = sphere_shading[k];
                        hit = true;
                    }
                    else
                    {
                        //this means there's already a hit from another sphere
                        //If that t is less than stored t, then replace stored t with the new t
                        if(t_new<t1)
                        {
                            t1 = t_new;
                            hit_sphere_center = sphere_center;
                            hit_sphere_radius = sphere_radius;
                            hit_sphere_color = Vector3d(sphere_colors.row(k));
                            hit_sphere_shading = sphere_shading[k];
                        }
                    }
                }

            }

            //if hit is true, then do lighting calculations
            if (hit == true)
            {
                Vector3d ray_intersection = ray_origin + t1 * ray_direction;
                Vector3d ray_normal = (ray_intersection-hit_sphere_center)/hit_sphere_radius;
                double color_intensity=0.0;

                //for each light source
                for(unsigned m=0;m<light_positions.rows();m++)
                {
                    Vector3d light_vector = (Vector3d(light_positions.row(m))-ray_intersection).normalized();
                    
                    // Simple diffuse model assuming light intensity I=1
                    color_intensity = light_vector.transpose() * ray_normal;
                    // Clamp to zero
                    color_intensity = max(color_intensity,0.);

                    //Specular shading, assuming Phong exponent p=100
                    if(hit_sphere_shading=='s')
                    {
                        Vector3d half_vector=(-ray_direction+light_vector).normalized();
                        color_intensity = color_intensity + pow(max(ray_normal.dot(half_vector),0.),100);
                    }
                    R(i,j) = R(i,j)+color_intensity;
                    G(i,j) = G(i,j)+color_intensity;
                    B(i,j) = B(i,j)+color_intensity;
                }

                //adjust RGB based on the hit sphere's color
                //assume kd and ks are same, based on RGB values provided in arguments
                R(i,j)=R(i,j)*hit_sphere_color.coeff(0)/255;
                G(i,j)=G(i,j)*hit_sphere_color.coeff(1)/255;
                B(i,j)=B(i,j)*hit_sphere_color.coeff(2)/255;

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }

        }
    }
    write_matrix_to_png(R,G,B,A,filename);
}

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

tri_mesh load_mesh(string off_filepath)
{
    // returns a triangule mesh object at the file path for the OFF file
    tri_mesh mesh_structure;

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
        mesh_structure.V.resize(num_vertices,3);
        mesh_structure.F.resize(num_faces,3);

        for (unsigned i=0;i<num_vertices;i++)
        {
            //load matrix V with vertices
            in.getline(str,255);
            substrings = space_sep_string(string(str));
            mesh_structure.V.row(i) << stof(substrings[0])*2,stof(substrings[1])*2,stof(substrings[2])*2;
        }

        for (unsigned i=0;i<num_faces;i++)
        {
            in.getline(str,255);
            substrings = space_sep_string(string(str));
            mesh_structure.F.row(i) << stoi(substrings[1]),stoi(substrings[2]),stoi(substrings[3]);
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

double tri_intersection(ray r, Vector3d v1, Vector3d v2, Vector3d v3)
{
    //Cramer's rule
    double a = v1[0] - v2[0];
    double b = v1[1] - v2[1];
    double c = v1[2] - v2[2];
    double d = v1[0] - v3[0];
    double e = v1[1] - v3[1];
    double f = v1[2] - v3[2];
    double g = r.d[0];
    double h = r.d[1];
    double i = r.d[2];
    double j = v1[0]-r.e[0];
    double k = v1[1]-r.e[1];
    double l = v1[2]-r.e[2];

    double ei_minus_hf = e*i - h*f;
    double gf_minus_di = g*f - d*i;
    double dh_minus_eg = d*h - e*g;

    double M = a*(ei_minus_hf) + b*(gf_minus_di) + c*(dh_minus_eg);
    double t = -(f*(a*k-j*b)+e*(j*c-a*l)+d*(b*l-k*c))/M;
    if(t<0)
    {
        return -1;
    }
    double gamma = (i*(a*k-j*b)+h*(j*c-a*l)+g*(b*l-k*c))/M;
    if(gamma < 0 || gamma >1){
        return -1;
    }
    double beta = (j*ei_minus_hf+k*gf_minus_di+l*dh_minus_eg)/M;
    if(beta < 0  || beta > 1-gamma){
        return -1;
    }
    return t;
}

class model
{
    public:
        Vector3d position;
        Vector3d color;
        Vector3d getPosition(){
            return position;
        }
};

class sphere: public model
{
    public:
        double radius;
};

class triangle: public model
{
    public:
        Vector3d v1;
        Vector3d v2;
        Vector3d v3;
};

class hit_object
{
    public:
        bool hit = false;
        double t;
        model hit_object;
};

// hit_object ray_hit(ray r, sphere m, string sphere_or_triangle)
// {
//     hit_object test;
//     if(sphere_or_triangle.compare("sphere")==0)
//     {
//         const double sphere_radius = sphere(m).radius;
//         Vector3d sphere_center = m.getPosition();

//         //Calculate discriminant
//         Vector3d e_c = r.e - sphere_center;
//         //quadratic equation
//         double a = ray_direction.dot(r.d);
//         double b = ray_direction.dot(e_c);
//         double c = e_c.dot(e_c) - sphere_radius * sphere_radius;

//         //discriminant
//         double discriminant = b*b - a*c;

//         if(discriminant>=0)
//         {
//             //if discriminant is greater than or equal to 0
//             //Calculate sphere intersection parameter, t_new
//             if(discriminant ==0)
//             {
//                 t_new = -b/a;
//             }
//             else
//             {
//                 t_new = (-b - sqrt(discriminant)) / a;
//                 if(t_new < 0)
//                 {
//                     t_new = (-b + sqrt(discriminant)) / a;
//                 }
//                 else
//                 {
//                     double t2 = (-b + sqrt(discriminant)) / a;
//                     if (t2>0 && t2<t_new)
//                     {
//                         t_new = t2;
//                     }
//                 }
//             }
//         }
//     }
//     return test;
// }

void task1_4(tri_mesh mesh_struct)
{
    //TODO: task1_4 needs to take an array of mesh structures
    std::cout << "Task 1.4: Triangle mesh" << std::endl;
    //camera coordinates
    double l = -1.0;
    double r = 1.0;
    double t = 1.0;
    double b = -1.0;

    const std::string filename("task1_4.png");
    MatrixXd R = MatrixXd::Zero(500,500); // Store the red
    MatrixXd G = MatrixXd::Zero(500,500); // green
    MatrixXd B = MatrixXd::Zero(500,500); // blue
    MatrixXd A = MatrixXd::Zero(500,500); // Store the alpha mask

    // Perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d cam_origin(0,0,0.4);
    Vector3d cam_d = Vector3d(0,0,-0.35);
    Vector3d x_displacement;
    Vector3d y_displacement;

    // Multiple Light Sources
    MatrixXd light_positions(2,3);
    light_positions <<  -0.8,  1.0,  1.2,
                         1.0, -1.0, -1.0;

    // Ray intersection parameter
    double t1 = 0;
    double t_new = 0;
    bool hit = false;

    //hit triangle properties
    int hit_mesh;
    int hit_face;
    Vector3d v1;
    Vector3d v2;
    Vector3d v3;
    vector<int> hit_color;

    // For each ray
    for (unsigned i=0;i<R.cols();i++)
    {
        for (unsigned j=0;j<R.rows();j++)
        {
            // Prepare the ray, perspective camera
            //std::cout << "Ray " << i << "," << j << std::endl;
            x_displacement = Vector3d(l+(r-l)*(i)/R.cols(),0,0);
            y_displacement = Vector3d(0,t-(t-b)*(j)/R.rows(),0);
            ray new_ray;
            new_ray.e = cam_origin;
            new_ray.d = cam_d+x_displacement+y_displacement;

            hit = false;

            // For each mesh
            // For each mesh triangle face
            for (unsigned k=0;k<mesh_struct.F.rows();k++)
            {
                v1 = Vector3d(mesh_struct.V.row(mesh_struct.F.coeffRef(k,0)).coeff(0),
                        mesh_struct.V.row(mesh_struct.F.coeffRef(k,0)).coeff(1),
                        mesh_struct.V.row(mesh_struct.F.coeffRef(k,0)).coeff(2));
                v2 = Vector3d(mesh_struct.V.row(mesh_struct.F.coeffRef(k,1)).coeff(0),
                        mesh_struct.V.row(mesh_struct.F.coeffRef(k,1)).coeff(1),
                        mesh_struct.V.row(mesh_struct.F.coeffRef(k,1)).coeff(2));
                v3 = Vector3d(mesh_struct.V.row(mesh_struct.F.coeffRef(k,2)).coeff(0),
                        mesh_struct.V.row(mesh_struct.F.coeffRef(k,2)).coeff(1),
                        mesh_struct.V.row(mesh_struct.F.coeffRef(k,2)).coeff(2));
                t_new = tri_intersection(new_ray,v1,v2,v3);
                if (t_new>=0)
                {
                    if(hit == false)
                    {
                        //this must be the first ray hit
                        t1 = t_new;
                        //store hit_mesh
                        hit_face = k;
                        //store color of hit_face to hit_color
                        hit = true;
                    }
                    else
                    {
                        //this means there's already a hit on another triangle face
                        //If that t is less than stored t, then replace stored t with the new t
                        if(t_new<t1)
                        {
                            t1 = t_new;
                            //store hit_mesh
                            hit_face = k;
                            //store color of hit_face to hit_color
                        }
                    }
                }

            }

            //if hit is true, then do lighting calculations
            if (hit)
            {
                Vector3d ray_intersection = new_ray.e + t1 * new_ray.d;

                Vector3d h1 = Vector3d(mesh_struct.V.row(mesh_struct.F.coeffRef(hit_face,0)).coeff(0),
                mesh_struct.V.row(mesh_struct.F.coeffRef(hit_face,0)).coeff(1),
                mesh_struct.V.row(mesh_struct.F.coeffRef(hit_face,0)).coeff(2));
                Vector3d h2 = Vector3d(mesh_struct.V.row(mesh_struct.F.coeffRef(hit_face,1)).coeff(0),
                mesh_struct.V.row(mesh_struct.F.coeffRef(hit_face,1)).coeff(1),
                mesh_struct.V.row(mesh_struct.F.coeffRef(hit_face,1)).coeff(2));
                Vector3d h3 = Vector3d(mesh_struct.V.row(mesh_struct.F.coeffRef(hit_face,2)).coeff(0),
                mesh_struct.V.row(mesh_struct.F.coeffRef(hit_face,2)).coeff(1),
                mesh_struct.V.row(mesh_struct.F.coeffRef(hit_face,2)).coeff(2));
                Vector3d ray_normal = ((h2-h1).cross(h3-h1)).normalized();
                double color_intensity=0.0;

                //for each light source
                for(unsigned m=0;m<light_positions.rows();m++)
                {
                    // Simple diffuse model assuming light intensity I=1
                    Vector3d light_vector = (Vector3d(light_positions.row(m))-ray_intersection).normalized();
                    color_intensity = light_vector.transpose() * ray_normal;
                    // Clamp to zero
                    color_intensity = max(color_intensity,0.);

                    R(i,j) = R(i,j)+color_intensity;
                    G(i,j) = G(i,j)+color_intensity;
                    B(i,j) = B(i,j)+color_intensity;
                }

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }

        }
    }
    write_matrix_to_png(R,G,B,A,filename);
}

bool surface_light(ray r, vector<float> sphere_radii, MatrixXd sphere_centers)
{
    double t1 = 0;
    double epsilon = 0.001;
    for (unsigned k=0;k<sphere_radii.size();k++)
    {
        const double sphere_radius = sphere_radii[k];
        Vector3d sphere_center(sphere_centers.coeffRef(k, 0),sphere_centers.coeffRef(k, 1),sphere_centers.coeffRef(k, 2));

        //Calculate discriminant
        Vector3d e_c = r.e - sphere_center;
        //quadratic equation
        double a = (r.d).dot(r.d);
        double b = (r.d).dot(e_c);
        double c = e_c.dot(e_c) - sphere_radius * sphere_radius;

        //discriminant
        double discriminant = b*b - a*c;
        
        if(discriminant>=0)
        {
            if(discriminant ==0)
            {
                t1 = -b/a;
            }
            else
            {
                t1 = (-b - sqrt(discriminant)) / a;
                if(t1 < 0)
                {
                    t1 = (-b + sqrt(discriminant)) / a;
                }
                else
                {
                    double t2 = (-b + sqrt(discriminant)) / a;
                    if (t2>0 && t2<t1)
                    {
                        t1 = t2;
                    }
                }
            }
            if (t1 > epsilon){
                return true;
            }
        }
    }
    //std::cout << "No hits to spheres" << std::endl;
    return false;
}

void task1_5(vector<float> sphere_radii, MatrixXd sphere_centers, MatrixXd sphere_colors, vector<char> sphere_shading)
{
    std::cout << "Task 1.5: Shadows" << std::endl;
    //camera coordinates
    double l = -1.0;
    double r = 1.0;
    double t = 1.0;
    double b = -1.0;

    const std::string filename("task1_5.png");
    MatrixXd R = MatrixXd::Zero(800,800); // Store the red
    MatrixXd G = MatrixXd::Zero(800,800); // green
    MatrixXd B = MatrixXd::Zero(800,800); // blue
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // Perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d cam_origin(0,0,2);
    Vector3d cam_d = Vector3d(0,0,-1.5);
    Vector3d x_displacement;
    Vector3d y_displacement;

    // Multiple Light Sources
    MatrixXd light_positions(2,3);
    light_positions <<  -2.0,  2.0,  2.0,
                         1.0, -1.0, -1.0;

    // Ray intersection parameter
    double t1 = 0;
    double t_new = 0;
    bool hit = false;

    //hit sphere properties
    double hit_sphere_radius;
    Vector3d hit_sphere_center;
    Vector3d hit_sphere_color;
    char hit_sphere_shading;

    // For each ray with direction -w
    for (unsigned i=0;i<R.cols();i++)
    {
        for (unsigned j=0;j<R.rows();j++)
        {
            // Prepare the ray, perspective camera
            x_displacement = Vector3d(l+(r-l)*(i)/R.cols(),0,0);
            y_displacement = Vector3d(0,t-(t-b)*(j)/R.rows(),0);
            Vector3d ray_origin = cam_origin;
            Vector3d ray_direction = cam_d+x_displacement+y_displacement;
            hit = false;

            // For each sphere
            for (unsigned k=0;k<sphere_radii.size();k++)
            {
                const double sphere_radius = sphere_radii[k];
                Vector3d sphere_center(sphere_centers.coeffRef(k, 0),sphere_centers.coeffRef(k, 1),sphere_centers.coeffRef(k, 2));

                //Calculate discriminant
                Vector3d e_c = ray_origin - sphere_center;
                //quadratic equation
                double a = ray_direction.dot(ray_direction);
                double b = ray_direction.dot(e_c);
                double c = e_c.dot(e_c) - sphere_radius * sphere_radius;

                //discriminant
                double discriminant = b*b - a*c;

                if(discriminant>=0)
                {
                    //if discriminant is greater than or equal to 0
                    //Calculate sphere intersection parameter, t_new
                    if(discriminant ==0)
                    {
                        t_new = -b/a;
                    }
                    else
                    {
                        t_new = (-b - sqrt(discriminant)) / a;
                        if(t_new < 0)
                        {
                            t_new = (-b + sqrt(discriminant)) / a;
                        }
                        else
                        {
                            double t2 = (-b + sqrt(discriminant)) / a;
                            if (t2>0 && t2<t_new)
                            {
                                t_new = t2;
                            }
                        }
                    }
                    if(hit == false)
                    {
                        //this must be the first ray hit
                        t1 = t_new;
                        hit_sphere_center = sphere_center;
                        hit_sphere_radius = sphere_radius;
                        hit_sphere_color = Vector3d(sphere_colors.row(k));
                        hit_sphere_shading = sphere_shading[k];
                        hit = true;
                    }
                    else
                    {
                        //this means there's already a hit from another sphere
                        //If that t is less than stored t, then replace stored t with the new t
                        if(t_new<t1)
                        {
                            t1 = t_new;
                            hit_sphere_center = sphere_center;
                            hit_sphere_radius = sphere_radius;
                            hit_sphere_color = Vector3d(sphere_colors.row(k));
                            hit_sphere_shading = sphere_shading[k];
                        }
                    }
                }

            }

            //if hit is true, then do lighting calculations
            if (hit == true)
            {
                Vector3d ray_intersection = ray_origin + t1 * ray_direction;
                Vector3d ray_normal = (ray_intersection-hit_sphere_center)/hit_sphere_radius;

                double color_intensity=0.0;

                //for each light source
                for(unsigned m=0;m<light_positions.rows();m++)
                {
                    Vector3d light_vector = (Vector3d(light_positions.row(m))-ray_intersection);
                    ray shadow_ray;
                    shadow_ray.e = ray_intersection;
                    shadow_ray.d = light_vector;

                    if(surface_light(shadow_ray,sphere_radii,sphere_centers)==false)
                    {
                        light_vector = light_vector.normalized();
                        // Simple diffuse model assuming light intensity I=1
                        color_intensity = light_vector.transpose() * ray_normal;
                        // Clamp to zero
                        color_intensity = max(color_intensity,0.);

                        //Specular shading, assuming Phong exponent p=100
                        if(hit_sphere_shading=='s')
                        {
                            Vector3d half_vector=(-ray_direction+light_vector).normalized();
                            color_intensity = color_intensity + pow(max(ray_normal.dot(half_vector),0.),100);
                        }
                        R(i,j) = R(i,j)+color_intensity;
                        G(i,j) = G(i,j)+color_intensity;
                        B(i,j) = B(i,j)+color_intensity;
                    }
                }

                //adjust RGB based on the hit sphere's color
                //assume kd and ks are same, based on RGB values provided in arguments
                R(i,j)=R(i,j)*hit_sphere_color.coeff(0)/255;
                G(i,j)=G(i,j)*hit_sphere_color.coeff(1)/255;
                B(i,j)=B(i,j)*hit_sphere_color.coeff(2)/255;

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }

        }
    }
    write_matrix_to_png(R,G,B,A,filename);
}

int main()
{
    //part1();
    //part2();
    vector<float> spheresRadii;
    MatrixXd spheresCenters(2,3);
    spheresRadii.push_back(0.5);
    spheresRadii.push_back(0.1);
    spheresCenters << 0, 0, 0,
                      -0.6, 0.6, 0.6;
    task1_1(spheresRadii,spheresCenters);
    MatrixXd spheres_colors(2,3); //in RGB
    spheres_colors << 66, 135, 245,
                     0, 255, 166;
    vector<char> sphere_shading;
    sphere_shading.push_back('d');
    sphere_shading.push_back('s');
    //task1_2(spheresRadii,spheresCenters,spheres_colors,sphere_shading);
    //task1_3(spheresRadii,spheresCenters,spheres_colors,sphere_shading);

    tri_mesh ground_plane;
    //TODO: create a simple ground plane...in a separate OFF file.

    //tri_mesh test_mesh = load_mesh("../data/bunny.off");

    //add mesh to mesh array

    // task1_4(test_mesh);
    task1_5(spheresRadii,spheresCenters,spheres_colors,sphere_shading);
    return 0;
}
