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

class model
{
    public:
        Vector3d position;
        Vector3d color;
        char shader_type;
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
        bool mirror = false;
};

class hit_sphere
{
    public:
        bool hit = false;
        double t;
        sphere hit_obj;
};

class hit_triangle
{
    public:
        bool hit = false;
        double t;
        triangle hit_obj;
};


class tri_mesh
{
    //basic triangular mesh
    public:
        //V is a matrix (dimension #V × 3 where #V is the number of vertices) that contains the positions of the vertices of the mesh
        MatrixXd V;
        //F is an matrix (dimension #faces × 3 where #F is the number of faces) which contains the descriptions of the triangles in the mesh. 
        MatrixXd F;
        Vector3d color;
        char shader_type;
        bool mirror=false;
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


hit_sphere ray_hit_sphere(ray r, double t_lower, sphere test_sphere)
{
    hit_sphere test;
    const double sphere_radius = test_sphere.radius;
    Vector3d sphere_center = test_sphere.position;

    //Calculate discriminant
    Vector3d e_c = r.e - sphere_center;
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

hit_triangle ray_hit_triangle(ray r, double t_lower, triangle test_triangle)
{
    hit_triangle test;
    Vector3d v1 = test_triangle.v1;
    Vector3d v2 = test_triangle.v2;
    Vector3d v3 = test_triangle.v3;
    //Cramer's rule
    double a = v1.coeff(0) - v2.coeff(0);
    double b = v1.coeff(1) - v2.coeff(1);
    double c = v1.coeff(2) - v2.coeff(2);
    double d = v1.coeff(0) - v3.coeff(0);
    double e = v1.coeff(1) - v3.coeff(1);
    double f = v1.coeff(2) - v3.coeff(2);
    double g = r.d.coeff(0);
    double h = r.d.coeff(1);
    double i = r.d.coeff(2);
    double j = v1.coeff(0)-r.e.coeff(0);
    double k = v1.coeff(1)-r.e.coeff(1);
    double l = v1.coeff(2)-r.e.coeff(2);

    double ei_minus_hf = e*i - h*f;
    double gf_minus_di = g*f - d*i;
    double dh_minus_eg = d*h - e*g;

    double M = a*(ei_minus_hf) + b*(gf_minus_di) + c*(dh_minus_eg);
    double t = -(f*(a*k-j*b)+e*(j*c-a*l)+d*(b*l-k*c))/M;
    if(t<t_lower)
    {
        return test;
    }
    double gamma = (i*(a*k-j*b)+h*(j*c-a*l)+g*(b*l-k*c))/M;
    if(gamma < 0 || gamma >1){
        return test;
    }
    double beta = (j*ei_minus_hf+k*gf_minus_di+l*dh_minus_eg)/M;
    if(beta < 0  || beta > 1-gamma){
        return test;
    }
    test.hit = true;
    test.hit_obj = test_triangle;
    test.t = t;
    return test;
}

void task1_1(vector<sphere> spheres)
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

    // For each ray with direction -w
    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {
            // Prepare the ray
            x_displacement = Vector3d(l+(r-l)*(i)/C.cols(),0,0);
            y_displacement = Vector3d(0,t-(t-b)*(j)/C.rows(),0);
            
            ray camera_ray;
            camera_ray.e = cam_origin + x_displacement + y_displacement;
            camera_ray.d = Vector3d(0,0,-1);

            hit_sphere hit_record;

            // For each sphere
            for (unsigned k=0;k<spheres.size();k++)
            {
                sphere test_sphere;
                test_sphere.radius = spheres[k].radius;
                test_sphere.position = spheres[k].position;
                test_sphere.color = spheres[k].color;
                test_sphere.shader_type = spheres[k].shader_type;
                hit_sphere hit_test = ray_hit_sphere(camera_ray,0,test_sphere);
                if(hit_test.hit)
                {
                    if(hit_record.hit==false)
                    {
                        hit_record = hit_test;
                    }
                    else
                    {
                        if(hit_test.t<hit_record.t)
                        {
                            hit_record = hit_test;
                        }
                    }

                }
            }
            if (hit_record.hit)
            {
                string hit_type;
                double t_int = hit_record.t;

                Vector3d ray_intersection = camera_ray.e + t_int * camera_ray.d;
                
                Vector3d ray_normal;
                if(hit_type.compare("sphere")==0)
                {
                    ray_normal = (ray_intersection-hit_record.hit_obj.position)/hit_record.hit_obj.radius;
                }

                double color_intensity=0.0;

                    Vector3d light_vector = light_position-ray_intersection;
                    ray shadow_ray;
                    shadow_ray.e = ray_intersection;
                    shadow_ray.d = light_vector;

                    light_vector = light_vector.normalized();
                    // Simple diffuse model assuming light intensity I=1
                    // Clamp to zero
                    color_intensity = light_vector.transpose() * ray_normal;
                    color_intensity = max(color_intensity,0.);

                    //Specular shading, assuming Phong exponent p=100
                    if(hit_record.hit_obj.shader_type=='s')
                    {
                        Vector3d half_vector=(-camera_ray.d+light_vector).normalized();
                        color_intensity = color_intensity + pow(max(ray_normal.dot(half_vector),0.),100);
                    }
                
                    C(i,j) = C(i,j)+color_intensity;


                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }

        }
    }
    write_matrix_to_png(C,C,C,A,filename);
}

void task1_2(vector<sphere> spheres)
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
            for (unsigned k=0;k<spheres.size();k++)
            {
                const double sphere_radius = spheres[k].radius;
                Vector3d sphere_center = spheres[k].position;

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
                        hit_sphere_color = spheres[k].color;
                        hit_sphere_shading = spheres[k].shader_type;
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
                            hit_sphere_color = spheres[k].color;
                            hit_sphere_shading = spheres[k].shader_type;
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

void task1_3(vector<sphere> spheres)
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
            for (unsigned k=0;k<spheres.size();k++)
            {
                const double sphere_radius = spheres[k].radius;
                Vector3d sphere_center = spheres[k].position;

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
                        hit_sphere_color = spheres[k].color;
                        hit_sphere_shading = spheres[k].shader_type;
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
                            hit_sphere_color = spheres[k].color;
                            hit_sphere_shading = spheres[k].shader_type;
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

tri_mesh load_mesh(string off_filepath, Vector3d position, double scale)
{
    // returns a triangule mesh object at the file path for the OFF file
    tri_mesh mesh_structure;
    mesh_structure.color = Vector3d(255,255,255);
    mesh_structure.shader_type = 'd';

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
            mesh_structure.V.row(i) << stod(substrings[0])*scale+position.coeff(0),
            stod(substrings[1])*scale+position.coeff(1),
            stod(substrings[2])*scale+position.coeff(2);
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
    Vector3d cam_origin(0,0,2);
    Vector3d cam_d = Vector3d(0,0,-1.5);
    Vector3d x_displacement;
    Vector3d y_displacement;

    // Multiple Light Sources
    MatrixXd light_positions(2,3);
    light_positions <<  5.0,  5.0,  5.0,
        -6.0, -6.0, 5.0;

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

bool surface_light(ray r, vector<sphere> spheres, vector<tri_mesh> tri_meshes)
{
    double t1 = 0;
    double epsilon = 0.001;
    for (unsigned k=0;k<spheres.size();k++)
    {
        const double sphere_radius = spheres[k].radius;
        Vector3d sphere_center = spheres[k].position;
        //TODO: use ray_hit_sphere() instead
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

    for (unsigned k=0;k<tri_meshes.size();k++)
    {
        for (unsigned a=0;a<tri_meshes[k].F.rows();a++) 
        {
            triangle test_tri;
            test_tri.v1 = Vector3d(tri_meshes[k].V.row(tri_meshes[k].F.coeff(a,0)));
            test_tri.v2 = Vector3d(tri_meshes[k].V.row(tri_meshes[k].F.coeff(a,1)));
            test_tri.v3 = Vector3d(tri_meshes[k].V.row(tri_meshes[k].F.coeff(a,2)));
            hit_triangle hit_test2 = ray_hit_triangle(r,epsilon,test_tri);
            if(hit_test2.hit)
            {
                return true;
            }
        }
    }
    return false;
}

void task1_5(vector<sphere> spheres, vector<tri_mesh> tri_meshes)
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

    // For each ray with direction -w
    for (unsigned i=0;i<R.cols();i++)
    {
        for (unsigned j=0;j<R.rows();j++)
        {
            // Prepare the ray, perspective camera
            x_displacement = Vector3d(l+(r-l)*(i)/R.cols(),0,0);
            y_displacement = Vector3d(0,t-(t-b)*(j)/R.rows(),0);
            ray camera_ray;
            camera_ray.e = cam_origin;
            camera_ray.d = cam_d+x_displacement+y_displacement;

            hit_sphere hit_record;
            hit_triangle hit_record2;

            // For each sphere
            for (unsigned k=0;k<spheres.size();k++)
            {
                sphere test_sphere;
                test_sphere.radius = spheres[k].radius;
                test_sphere.position = spheres[k].position;
                test_sphere.color = spheres[k].color;
                test_sphere.shader_type = spheres[k].shader_type;
                hit_sphere hit_test = ray_hit_sphere(camera_ray,0,test_sphere);
                if(hit_test.hit)
                {
                    if(hit_record.hit==false)
                    {
                        hit_record = hit_test;
                    }
                    else
                    {
                        if(hit_test.t<hit_record.t)
                        {
                            hit_record = hit_test;
                        }
                    }

                }
            }

            //For each mesh
            for (unsigned k=0;k<tri_meshes.size();k++)
            {
                for (unsigned a=0;a<tri_meshes[k].F.rows();a++) 
                {
                    
                    triangle test_tri;
                    test_tri.v1 = Vector3d(tri_meshes[k].V.row(tri_meshes[k].F.coeff(a,0)));
                    test_tri.v2 = Vector3d(tri_meshes[k].V.row(tri_meshes[k].F.coeff(a,1)));
                    test_tri.v3 = Vector3d(tri_meshes[k].V.row(tri_meshes[k].F.coeff(a,2)));
                    test_tri.color = tri_meshes[k].color;
                    test_tri.shader_type = tri_meshes[k].shader_type;
                    hit_triangle hit_test2 = ray_hit_triangle(camera_ray,0,test_tri);
                    if(hit_test2.hit)
                    {
                        if(hit_record2.hit == false)
                        {
                            hit_record2 = hit_test2;
                        }
                        else
                        {
                            if(hit_test2.t<hit_record2.t)
                            {
                                hit_record2 = hit_test2;
                            }
                        }
                    }
                }
            }

            //if hit is true, then do lighting calculations
            if (hit_record.hit || hit_record2.hit)
            {
                string hit_type;
                double t_int;
                if(hit_record.hit && hit_record2.hit)
                {
                    if(hit_record.t<=hit_record2.t)
                    {
                        hit_type = "sphere";
                        t_int = hit_record.t;
                    }
                    else
                    {
                        hit_type = "triangle";
                        t_int = hit_record2.t;
                    }
                }
                else
                {
                    if(hit_record.hit)
                    {
                        hit_type = "sphere";
                        t_int = hit_record.t;
                    }
                    if(hit_record2.hit)
                    {
                        hit_type = "triangle";
                        t_int = hit_record2.t;
                    }
                }
                
                
                Vector3d ray_intersection = camera_ray.e + t_int * camera_ray.d;
                
                Vector3d ray_normal;
                if(hit_type.compare("sphere")==0)
                {
                    ray_normal = (ray_intersection-hit_record.hit_obj.position)/hit_record.hit_obj.radius;
                }
                if(hit_type.compare("triangle")==0)
                {
                    Vector3d h1 = hit_record2.hit_obj.v1;
                    Vector3d h2 = hit_record2.hit_obj.v2;
                    Vector3d h3 = hit_record2.hit_obj.v3;
                    ray_normal = ((h2-h1).cross(h3-h1)).normalized();
                }

                double color_intensity=0.0;

                //for each light source
                for(unsigned m=0;m<light_positions.rows();m++)
                {
                    Vector3d light_vector = (Vector3d(light_positions.row(m))-ray_intersection);
                    ray shadow_ray;
                    shadow_ray.e = ray_intersection;
                    shadow_ray.d = light_vector;

                    if(surface_light(shadow_ray, spheres, tri_meshes)==false)
                    {
                        light_vector = light_vector.normalized();
                        // Simple diffuse model assuming light intensity I=1
                        // Clamp to zero
                        color_intensity = light_vector.transpose() * ray_normal;
                        color_intensity = max(color_intensity,0.);
                        if(hit_type.compare("sphere")==0)
                        {
                            //Specular shading, assuming Phong exponent p=100
                            if(hit_record.hit_obj.shader_type=='s')
                            {
                                Vector3d half_vector=(-camera_ray.d+light_vector).normalized();
                                color_intensity = color_intensity + pow(max(ray_normal.dot(half_vector),0.),100);
                            }
                        }
                        if(hit_type.compare("triangle")==0)
                        {
                            //Specular shading, assuming Phong exponent p=100
                            if(hit_record2.hit_obj.shader_type=='s')
                            {
                                Vector3d half_vector=(-camera_ray.d+light_vector).normalized();
                                color_intensity = color_intensity + pow(max(ray_normal.dot(half_vector),0.),100);
                            }
                        }
                        R(i,j) = R(i,j)+color_intensity;
                        G(i,j) = G(i,j)+color_intensity;
                        B(i,j) = B(i,j)+color_intensity;
                    }
                }

                //adjust RGB based on the hit sphere's color
                //assume kd and ks are same, based on RGB values provided in arguments
                if(hit_type.compare("sphere")==0)
                {
                    R(i,j)=R(i,j)*hit_record.hit_obj.color.coeff(0)/255;
                    G(i,j)=G(i,j)*hit_record.hit_obj.color.coeff(1)/255;
                    B(i,j)=B(i,j)*hit_record.hit_obj.color.coeff(2)/255;
                }
                if(hit_type.compare("triangle")==0)
                {
                    R(i,j)=R(i,j)*hit_record2.hit_obj.color.coeff(0)/255;
                    G(i,j)=G(i,j)*hit_record2.hit_obj.color.coeff(1)/255;
                    B(i,j)=B(i,j)*hit_record2.hit_obj.color.coeff(2)/255;
                }

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }

        }
    }
    write_matrix_to_png(R,G,B,A,filename);
}

Vector4d raycolor(ray r, double t_lower, vector<sphere> spheres, vector<tri_mesh> tri_meshes, MatrixXd light_positions, int recur_count)
{
    Vector4d rgba (0,0,0,0);
    if(recur_count<2)
    {
        double t1 = 0;
        double epsilon = 0.001;
        hit_sphere hit_record;
        hit_triangle hit_record2;

        // For each sphere
        for (unsigned k=0;k<spheres.size();k++)
        {
            sphere test_sphere;
            test_sphere.radius = spheres[k].radius;
            test_sphere.position = spheres[k].position;
            test_sphere.color = spheres[k].color;
            test_sphere.shader_type = spheres[k].shader_type;
            hit_sphere hit_test = ray_hit_sphere(r,t_lower,test_sphere);
            if(hit_test.hit)
            {
                if(hit_record.hit==false)
                {
                    hit_record = hit_test;
                }
                else
                {
                    if(hit_test.t<hit_record.t)
                    {
                        hit_record = hit_test;
                    }
                }

            }
        }

        //For each mesh
        for (unsigned k=0;k<tri_meshes.size();k++)
        {
            for (unsigned a=0;a<tri_meshes[k].F.rows();a++) 
            {
                triangle test_tri;
                test_tri.v1 = Vector3d(tri_meshes[k].V.row(tri_meshes[k].F.coeff(a,0)));
                test_tri.v2 = Vector3d(tri_meshes[k].V.row(tri_meshes[k].F.coeff(a,1)));
                test_tri.v3 = Vector3d(tri_meshes[k].V.row(tri_meshes[k].F.coeff(a,2)));
                test_tri.color = tri_meshes[k].color;
                test_tri.shader_type = tri_meshes[k].shader_type;
                test_tri.mirror = tri_meshes[k].mirror;
                hit_triangle hit_test2 = ray_hit_triangle(r,t_lower,test_tri);
                if(hit_test2.hit)
                {
                    if(hit_record2.hit == false)
                    {
                        hit_record2 = hit_test2;
                    }
                    else
                    {
                        if(hit_test2.t<hit_record2.t)
                        {
                            hit_record2 = hit_test2;
                        }
                    }
                }
            }
        }

        //if hit is true, then do lighting calculations
        if (hit_record.hit || hit_record2.hit)
        {
            string hit_type;
            double t_int;
            if(hit_record.hit && hit_record2.hit)
            {
                if(hit_record.t<=hit_record2.t)
                {
                    hit_type = "sphere";
                    t_int = hit_record.t;
                }
                else
                {
                    hit_type = "triangle";
                    t_int = hit_record2.t;
                }
            }
            else
            {
                if(hit_record.hit)
                {
                    hit_type = "sphere";
                    t_int = hit_record.t;
                }
                if(hit_record2.hit)
                {
                    hit_type = "triangle";
                    t_int = hit_record2.t;
                }
            }
            
            Vector3d ray_intersection = r.e + t_int * r.d;

            Vector3d ray_normal;
            ray reflect_ray;

            if(hit_type.compare("sphere")==0)
            {
                ray_normal = (ray_intersection-hit_record.hit_obj.position)/hit_record.hit_obj.radius;
            }
            if(hit_type.compare("triangle")==0)
            {
                Vector3d h1 = hit_record2.hit_obj.v1;
                Vector3d h2 = hit_record2.hit_obj.v2;
                Vector3d h3 = hit_record2.hit_obj.v3;
                ray_normal = ((h2-h1).cross(h3-h1)).normalized();

                //if hit type is a mirror triangle
                if(hit_record2.hit_obj.mirror)
                {
                    //construct a new ray with origin as ray intersection and direction based on ray normal
                    reflect_ray.e = ray_intersection;
                    reflect_ray.d = r.d - 2*(r.d).dot(ray_normal)*ray_normal;
                }
            }

            double color_intensity=0.0;

            //for each light source
            for(unsigned m=0;m<light_positions.rows();m++)
            {
                Vector3d light_vector = (Vector3d(light_positions.row(m))-ray_intersection);
                ray shadow_ray;
                shadow_ray.e = ray_intersection;
                shadow_ray.d = light_vector;

                if(surface_light(shadow_ray, spheres, tri_meshes)==false)
                {
                    //if there is a detected reflection, 
                    //and the point is lit by at least one light source
                    //then no further lighting calcs are needed because ideally that color is the reflected point's color
                    if(hit_type == "triangle" && hit_record2.hit_obj.mirror)
                    {
                        rgba = raycolor(reflect_ray,epsilon,spheres,tri_meshes,light_positions,recur_count+1);
                    }
                    else
                    {
                        light_vector = light_vector.normalized();
                        // Simple diffuse model assuming light intensity I=1
                        // Clamp to zero
                        color_intensity = light_vector.transpose() * ray_normal;
                        color_intensity = max(color_intensity,0.);
                        if(hit_type.compare("triangle")==0)
                        {
                            //Specular shading, assuming Phong exponent p=100
                            if(hit_record2.hit_obj.shader_type=='s')
                            {
                                Vector3d half_vector=(-r.d+light_vector).normalized();
                                color_intensity = color_intensity + pow(max(ray_normal.dot(half_vector),0.),100);
                            }
                            rgba = rgba + color_intensity * Vector4d::Ones();
                            rgba(0) = rgba.coeff(0)*hit_record2.hit_obj.color.coeff(0)/255;
                            rgba(1) = rgba.coeff(1)*hit_record2.hit_obj.color.coeff(1)/255;
                            rgba(2) = rgba.coeff(2)*hit_record2.hit_obj.color.coeff(2)/255;
                        }
                        else
                        {
                            if(hit_type.compare("sphere")==0)
                            {
                                //Specular shading, assuming Phong exponent p=100
                                if(hit_record.hit_obj.shader_type=='s')
                                {
                                    Vector3d half_vector=(-r.d+light_vector).normalized();
                                    color_intensity = color_intensity + pow(max(ray_normal.dot(half_vector),0.),100);
                                }
                            }
                            rgba = rgba + color_intensity * Vector4d::Ones();
                            rgba(0) = rgba.coeff(0)*hit_record.hit_obj.color.coeff(0)/255;
                            rgba(1) = rgba.coeff(1)*hit_record.hit_obj.color.coeff(1)/255;
                            rgba(2) = rgba.coeff(2)*hit_record.hit_obj.color.coeff(2)/255;
                        }
                        //TODO: think about where this should go
                        rgba(3) = 1;
                    }
                }
                else
                {
                    rgba(3) = 1;
                }
            }
        }
    }
    return rgba;
}

void task1_6(vector<sphere> spheres, vector<tri_mesh> tri_meshes)
{
    std::cout << "Task 1.6: Reflection" << std::endl;
    //camera coordinates
    double l = -1.0;
    double r = 1.0;
    double t = 1.0;
    double b = -1.0;

    const std::string filename("task1_6.png");
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

    // For each ray 
    for (unsigned i=0;i<R.cols();i++)
    {
        for (unsigned j=0;j<R.rows();j++)
        {
            // Prepare the ray, perspective camera
            x_displacement = Vector3d(l+(r-l)*(i)/R.cols(),0,0);
            y_displacement = Vector3d(0,t-(t-b)*(j)/R.rows(),0);
            ray camera_ray;
            camera_ray.e = cam_origin;
            camera_ray.d = cam_d+x_displacement+y_displacement;

            Vector4d rgba = raycolor(camera_ray, 0, spheres, tri_meshes, light_positions,0);

            R(i,j)=rgba.coeff(0);
            G(i,j)=rgba.coeff(1);
            B(i,j)=rgba.coeff(2);
            A(i,j)=rgba.coeff(3);

        }
    }
    write_matrix_to_png(R,G,B,A,filename);
}

int main()
{
    //part1();
    //part2();

    sphere sphere1;
    sphere sphere2;
    sphere sphere3;

    sphere1.radius = 0.5;
    sphere1.position = Vector3d(0,0.5,0);
    sphere1.color = Vector3d(66, 135, 245);
    sphere1.shader_type = 'd';

    sphere2.radius = 0.1;
    sphere2.position = Vector3d(-0.6, 0.6, 0.6);
    sphere2.color = Vector3d(0, 255, 166);
    sphere2.shader_type = 's';

    sphere3.radius = 0.3;
    sphere3.position = Vector3d(0.3, 0.3, 0.3);
    sphere3.color = Vector3d(255, 0, 255);
    sphere3.shader_type = 's';

    vector<sphere> spheres;
    spheres.push_back(sphere1);
    spheres.push_back(sphere2);
    spheres.push_back(sphere3);

    task1_1(spheres);
    task1_2(spheres);
    task1_3(spheres);

    tri_mesh ground_plane = load_mesh("../data/ground_plane.off",Vector3d(0,0,0),1.0);;

    //tri_mesh test_mesh = load_mesh("../data/bumpy_cube.off",Vector3d(0,0,0),1/4.5);
    tri_mesh test_mesh = load_mesh("../data/bunny.off",Vector3d(-0.3,0,0.5),3);
    test_mesh.color = Vector3d(245, 226, 20);
    //add mesh to mesh array

    //task1_4(test_mesh);
    
    vector<tri_mesh> tri_meshes;
    tri_meshes.push_back(ground_plane);
    task1_5(spheres, tri_meshes);

    ground_plane.mirror=true;
    tri_meshes.pop_back();
    tri_meshes.push_back(ground_plane);

    MatrixXd light_positions(2,3);
    light_positions <<  -2.0,  2.0,  2.0,
                         1.0, -1.0, -1.0;

    task1_6(spheres,tri_meshes);
    return 0;
}
