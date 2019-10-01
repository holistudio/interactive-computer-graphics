// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <math.h>

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"

// Shortcut to avoid Eigen:: and std:: everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;

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
    const Vector3d light_position(-0.8,1,1.2);

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

int main()
{
    //part1();
    //part2();
    vector<float> spheresRadii;
    MatrixXd spheresCenters(3,3);
    spheresRadii.push_back(0.9);
    spheresRadii.push_back(0.3);
    spheresRadii.push_back(0.4);
    spheresCenters << 0, 0, 0,
                      0.2, -0.2, 1.0,
                      0.5,0.5,0.5;
    //task1_1(spheresRadii,spheresCenters);
    MatrixXd spheres_colors(3,3); //in RGB
    spheres_colors << 66, 135, 245,
                     0, 255, 166,
                     237, 0, 206;
    vector<char> sphere_shading;
    sphere_shading.push_back('d');
    sphere_shading.push_back('s');
    sphere_shading.push_back('s');
    //task1_2(spheresRadii,spheresCenters,spheres_colors,sphere_shading);
    task1_3(spheresRadii,spheresCenters,spheres_colors,sphere_shading);
    return 0;
}
