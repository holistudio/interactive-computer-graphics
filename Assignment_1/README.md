# Basic Ray Tracer

A basic C++ raytracing program using the Eigen library (basically, if you didn't have OpenGL in your life). Spheres are added in the scene manually in the code. Triangular meshes are loaded as .OFF files from the data folder, then scaled and moved based on input parameters (see load_mesh function in the source code).

Each task involves an additional raytracing features, such as different lighting and shading techniques. For the final task, the ground plane is converted into a mirror.

All images below are from the same camera position and focal length (for perspective projection only).

## Shading Spheres - Orthographic Projection

<img src="img/task1_1.png" width="500">

## Ambient, Diffuse, and Specular Lighting - Orthographic Projection

<img src="img/task1_2.png" width="500">

## Perspective Projection

<img src="img/task1_3.png" width="500">

## Triangular Meshes

<img src="img/task1_4.png" width="500">

## Shadows

<img src="img/task1_5.png" width="500">

## Ground Plane Mirror

<img src="img/task1_6.png" width="500">

# Challenges

## Epsilon Criteria for Shadows

<img src="img/task1_5_epsilon_01.png" width="500">

Above, with epsilon variable set to 0.01. At first I thought I needed to set a smaller epsilon for shadows.

<img src="img/task1_5_epsilon_001.png" width="500">

Above, with epsilon variable set to 0.001. Then I realized I didn't apply epsilon criteria to triangle meshes.

<img src="img/task1_5_epsilon-criteria-applied-to-meshes.png" width="500">

## Eigen C++ Library

Eigen doesn't seem to have a clear way of casting a matrix of a row directly into a VectorXd.

Eigen's cross product function `x.cross(y)` doesn't work unless you `#include <Eigen/Geometry>`. The error message for not doing just mentions that you are using the wrong data type, which doesn't exactly help beginners. After hours of Googling, Eigen's common pitfalls webpage came up (https://eigen.tuxfamily.org/dox/TopicPitfalls.html#title4).

## Basic compiling instructions

If this is the first time you cloned or downloaded this repository, you'll need to set up with CMake (so install CMake before continuing).

Then, with Assignment_1 as the current directory in terminal:

```
mkdir build
cd build
cmake ../
```

To compile the code:
```
make
./Assignment_1
```

Output image files should be in the build folder.

If only changes are made to the files, only the last two lines of terminal commands are needed. If you ever add files, of course you'll need to `cmake ../` again.

# References

Fundamentals of Computer Graphics
