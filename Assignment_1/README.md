# Shading Spheres - Orthographic Projection
<img src="img/task1_1.png" width="500">

# Ambient, Diffuse, and Specular Lighting - Orthographic Projection
<img src="img/task1_2.png" width="500">

# Perspective Projection
<img src="img/task1_3.png" width="500">

# Triangular Meshes
<img src="img/task1_4.png" width="500">

# Shadows
<img src="img/task1_5.png" width="500">

# Ground Plane Mirror
<img src="img/task1_6.png" width="500">

# Challenges

Eigen doesn't seem to have a clear way of casting a matrix of a row directly into a VectorXd.

Eigen's cross product function `x.cross(y)` doesn't work unless you `#include <Eigen/Geometry>`. The error message for not doing just mentions that you are using the wrong data type, which doesn't exactly help beginners. After hours of Googling, Eigen's common pitfalls webpage came up (https://eigen.tuxfamily.org/dox/TopicPitfalls.html#title4).

Epsilon Criteria for Shadows
<img src="img/task1_5_epsilon_01.png" width="500">

Above, with epsilon variable set to 0.01. At first I thought I needed to set a smaller epsilon for shadows.

<img src="img/task1_5_epsilon_001.png" width="500">

Above, with epsilon variable set to 0.001. Then I realized I didn't apply epsilon criteria to triangle meshes.

<img src="img/task1_5_epsilon-criteria-applied-to-meshes.png" width="500">
