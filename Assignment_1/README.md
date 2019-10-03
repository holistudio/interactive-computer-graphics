# Issues

Eigen's cross product function `x.cross(y)` doesn't work unless you `#include <Eigen/Geometry>`. The error message for not doing just mentions that you are using the wrong data type, which doesn't exactly help beginners. After hours of Googling, Eigen's common pitfalls webpage came up (https://eigen.tuxfamily.org/dox/TopicPitfalls.html#title4).