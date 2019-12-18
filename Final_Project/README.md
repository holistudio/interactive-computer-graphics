# Rendering Meshes with OpenGL

Using C++ and OpenGL, this program allows for a user to insert and animate mesh model based on 3D human pose data. This project builds off of earlier human pose estimation work [here](https://github.com/holistudio/nomocap). The full "pipeline" of data and visualization can be illustrated as follows:

<img src="img/pipeline-01.png" width="500">

Full results visualizations can be found [here](https://vimeo.com/380161358)

## Demo / Getting Started

### CMake

If this is the first time you cloned or downloaded this repository, you'll need to set up with CMake (so install CMake before continuing).

Then, with Final_Project as the current directory in terminal:

```
mkdir build
cd build
cmake ../
```

To compile the code:
```
make
./Final_Project_bin
```

If only changes are made to the files, only the last two lines of terminal commands are needed. If you ever add files, of course you'll need to `cmake ../` again.

### OK it's running...

Once the build is running, you should see a window titled "Wushu!"

Press the '1' key to load the pose data. Specifically `data/vertices.csv` contains the vertex coordinates for all poses detected in a video sequence that follows the [Human3.6M](http://vision.imar.ro/human3.6m/description.php) dataset conventions. To visualize other human poses, you'll need to replace the csv file with your own following the same conventions for pose joints.

Once you press the '1' key, the console should display:
```
Pose coordinates loaded
All 65 keyframe poses loaded
```
The human pose mesh should also appear in the window.

Press the '/' key to begin the animation.

### Camera Control

As the animation is running (or even before you start the animation), you can adjust the camera as follows:

 - 'o' / 'p' keys switch between orthographic and perspective
 - 'a' / 'd' moves the camera on the x-axis directions
 - 'q' / 'z' moves the camera on the y-axis directions
 - 'w' / 's' moves the camera on the z-axis directions

Currently the camera gaze direction always points towards the world origin.

## How the Code Works

The current version of the code generates meshes around the 3D pose skeletons in `data/vertices.csv` using the unit cube mesh (`data/cube.off`).

<img src="img/mesh_make-01.png" width="500">

Meshes are then visualized in OpenGL. Future development of the code will look at fancier means of generating mesh around the 3D pose skeletons like MD5 and FBX files.
