# IRON
IRON: A Fast Interest Point Descriptor for Robust NDT-Map Matching - Reference Implementation

![alt tag](https://raw.github.com/thoschm/IRON/master/img/img1.png)

**Features:**
- IRON provides a simple interface for the alignment of 3D point clouds using 6 degrees of freedom (x, y, z, roll, pitch, yaw)
- Almost any given point cloud format can be directly used without expensive point cloud type conversion beforehand (see src/demo.cpp for an example)
- It's robust and insensitive to a large initial displacement between two point clouds (e.g. when they were captured from different locations)
- IRON uses NDT-maps as an efficient 3D data structure and is able to create them with sufficient accuracy in about 3 ms (average conversion time for a typical Kinect 1 / Asus Xtion point cloud)
- The complete alignment process takes only 10-20 ms on a single core of a modern i7 CPU (depending on registration parameters and map sizes) and is therefore well suited for real-time applications

**For algorithmic details and benchmarking results please refer to the original paper:**
*************************************************************************************************************
Schmiedel, Th., Einhorn, E., Gross, H.-M.
IRON: A Fast Interest Point Descriptor for Robust NDT-Map Matching and Its Application to Robot Localization.
IEEE/RSJ Int. Conf. on Intelligent Robots and Systems (IROS), Hamburg, Germany, 2015
*************************************************************************************************************

You may also visit my [research page](http://research.thomas-schmiedel.com/?page_id=35) for a thorough introduction.
 
 
**Required packages** (Ubuntu 14.04; please install equivalent packages for your Linux distribution):
- build-essential
- libeigen3-dev 
- gnuplot-x11 


**How to run demo?**
- cd IRON
- make
- cd examples
- ./demo example1_cloud1.pcd example1_cloud2.pcd
- ./plot_means_keypoints.sh
- ./plot_result.sh


**How to use the IRON registration pipeline?**
- just include IRON.h into your project
- have a look at src/demo.cpp for a detailed explanation


**Please contact me if you have any questions!**

thomas.schmiedel@tu-ilmenau.de
