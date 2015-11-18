# IRON
IRON: A Fast Interest Point Descriptor for Robust NDT-Map Matching and Its Application to Robot Localization

IRON Reference Implementation:
*************************************************************************************************************
Schmiedel, Th., Einhorn, E., Gross, H.-M.
IRON: A Fast Interest Point Descriptor for Robust NDT-Map Matching and Its Application to Robot Localization.
IEEE/RSJ Int. Conf. on Intelligent Robots and Systems (IROS), Hamburg, Germany, 2015
*************************************************************************************************************
Copyright (C) 2015, Thomas Schmiedel (thomas.schmiedel@tu-ilmenau.de)
 
 
Required packages (Linux Mint 17.1; please install equivalent packages for your Linux distribution):
- build-essential
- libeigen3-dev 
- gnuplot-x11 


How to run demo?
- cd IRON
- make
- cd examples
- ./demo example1_cloud1.pcd example1_cloud2.pcd
- ./plot_means_keypoints.sh
- ./plot_result.sh


How to use IRON registration pipeline?
- just include IRON.h into your project
- have a look at src/demo.cpp for a detailed explanation


Please contact me if you have any questions!

