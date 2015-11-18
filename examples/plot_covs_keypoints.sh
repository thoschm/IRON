#!/usr/bin/gnuplot

splot "covs1.txt", "covs_keypoints1.txt", "covs2.txt", "covs_keypoints2.txt", "inlierset.txt" w l
pause -1
