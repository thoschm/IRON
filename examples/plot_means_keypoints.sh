#!/usr/bin/gnuplot

splot "means1.txt", "means_keypoints1.txt",  "means2.txt",  "means_keypoints2.txt", "inlierset.txt" w l
pause -1
