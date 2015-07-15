GPS-INS program - Written by Duncan Greer (c) 2006
==================================================
Revision Information: $Id: Readme.txt 340 2007-10-08 11:47:09Z greerd $

This program implements a GPS-INS estimator using GSL and RTAI routines in C.  GPS and IMU data is collected from the serial ports and stored in shared memory.  A separate program then reads the measurements from the shared memory and calculates the position solution.


Files
=====

gpsins.c - Contains the main function and control thread.


SVN
===

The subversion repository for this project is http://192.168.0.35:81/svn-duncan/gpsins/trunk/

