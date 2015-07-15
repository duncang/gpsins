/** 
 * \file inslib.h Implements INS helper functions
 *   $Id: inslib.h 355 2007-10-15 21:59:57Z greerd $
 */


#ifndef _INSLIB_H
#define _INSLIB_H

/**
 * Estimate the attitude quaternion based on the Accelerometer and Magnetometer 
 * measurements.
 *
 * \param dQuatResult - Storage for the 4-element quaternion vector result
 * \param dAccel - The 3-axis accelerometer measurements: Abx, Aby, Abz
 * \param dMag - The 3-axis magnetometer measurements: Mbx, Mby, Mbz
 * \return -1 if attitude is indeterminate
 */
int QuaternionFromAccel(double dQuatResult[4], double dAccel[3], double dMag[3]);

/** find the sign (positive or negative) of a double */
int sign(double dValue);

/** get the Earth Frame to Navigation Frame (NED) transformation matrix
 * 
 * \param dCEN[3][3]	The destination matrix
 * \param dLongtitude 	The longitude in radians
 * \param dLatitude	The geodetic latitude in radians
 * \return Always 0
 */

int GenerateCEN(double dCEN[3][3],double dLongitude, double dLatitude);

#endif

