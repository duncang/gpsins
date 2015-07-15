/** 
 * \file inslib.c Implements INS helper functions
 *   $Id: inslib.c 355 2007-10-15 21:59:57Z greerd $
 */


#include "inslib.h"
#include <math.h>
#include "quat.h"




/**
 * Perform coarse levelling of the Strapdown IMU.  This basically involves calculating 
 * the attitude of the platform from the accelerometer outputs assuming that the 
 * system is stationary.
 *
 *
 */
void CoarseLevelling()
{

	



}


/**
 * Estimate the attitude quaternion based on the Accelerometer and Magnetometer 
 * measurements.
 *
 * \param dQuatResult - Storage for the 4-element quaternion vector result
 * \param dAccel - The 3-axis accelerometer measurements: Abx, Aby, Abz
 * \param dMag - The 3-axis magnetometer measurements: Mbx, Mby, Mbz
 * \return -1 if attitude is indeterminate
 */
int QuaternionFromAccel(double dQuatResult[4], double dAccel[3], double dMag[3])
{
	double cos_theta, cos_theta_2, sin_theta_2, cos_phi, cos_phi_2, sin_phi_2;	
	double dMag2[4];
	double dq[4], dqinv[4];
	double dq1[4];
	double dq2[4];
	double dq3[4];

	//double T_nb[3][3], T_bn[3][3];
	//double dMag_n[3];

	cos_theta = 1 / (sqrt(1 + pow((dAccel[0] / - dAccel[2]),2)));
	cos_theta_2 = sqrt(0.5 * (1 + cos_theta));
	sin_theta_2 = sqrt(0.5 * (1 - cos_theta)) * (double)sign(dAccel[0]);

	cos_phi = 1 / (sqrt(1 + pow((-dAccel[1] / sqrt(pow(dAccel[0],2) + pow(dAccel[2],2))),2)));
	cos_phi_2 = sqrt(0.5 * (1 + cos_phi));
	sin_phi_2 = sqrt(0.5 * (1 - cos_phi)) * (double)sign(-dAccel[1]);	

	dq[0] = cos_theta_2 * cos_phi_2;
	dq[1] = cos_theta_2 * sin_phi_2;
	dq[2] = sin_theta_2 * cos_phi_2;
	dq[3] = -sin_theta_2 * sin_phi_2;

	/* renormalise */
	normq(dq);

	/* find the heading from the magnetometer */

	dMag2[0] = 0;
	dMag2[1] = dMag[0];
	dMag2[2] = dMag[1];
	dMag2[3] = dMag[2];
	normq(dMag2);

	quatinv(dqinv,dq);
	
	quatmul(dq1,dqinv,dMag2);
	quatmul(dq2,dq1,dq);
	
	normq(dq2);

	dq3[0] = -atan2(dq2[2],dq2[1]);
	dq3[1] = 0;
	dq3[2] = 0;
	dq3[3] = 1;
	
	quatmul(dQuatResult,dq,dq3);


	/* renormalise */
	normq(dQuatResult);


	return 0;
}

/** find the sign (positive or negative) of a double */
int sign(double dValue)
{
	if(dValue < 0.0)
	{
		return -1;
	} else
	{
		return 1;
	}
}

int GenerateCEN(double dCEN[3][3],double dLongitude, double dLatitude)
{
	double cLong,cLat,sLong,sLat;

	cLong = cos(dLongitude);
	cLat = cos(dLatitude);
	sLong = sin(dLongitude);
	sLat = sin(dLatitude);

	dCEN[0][0] = -cLong*sLat;
	dCEN[0][1] = -sLong*sLat;
	dCEN[0][2] = cLat;
	dCEN[1][0] = -sLong;
	dCEN[1][1] = cLong;
	dCEN[1][2] = 0.0;
	dCEN[2][0] = -cLong*cLat;
	dCEN[2][1] = -sLong*cLat;
	dCEN[2][2] = -sLat;
	
	return 0;

}

