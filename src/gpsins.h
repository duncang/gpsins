/**
 * \file gpsins.h
 *
 * \brief Contains headers for GPS-Inertial system
 */

#ifndef _gpsins_h
#define _gpsins_h

#include "../../novatelstore/src/novatel.h"

/** \brief IMUData structure to store IMU measurements
*/
typedef struct _IMUData
{
	double ddt; /* time step */
	double dAcc[3];  /* acceleration measurements - body axes */
	double dRate[3]; /* rate measurements - body axes */
	double dMag[3];  /* magnetometer measurements */
} IMUData;



/**
 *	_message: A generic message type to transmit data between threads and processes
 */
typedef struct _message
{
	int iMessageType;
	char azucMessageData[100];
} message;


/* KF stuff  */
#define NUMBER_STATES		17

/* shared memory keys */
#define GPSDATA_KEY	1234
#define IMUDATA_KEY	GPSDATA_KEY + 1


#define WGS84_RE 6378136.0
#define LOCAL_GRAVITY 9.81

#define IMU_RATE 100.0
#define INS_DT 1.0/IMU_RATE
#define GPS_RATE 1.0
#define GPS_DT 1.0/GPS_RATE

#define OMEGA_e 7.292115e-5



#endif /* _gpsins_h  */

