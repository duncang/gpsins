/** 
 * \file gps.h Implements various algorithms used in GPS
 * \author Duncan Greer
 * $Id: gps.h 2834 2009-08-25 07:12:42Z greerd $
 */

#ifndef _GPS_H
#define _GPS_H

#include "novatel.h"


#define MAX_SV 35
#define GPS_SPEEDOFLIGHT 2.99792458e8
#define GPS_OMEGAEDOT 7.2921151467e-5
#define GPS_PI 3.1415926535898 // ICD value of PI
#define GPS_MU 3.986005e14 // WGS-84 valeu of the earths universal gravitational parameter in m^3/s^2
#define GPS_EARTHRADIUS 6378136  // m
#define GPS_F -4.442807633e-10 //a random number from the ICD page 88
#define GPS_L1_f 1575.42e6  // Hz
#define GPS_L2_f 1227.6e6  // Hz
#define GPS_GAMMA (L1_f/L2_f)*(L1_f/L2_f) //  unitless
#define GPS_L1_Wavelength GPS_SPEEDOFLIGHT/GPS_L1_f  //  Metres

#define RAD2DEG 180.0/GPS_PI
#define DEG2RAD GPS_PI/180.0

// WGS-84 ellipsoid parameters
#define WGS84_a 6378137.0 	// semi-major axi
#define WGS84_b 6356752.3142    // semi-minor axis
#define WGS84_Esq (1.0 - (WGS84_b*WGS84_b) / (WGS84_a*WGS84_a))
#define WGS84_EPsq ((WGS84_a*WGS84_a)/(WGS84_b*WGS84_b) - 1.0)




int GPSOrbitPropagator(double dTransmissionTime, 
			unsigned short usPRN, 
			const GPSEPHEM *GPSEphemData, 
			double dEphemValidityTime, 
			double *SVPos);
						
int GPSOrbitPropagatorVelocities(double dTransmissionTime, 
			unsigned short usPRN, 
			const GPSEPHEM *GPSEphemData, 
			double dEphemValidityTime, 
			double *SVVel,
			double *SVAcc);

double KepplerSolver(double Mk, double e);

int WGS84_ECEF2LLH(const double *ecef, double *llh);

int ICD200_IonoModel(GPSTime gtTimeStamp, double *Xu, double *Su, double *alpha, double *beta,
	double *IonoDelay, double *Azimuth, double *Elevation);



void T_ECEF2ENU(double Latitude, double Longitude, double **Result);

/**
 *	\fn TropoDelay
 *   Calculates the ionospheric delay for a satellite pseudorange from a simplified tropo model.
 *	\author Duncan Greer
 *
 */
double TropoDelayModel(double Elevation,double UserHeight);



void T_Body2NED(double PHI, double THETA, double PSI, double **Result);

void T_ECEF2NED(double Latitude, double Longitude, double **Result);
void T_NED2ECEF(double Latitude, double Longitude, double **Result);

int WGS84_LLH2ECEF(const double *llh,  double *ecef);

void WGS84_calcRnRe(double lat,  double *Rn, double *Re);



#endif
