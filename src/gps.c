/** 
 * \file gps.c Implements various algorithms used in GPS
 * \author Duncan Greer
 * $Id: gps.c 452 2007-11-23 22:19:47Z n2523710 $
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "gps.h"
#include "matrix.h"

int GPSOrbitPropagator(double dTransmissionTime, 
			unsigned short usPRN, 
			const GPSEPHEM *GPSEphemData, 
			double dEphemValidityTime, 
			double *SVPos)
{

	double ts,delta_ts,toc,TGD,A,n0,n;
	double toe,af0,af1,af2,vk,delta_tr;
	double tk,e,Ek,Mk;
	double uk,rk,ik,xk,yk,zk;
	int ValidData;
	double OMEGA0, OMEGAdot, OMEGAk;
	double Cis,Cus,Cic,Cuc,Crs,Crc;
	double IDOT,i0,xk_dash,yk_dash;
	double delta_ik,delta_uk,delta_rk;
	double PHIk;	
	double SVPosVector[4];
	double calcGPSTime;
	
        ts = dTransmissionTime;

        // make initial estimate of Ek - this is revised later with the
        // corrected tk

        toe = GPSEphemData->dTOE;
        A = GPSEphemData->dA;
        n0 = sqrt(GPS_MU/(A*A*A));
        tk = ts - toe;

        // check for end of week roll-over
        if (tk > 302400)
		{
            tk = tk - 604800;
		} else if (tk < -302400)
		{
            tk = tk + 604800;
        }

        

        n = n0 + GPSEphemData->dDeltaN;
        Mk = GPSEphemData->dM0 + n*tk;
        e = GPSEphemData->dEccentricity;

        Ek = KepplerSolver(Mk,e);


        // get the clock model coefficients
        af0 = GPSEphemData->dA_f0;
        af1 = GPSEphemData->dA_f1;
        af2 = GPSEphemData->dA_f2;

        // calculate the relativistic correction
        delta_tr = GPS_F * e * sqrt(A) * sin(Ek);

        toc = GPSEphemData->dTOC;

        // calculate the satellite time offset
        delta_ts = af0 + af1*(ts - toc) + af2*(ts-toc)*(ts-toc) + delta_tr;

        // correct delta_ts for iono group delay
        TGD = GPSEphemData->dTGD;
        delta_ts = delta_ts - TGD;

        
        // correct time for satellite clock offset
        ts = dTransmissionTime - delta_ts;
       
        // make revised delta_ts

        toe = GPSEphemData->dTOE;
        A = GPSEphemData->dA;
        n0 = sqrt(GPS_MU/(A*A*A));
        tk = ts - toe;

        // check for end of week roll-over
        if (tk > 302400)
		{
            tk = tk - 604800;
		} else if (tk < -302400)
		{
            tk = tk + 604800;
        }

        

        n = n0 + GPSEphemData->dDeltaN;
        Mk = GPSEphemData->dM0 + n*tk;
        e = GPSEphemData->dEccentricity;

        Ek = KepplerSolver(Mk,e);


        // get the clock model coefficients
        af0 = GPSEphemData->dA_f0;
        af1 = GPSEphemData->dA_f1;
        af2 = GPSEphemData->dA_f2;

        // calculate the relativistic correction
        delta_tr = GPS_F * e * sqrt(A) * sin(Ek);

        toc = GPSEphemData->dTOC;

        // calculate the satellite time offset
        delta_ts = af0 + af1*(ts-toc) + af2*(ts-toc)*(ts-toc) + delta_tr;

        // correct delta_ts for iono group delay
        TGD = GPSEphemData->dTGD;
        delta_ts = delta_ts - TGD;

        // calculate the GPS system time
        calcGPSTime = ts - delta_ts;                   

        // calculate the orbital coeficients or extract from navigation
        // message

        // revise estimate of Ek based on updated tk
        tk = calcGPSTime - toe;
        
        // check for end of week roll-over
        if (tk > 302400)
		{
            tk = tk - 604800;
		} else if (tk < -302400)
		{
            tk = tk + 604800;
        }
        
        
        n = n0 + GPSEphemData->dDeltaN;
        Mk = GPSEphemData->dM0 + n*tk;
        e = GPSEphemData->dEccentricity;

        Ek = KepplerSolver(Mk,e);

        vk = atan2((sqrt(1 - e*e)*sin(Ek)/(1-e*cos(Ek))),((cos(Ek) - e)/(1 - e*cos(Ek))));
        PHIk = vk + GPSEphemData->dOmega;

        // orbital pertubations
        Cus = GPSEphemData->dcus;
        Cuc = GPSEphemData->dcuc;
        Crs = GPSEphemData->dcrs;
        Crc = GPSEphemData->dcrc;
        Cis = GPSEphemData->dcis;
        Cic = GPSEphemData->dcic;


        delta_uk = Cus * sin(2*PHIk) + Cuc * cos(2*PHIk);
        delta_rk = Crs * sin(2*PHIk) + Crc * cos(2*PHIk);
        delta_ik = Cis * sin(2*PHIk) + Cic * cos(2*PHIk);

        IDOT = GPSEphemData->dInclination_dot;
        i0 = GPSEphemData->dInclination0;

        uk = PHIk + delta_uk;
        rk = A*(1 - e*cos(Ek)) + delta_rk;
        ik = i0 + delta_ik + IDOT * tk;

        xk_dash = rk * cos(uk);
        yk_dash = rk * sin(uk);

        OMEGA0 = GPSEphemData->dOmega0;
        OMEGAdot = GPSEphemData->dOmega_dot;

        OMEGAk = OMEGA0 + (OMEGAdot - GPS_OMEGAEDOT) * tk - GPS_OMEGAEDOT*toe;

        // final satellite positions in ecef
        xk = xk_dash * cos(OMEGAk) - yk_dash * cos(ik)*sin(OMEGAk);
        yk = xk_dash * sin(OMEGAk) + yk_dash * cos(ik)*cos(OMEGAk);
        zk = yk_dash * sin(ik);

       
		// Troy - testing 
        SVPosVector[0] = xk;
        SVPosVector[1] = yk;
        SVPosVector[2] = zk;
        SVPosVector[3] = delta_ts;
        ValidData = 1;

        //NOTE: Don't access SVPos directly like this below: because it 
        //made some funny problems when running the code
        //e.g. the output SVPos[0] was 0.0000
		//SVPos[0] = xk;
        //SVPos[1] = yk;
        //SVPos[2] = zk;
        //SVPos[3] = delta_ts;
        //ValidData = 1;
       
		memcpy(SVPos,SVPosVector,sizeof(SVPosVector));

		return ValidData;
}

int GPSOrbitPropagatorVelocities(double dTransmissionTime, 
			unsigned short usPRN, 
			const GPSEPHEM *GPSEphemData, 
			double dEphemValidityTime, 
			double *SVVel,
			double *SVAcc)
{

	double ts,delta_ts,toc,TGD,A,n0,n;
	double toe,af0,af1,af2,vk,delta_tr;
	double tk,e,Ek,Mk;
	double uk,rk,ik,xk,yk,zk;
//	int ValidData;
	double OMEGA0, OMEGAdot, OMEGAk;
	double Cis,Cus,Cic,Cuc,Crs,Crc;
	double IDOT,i0,xk_dash,yk_dash;
	double delta_ik,delta_uk,delta_rk;
	double PHIk;	
//	double SVPosVector[4];
	double calcGPSTime;
	double Ekdot,Ekdotdot,Mkdot;
	double tmp,tmp2;
	double xkacc,ykacc,zkacc;
	double xkvel,ykvel,zkvel;
	double OMEGAkdot;
	double mu_div_rk3;
	double xkdot_dash,ykdot_dash;
	double delta_tr_dot,delta_tr_dotdot,delta_tsdot,delta_tsdotdot;
	double PHIkdot,ukdot,rkdot,ikdot;

	double SVVelVector[4];
	double SVAccVector[4];
	
        ts = dTransmissionTime;

        // make initial estimate of Ek - this is revised later with the
        // corrected tk

        toe = GPSEphemData->dTOE;
        A = GPSEphemData->dA;
        n0 = sqrt(GPS_MU/(A*A*A));
        tk = ts - toe;

        // check for end of week roll-over
        if (tk > 302400)
		{
            tk = tk - 604800;
		} else if (tk < -302400)
		{
            tk = tk + 604800;
        }

        

        n = n0 + GPSEphemData->dDeltaN;
        Mk = GPSEphemData->dM0 + n*tk;

        Mkdot = sqrt(GPS_MU/(A*A*A)) + GPSEphemData->dDeltaN; //derivative of Mk in rad/s

        e = GPSEphemData->dEccentricity;
        Ek = KepplerSolver(Mk,e);

		Ekdot = Mkdot/(1.0-e*cos(Ek));

        Ekdotdot = -(Mkdot*e*sin(Ek)*Ekdot)/((1.0-e*cos(Ek))*(1.0-e*cos(Ek)));


        // get the clock model coefficients
        af0 = GPSEphemData->dA_f0;
        af1 = GPSEphemData->dA_f1;
        af2 = GPSEphemData->dA_f2;

        // calculate the relativistic correction
        delta_tr = GPS_F * e * sqrt(A) * sin(Ek);



        delta_tr_dot = GPS_F * e * sqrt(A) * cos(Ek)*Ekdot; // I don't know if this is correct or not but think it is TB 1/9/05
        delta_tr_dotdot = GPS_F * e * sqrt(A) * cos(Ek)*Ekdotdot - Ekdot*GPS_F*e*sqrt(A)*sin(Ek)*Ekdot; // I don't know if this is correct or not but seems to be, in the same order of magnitude if you just
        // subtract the velocities TB 1/9/05..about 3e-16 seconds which is negligable

        toc = GPSEphemData->dTOC;

        // calculate the satellite time offset
        delta_ts = af0 + af1*(ts - toc) + af2*(ts-toc)*(ts-toc) + delta_tr;

        // calculate the satellite time drift
        delta_tsdot = af1 + 2.0*af2*(ts - toc) + delta_tr_dot;

        // calculate the satellite time drift rate (acceleration)
        delta_tsdotdot =  2.0*af2 + delta_tr_dotdot;

        // correct delta_ts for iono group delay
        TGD = GPSEphemData->dTGD;
        delta_ts = delta_ts - TGD;




        // calculate the GPS system time
        calcGPSTime = ts - delta_ts; 


        // revise estimate of Ek based on updated tk
        tk = calcGPSTime - toe;

        // check for end of week roll-over
        if (tk > 302400)
		{
            tk = tk - 604800;
		} else if (tk < -302400)
		{
            tk = tk + 604800;
        }
        
        n = n0 + GPSEphemData->dDeltaN;
        Mk = GPSEphemData->dM0 + n*tk;

        Mkdot = sqrt(GPS_MU/(A*A*A)) + GPSEphemData->dDeltaN; //derivative of Mk in rad/s

        e = GPSEphemData->dEccentricity;
        Ek = KepplerSolver(Mk,e);
		Ekdot = Mkdot/(1.0-e*cos(Ek));
        Ekdotdot = -(Mkdot*e*sin(Ek)*Ekdot)/((1.0-e*cos(Ek))*(1.0-e*cos(Ek)));


        vk = atan2((sqrt(1.0 - e*e)*sin(Ek)),(cos(Ek) - e));
        

        PHIk = vk + GPSEphemData->dOmega;


		PHIkdot = sqrt(1.0 - e*e)*Ekdot/(1.0-e*cos(Ek));

        // orbital pertubations
        Cus = GPSEphemData->dcus;
        Cuc = GPSEphemData->dcuc;
        Crs = GPSEphemData->dcrs;
        Crc = GPSEphemData->dcrc;
        Cis = GPSEphemData->dcis;
        Cic = GPSEphemData->dcic;


        delta_uk = Cus * sin(2*PHIk) + Cuc * cos(2*PHIk);
        delta_rk = Crs * sin(2*PHIk) + Crc * cos(2*PHIk);
        delta_ik = Cis * sin(2*PHIk) + Cic * cos(2*PHIk);

        IDOT = GPSEphemData->dInclination_dot;
        i0 = GPSEphemData->dInclination0;

        uk = PHIk + delta_uk;
        ukdot = PHIkdot*(1.0 + 2.0*(Cus*cos(2.0*PHIk) - Cuc*sin(2.0*PHIk)));




        rk = A*(1 - e*cos(Ek)) + delta_rk;
        rkdot = A*e*sin(Ek)*Ekdot + 2.0*PHIkdot*(Crs*cos(2.0*PHIk) - Crc*sin(2.0*PHIk));


        ik = i0 + delta_ik + IDOT * tk;
        ikdot = IDOT + 2.0*PHIkdot*(Cis*cos(2.0*PHIk) - Cic*sin(2.0*PHIk));


        xk_dash = rk * cos(uk);
        yk_dash = rk * sin(uk);


        xkdot_dash = rkdot*cos(uk) - yk_dash*ukdot;
        ykdot_dash = rkdot*sin(uk) + xk_dash*ukdot;


        OMEGA0 = GPSEphemData->dOmega0;
        OMEGAdot = GPSEphemData->dOmega_dot;

        OMEGAk = OMEGA0 + (OMEGAdot - GPS_OMEGAEDOT) * tk - GPS_OMEGAEDOT*toe;


        OMEGAkdot = OMEGAdot - GPS_OMEGAEDOT;

        // final satellite positions in ecef
        xk = xk_dash * cos(OMEGAk) - yk_dash * cos(ik)*sin(OMEGAk);
        yk = xk_dash * sin(OMEGAk) + yk_dash * cos(ik)*cos(OMEGAk);
        zk = yk_dash * sin(ik);

        // final satellite velocities in ecef

        tmp = ykdot_dash*cos(ik) - yk_dash*sin(ik)*ikdot;

        xkvel = -OMEGAkdot*yk + xkdot_dash*cos(OMEGAk) - tmp*sin(OMEGAk);
        ykvel = OMEGAkdot*xk + xkdot_dash*sin(OMEGAk) + tmp*cos(OMEGAk);
        zkvel = yk_dash*cos(ik)*ikdot + ykdot_dash*sin(ik);

        mu_div_rk3 = - GPS_MU/(rk*rk*rk);
        tmp2 = mu_div_rk3 + GPS_OMEGAEDOT*GPS_OMEGAEDOT;

        xkacc = tmp2*xk + 2.0*ykvel*GPS_OMEGAEDOT;
        ykacc = tmp2*yk - 2.0*xkvel*GPS_OMEGAEDOT;
        zkacc = mu_div_rk3*zk;

		
		
        SVVelVector[0] = xkvel;
        SVVelVector[1] = ykvel;
        SVVelVector[2] = zkvel;
        SVVelVector[3] = delta_tsdot;
       // ValidData = 1;

	       

	
        SVAccVector[0] = xkacc;
        SVAccVector[1] = ykacc;
        SVAccVector[2] = zkacc;
        SVAccVector[3] = delta_tsdotdot;
        //ValidData = 1;

	       
		memcpy(SVVel,SVVelVector,sizeof(SVVelVector));
		memcpy(SVAcc,SVAccVector,sizeof(SVAccVector));

        return 1;
}


/** 
 * Keppler Solver iteraly solves for the eccentric anomoly using iterative
 * method described on page 164 of GPS theory and applciations volume 1
 * written by Duncan Greer (c) CRCSS 2005 Converted to C 9 OCT 2007
 */
double KepplerSolver(double Mk, double e)
{

	// M error threshold
	double Threshold = 1e-10;
	double E0,M,M_check;
	int n;
	double E_new,E;
	double M_error;
	
	// initial guess
	M=Mk;
	E0 = M + e*sin(M) / (1 - sin(M+e) + sin(M));

	// do iterations until the error is less than the thresshold
	n = 1;
	E = E0;

	while (1)
	{
		E_new = E - (E - e*sin(E) - M) / (1 - e*cos(E));
		
		// break loop when M can be solved from E
		M_check = E_new - e * sin(E_new);
		
		M_error = M - M_check;
		if(fabs(M_error) < Threshold)
		{
			break;
		}
		n = n + 1;
		M = M_check;
		
		if (n > 10000)
		{
			fprintf(stderr,"WARNING: KepplerSolver: n=%d, M=%f\n",
				n, M);
			return E_new;
		}
	}

	return E_new;

}


/**
 * \function WGS84_ECEF2LLH
 * Implements ellipsoidal conversion from ECEF to Lat Long Height using WGS84 paramaters
 *
 * \param double *ecef 3x1 vector containing ecef coordinates to be converted
 * \param double *llh 3x1 vector where result is to be stored
 * \return 0 on success, 1 on failure
 */
int WGS84_ECEF2LLH(const double *ecef, double *llh)
{
	double X = ecef[0];
	double Y = ecef[1];
	double Z = ecef[2];


	double P,Theta;
	double Latitude,Longitude,Height;
	double n;

	double cosTheta3,sinTheta3;

	P = sqrt(X*X + Y*Y);

	Theta = atan((Z * WGS84_a) / (P * WGS84_b));

	sinTheta3 = sin(Theta)*sin(Theta)*sin(Theta);
	cosTheta3 = cos(Theta)*cos(Theta)*cos(Theta);

	Latitude = atan((Z + WGS84_EPsq * WGS84_b * sinTheta3) / 
			(P - WGS84_Esq  * WGS84_a * cosTheta3));
	
	
	Longitude = atan2(Y,X);
	n = WGS84_a*WGS84_a / sqrt(WGS84_a*WGS84_a * cos(Latitude)*cos(Latitude) + 
		WGS84_b*WGS84_b*sin(Latitude)*sin(Latitude));
	Height = P / cos(Latitude) - n;

	llh[0] = Latitude;
	llh[1] = Longitude;
	llh[2] = Height;

	return 0;
}



int ICD200_IonoModel(GPSTime gtTimeStamp, double *Xu, double *Su, double *alpha, double *beta,
	double *IonoDelay, double *Azimuth, double *Elevation)
{

	double phi_u,lambda_u,h_u;
	double LLHu[3];
	double SV_LOS_ecef[3];
	double TMatrix_ECEF2ENU[3][3];
	double SV_LOS_enu[3];
	double x,AMP,D_IONO,T_IONO,F;
	double t_local;
	double GPSTimeSecs;
	double PER;
	double phi_i,lambda_i,phi_m,psi,A,E;
	
	GPSTimeSecs = (double)gtTimeStamp.fGPSSecondOfWeek;

	// get user geodetic latitude and longitude
	WGS84_ECEF2LLH(Xu,LLHu);

	phi_u = LLHu[0];
	lambda_u = LLHu[1];
	h_u = LLHu[2];

//	printf("phi: %f lambda: %f h: %f ",phi_u,lambda_u,h_u);

	// calculate satellite line of sight vector
	SV_LOS_ecef[0] = Su[0] - Xu[0];
	SV_LOS_ecef[1] = Su[1] - Xu[1];
	SV_LOS_ecef[2] = Su[2] - Xu[2];


	// generate rotation matrix
	T_ECEF2ENU(phi_u,lambda_u,(double **)TMatrix_ECEF2ENU);

	// convert to LTP coordinates
	MatMul331(TMatrix_ECEF2ENU,SV_LOS_ecef,SV_LOS_enu);

	// calculate azimuth and elevation
	A = atan2(SV_LOS_enu[0],SV_LOS_enu[1]) / GPS_PI; // semi-circles
	
	E = atan2(SV_LOS_enu[2],sqrt((SV_LOS_enu[0]*SV_LOS_enu[0]) + (SV_LOS_enu[1]*SV_LOS_enu[1]))) / GPS_PI; // semi-circles


	// if the elevation is less than 0, return a delay of 0
	if (E < 0.0)
	{
		D_IONO = 0.0;

		*IonoDelay = D_IONO;
		*Azimuth = A*GPS_PI;// note that A and E are in semicircles, not radians
		*Elevation = E*GPS_PI;

		return 1;
	}

	// calculate the Earth-centred angle, psi
	psi = (0.0137 / (E + 0.11)) - 0.022; // semi-circles

	// calculate the subionospheric latitude, phi_i
	phi_i = phi_u/GPS_PI + psi * (cos(A*GPS_PI));

	// limit phi_i
	if (phi_i > 0.416)
	{
		phi_i = 0.416;
	} else 
	if (phi_i < -0.416)
	{
		phi_i = -0.416;
	}


	// geodetic longitude of the earth projection of hte ionopsheric
	// intersection point (subionospheric longitude)
	lambda_i = lambda_u/GPS_PI + psi * (sin(A*GPS_PI)/cos(phi_i*GPS_PI)); // semi-circles

	// geomagnetic latitude of the earth projection of the ionophseric intersection point
	phi_m = phi_i + 0.064 * cos((lambda_i - 1.617)*GPS_PI); // semi-circles

	// calculate local time
	t_local = 4.32e4 * lambda_i + GPSTimeSecs; // seconds

	// limit local time to +- 86400 seconds (24 hours)
	if (t_local < -86400.0)
	{
		t_local = fmod(t_local,-86400.0);
	} else 
	if (t_local >= 86400.0)
	{
		t_local = fmod(t_local,86400.0);
	}

	// calculate the slant factor, F
	F = 1.0 + 16.0 * pow((0.53 - E),3);


	// calculate the amplitude and phase of the delay
	AMP = alpha[0] + 
	      alpha[1] * phi_m + 
	      alpha[2] * phi_m*phi_m + 
	      alpha[3] * phi_m*phi_m*phi_m;

	if (AMP < 0.0)
	{
		AMP = 0.0;
	}

	PER = beta[0] + 
	      beta[1] * phi_m + 
	      beta[2] * phi_m*phi_m + 
	      beta[3] * phi_m*phi_m*phi_m;

	if (PER < 72000.0)
	{
		PER = 72000.0;
	}

	// calculate x - what is x?!??
	x = 2.0 * GPS_PI * (t_local - 50400.0) / PER;

	// calculate the delay in seconds using x
	if (fabs(x) < 1.57)
	{
		T_IONO = F * (5e-9 + AMP*(1.0 - (x*x)/2.0 + (x*x*x*x)/24.0));
	} else
	{
		T_IONO = F * 5e-9;
	}

	// convert seconds to metres by multiplying by speed of light.
	D_IONO = T_IONO * GPS_SPEEDOFLIGHT;

	*IonoDelay = D_IONO;
	*Azimuth = A*GPS_PI;// note that A and E are in semicircles, not radians
	*Elevation = E*GPS_PI;
	return 0;
}



void T_ECEF2ENU(double Latitude, double Longitude, double **Result)
{
	double clong,clat;
	double slong,slat;

	double LocalResult[3][3];
	
	clong = cos(Longitude);
	clat = cos(Latitude);
	slong = sin(Longitude);
	slat = sin(Latitude);
	
	LocalResult[0][0] = -slong;
	LocalResult[0][1] = clong;
	LocalResult[0][2] = 0;
	
	LocalResult[1][0] = -clong*slat;
	LocalResult[1][1] = -slong*slat;
	LocalResult[1][2] = clat;
	
	LocalResult[2][0] = clong*clat;
	LocalResult[2][1] = slong*clat;
	LocalResult[2][2] = slat;
/*	
	LocalResult[0][0] = -slong*clat;
	LocalResult[0][1] = clong*clat;
	LocalResult[0][2] = slat;
	
	LocalResult[1][0] = -clong*slat;
	LocalResult[1][1] = -slong*slat;
	LocalResult[1][2] = clat;
	
	LocalResult[2][0] = clong*clat;
	LocalResult[2][1] = slong*clat;
	LocalResult[2][2] = slat;
*/
	memcpy(Result,LocalResult,sizeof(double)*3*3);

}

/**
 *	\fn TropoDelay
 *   Calculates the ionospheric delay for a satellite pseudorange from a simplified tropo model.
 *	\author Duncan Greer
 *
 */
double TropoDelayModel(double Elevation,double UserHeight)
{


	return 2.47 * exp(-0.133*UserHeight/1000.0) / (sin(Elevation) + 0.0121);

}



void T_Body2NED(double PHI, double THETA, double PSI, double **Result)
{

/*% generate  body frame to Fixed frame (NED) transformation matrix based on Euler
% angles. 
%by Troy Bruggemann 14 June. Taken from p 102 Nelson.  
%this has been verified with other sources to be correct (seeing Nelson has many mistakes in it)
%note: the transpose of this matrix can be used for NED2Body
%transformation. */



	double Cphi,Ctheta,Cpsi;
	double Sphi,Stheta,Spsi;
	double LocalResult[3][3];
	
	
	Cpsi = cos(PSI);
	Ctheta = cos(THETA);
	Cphi =  cos(PHI);

	Spsi = sin(PSI);
	Stheta = sin(THETA);
	Sphi = sin(PHI);
		

	
	LocalResult[0][0] = Ctheta*Cpsi;
	LocalResult[0][1] = Sphi*Stheta*Cpsi - Cphi*Spsi;
	LocalResult[0][2] = Cphi*Stheta*Cpsi + Sphi*Spsi;
	
	LocalResult[1][0] = Ctheta*Spsi;
	LocalResult[1][1] = Sphi*Stheta*Spsi + Cphi*Cpsi;
	LocalResult[1][2] = Cphi*Stheta*Spsi - Sphi*Cpsi;
	
	LocalResult[2][0] = -Stheta;
	LocalResult[2][1] = Sphi*Ctheta;
	LocalResult[2][2] = Cphi*Ctheta;



	memcpy(Result,LocalResult,sizeof(double)*3*3);

}




void T_ECEF2NED(double Latitude, double Longitude, double **Result)
{


	double clong,clat;
	double slong,slat;
	double LocalResult[3][3];

		
	clong = cos(Longitude);
	clat = cos(Latitude);
	slong = sin(Longitude);
	slat = sin(Latitude);

	LocalResult[0][0] = -clong*slat;
	LocalResult[0][1] = -slong*slat;
	LocalResult[0][2] =  clat;

	LocalResult[1][0] = -slong;
	LocalResult[1][1] = clong;
	LocalResult[1][2] = 0.0;

	LocalResult[2][0] = -clong*clat;
	LocalResult[2][1] = -slong*clat;
	LocalResult[2][2] = -slat;

	memcpy(Result,LocalResult,sizeof(double)*3*3);

}

void T_NED2ECEF(double Latitude, double Longitude, double **Result)
{


	double clong,clat;
	double slong,slat;
	double LocalResult[3][3];

		
	clong = cos(Longitude);
	clat = cos(Latitude);
	slong = sin(Longitude);
	slat = sin(Latitude);

	LocalResult[0][0] = -clong*slat;
	LocalResult[1][0] = -slong*slat;
	LocalResult[2][0] =  clat;

	LocalResult[0][1] = -slong;
	LocalResult[1][1] = clong;
	LocalResult[2][1] = 0.0;

	LocalResult[0][2] = -clong*clat;
	LocalResult[1][2] = -slong*clat;
	LocalResult[2][2] = -slat;

	memcpy(Result,LocalResult,sizeof(double)*3*3);

}










int WGS84_LLH2ECEF(const double *llh,  double *ecef)
{
	double Latitude = llh[0];
	double Longitude = llh[1];
	double Height = llh[2];

	double N, X, Y, Z;

	double e;

	//% calculate the first eccentricity

	e = sqrt(WGS84_a*WGS84_a - WGS84_b*WGS84_b) / WGS84_a;


	N = WGS84_a*WGS84_a / sqrt(WGS84_a*WGS84_a * cos(Latitude)*cos(Latitude) + WGS84_b*WGS84_b * sin(Latitude)*sin(Latitude));

	X = (N + Height) * cos(Latitude) * cos (Longitude);
	Y = (N + Height) * cos(Latitude) * sin (Longitude);
	Z = ((1 - e*e)*N + Height) * sin(Latitude);

	ecef[0] = X;
	ecef[1] = Y;
	ecef[2] = Z;


	return 0;
}


void WGS84_calcRnRe(double lat,  double *Rn, double *Re)
{
	//intput lat:
	//output Rn = Rnorth = Rmeridian (north -south)
	// Re = Reast = Rprime (east-west)
	//a = 6378137.0;   //% semi-major axis (metres)
	//f = 1/298.2572;
	double f = WGS84_Esq; //% flattening	
	double e2 = f * (2-f); //% eccentricity squared
	//double e = sqrt(e2);   //% first eccentricity


	//Rmeridian (Rnorth-south, called RM in WGS_84
	*Rn = WGS84_a * (1 - e2) / ((sqrt(1 - e2 * sin(lat)*sin(lat)))*(sqrt(1 - e2 * sin(lat)*sin(lat)))*(sqrt(1 - e2 * sin(lat)*sin(lat))));

	//R prime (Reast-west), called RN in WGS_84 //%find the normal radius of curvature
	*Re = WGS84_a / sqrt(1 - e2 * sin(lat)*sin(lat));


}

