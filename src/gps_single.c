/** 
 * \file gps_single.c Implements GPS single point solution using Least Squares
 * \author Duncan Greer
 * $Id: gps_single.c 3579 2010-06-26 10:57:31Z greerd $
 */



#include <stdio.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include "../../novatelstore/src/novatel.h"
#include "gps.h"
#include "../../novatelstore/src/shmdef.h"

#include "matrix.h"

/** Controls whether or not the program should keep running. */
int iShutdown = 0;

#define MODULE_NAME "[GPS_SINGLE]"

#define log(string) 	fprintf(stdout,"%s %s\n", MODULE_NAME, string)
#define err(string)	fprintf(stderr,"%s ERROR: %s\n", MODULE_NAME, string)

#define ELEVATION_CUTOFF_RAD (10.0 * DEG2RAD)


/** 
	\function CalculateDOP
	Calculates the Dilution of Precision based on the provided G Matrix.
	\author Duncan Greer and Troy Bruggemann
	
	Note that the G-Matrix must be provided already transformed to 
	the horizontal plane.
*/
void CalculateDOP(gsl_matrix *G, double *DOP);


/**
 *	\fn HandleSignal
 *	\brief Handle caught signals
 *	
 * 	This function sets a flag (iShutdown) to shut the program down.
 */
void HandleSignal(int signal) 
{ 
	fprintf(stdout,"Got Signal: %d\n",signal);
	iShutdown = 1; 
}

/**
 * 	\fn Usage
 *	\brief Print Usage Information
 */
void Usage(char *arg)
{

	/* options */
	/*
	    b = baudrate
	    p = port
	    l = logging 
	    f = log filename
	    t = enable tcp server
	*/
	fprintf(stdout,"Usage: %s <options>\n",arg);
	fprintf(stdout,"    -l Log results to file (OFF)\n");
	fprintf(stdout,"    -f log filename\n");
	fprintf(stdout,"    -L Use specified log file path (/data/gps)\n");

}

int main (int argc, char **argv)
{
	gsl_matrix *G;
	gsl_matrix *G_LTP;
//	gsl_matrix *GT;
//	gsl_matrix *tempNN;
//	gsl_matrix *V;
//	gsl_vector *S;
//	gsl_matrix *U;
	gsl_vector *W;
	gsl_vector *y;
	gsl_vector *delta_y;
	gsl_vector *x;
	gsl_vector *delta_x;
	gsl_vector *x_error;
	gsl_vector *x_true;
	gsl_vector *temp1;
	
	gsl_multifit_linear_workspace *gsl_workspace;
	gsl_vector *SVPos[MAX_SV];
	double SVPos_calc[4];
	double dSSE;
	gsl_matrix *SolutionCovariance;
	int iReturnValue = 0;

	int iNumberStates = 4;
	int iNumberMeasurements = 9;
	int iNumberValidEphemeris = 0;
	int iNumberValidObservations = 0;
	int iIndex;
	GPSTime gtObservationTime;
	unsigned short usPRN;
	int iIteration;

	int iUseSV[NOVATEL_MAXPRNS];

	double user_x,user_y,user_z,user_t = 0.0;
	double sat_x,sat_y,sat_z,sat_t;
//	double sat_x_calc,sat_y_calc,sat_z_calc,sat_t_calc;
	double dTemp;
	
	double range;
	double delta_pr_omegaedot;

	double ecef[3];
	double llh[3];
	double TECEF2ENU[3][3];
	gsl_matrix *T_ECEF2NED_44;
	double NEDError[3];
	double ECEFError[3];
	double UserPos[4];
	double ION_Alpha[4];
	double ION_Beta[4];
	double IonoDelay[NOVATEL_MAXCHANNELS];
	double TropoDelay[NOVATEL_MAXCHANNELS];
	double SV_Azimuth[NOVATEL_MAXCHANNELS];
	double SV_Elevation[NOVATEL_MAXCHANNELS];

	double dTransmissionTime;
	
	GPSEPHEM  GPSEphemData[NOVATEL_MAXPRNS];
	RANGE GPSRangeData;

	double GDOP,PDOP;
	double DOPS[5];
	double limit;

	int iUseSHM = 1;
	int iSHM_FID;

	int iRow,iCol;
	
	shm_struct *shm = NULL;
	
	unsigned long ulSignalType = 0;	
	
	char copt;

	/* log files */
	FILE *fpLogFile = NULL;    /* pointer to log file stream */
	char azcLogFileName[256];  /* log file name */
	char azcLogFilePath[255];
	char azcLogFileNameWithPath[1024];

	struct tm *timestruct;
	time_t tCurrentTime = 0;
	time_t tLogFileStart = 0;

	/* options flags */
	struct
	{
		int iLogResults;
		int iCustomLogFile;
	} options;

	// default options
	options.iLogResults = 0;
	options.iCustomLogFile = 0;

	
	
	/* setup signal handlers */
	signal(SIGHUP, HandleSignal);
	signal(SIGTERM, HandleSignal);
	signal(SIGKILL, HandleSignal);
	signal(SIGALRM, HandleSignal);
	signal(SIGINT, HandleSignal);
	
	/* get the current time to use as a timestamp in the log file name */
	tCurrentTime = time(NULL);
	timestruct = gmtime(&tCurrentTime);

	/* setting default log file path */
	sprintf(azcLogFilePath,"%s","/data/gps");

	/* read input parameters */
	while((copt = getopt(argc,argv,"lf:L:")) && copt != -1)
	{
	  switch(copt)
	  {
	   case 'l':
		/* log results to file */
		options.iLogResults = 1;
		log("Logging enabled");
		strftime(azcLogFileName,sizeof(azcLogFileName),"gps_single_%g%m%d%H%M%S.csv",timestruct);
		break;

	  case 'f':
		/* set the log file name */
		sprintf(azcLogFileName,"%s",optarg);
		options.iCustomLogFile = 1; /* since we're using a user specified log file name, we dont want to recycle it every hour */
		break;

	  case 'L':
	  	/* setting log file path */
	  	sprintf(azcLogFilePath,"%s",optarg);

	  	fprintf(stdout,"%s Setting log file path to %s",MODULE_NAME,azcLogFilePath);
	  	break;
	  case '?':
		Usage(argv[0]);
		return 0;
		break;
	  }
	}



	
	for(iIndex=0;iIndex<NOVATEL_MAXPRNS;iIndex++)
	{
		iUseSV[iIndex] = 0;
	}
	
	fprintf(stdout,"Allocating Matrices\n");
	x = gsl_vector_alloc(iNumberStates);
	delta_x = gsl_vector_alloc(iNumberStates);
	x_true = gsl_vector_alloc(iNumberStates);
	x_error = gsl_vector_alloc(iNumberStates);
	SolutionCovariance = gsl_matrix_alloc(iNumberStates,iNumberStates);



	temp1 = gsl_vector_alloc(3);

	
	for(iIndex=0;iIndex<MAX_SV;iIndex++)
	{
		SVPos[iIndex] = gsl_vector_alloc(4);
	}




	T_ECEF2NED_44 = gsl_matrix_alloc(4,4);

       	/* setup shared memory */
	if(iUseSHM == 1)
	{
	  iSHM_FID = shm_open("/novatel0",O_RDWR,0777);
	  if(iSHM_FID == -1)
	  {
		perror("shm_open: ");
		iUseSHM = 0;
 	  } else
	  {
	     	/* memory map the shm object */
	     	shm = mmap(0,sizeof(*shm),PROT_READ,MAP_SHARED,iSHM_FID,0);
	     	if(shm == MAP_FAILED)
	     	{
		   perror("shm mmap: ");
		   iUseSHM = 0;
	     	} else
		{
		  fprintf(stdout,"%s SHM Address Is: 0x%08x\n",MODULE_NAME,(unsigned int)shm);
		}
	  }
	}



	/* setup time */
	gtObservationTime.usGPSWeek = 1416;
	gtObservationTime.fGPSSecondOfWeek = 135900.0;

	/* setup measurements */
/*	GPSRangeData.usPRN[0] = 10;
	GPSRangeData.usPRN[1] = 12;
	GPSRangeData.usPRN[2] = 30;
	GPSRangeData.usPRN[3] = 29;
	GPSRangeData.usPRN[4] = 7;
	GPSRangeData.usPRN[5] = 26;
	GPSRangeData.usPRN[6] = 22;
	GPSRangeData.usPRN[7] = 24;
	GPSRangeData.dPseudorange[0] = 23677691.85548;
	GPSRangeData.dPseudorange[1] = 24609189.10747;
	GPSRangeData.dPseudorange[2] = 22288602.64148;
	GPSRangeData.dPseudorange[3] = 22641003.92448;
	GPSRangeData.dPseudorange[4] = 20042814.45649;
	GPSRangeData.dPseudorange[5] = 23024522.23148;
	GPSRangeData.dPseudorange[6] = 24837761.16146;
	GPSRangeData.dPseudorange[7] = 21114887.48349;
*/
	/* setup sat positions */
	
	/* estimated position - starting point */
	gsl_vector_set(x,0,-5046773.0);
	gsl_vector_set(x,1, 2568446.0);
	gsl_vector_set(x,2,-2925289.0);
	gsl_vector_set(x,3, 0.0);
	
	llh[0] = -0.47;
	llh[1] = 2.53;
	llh[2] = 0.0;

	//gsl_vector_fprintf(stdout,x,"%10.3f");

	
	/* true position */
	gsl_vector_set(x_true,0,-5053034.84 );
	gsl_vector_set(x_true,1, 2562951.08);
	gsl_vector_set(x_true,2,-2919245.94);
	gsl_vector_set(x_true,3,0.0);
	
/*
	// PRN 10
	gsl_vector_set(SVPos[10],0,-11106779.576);
	gsl_vector_set(SVPos[10],1,-10162086.179);
	gsl_vector_set(SVPos[10],2,-21984106.076);
	gsl_vector_set(SVPos[10],3, 98.034009e-6*GPS_SPEEDOFLIGHT);
	
	// PRN 12
	gsl_vector_set(SVPos[12],0,-20114520.008);
	gsl_vector_set(SVPos[12],1,  6748657.642);
	gsl_vector_set(SVPos[12],2, 16058860.089);
	gsl_vector_set(SVPos[12],3,-47.814630e-6*GPS_SPEEDOFLIGHT);

	// PRN 30
	gsl_vector_set(SVPos[30],0, -19111202.845);
	gsl_vector_set(SVPos[30],1, 17011909.480);
	gsl_vector_set(SVPos[30],2, 6599201.305);
	gsl_vector_set(SVPos[30],3, 26.513082e-6*GPS_SPEEDOFLIGHT);
	
	
	// PRN 29
	gsl_vector_set(SVPos[29],0, -22931680.210);
	gsl_vector_set(SVPos[29],1, -10007127.790);
	gsl_vector_set(SVPos[29],2,  -9208106.350);
	gsl_vector_set(SVPos[29],3, 358.103982e-6*GPS_SPEEDOFLIGHT);
	
	// PRN 7
	gsl_vector_set(SVPos[7],0, -18000552.275);
	gsl_vector_set(SVPos[7],1,  10196913.772);
	gsl_vector_set(SVPos[7],2, -16363130.228);
	gsl_vector_set(SVPos[7],3, 402.849928e-6*GPS_SPEEDOFLIGHT);
	
	// PRN 26
	gsl_vector_set(SVPos[26],0, -24783876.733);
	gsl_vector_set(SVPos[26],1,  -8910991.974);
	gsl_vector_set(SVPos[26],2,  -5742918.813);
	gsl_vector_set(SVPos[26],3, -62.236139e-6*GPS_SPEEDOFLIGHT);
	

	// PRN 22
	gsl_vector_set(SVPos[22],0, -8616964.914);
	gsl_vector_set(SVPos[22],1, 22105308.026);
	gsl_vector_set(SVPos[22],2, 12074730.766);
	gsl_vector_set(SVPos[22],3, 169.425839e-6*GPS_SPEEDOFLIGHT);

	// PRN 24
	gsl_vector_set(SVPos[24],0, -15110809.416);
	gsl_vector_set(SVPos[24],1,   2514877.198);
	gsl_vector_set(SVPos[24],2, -21514889.442);
	gsl_vector_set(SVPos[24],3, 80.656626e-6*GPS_SPEEDOFLIGHT);
	
	*/
	
	/* open log file */
	if(options.iLogResults == 1)
	{
		sprintf(azcLogFileNameWithPath,"%s/%s",azcLogFilePath,azcLogFileName);
	
		if((fpLogFile = fopen(azcLogFileNameWithPath,"a")) == NULL)
		{	
			fprintf(stderr,"Error opening %s for writing: %s\n",azcLogFileNameWithPath,strerror(errno));
			
			return -1;
		} else
		{
			tLogFileStart = tCurrentTime;
			fprintf(stdout,"Using log file: %s\n",azcLogFileNameWithPath);
		}

		
	}


  while(iShutdown == 0)
  {	

	/* if log file has been opend for > 1 hour, start a new one */
	if(options.iLogResults == 1 && options.iCustomLogFile == 0)
	{
		tCurrentTime = time(NULL);
		if(difftime(tCurrentTime,tLogFileStart) > 3600.0)
		{
			/* start a new log file */
			fclose(fpLogFile);
					
			timestruct = gmtime(&tCurrentTime);
			strftime(azcLogFileName,sizeof(azcLogFileName),"gps_single_%g%m%d%H%M%S.csv",timestruct);
		
			sprintf(azcLogFileNameWithPath,"%s/%s",azcLogFilePath,azcLogFileName);

			if((fpLogFile = fopen(azcLogFileNameWithPath,"a")) == NULL)
			{	
				fprintf(stderr,"Error opening %s for writing: %s\n",azcLogFileNameWithPath,strerror(errno));
				options.iLogResults = 0;
			} else
			{
				tLogFileStart = tCurrentTime;
				fprintf(stdout,"%s Started New Logfile: %s\n",MODULE_NAME,azcLogFileNameWithPath);
				
			}
		}

	}

	/* TODO:  replace this sleep with a semaphore to wait for new data */
  	sleep(1);
	
	fprintf(stdout,"====== Single-Point Solution ======\n");
	/* retrive the current data from shm */
	if(iUseSHM == 1)
	{
	  memcpy(&GPSEphemData,&shm->CurrentGPSEPHEM,sizeof(GPSEPHEM)*NOVATEL_MAXPRNS);
	  
	  
	  // get the current observation time to compare with ephemeris
	  memcpy(&gtObservationTime,&shm->CurrentRANGE.gtTimeStamp,sizeof(GPSTime));


	  /* ensure we have ephemeris for all measurements, otherwise exclude */
	  iNumberValidEphemeris = 0;
	  for(iIndex=0;iIndex<NOVATEL_MAXPRNS;iIndex++)
	  {
	  	if((GPSEphemData[iIndex].ulPRN != 0) &&
		   (GPSEphemData[iIndex].ulHealth == 0) &&
		   (GPSEphemData[iIndex].dTOE < (gtObservationTime.fGPSSecondOfWeek+7500.0)) &&
		   (GPSEphemData[iIndex].dTOE > (gtObservationTime.fGPSSecondOfWeek-7500.0))
	  	)
	  	{
	  	  iUseSV[iIndex] = 1;
	  	  iNumberValidEphemeris++;
//	  	  fprintf(stdout,"SV%d ",iIndex);
	  	} else
	  	{
	  	  iUseSV[iIndex] = 0;
	  	}
	  }
	  
	  /* loop through measurements and see 
	     if we have valid ephem for each */
	  
	  memcpy(&GPSRangeData.gtTimeStamp,&shm->CurrentRANGE.gtTimeStamp,sizeof(GPSTime));
	  iNumberValidObservations = 0;
	  
	  for(iIndex=0;iIndex<shm->CurrentRANGE.lNumberObservations;iIndex++)
	  {
	  	if(iUseSV[shm->CurrentRANGE.usPRN[iIndex]] == 1)
	  	{
	  		ulSignalType = shm->CurrentRANGE.ulTrackingStatus[iIndex] & 0x3E00000;
	  		
	  		//  check that we're only getting hte L1 signal
	  		if (ulSignalType == 0)
	  		{
	  			// use current range data for now
	  			memcpy(&GPSRangeData.usPRN[iNumberValidObservations],
				        &shm->CurrentRANGE.usPRN[iIndex],sizeof(unsigned short));
				memcpy(&GPSRangeData.dPseudorange[iNumberValidObservations],
					&shm->CurrentRANGE.dPseudorange[iIndex], sizeof(double));
			
			
					iNumberValidObservations++;
	  		
	  		} else
	  		{
	  			

				// not an L1 obs
				//fprintf(stdout,"Value of Signal Type for PRN%02d was 0x%0x\n",
				//	shm->CurrentRANGE.usPRN[iIndex], 
				//	ulSignalType);
	  		}
	  	} else
	  	{
	  	  fprintf(stdout,"NE: %d; ",shm->CurrentRANGE.usPRN[iIndex]);
	  	}
	  }

	  fprintf(stdout,"Obs: %d; \n",iNumberValidObservations);	  
	  GPSRangeData.lNumberObservations = iNumberValidObservations;
	  
	  iNumberMeasurements = GPSRangeData.lNumberObservations;
	  
	  if(iNumberMeasurements < 4)
	  {
		fprintf(stdout,"Less than 4 measurements\n");
		continue;
	  }

	/*
		  if(iNumberMeasurements > 9)
		  {
			 iNumberMeasurements = 9; 
		  }
		  
		  if(iNumberMeasurements < 9)
		  {
			err("Less than 9 measurements!");
		  }
	*/
		G = gsl_matrix_alloc(iNumberMeasurements,iNumberStates);
		G_LTP = gsl_matrix_alloc(iNumberMeasurements,iNumberStates);
		y = gsl_vector_alloc(iNumberMeasurements);
		delta_y = gsl_vector_alloc(iNumberMeasurements);
		W = gsl_vector_alloc(iNumberMeasurements);
		gsl_workspace = gsl_multifit_linear_alloc (iNumberMeasurements,iNumberStates);

		/* calculate G-matrix */
		gsl_matrix_set_zero (G);
		//gsl_matrix_fprintf(stdout,G,"%5.3f");

		/* calculate W-matrix */
		gsl_vector_set_all(W,1.0);




	} else
	{
		log("No Shared Memory");



	}

		for(iIteration=0;iIteration<20;iIteration++)
		{
		//	fprintf(stdout,"\n===== ITERATION %d =====\n",iIteration+1);

			user_x = gsl_vector_get(x,0);
			user_y = gsl_vector_get(x,1);
			user_z = gsl_vector_get(x,2);
			user_t = gsl_vector_get(x,3);
			
			UserPos[0] = user_x;
			UserPos[1] = user_y;
			UserPos[2] = user_z;
			UserPos[3] = user_t;

			ION_Alpha[0] = shm->CurrentIONUTC.a0;
			ION_Alpha[1] = shm->CurrentIONUTC.a1;
			ION_Alpha[2] = shm->CurrentIONUTC.a2;
			ION_Alpha[3] = shm->CurrentIONUTC.a3;
			
			ION_Beta[0] = shm->CurrentIONUTC.b0;
			ION_Beta[1] = shm->CurrentIONUTC.b1;
			ION_Beta[2] = shm->CurrentIONUTC.b2;
			ION_Beta[3] = shm->CurrentIONUTC.b3;
			
			
		
	//		fprintf(stdout,"User: [%+10.3f\t%+10.3f\t%+10.3f\t%+10.3f]\n",user_x,user_y,user_z,user_t);
			
			for (iIndex=0;iIndex<iNumberMeasurements;iIndex++)
			{
		
				dTransmissionTime = (double)gtObservationTime.fGPSSecondOfWeek - GPSRangeData.dPseudorange[iIndex]/GPS_SPEEDOFLIGHT;
				
				GPSOrbitPropagator(dTransmissionTime, 
						GPSRangeData.usPRN[iIndex],
						&GPSEphemData[GPSRangeData.usPRN[iIndex]], 
						7500.0, 
						SVPos_calc);

				
				
				usPRN = GPSRangeData.usPRN[iIndex];
				
				sat_x = SVPos_calc[0];
				sat_y = SVPos_calc[1];
				sat_z = SVPos_calc[2];
				sat_t = SVPos_calc[3]*GPS_SPEEDOFLIGHT;
				
				// find the range
				//gsl_vector_memcpy(temp1,SVPos[usPRN]);
				//gsl_vector_sub(temp1,x);
				gsl_vector_set(temp1,0,sat_x - user_x);
				gsl_vector_set(temp1,1,sat_y - user_y);
				gsl_vector_set(temp1,2,sat_z - user_z);
				
				range = gsl_blas_dnrm2(temp1);
		
				gsl_matrix_set(G,iIndex,0,-(sat_x - user_x)/range);
				gsl_matrix_set(G,iIndex,1,-(sat_y - user_y)/range);
				gsl_matrix_set(G,iIndex,2,-(sat_z - user_z)/range);
				gsl_matrix_set(G,iIndex,3,1.0);
		
				gsl_vector_set(y,iIndex,GPSRangeData.dPseudorange[iIndex]);

				// calculate earth rotation correction
				delta_pr_omegaedot = -(GPS_OMEGAEDOT/GPS_SPEEDOFLIGHT) * (sat_x*user_y - sat_y*user_x);

				// calculate iono correction
				ICD200_IonoModel(gtObservationTime, 
						UserPos, 
						SVPos_calc, 
						ION_Alpha, 
						ION_Beta,
						&IonoDelay[iIndex],
						&SV_Azimuth[iIndex],
						&SV_Elevation[iIndex]);

				// examine sv elevations - if less than cutoff, set corresponding row of 'G' to zeros
				if (SV_Elevation[iIndex] < ELEVATION_CUTOFF_RAD)
				{
					gsl_matrix_set(G,iIndex,0,0.0);
					gsl_matrix_set(G,iIndex,1,0.0);
					gsl_matrix_set(G,iIndex,2,0.0);
					gsl_matrix_set(G,iIndex,3,0.0);
					fprintf(stdout,"Ignoring PRN%02d due low elevation (%3.1fdeg)\n",
						GPSRangeData.usPRN[iIndex],SV_Elevation[iIndex]*RAD2DEG);
					
					TropoDelay[iIndex] = 0.0;
					gsl_vector_set(delta_y,iIndex,0.0);					
				} else
				{

					TropoDelay[iIndex] = TropoDelayModel(SV_Elevation[iIndex],llh[2]);
				
				
					gsl_vector_set(delta_y,iIndex,GPSRangeData.dPseudorange[iIndex]-range-(user_t-sat_t)+delta_pr_omegaedot-IonoDelay[iIndex]-TropoDelay[iIndex]);
				}		
	//			fprintf(stdout,"Sat %d:\t[%+10.3f\t%+10.3f\t%+10.3f\t%+10.3f] range=%f, delta_pr=%f\n",usPRN,sat_x,sat_y,sat_z,sat_t,range,delta_pr_omegaedot);
	//			fprintf(stdout,"Sat %d:\t[%+10.3f\t%+10.3f\t%+10.3f\t%+10.3f]\n",usPRN,sat_x,sat_y,sat_z,sat_t);
				
				
			}

			
	//		fprintf(stdout,"delta_y: \n");
	//		gsl_vector_fprintf(stdout,delta_y,"%10.3f");
			
//			fprintf(stdout,"G-Matrix: \n");
//			for (iIndex=0;iIndex<iNumberMeasurements;iIndex++)
//			{
//				fprintf(stdout,"[%+10.3f\t%+10.3f\t%+10.3f\t%+10.3f]\n",
//					gsl_matrix_get(G,iIndex,0),
//					gsl_matrix_get(G,iIndex,1),
//					gsl_matrix_get(G,iIndex,2),
//					gsl_matrix_get(G,iIndex,3));
//					
//			}


			/* weight measurements according to elevation */
//			for(iIndex=0;iIndex<iNumberMeasurements;iIndex++)
//			{
//			  gsl_vector_set(W,iIndex,SV_Elevation[iIndex]);
//			}
			
			//gsl_vector_fprintf(stdout,W,"%f");

			iReturnValue = gsl_multifit_wlinear(G, 
							W, 
							delta_y, 
							delta_x, 
							SolutionCovariance, 
							&dSSE, 
							gsl_workspace);
			
			/*iReturnValue = gsl_multifit_linear(G, 
							delta_y, 
							delta_x, 
							SolutionCovariance, 
							&dSSE, 
							gsl_workspace);
			*/
			gsl_vector_add(x,delta_x);


		//	fprintf(stdout,"Least-squares fit returned value: %d; RSS=%f\n",iReturnValue,sqrt(dSSE));
		//	fprintf(stdout,"delta_x: \n");
		//	gsl_vector_fprintf(stdout,delta_x,"%20.12f");

			gsl_vector_memcpy(x_error,x);
			gsl_vector_sub(x_error,x_true);


			
			/* test if we should break out */
			limit = gsl_blas_dnrm2(delta_x);
			if (limit < 1e-8)
			{	
				fprintf(stdout,"%f:Least-Squares solved in %d iterations; RSS=%f\n",
					gtObservationTime.fGPSSecondOfWeek,iIteration+1,sqrt(dSSE));
				
				/* convert to LLH */
				ecef[0] = gsl_vector_get(x,0);
				ecef[1] = gsl_vector_get(x,1);
				ecef[2] = gsl_vector_get(x,2);
				
				WGS84_ECEF2LLH(ecef, llh);
				T_ECEF2ENU(llh[0], llh[1], (double **)TECEF2ENU);

	//			fprintf(stdout,"[%f %f %f]\n",TECEF2ENU[0][0],TECEF2ENU[0][1],TECEF2ENU[0][2]);
	//			fprintf(stdout,"[%f %f %f]\n",TECEF2ENU[1][0],TECEF2ENU[1][1],TECEF2ENU[1][2]);
	//			fprintf(stdout,"[%f %f %f]\n",TECEF2ENU[2][0],TECEF2ENU[2][1],TECEF2ENU[2][2]);
				

				ECEFError[0] = gsl_vector_get(x_error,0);
				ECEFError[1] = gsl_vector_get(x_error,1);
				ECEFError[2] = gsl_vector_get(x_error,2);
				
				MatMul331(TECEF2ENU, ECEFError, NEDError);

				/* must swap ENU to NED */
				dTemp = NEDError[0];
				NEDError[0] = NEDError[1];
				NEDError[1] = dTemp;
				NEDError[2] = -NEDError[2];

				//fprintf(stdout,"X error: \n");
				//gsl_vector_fprintf(stdout,x_error,"%10.5f m");
				
				//fprintf(stdout,"X: \n");
				//gsl_vector_fprintf(stdout,x,"%10.5f m");

				/* calculate DOPS */
				for(iRow=0;iRow<3;iRow++)
				{
				  for(iCol=0;iCol<3;iCol++)
				  {
					// note that TECEF2ENU must be transposed
					gsl_matrix_set(T_ECEF2NED_44,
						iRow,iCol,TECEF2ENU[iCol][iRow]); 
				  }
				  gsl_matrix_set(T_ECEF2NED_44,iRow,3,0.0);
				}
				for(iCol=0;iCol<3;iCol++)
				{
				  gsl_matrix_set(T_ECEF2NED_44,3,iCol,0.0);
				}  
				gsl_matrix_set(T_ECEF2NED_44,3,3,1.0);
				
				gsl_matrix_set_zero(G_LTP);
				gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 
					1.0, G, T_ECEF2NED_44, 1.0, G_LTP);
				
				CalculateDOP(G_LTP,DOPS);
				GDOP = DOPS[0];
				PDOP = DOPS[1];


				/* write results to shared memory */
//				memcpy(&shm->dDOPS[0],&DOPS[0],sizeof(double));
//				memcpy(shm->dLLH_LSQ,llh,sizeof(double)*3);
//				memcpy(shm->dECEF_LSQ,ecef,sizeof(double)*3);
//				shm->dRxClock_LSQ = gsl_vector_get(x,3);
				

				fprintf(stdout,"Lat: %f, Long: %f, Height: %f\n",
					llh[0]*RAD2DEG,llh[1]*RAD2DEG,llh[2]);
				fprintf(stdout,"GDOP=%f, PDOP=%f, HDOP=%f, VDOP=%f, TDOP=%f\n",
					GDOP,PDOP,DOPS[2],DOPS[3],DOPS[4]);
				fprintf(stdout,"NED Error: %fm %fm %fm\n",NEDError[0],
					NEDError[1],
					NEDError[2]);

				fprintf(stdout,"SV: ");
				for(iIndex=0;iIndex<iNumberMeasurements;iIndex++)
				{
				  fprintf(stdout,"%d ",GPSRangeData.usPRN[iIndex]);
				}
				fprintf(stdout,"\n");

				fprintf(stdout,"Iono: ");
				for(iIndex=0;iIndex<iNumberMeasurements;iIndex++)
				{
				  fprintf(stdout,"%3.1fm ",IonoDelay[iIndex]);
				}
				fprintf(stdout,"\n");

				fprintf(stdout,"Tropo: ");
				for(iIndex=0;iIndex<iNumberMeasurements;iIndex++)
				{
				  fprintf(stdout,"%3.1fm ",TropoDelay[iIndex]);
				}
				fprintf(stdout,"\n");

				fprintf(stdout,"Elevation: ");
				for(iIndex=0;iIndex<iNumberMeasurements;iIndex++)
				{
				  fprintf(stdout,"%3.1f ",SV_Elevation[iIndex]*RAD2DEG);
				}
				fprintf(stdout,"\n");
				fprintf(stdout,"Azimuth: ");
				for(iIndex=0;iIndex<iNumberMeasurements;iIndex++)
				{
				  fprintf(stdout,"%3.1f ",SV_Azimuth[iIndex]*RAD2DEG);
				}
				fprintf(stdout,"\n");
		
				// check if we should log to disk
				if(options.iLogResults == 1 && fpLogFile != NULL)
				{
					
					fprintf(fpLogFile, "%d,%f,%d,%lf,%d,%lf,%lf,%lf,%.10lf,%.10lf,%lf,%f,%f,%f,%f,%f\n",
						gtObservationTime.usGPSWeek,
						gtObservationTime.fGPSSecondOfWeek,
						iIteration+1,
						sqrt(dSSE),
						iNumberMeasurements,
						ecef[0],
						ecef[1],
						ecef[2],
						llh[0],
						llh[1],
						llh[2],
						DOPS[0],
						DOPS[1],
						DOPS[2],
						DOPS[3],
						DOPS[4]
					);
					fflush(fpLogFile);
				}


				break;
				
				

			}

			
		}

		if (iIteration == 20)
		{
			fprintf(stderr,"Iterated 20 times with no solution (limit=%f)", limit);
		}

		gsl_matrix_free(G);
		gsl_matrix_free(G_LTP);
		gsl_vector_free(W);
		gsl_vector_free(y);
		gsl_vector_free(delta_y);
		gsl_multifit_linear_free(gsl_workspace);
	}


	fclose(fpLogFile);

	fprintf(stdout,"Freeing Matrices\n");


	gsl_vector_free(x);
	gsl_vector_free(delta_x);
	gsl_vector_free(x_true);
	gsl_vector_free(x_error);

	gsl_vector_free(temp1);
	gsl_matrix_free(SolutionCovariance);

	
	gsl_matrix_free(T_ECEF2NED_44);
	
	for(iIndex=0;iIndex<MAX_SV;iIndex++)
	{
		gsl_vector_free(SVPos[iIndex]);
	}

	return 0;
}


void CalculateDOP(gsl_matrix *G, double *DOP)
{
	gsl_matrix *GT;
	gsl_matrix *tempNN;
	gsl_matrix *V;
	gsl_vector *S;
	gsl_vector *Sinv;
	gsl_matrix *SinvMat;
	gsl_matrix *A;  // = inv(GT*G);
	gsl_vector *SVDWork;
	gsl_matrix *U;
	
	int iReturnValue,iIndex;
	int iNumberMeasurements,iNumberStates;
	
	double var_x,var_y,var_z,var_dt;
	
	iNumberMeasurements = G->size1;
	iNumberStates = G->size2;

	GT = gsl_matrix_alloc(iNumberStates,iNumberMeasurements);
	tempNN = gsl_matrix_alloc(iNumberStates,iNumberStates);
	V = gsl_matrix_alloc(iNumberStates,iNumberStates);
	S = gsl_vector_alloc(iNumberStates);
	Sinv = gsl_vector_alloc(iNumberStates);
	SinvMat = gsl_matrix_alloc(iNumberStates,iNumberStates);
	A = gsl_matrix_alloc(iNumberStates,iNumberStates);
	SVDWork = gsl_vector_alloc(iNumberStates);
	U = gsl_matrix_alloc(iNumberStates,iNumberStates);
	
	
	// calculate DOP
	gsl_matrix_transpose_memcpy(GT,G);
	gsl_matrix_set_zero(tempNN);
	gsl_blas_dgemm (CblasNoTrans,CblasNoTrans, 1.0, GT, G, 1.0, tempNN);
	
	// take the SVD of GT*T (stored in TempNN) to find inverse - tempNN now contains 'U'
	iReturnValue =  gsl_linalg_SV_decomp (tempNN,V,S,SVDWork);
		
	// get element wise inverse of S
	gsl_vector_set_all(Sinv,1.0);
	gsl_vector_div(Sinv, S);

	gsl_matrix_set_zero (SinvMat);
	for(iIndex=0;iIndex<iNumberStates;iIndex++)
	{
		gsl_matrix_set(SinvMat,iIndex,iIndex,gsl_vector_get(Sinv,iIndex));
	}
		
	// inv(S)*UT
	gsl_matrix_memcpy(U,tempNN);
	gsl_matrix_set_zero (tempNN);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, SinvMat, U, 1.0, tempNN);
	
	// A = V * inv(S) * UT
	gsl_matrix_set_zero(A);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, V, tempNN, 1.0, A);


	var_x = gsl_matrix_get(A,0,0);
	var_y = gsl_matrix_get(A,1,1);
	var_z = gsl_matrix_get(A,2,2);
	var_dt = gsl_matrix_get(A,3,3);

	DOP[0] = sqrt(var_x+var_y+var_z+var_dt);
	DOP[1] = sqrt(var_x+var_y+var_z);
	DOP[2] = sqrt(var_x+var_y);
	DOP[3] = sqrt(var_z);
	DOP[4] = sqrt(var_dt);


	gsl_matrix_free(GT);
	gsl_matrix_free(tempNN);

	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(Sinv);
	gsl_matrix_free(U);
	gsl_vector_free(SVDWork);
	gsl_matrix_free(SinvMat);
	
}

