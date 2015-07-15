/**
 *	
 *	\file gpsins.c - Implements an integrated GPS-Inertial navigation system.
 *
 * 	This program implements an integrated GPS-Inertial navigation system.  Data is 
 *	retrieved from shared memory and processed in two threads.  The high speed thread
 *	implements the inertial navigation mechanisation.  The low speed thread combines
 *	the inertial navigation solution with GPS measurements in a tightly coupled
 *	architecture.  The output is stored back into shared memory.
 *
 * 	This program implements a 17 state tightly coupled EKF.  Note 
 * 	that the filter states are organised such that the first 15 states
 *	are identical to those used by a 15 state loosely coupled EKF.  
 * 
 *	The filter uses an error state formulation.  The EKF_x_hat vector contains the state
 *	error rather than the full state.  The full state is stored in regular arrays of
 *	Pos_LLH, Vel_NED and C_BN representing the Latitude, Longitude, Height, North, East
 *	and Down Velocity, and attitude direction cosine matrix (DCM).  
 *
 *      \code
 *	Filter States
 *	=============
 *	x1  - Latitude Error (radians)
 * 	x2  - Longitude Error (radians)
 * 	x3  - Height Error (metres)
 * 	x4  - North Velocity Error (m/s)
 * 	x5  - East Velocity Error (m/s)
 *	x6  - Down Velocity Error (m/s)
 * 	x7  - North Tilt Error (radians)
 *	x8  - East Tilt Error (radians)
 *	x9  - Down Tilt Error (radians) 
 * 	x10 - X Accelerometer Bias Error (m/s/s)
 *	x11 - Y Accelerometer Bias Error (m/s/s)
 * 	x12 - Z Accelerometer Bias Error (m/s/s)
 * 	x13 - X Gyroscope Bias Error (rad/s)
 * 	x14 - Y Gyroscope Bias Error (rad/s)
 * 	x15 - Z Gyroscope Bias Error (rad/s)
 * 	x16 - Receiver Clock Bias Error (metres)
 * 	x17 - Receiver Clock Frequency Error (metres/sec)
 *
 *	Measurement States - assuming 'n' satellite measurements available
 * 	==================
 * 	z1 - zn - GPS Pseudorange Measurements
 * 	zn+1 - 2zn - GPS Pseudorange Rate (doppler) measurements
 *
 *	These states are not yet implemented...
 *	2nz+1 - Barometric Altitude Error (metres)
 *	2zn+2 - Compass Heading Error (radians)
 *	\endcode
 *
 *	\author Duncan Greer
 *	\version "$Id: gpsins.c 955 2007-12-06 06:38:04Z greerd $"
 */


/* this is an additonal comment */

/* if we are using rtai, declare it here */
/* #define USE_RTAI */
#undef USE_RTAI

#include <stdio.h>
#include <termios.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <strings.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>
#include <stdlib.h>
#include <signal.h>

#ifndef _SVID_SOURCE
#define _SVID_SOURCE
#endif

#include <sys/shm.h>
#include <sys/ipc.h>

#include <mqueue.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <spawn.h>


#include <sys/mman.h>

#ifdef USE_RTAI
#include <rtai_lxrt.h>
#include <rtai_bits.h>
#else
#include <sys/ioctl.h>
#include <linux/rtc.h>
#endif

#include <math.h>
#include "gpsins.h"
#include "gps.h"
#include "../../novatelstore/src/shmdef.h"
#include "matrix.h"

#define MODULE_NAME 	"[GPSINS]"


#define log(string)	fprintf(stdout,"%s %s\n", MODULE_NAME, string)
#define err(string)	fprintf(stderr,"%s ERROR: %s\n",MODULE_NAME, string)


#define CONTROL_MQUEUE_NAME 	"/control-queue"
#define GPSINS_MQUEUE_NAME 	"/gpsins-queue"
#define	INS_MQUEUE_NAME		"/ins-queue"


#define CMD_SHUTDOWN 1
#define CMD_RUNLOOP 2


#define MAX_DIR_LENGTH 1024

//int Shutdown = 0;

/** thread to control execution of the program */
void *ControlThread(void *ThreadArgs);

/** thread to perform INS solution - 100Hz */
void *INSThread(void *ThreadArgs);

/** thread to perform GPS-INS solution - 1 Hz */
void *GPSINSThread(void *ThreadArgs);


/* INS globals */
double C_BN[3][3];
double Pos_LLH[3];
double Vel_NED[3];
double Acc_NED[3];
double UserClock[2];
double GyroBias[3];
double AccBias[3];

int iShutdown = 0;

/* state variables */
gsl_matrix *EKF_P;
gsl_vector *EKF_x_hat;


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
* main
*/
int main(int argc, char *argv[])
{
	/* generic return value */
	int iReturnValue;

#ifndef USE_RTAI
	int iRTCDevice;
	long lRTCBuffer;
#endif

	/*
	int iMessage = 0;
	char cMessage = 0;
	*/

	message mMessage;

	/* thread paramter/control objects */
	pthread_t ptControlThread;
	pthread_t ptINSThread;
	pthread_t ptGPSINSThread;	


	/* pointers to GPS and IMU data respectively */
//	GPSData *ptGPSData;
//	IMUData *ptIMUData;

	/* shared memory objects */
//	int shmGPSId, shmIMUId;
//	key_t shmGPSKey, shmIMUKey;
	
	/* message queue */
	mqd_t	mqdControl;
	
	struct mq_attr mqdControlAttr;

	/** GPS INS Thread Message Queue (POSIX) */
	mqd_t	mqdGPSINS;

	/** INS Thread Message Queue (POSIX) */
	mqd_t	mqdINS;

	/** storage for gps store pid */
//	pid_t pidGPSStore;
	
	/* storage for current working directory */
//	char azcCurrentWorkingDirectory[MAX_DIR_LENGTH];
//	char azcLaunchPath[MAX_DIR_LENGTH+20];
//	char *ptr;

	/* storage for GPSStore Arguments */
//	char *azcGPSStoreArgs[4] = {"gpsstore","/dev/ttyUSB0","230400",NULL};


	sigset_t ssSigSet;

	log("Initialising Program");

	/* setup signal handlers */
	signal(SIGHUP, HandleSignal);
	signal(SIGTERM, HandleSignal);
	signal(SIGKILL, HandleSignal);
	signal(SIGALRM, HandleSignal);
	signal(SIGINT, HandleSignal);

	sigemptyset(&ssSigSet);
	sigaddset(&ssSigSet,SIGRTMIN);
	pthread_sigmask(SIG_BLOCK, &ssSigSet, NULL); 

	mqdControlAttr.mq_maxmsg = 10;
	mqdControlAttr.mq_msgsize = sizeof(message);
	mqdControlAttr.mq_flags = 0;

	/*  initialise message	*/
	memset(&mMessage,0,sizeof(message));

	/*  create a control message queue */
	mqdControl = mq_open(CONTROL_MQUEUE_NAME,O_CREAT|O_RDWR|O_EXCL|O_NONBLOCK, (mode_t)0666, &mqdControlAttr);
	if(mqdControl == -1)
	{
		perror("mq_open(): CONTROL_MQUEUE_NAME");
		exit(-1);
	}

	/* setup gps ins message queue */
#ifdef USE_RTAI
	mqdGPSINS = mq_open(GPSINS_MQUEUE_NAME,O_CREAT|O_RDWR|O_EXCL|O_NONBLOCK, (mode_t)0666, &mqdControlAttr);
#else
	mqdGPSINS = mq_open(GPSINS_MQUEUE_NAME,O_CREAT|O_RDWR|O_EXCL, (mode_t)0666, &mqdControlAttr);
#endif
	if(mqdGPSINS == -1)
	{
		perror("mq_open(): GPSINS_MQUEUE_NAME");
		exit(-1);
	}
	/* setup ins message queue */
	mqdINS = mq_open(INS_MQUEUE_NAME,O_CREAT|O_RDWR|O_EXCL|O_NONBLOCK, (mode_t)0666, &mqdControlAttr);
	if(mqdINS == -1)
	{
		perror("mq_open(): INS_MQUEUE_NAME");
		exit(-1);
	}
	/*
	ptmqdControl = (mqd_t *)malloc(sizeof(mqd_t));
	*ptmqdControl = mqdControl;
	*/

	/* start control thread */
	iReturnValue = pthread_create(&ptControlThread,NULL,ControlThread,(void *)&mqdControl);
	if(iReturnValue < 0)
	{
		perror("pthread_create():");
	}

	/* debug message */
	/*
	fprintf(stdout,"%s GPS Data: %d bytes, IMU Data: %d bytes\n",MODULE_NAME, sizeof(GPSData), sizeof(IMUData));
	*/
	
//	shmGPSKey = GPSDATA_KEY;
//	shmIMUKey = IMUDATA_KEY;

	/*. create shared memory space for storing the data */
//	if((shmGPSId = shmget(shmGPSKey, sizeof(GPSData), IPC_CREAT | 0666)) < 0) 
//	{
//      	perror("shmget(): allocating GPS data storage");
//   		exit(-1);
//	}

	/* attach */
//	if((ptGPSData = shmat(shmGPSId, NULL, 0)) == (GPSData *) -1) 
//	{
//       	perror("shmat(): attaching GPS Data Storage");
//        	exit(-1);
 //   	}

//	if((shmIMUId = shmget(shmIMUKey,sizeof(GPSData), IPC_CREAT | 0666)) < 0)
//	{
//		perror("shmget(): allocating IMU data storage");
//		exit(-1);
//	}
//	if((ptIMUData = shmat(shmIMUId, NULL, 0)) == (IMUData *)-1)
//	{
//		perror("shmat(): attaching IMU data storage");
//		exit(-1);
//	}

	/* intialise shared memory */
//	memset(ptGPSData,0,sizeof(GPSData));
//	memset(ptIMUData,0,sizeof(IMUData));

	/* get current working directory - used for spawn call */
//	ptr = getcwd(azcCurrentWorkingDirectory,(size_t) MAX_DIR_LENGTH);
//	strcat(azcLaunchPath,azcCurrentWorkingDirectory);
//	strcat(azcLaunchPath,"/");
//	strcat(azcLaunchPath,azcGPSStoreArgs[0]);

	/* start GPS and INS data stores */
//	fprintf(stdout,"%s Launching GPSStore: %s\n",MODULE_NAME,azcLaunchPath);
//	iReturnValue = posix_spawn(&pidGPSStore,		/* storage for child pid */
//		    azcCurrentWorkingDirectory, /* current working directory */
//		    NULL,			/* posix_spawn_file_actions_t */
//		    NULL,			/* posix_spawnattr_t */
//		    azcGPSStoreArgs,		/* GPS Store Args - serial port, baud */
//		    NULL);			/* environment */

//	if(iReturnValue != 0)
//	{
//		perror("posix_spawn(): GPSStore");
//	}

	/* start GPS and INS threads  */
	iReturnValue = pthread_create(&ptINSThread,NULL,INSThread,(void *)&mqdINS);
	if(iReturnValue < 0)
	{
		perror("pthread_create(): INS Thread");
	}

	iReturnValue = pthread_create(&ptGPSINSThread,NULL,GPSINSThread,(void *)&mqdGPSINS);
	if(iReturnValue < 0)
	{
		perror("pthread_create(): GPS-INS Thread");
	}
	
#ifndef USE_RTAI
	/* open real time clock if we are not using RTAI */
	iRTCDevice = open("/dev/rtc",O_RDONLY);
	if (iRTCDevice < 0)
	{
		perror("open: RTC: ");
	} else
	{
		ioctl(iRTCDevice,RTC_UIE_ON,NULL);
		log("Successfully opened RTC");
	}
#endif
	while(iShutdown == 0)
	{
		iReturnValue = mq_receive(mqdControl,(char *)&mMessage,sizeof(message),NULL);
		if(iReturnValue < 0)
		{
			/* 
				perror("mq_receive():"); 
			*/
		}
		/*
		log("Received Message");
		fprintf(stdout,"Rx: %d, %d\n",iReturnValue, mMessage.iMessageType);
		*/
		if(mMessage.iMessageType == 1)
		{
			// exit 
			log("Received Shutdown Message"); 
			iShutdown = 1; 
			break;
		}
		
		/* this read will return at 1 second intervals */
		read(iRTCDevice,&lRTCBuffer,sizeof(long));

		/* send message to run threads */
		mMessage.iMessageType = CMD_RUNLOOP;
		iReturnValue = mq_send(mqdGPSINS,(const char *)&mMessage,sizeof(message),1);
		


	}

	log("Shutting Down");

	/* send message to shutdown threads */
	mMessage.iMessageType = CMD_SHUTDOWN;
	iReturnValue = mq_send(mqdGPSINS,(const char *)&mMessage,sizeof(message),1);
	if(iReturnValue < 0)
	{
		perror("mq_send(): Thread Shutdown");
	}
	iReturnValue = mq_send(mqdINS,(const char *)&mMessage,sizeof(message),1);
	if(iReturnValue < 0)
	{
		perror("mq_send(): Thread Shutdown");
	}

	iReturnValue = mq_send(mqdControl,(const char *)&mMessage,sizeof(message),1);
	if(iReturnValue < 0)
	{
		perror("mq_send(): Thread Shutdown");
	}

	/* wait for threads to shutdown */
	//sleep(2);
	fprintf(stdout,"%s Waiting for INS thread to join...", MODULE_NAME);
	pthread_join(ptINSThread,NULL); 
	fprintf(stdout,"DONE\n");
	
	fprintf(stdout,"%s Waiting for GPS-INS thread to join...", MODULE_NAME);
	pthread_join(ptGPSINSThread,NULL);
	fprintf(stdout,"DONE\n");
	
	
	fprintf(stdout,"%s Waiting for Control thread to join...", MODULE_NAME);
	pthread_join(ptControlThread,NULL);
	fprintf(stdout,"DONE\n");
	
#ifndef USE_RTAI
	/* close rtc */
	close(iRTCDevice);
#endif

	/* close and unlink message queues */
	fprintf(stdout,"%s Unlinking Message Queues...",MODULE_NAME);
	mq_close(mqdControl);
	mq_unlink(CONTROL_MQUEUE_NAME);
	mq_close(mqdGPSINS);
	mq_unlink(GPSINS_MQUEUE_NAME);
	mq_close(mqdINS);
	mq_unlink(INS_MQUEUE_NAME);

	fprintf(stdout,"DONE\n");
	
	log("Exiting");
	pthread_exit(NULL);
	return 0;
} /* end main */
 


/*
* ControlThread - controls execution of the program
*/
void *ControlThread(void *ThreadArgs)
{
	char c[100];
	mqd_t *ptmqdControl;
	int iReturnValue;

	message mMessage;

	fprintf(stdout,"%s Control Thread Initialised - press 'q <enter>' to exit!\n",MODULE_NAME);

	/* initialise message	 */
	memset(&mMessage,0,sizeof(message));

	/* read arguments */
	ptmqdControl = (mqd_t *)ThreadArgs;

	while(scanf("%s",c))
	{
		/* check for message (non blocking) */
		iReturnValue = mq_receive(*ptmqdControl,(char *)&mMessage,sizeof(message),NULL);
		if(iReturnValue < 0)
		{
			perror("mq_recieve(): Control Thread");
		}

		if(mMessage.iMessageType == CMD_SHUTDOWN)
		{
			log("Control Thread Received SHUTDOWN command");
			//pthread_exit(NULL);
			break;
		}


		switch(strlen(c))
		{
			case 1:
				/* command character */
				switch(c[0])
				{
					case 'q':
						mMessage.iMessageType = CMD_SHUTDOWN;

						/* send shutdown message */
						iReturnValue = mq_send(*ptmqdControl,
								(const char *)&mMessage,
								sizeof(message),
								1);
						if(iReturnValue < 0)
						{
							perror("mq_send():");
						}
						break;

					default:
						fprintf(stderr,"%s Unknown command: %c\n",MODULE_NAME, c[0]);
						break;
				}
			default:
				/* do nothing */
				break;
		}

	}

	fprintf(stdout,"%s Control Thread Exiting\n", MODULE_NAME);

	/* never reached */
	pthread_exit(NULL);
}

int iINSThreadCounter;

void *INSThread(void *ThreadArgs)
{
	/* message storage */
	message mMessage;
	


	/* generic return storage */
	int iReturnValue;

	/* message queue pointer */
	mqd_t *ptmqdINS;	

	/* timer variables */
	int iSigNumber;
	sigset_t ssSigSet;
	struct sigevent seSigEvent;
	timer_t tPeriodicTimer;
	struct timespec first,period;
	struct itimerspec required,old;
	long lThreadPeriod = 10000000; /* 10ms = 100 hz */
	
	/* shared memory variables */
	int iUseSHM = 1;
	shm_struct *shm = NULL;
	int iSHM_FID;

	double EarthRate[3];
	double Tecef2ned[3][3];
	double RM,RP;
	double Lat,Long;
	double LatDot,LongDot;
	double sLat,cLat;
	double OMEGA_b[3][3];
	double Coriolis[3];
	double Acc_XYZ[3];
	double omega_x,omega_y,omega_z;
	double OMEGA_e_n[3];
	double OMEGA_in[3][3];
	
	fprintf(stdout,"%s INS Thread Started\n",MODULE_NAME);

	/* initialise message	*/
	memset(&mMessage,0,sizeof(message));
	
	/* setup message queue */
	ptmqdINS = (mqd_t *)ThreadArgs;

	/* setup thread timer */
	seSigEvent.sigev_notify = SIGEV_SIGNAL;
	seSigEvent.sigev_signo = SIGRTMIN;
	
	clock_gettime(CLOCK_REALTIME, &first);
	first.tv_sec = first.tv_sec + 1;
	
	fprintf(stdout,"INS: First loop runs at %ld.%ld\n", (long)first.tv_sec,(long)first.tv_nsec);

	period.tv_sec = 0;
	period.tv_nsec = lThreadPeriod;

	required.it_value = first;
	required.it_interval = period;

	iReturnValue = timer_create(CLOCK_REALTIME,&seSigEvent,&tPeriodicTimer);
	if(iReturnValue < 0)
	{
		perror("INS: timer_create():");
	}

	sigemptyset(&ssSigSet);
	sigaddset(&ssSigSet,SIGRTMIN);
	//pthread_sigmask(SIG_UNBLOCK, &ssSigSet, NULL);
	

	iReturnValue = timer_settime(tPeriodicTimer,TIMER_ABSTIME,&required, &old);
	if(iReturnValue < 0)
	{
		perror("INS: timer_settime():");
	}

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



	/* initialise INS parameters */
	C_BN[0][0] = 1.0;
	C_BN[1][1] = 1.0;
	C_BN[2][2] = 1.0;

	GyroBias[0] = 0.0;
	GyroBias[1] = 0.0;
	GyroBias[2] = 0.0;
	
	AccBias[0] = 0.0;
	AccBias[1] = 0.0;
	AccBias[2] = 0.0;
	

	/* position */
	Pos_LLH[0] = -27.5*DEG2RAD;
	Pos_LLH[1] = 153.0*DEG2RAD;
	Pos_LLH[2] = 90.0;

	/* velocity */
	Vel_NED[0] = 0.0;
	Vel_NED[1] = 0.0;
	Vel_NED[2] = 0.0;
	
	/* acceleration */
	Acc_NED[0] = 0.0;
	Acc_NED[1] = 0.0;
	Acc_NED[2] = -LOCAL_GRAVITY;
	
	Acc_XYZ[0] = 0.0;
	Acc_XYZ[1] = 0.0;
	Acc_XYZ[2] = -LOCAL_GRAVITY;
	

	EarthRate[0] = 0.0;
	EarthRate[1] = 0.0;
	EarthRate[2] = OMEGA_e;

	Coriolis[0] = 0.0;
	Coriolis[1] = 0.0;
	Coriolis[2] = 0.0;
	
	RM = WGS84_a;
	RP = WGS84_a;
	WGS84_calcRnRe(Pos_LLH[0],  &RM, &RP);
		
	OMEGA_e_n[0] = 0.0;
	OMEGA_e_n[1] = 0.0;
	OMEGA_e_n[2] = 0.0;
	

	log("INS Thread Commencing Periodic");

	while(1)
	{

		//fprintf(stdout,"Waiting...");
		/* wait for signal */
		sigwait(&ssSigSet, &iSigNumber);
		
		iINSThreadCounter++;


		/* check for messages */
		iReturnValue = mq_receive(*ptmqdINS,(char *)&mMessage,sizeof(message),NULL);
		if(iReturnValue < 0 && errno != EAGAIN)
		{
			perror("INS: mq_recieve(): ");
			//fprintf(stdout,"%d",iReturnValue);
		}
		
		if(mMessage.iMessageType == CMD_SHUTDOWN)
		{
			break;
		}

		
		/* get current lat and long */
		Lat = Pos_LLH[0];
		Long = Pos_LLH[1];
		cLat = cos(Lat);
		sLat = sin(Lat);
		
		T_ECEF2NED(Lat, Long, (double **)Tecef2ned);
		
		/* attitude update */
		omega_x = shm->CurrentIMUData.dRate[0] - GyroBias[0];
		omega_y = shm->CurrentIMUData.dRate[1]  - GyroBias[1];
		omega_z = shm->CurrentIMUData.dRate[2]  - GyroBias[2];
    
		/* form the skew symmetric form */
		OMEGA_b[0][0] = 0.0;
		OMEGA_b[0][1] = -omega_z;
		OMEGA_b[0][2] =  omega_y;
		OMEGA_b[1][0] =  omega_z;
		OMEGA_b[1][1] = 0.0;
		OMEGA_b[1][2] = -omega_x;
		OMEGA_b[2][0] = -omega_y;
		OMEGA_b[2][1] = omega_x;
		OMEGA_b[2][2] = 0.0;

		/* navigation frame rotation correction */
		OMEGA_in[0][0] = 0.0;
		OMEGA_in[0][1] = -OMEGA_e * sLat;
		OMEGA_in[0][2] = 0.0;
		OMEGA_in[1][0] = OMEGA_e * sLat;
		OMEGA_in[1][1] = 0.0;
		OMEGA_in[1][2] = -OMEGA_e * cLat;
		OMEGA_in[2][0] = 0.0;
		OMEGA_in[2][1] = OMEGA_e * cLat;
		OMEGA_in[2][2] = 0.0;
		
		/* DCM Update */
		/* C_BN = C_BN + (C_BN * OMEGA_b - OMEGA_in * C_BN)*INS_DT */

		/* acceleration calculation */
		/* Acc_NED = C_BN * Acc_XYZ */
		Acc_XYZ[0] = shm->CurrentIMUData.dAcc[0] - AccBias[0];
		Acc_XYZ[1] = shm->CurrentIMUData.dAcc[1] - AccBias[1];
		Acc_XYZ[2] = shm->CurrentIMUData.dAcc[2] - AccBias[2] -LOCAL_GRAVITY;

		MatMul331(C_BN, Acc_XYZ, Acc_NED);

		/* correct vertical channel for gravity */
		Acc_NED[2] = Acc_NED[2] + LOCAL_GRAVITY;

		/* coriolis correction */
		MatMul331(Tecef2ned, EarthRate, OMEGA_e_n);
		
		/* Coriolis = cross(2 * OMEGA_e_n,[0,0,0]); */

		Acc_NED[0] = Acc_NED[0] - Coriolis[0];
		Acc_NED[1] = Acc_NED[1] - Coriolis[1];
		Acc_NED[2] = Acc_NED[2] - Coriolis[2];
		

		/* velocity update */
		Vel_NED[0] = Vel_NED[0] + Acc_NED[0] * INS_DT;
		Vel_NED[1] = Vel_NED[1] + Acc_NED[1] * INS_DT;
		Vel_NED[2] = Vel_NED[2] + Acc_NED[2] * INS_DT;
		
		WGS84_calcRnRe(Pos_LLH[0],  &RM, &RP);
		LatDot = Vel_NED[0] / (RM + Pos_LLH[2]);
		LongDot = Vel_NED[1] / (cos(Pos_LLH[0])*(RP + Pos_LLH[2]));

		/* position update */
		Pos_LLH[0] = Pos_LLH[0] + LatDot * INS_DT;
		Pos_LLH[1] = Pos_LLH[1] + LongDot * INS_DT;
		Pos_LLH[2] = Pos_LLH[2] + -Vel_NED[2] * INS_DT;
		
		

		/* covariance update */
		/* P = PHI * P * PHI' + Qd */

	}
	/* exit */
	fprintf(stdout,"%s INS Thread Exiting\n",MODULE_NAME);
	pthread_exit(NULL);
}





void *GPSINSThread(void *ThreadArgs)
{
	int iRunThread = 1;

	/* pointer to message queue designator passed from arguments */
	mqd_t *ptmqdGPSINS;

	/* generic return value storage */
	int iReturnValue;

	/* generic message storage */
	message mMessage;

#ifdef USE_RTAI
	/* real time task */
	static RT_TASK *rttGPSINS;
#endif

	/* KF stuff */
	int iNumberStates,iNumberMeasurements;

	/* matrices */
	gsl_matrix *EKF_Q;
	gsl_matrix *EKF_PHI; /* phi ~= I + Fdt */
	gsl_matrix *EKF_F;
	gsl_matrix *EKF_P_minus;
	gsl_matrix *tempP;
		
	gsl_matrix *EKF_H;
	gsl_matrix *EKF_R;
	gsl_matrix *EKF_K;
	gsl_matrix *EKF_V;
	gsl_matrix *EKF_Vinv;
	gsl_matrix *tempV;
	gsl_matrix *tempK;

	gsl_matrix *KH;
	gsl_matrix *I_17by17;
	gsl_matrix *IminusKH;
	gsl_matrix *KHPminus;
	gsl_matrix *KRKT;
	gsl_matrix *tempKRKT;
	
	/* svd stuff */
	gsl_vector *SVD_S;
	gsl_vector *SVD_Sinv;
	gsl_matrix *SVD_SinvMat;
	gsl_vector *SVDWork;
	gsl_matrix *SVD_V;
	gsl_matrix *SVD_U;
	gsl_matrix *tempNN;




	/* vectors */
	gsl_vector *EKF_z;
	gsl_vector *EKF_deltaz;
	gsl_vector *EKF_x_hat_out;
	

	gsl_vector *temp1;

	double x_gyro_beta;
	double y_gyro_beta;
	double z_gyro_beta;
	double x_accel_beta;
	double y_accel_beta;
	double z_accel_beta;

	double user_x,user_y,user_z,user_t;
	double user_xdot,user_ydot,user_zdot,user_tdot;
	double sat_x,sat_y,sat_z,sat_t;
	double sat_xdot,sat_ydot,sat_zdot,sat_tdot;
	double ION_Alpha[4];
	double ION_Beta[4];
	unsigned short usPRN;
	double range,rangerate;
	double geo_vel_to_sat;
	double delta_pr_omegaedot;
	
	double Rm,Rp;
	double Rmh,Rph,Rg;
	
	double IonoDelay[NOVATEL_MAXCHANNELS],TropoDelay[NOVATEL_MAXCHANNELS];
	double SV_Azimuth[NOVATEL_MAXCHANNELS];
	double SV_Elevation[NOVATEL_MAXCHANNELS];
//	double llh[3];
//	double ecef[3];
	double dTransmissionTime;
	
	GPSTime gtObservationTime;
	RANGE	GPSRangeData;
	GPSEPHEM GPSEphemData[NOVATEL_MAXPRNS];
	int iIndex;
	int iNumberValidObservations;
	int iNumberValidEphemeris;
	int iUseSV[NOVATEL_MAXCHANNELS];
	
	double UserPos[3];
	double UserVel[3];
	
	double Tecef2ned[3][3];
	double Tned2ecef[3][3];
	
	double SVPos_calc[4];
	double SVVel_calc[4];
	double SVAcc_calc[4];
	
	double C_BN_new[3][3];

	double EKF_H_E[3][3];
	double EKF_H_LTP[3][3];

	double del_alpha,del_beta,del_gamma;
	double del_att_skew[3][3];
	
	/* matrix counter variables */
	int iRowIndex,iColIndex;

	/* shared memory variables */
	int iUseSHM = 1;
	shm_struct *shm = NULL;
	int iSHM_FID;
	
	/* log message */
	fprintf(stdout,"%s GPS-INS Thread Started\n",MODULE_NAME);

	/* initialise message	 */
	memset(&mMessage,0,sizeof(message));

	/* setup message queue */
	ptmqdGPSINS = (mqd_t *)ThreadArgs;


	/* setup EKF */
	iNumberStates = NUMBER_STATES;
	iNumberMeasurements = 0;
	
	EKF_x_hat = gsl_vector_alloc(iNumberStates);
	EKF_x_hat_out = gsl_vector_alloc(iNumberStates);
	
	EKF_Q = gsl_matrix_alloc(iNumberStates,iNumberStates);
	EKF_P = gsl_matrix_alloc(iNumberStates,iNumberStates);
	EKF_P_minus = gsl_matrix_alloc(iNumberStates,iNumberStates);
	tempP = gsl_matrix_alloc(iNumberStates,iNumberStates);
	
	EKF_PHI = gsl_matrix_alloc(iNumberStates,iNumberStates);
	EKF_F = gsl_matrix_alloc(iNumberStates,iNumberStates);
	
	I_17by17 = gsl_matrix_alloc(iNumberStates,iNumberStates);
	
	KH = gsl_matrix_alloc(iNumberStates,iNumberStates);
	IminusKH = gsl_matrix_alloc(iNumberStates,iNumberStates);
	KHPminus = gsl_matrix_alloc(iNumberStates,iNumberStates);

	gsl_matrix_set_identity(EKF_P);
	gsl_matrix_scale(EKF_P,10.0);
	
	temp1 = gsl_vector_alloc(3);

	x_gyro_beta = 1/100.0;
	y_gyro_beta = 1/100.0;
	z_gyro_beta = 1/100.0;
	
	x_accel_beta = 1/100.0;
	y_accel_beta = 1/100.0;
	z_accel_beta = 1/100.0;
	
	UserClock[0] = 0.0;
	UserClock[1] = 0.0;
	
#ifdef USE_RTAI
	/* make this thread real time */
	rt_set_oneshot_mode();
	start_rt_timer(0);
	rttGPSINS = rt_task_init_schmod(nam2num("GPSINS"),
				10,
				0,
				0,
				SCHED_FIFO,
				0);
 	if (rttGPSINS == NULL) 
	{
		perror("rt_task_init(): GPSINS");
		exit(-1);
	}


	rt_make_hard_real_time();
	rt_task_make_periodic(rttGPSINS,
		rt_get_time()+nano2count(1000000000)*2, 
		nano2count(1000000000));

#endif	

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

	/* reset INS thread counter */
	iINSThreadCounter = 0;

	while(iRunThread == 1)
	{
		/* wait for message */ /* non blocking  if using RTAI*/
		iReturnValue = mq_receive(*ptmqdGPSINS,(char *)&mMessage,sizeof(message),NULL); 
		if(iReturnValue < 0)
		{
			/* check value of errno */
			switch(errno)
			{	
				case EAGAIN:
					/* return with no message - do nothing */
					break;
				default:
					perror("mq_recieve(): GPSINSThread");
					break;
			}
			
		}
		
// 		if(mMessage.iMessageType == CMD_SHUTDOWN)
// 		{
// 			break;
// 		}

		switch(mMessage.iMessageType)
		{
		  case CMD_RUNLOOP:
			break;

		  case CMD_SHUTDOWN:
			iRunThread = 0;
			break;
	
		  default:
			err("GPSINSThread Unknown Message");
			break;
		}

#ifdef USE_RTAI
		/* wait until its time to run */
		rt_task_wait_period();
#endif
//		fprintf(stdout,".\n");
		fprintf(stdout,"%s: GPSINS Counter: %d\n",MODULE_NAME,iINSThreadCounter);
		
		/* RUN PERIODIC TASKS HERE */

		/* get INS update from store */
		
		T_ECEF2NED(Pos_LLH[0],Pos_LLH[1],(double **)Tecef2ned);

		/* get GPS measurements */
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
	  			// use current range data for now
	  			memcpy(&GPSRangeData.usPRN[iNumberValidObservations],
				       &shm->CurrentRANGE.usPRN[iIndex],sizeof(unsigned short));
				memcpy(&GPSRangeData.dPseudorange[iNumberValidObservations],
				       &shm->CurrentRANGE.dPseudorange[iIndex], sizeof(double));
			
				iNumberValidObservations++;
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
			log("GPSINS Less than 4 Measurements.... ");
			if(iRunThread == 1)
			{
				log("GPSINS Continuing...");
				continue;
			} else
			{
				log("GPSINS Exiting..");
				break;
			}
		}

		fprintf(stdout,"Allocating Matrices\n");
			
		EKF_z = gsl_vector_alloc(2*iNumberMeasurements);
		EKF_deltaz = gsl_vector_alloc(2*iNumberMeasurements);
		EKF_R = gsl_matrix_alloc(2*iNumberMeasurements,2*iNumberMeasurements);
		EKF_H = gsl_matrix_alloc(2*iNumberMeasurements,iNumberStates);
		tempV = gsl_matrix_alloc(2*iNumberMeasurements,iNumberStates);
		EKF_V = gsl_matrix_alloc(2*iNumberMeasurements,2*iNumberMeasurements);
		EKF_Vinv = gsl_matrix_alloc(2*iNumberMeasurements,2*iNumberMeasurements);
		EKF_K = gsl_matrix_alloc(iNumberStates,2*iNumberMeasurements);
		tempK = gsl_matrix_alloc(iNumberStates,2*iNumberMeasurements);
		tempKRKT = gsl_matrix_alloc(iNumberStates,2*iNumberMeasurements);
		KRKT = gsl_matrix_alloc(iNumberStates,iNumberStates);
		
		SVDWork = gsl_vector_alloc(2*iNumberMeasurements);	
		SVD_S = gsl_vector_alloc(2*iNumberMeasurements);
		SVD_U = gsl_matrix_alloc(2*iNumberMeasurements,2*iNumberMeasurements);
		SVD_SinvMat = gsl_matrix_alloc(2*iNumberMeasurements,2*iNumberMeasurements);
		tempNN = gsl_matrix_alloc(2*iNumberMeasurements,2*iNumberMeasurements);
		SVD_V = gsl_matrix_alloc(2*iNumberMeasurements,2*iNumberMeasurements);
		SVD_Sinv = gsl_vector_alloc(2*iNumberMeasurements);
		
		/* now that all our memory is allocated, disable paging */
//		mlockall(MCL_CURRENT | MCL_FUTURE);

//		user_x = gsl_vector_get(EKF_x_hat,0);
//		user_y = gsl_vector_get(EKF_x_hat,1);
//		user_z = gsl_vector_get(EKF_x_hat,2);
//		user_t = gsl_vector_get(EKF_x_hat,15);
		
//		UserPos[0] = user_x;
//		UserPos[1] = user_y;
//		UserPos[2] = user_z;
//		UserPos[3] = user_t;

		WGS84_LLH2ECEF(Pos_LLH, UserPos);
		user_x = UserPos[0];
		user_y = UserPos[1];
		user_z = UserPos[2];
		user_t = UserClock[0];
		
//		user_xdot = gsl_vector_get(EKF_x_hat,3);
//		user_ydot = gsl_vector_get(EKF_x_hat,4);
//		user_zdot = gsl_vector_get(EKF_x_hat,5);
//		user_tdot = gsl_vector_get(EKF_x_hat,16);
		
//		UserVel[0] = user_xdot;
//		UserVel[1] = user_ydot;
//		UserVel[2] = user_zdot;
//		UserVel[3] = user_tdot;

		T_NED2ECEF(Pos_LLH[0], Pos_LLH[1], (double **)Tned2ecef);

		MatMul331(Tned2ecef, Vel_NED, UserVel);
		user_xdot = UserVel[0];
		user_ydot = UserVel[1];
		user_zdot = UserVel[2];
		user_tdot = UserClock[1];

		ION_Alpha[0] = shm->CurrentIONUTC.a0;
		ION_Alpha[1] = shm->CurrentIONUTC.a1;
		ION_Alpha[2] = shm->CurrentIONUTC.a2;
		ION_Alpha[3] = shm->CurrentIONUTC.a3;
			
		ION_Beta[0] = shm->CurrentIONUTC.b0;
		ION_Beta[1] = shm->CurrentIONUTC.b1;
		ION_Beta[2] = shm->CurrentIONUTC.b2;
		ION_Beta[3] = shm->CurrentIONUTC.b3;
			
			
		
		//	fprintf(stdout,"User: [%+10.3f\t%+10.3f\t%+10.3f\t%+10.3f]\n",user_x,user_y,user_z,user_t);
		fprintf(stdout,"Organising Measurements\n");
			
		for (iIndex=0;iIndex<iNumberMeasurements;iIndex++)
		{
	
			dTransmissionTime = (double)gtObservationTime.fGPSSecondOfWeek 
				- GPSRangeData.dPseudorange[iIndex]/GPS_SPEEDOFLIGHT;
			
			GPSOrbitPropagator(dTransmissionTime, 
					GPSRangeData.usPRN[iIndex],
					&GPSEphemData[GPSRangeData.usPRN[iIndex]], 
					7500.0, 
					SVPos_calc);
			
			GPSOrbitPropagatorVelocities(dTransmissionTime, 
						GPSRangeData.usPRN[iIndex],
						&GPSEphemData[GPSRangeData.usPRN[iIndex]], 
						7500.0, 
						SVVel_calc,
						SVAcc_calc);
	
				
			usPRN = GPSRangeData.usPRN[iIndex];
				
			sat_x = SVPos_calc[0];
			sat_y = SVPos_calc[1];
			sat_z = SVPos_calc[2];
			sat_t = SVPos_calc[3]*GPS_SPEEDOFLIGHT;

			sat_xdot = SVVel_calc[0];
			sat_ydot = SVVel_calc[1];
			sat_zdot = SVVel_calc[2];
			sat_tdot = SVVel_calc[3];
			
				
			// find the range
			gsl_vector_set(temp1,0,sat_x - user_x);
			gsl_vector_set(temp1,1,sat_y - user_y);
			gsl_vector_set(temp1,2,sat_z - user_z);
				
			range = gsl_blas_dnrm2(temp1);
			
			/* calculate the range rate (relative velocity) to the SV */
			geo_vel_to_sat =  (sat_xdot - user_xdot)*(sat_x-user_x) 
				  	+ (sat_ydot - user_ydot)*(sat_y-user_y) 
					+ (sat_zdot - user_zdot)*(sat_z-user_z);

			rangerate = geo_vel_to_sat / range;
	
			EKF_H_E[iIndex][0] = -(sat_x - user_x)/range;
			EKF_H_E[iIndex][1] = -(sat_y - user_y)/range;
			EKF_H_E[iIndex][2] = -(sat_z - user_z)/range;
			
//			gsl_matrix_set(EKF_H,iIndex,0,-(sat_x - user_x)/range);
//			gsl_matrix_set(EKF_H,iIndex,1,-(sat_y - user_y)/range);
//			gsl_matrix_set(EKF_H,iIndex,2,-(sat_z - user_z)/range);
//			gsl_matrix_set(EKF_H,iIndex,15,1.0);

//			gsl_matrix_set(EKF_H,iIndex+iNumberMeasurements,3,-(sat_x - user_x)/range);
//			gsl_matrix_set(EKF_H,iIndex+iNumberMeasurements,4,-(sat_y - user_y)/range);
//			gsl_matrix_set(EKF_H,iIndex+iNumberMeasurements,5,-(sat_z - user_z)/range);
//			gsl_matrix_set(EKF_H,iIndex+iNumberMeasurements,16,1.0);
	
			gsl_vector_set(EKF_z,iIndex,GPSRangeData.dPseudorange[iIndex]);
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

			TropoDelay[iIndex] = TropoDelayModel(SV_Elevation[iIndex],Pos_LLH[3]);

			/* set PR measurement */
			gsl_vector_set(EKF_deltaz,iIndex,
				GPSRangeData.dPseudorange[iIndex] - range 
				- (user_t-sat_t)
				+ delta_pr_omegaedot
				- IonoDelay[iIndex]
				- TropoDelay[iIndex]);

//			fprintf(stdout,"%f %f %f %f %f %f\n",
//				range,
//				user_t,
//				sat_t,
//				delta_pr_omegaedot,
//				IonoDelay[iIndex],
//				TropoDelay[iIndex]);			

			/* set doppler measurement */
			gsl_vector_set(EKF_deltaz,iIndex+iNumberMeasurements,
				-(double)GPSRangeData.fDoppler[iIndex]*GPS_L1_Wavelength
				- rangerate
				- (user_tdot));

			
						
				
		}  // end for i=1:NumberMeasurements

		// rotate H matrix to local
		MatMul333(Tecef2ned,EKF_H_E,EKF_H_LTP);
		
		for(iIndex=0;iIndex<iNumberMeasurements;iIndex++)
		{
			gsl_matrix_set(EKF_H,iIndex,0,EKF_H_LTP[iIndex][0]);
			gsl_matrix_set(EKF_H,iIndex,1,EKF_H_LTP[iIndex][1]);
			gsl_matrix_set(EKF_H,iIndex,2,EKF_H_LTP[iIndex][2]);
			gsl_matrix_set(EKF_H,iIndex,15,1.0);
			
			gsl_matrix_set(EKF_H,iIndex+iNumberMeasurements,3,EKF_H_LTP[iIndex][0]);
			gsl_matrix_set(EKF_H,iIndex+iNumberMeasurements,4,EKF_H_LTP[iIndex][1]);
			gsl_matrix_set(EKF_H,iIndex+iNumberMeasurements,5,EKF_H_LTP[iIndex][2]);
			gsl_matrix_set(EKF_H,iIndex+iNumberMeasurements,16,1.0);
			
			
		}	

//		gsl_vector_fprintf(stdout,EKF_deltaz,"[%f]");
//		gsl_matrix_fprintf(stdout,EKF_H,"%f");
		
		fprintf(stdout,"Setting up PHI matrix\n");
		WGS84_calcRnRe(Pos_LLH[0],  &Rm, &Rp);
		Rmh = Rm+Pos_LLH[3];
		Rph = Rp+Pos_LLH[3];
		Rg = sqrt(Rmh*Rmh + Rph*Rph);
		
		/* kalman filter */
		gsl_matrix_set_zero(EKF_F);
		gsl_matrix_set(EKF_F,0,3,1.0/Rmh);	
		gsl_matrix_set(EKF_F,1,4,1.0/(Rph*cos(Pos_LLH[0])));
		gsl_matrix_set(EKF_F,2,5,-1.0);
		gsl_matrix_set(EKF_F,0,2,-Vel_NED[0]/(Rg*Rg));
		gsl_matrix_set(EKF_F,1,0,Vel_NED[1]*sin(Pos_LLH[0])/(Rg*cos(Pos_LLH[0])*cos(Pos_LLH[0])));
		gsl_matrix_set(EKF_F,1,2,-Vel_NED[1]/(Rg*Rg*cos(Pos_LLH[0])));
		gsl_matrix_set(EKF_F,1,4,1.0/(Rg*cos(Pos_LLH[0])));
		gsl_matrix_set(EKF_F,3,0,-Vel_NED[1]*2.0*OMEGA_e*cos(Pos_LLH[0]) - Vel_NED[1]*Vel_NED[1]/(Rph*cos(Pos_LLH[0])*cos(Pos_LLH[0])));
		gsl_matrix_set(EKF_F,3,2,(Vel_NED[1]*Vel_NED[1] * tan(Pos_LLH[0]) - Vel_NED[0]*Vel_NED[2])/(Rg*Rg));
		gsl_matrix_set(EKF_F,3,3,Vel_NED[2]/Rg);
		gsl_matrix_set(EKF_F,3,4,-2.0*(OMEGA_e*sin(Pos_LLH[0]) + Vel_NED[1]*tan(Pos_LLH[0])/Rg));
		gsl_matrix_set(EKF_F,3,5,Vel_NED[0]/Rg);
		gsl_matrix_set(EKF_F,4,0,2.0*OMEGA_e*(Vel_NED[0]*cos(Pos_LLH[0]) - Vel_NED[2]*sin(Pos_LLH[0])) + 
								(Vel_NED[0]*Vel_NED[1]/(Rg*cos(Pos_LLH[0])*cos(Pos_LLH[0]))));
		gsl_matrix_set(EKF_F,4,2,-Vel_NED[1]*(Vel_NED[0]*tan(Pos_LLH[0]) - Vel_NED[1]*Vel_NED[2])/(Rg*Rg));
		gsl_matrix_set(EKF_F,4,3,2.0*OMEGA_e*sin(Pos_LLH[0]) + Vel_NED[1]*tan(Pos_LLH[0])/Rg);
		gsl_matrix_set(EKF_F,4,4,(Vel_NED[0]*tan(Pos_LLH[0])+Vel_NED[2])/Rg);
		gsl_matrix_set(EKF_F,4,5,2.0*OMEGA_e * cos(Pos_LLH[0]) + Vel_NED[1]/Rg);

		// velocity error due to tilt error
		gsl_matrix_set(EKF_F,3,7,-(Acc_NED[2] - LOCAL_GRAVITY));  // -fd
		gsl_matrix_set(EKF_F,3,8,(Acc_NED[1])); // fe
		gsl_matrix_set(EKF_F,4,6,(Acc_NED[2] - LOCAL_GRAVITY));  // fd
		gsl_matrix_set(EKF_F,4,8,-(Acc_NED[0])); // -fn
		gsl_matrix_set(EKF_F,5,6,-(Acc_NED[1])); // -fe
		gsl_matrix_set(EKF_F,5,7,(Acc_NED[0])); // fn

		gsl_matrix_set(EKF_F,5,0,2.0*OMEGA_e*Vel_NED[1]*sin(Pos_LLH[0]));
		gsl_matrix_set(EKF_F,5,2,(Vel_NED[0]*Vel_NED[0] + Vel_NED[1]*Vel_NED[1])/(Rg*Rg));  // note: gravity correction term missing
		gsl_matrix_set(EKF_F,5,3,-2.0*Vel_NED[0]/Rg);
		gsl_matrix_set(EKF_F,5,4,-2.0*(OMEGA_e * cos(Pos_LLH[0]) + Vel_NED[1]/Rg));
		gsl_matrix_set(EKF_F,6,0,-OMEGA_e * sin(Pos_LLH[0]));
		gsl_matrix_set(EKF_F,6,2,-Vel_NED[1]/(Rg*Rg));
		gsl_matrix_set(EKF_F,7,2,Vel_NED[0] / (Rg*Rg));
		gsl_matrix_set(EKF_F,8,0,-OMEGA_e * cos(Pos_LLH[0]) - Vel_NED[1]/(Rg*cos(Pos_LLH[0])*cos(Pos_LLH[0])));
		gsl_matrix_set(EKF_F,8,2,Vel_NED[1]*tan(Pos_LLH[0])/(Rg*Rg));

		// tilt error due to velocity error
		gsl_matrix_set(EKF_F,6,4,1.0/Rg);
		gsl_matrix_set(EKF_F,7,3,-1.0/Rg);
		gsl_matrix_set(EKF_F,8,4,-tan(Pos_LLH[0])/Rg);


		// tilt error due to inertial rotation of the reference frame
		gsl_matrix_set(EKF_F,6,7,-OMEGA_e*sin(Pos_LLH[0]) - Vel_NED[1]*tan(Pos_LLH[0])/Rg);
		gsl_matrix_set(EKF_F,6,8,Vel_NED[0]/Rg);
		gsl_matrix_set(EKF_F,7,6,OMEGA_e*sin(Pos_LLH[0]) + Vel_NED[1]*tan(Pos_LLH[0])/Rg);
		gsl_matrix_set(EKF_F,8,6,-Vel_NED[0]/Rg);
		gsl_matrix_set(EKF_F,8,7,-OMEGA_e*cos(Pos_LLH[0]) - Vel_NED[1]/Rg);
		gsl_matrix_set(EKF_F,7,8,OMEGA_e*cos(Pos_LLH[0]) + Vel_NED[1]/Rg);


		for(iRowIndex=0;iRowIndex<3;iRowIndex++)
		{
			for(iColIndex=0;iColIndex<3;iColIndex++)
			{
				gsl_matrix_set(EKF_F,6+iRowIndex,9+iColIndex,-C_BN[iRowIndex][iColIndex]);
				gsl_matrix_set(EKF_F,3+iRowIndex,12+iColIndex,C_BN[iRowIndex][iColIndex]);
			}
		}

		gsl_matrix_set(EKF_F,9,9, -x_gyro_beta);
		gsl_matrix_set(EKF_F,10,10, -y_gyro_beta);
		gsl_matrix_set(EKF_F,11,11, -z_gyro_beta);
		gsl_matrix_set(EKF_F,12,12, -x_accel_beta);
		gsl_matrix_set(EKF_F,13,13, -y_accel_beta);
		gsl_matrix_set(EKF_F,14,14, -z_accel_beta);
		

		/*  PHI = I + F * dt */
		gsl_matrix_set_identity(EKF_PHI);
		gsl_matrix_scale (EKF_F, 1.0/GPS_RATE);
		gsl_matrix_add (EKF_PHI, EKF_F);

//		gsl_matrix_fprintf(stdout,EKF_PHI,"%f");

		fprintf(stdout,"Q-Matrix\n");
		/* Continuous Q-Matrix */
		gsl_matrix_set_identity(EKF_Q);
		
		/* Discrete Qd = PHI * G * Q * G' * PHI' */
		// TODO: Qd

		fprintf(stdout,"R-Matrix\n");
		
		/* R - Matrix */
		gsl_matrix_set_identity(EKF_R);
		
		for(iIndex=0;iIndex<iNumberMeasurements;iIndex++)
		{
			gsl_matrix_set(EKF_R,iIndex,iIndex,3.0);
			gsl_matrix_set(EKF_R,iIndex+iNumberMeasurements,iIndex+iNumberMeasurements,0.3);
		}

		/* prediction step */
		/* x- = PHI * x+ */
		/* z- = H * x- */
		/* P- = PHI * P+ * PHI' + Q */
		gsl_matrix_memcpy(EKF_P_minus,EKF_Q);
		gsl_matrix_set_zero(tempP);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,EKF_PHI,EKF_P,1.0, tempP);
		gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,tempP,EKF_PHI,1.0,EKF_P_minus);
			
		//V = H * P_minus * H' + R;
		fprintf(stdout,"HPH'+R\n");
		
		gsl_matrix_set_zero (tempV);
		gsl_matrix_set_zero(EKF_V);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, EKF_H, EKF_P_minus, 1.0, tempV);
		gsl_matrix_memcpy(EKF_V,EKF_R); 
		gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0,tempV, EKF_H, 1.0,EKF_V);	 //this adds V = V + A*B
	
	
		/* calculate inv(EKF_V) */
		/*
		  
		*/
		fprintf(stdout,"inv(V)\n");
		
		iReturnValue =  gsl_linalg_SV_decomp (EKF_V,SVD_V,SVD_S,SVDWork);
		gsl_matrix_memcpy(SVD_U,EKF_V);
		
		// get element wise inverse of S
		gsl_vector_set_all(SVD_Sinv,1.0);
		gsl_vector_div(SVD_Sinv, SVD_S);
	
		gsl_matrix_set_zero (SVD_SinvMat);
		for(iIndex=0;iIndex<2*iNumberMeasurements;iIndex++)
		{
			gsl_matrix_set(SVD_SinvMat,iIndex,iIndex,gsl_vector_get(SVD_Sinv,iIndex));
		}

		// inv(S)*UT
		gsl_matrix_set_zero (tempNN);
		gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, SVD_SinvMat, SVD_U, 1.0, tempNN);
		
		// A = V * inv(S) * UT
		gsl_matrix_set_zero(EKF_Vinv);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, SVD_V, tempNN, 1.0, EKF_Vinv);
	
		fprintf(stdout,"PH'*inv(V)\n");
		
		//K = P_minus * H' * inv(V);		
		gsl_matrix_set_zero (tempK);
		gsl_matrix_set_zero(EKF_K);
		gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, EKF_P_minus, EKF_H, 1.0, tempK);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0,tempK, EKF_Vinv, 1.0,EKF_K);	
	
	
		// x_hat_out = K * (z);
		gsl_vector_set_zero(EKF_x_hat_out);
		
		fprintf(stdout,"Kz\n");
		gsl_blas_dgemv (CblasNoTrans, 1.0, EKF_K, EKF_deltaz, 1.0, EKF_x_hat_out);
	
		// %symmetric form of the P-update equation to ward off divergence p347
		//%Brown and Hwang
		// P_out = (eye(NumberStates) - K*H)*P_minus*(eye(NumberStates)-K*H)' + K*R*K';
	
	
		//K*H
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0,EKF_K, EKF_H, 1.0,KH);	
	
	
		//(eye(NumberStates) - K*H)
		gsl_matrix_set_identity(I_17by17); 
		gsl_matrix_memcpy(IminusKH,I_17by17);
		gsl_matrix_sub(IminusKH, KH);
	
		//(eye(NumberStates) - K*H)*P_minus
		gsl_matrix_set_zero(KHPminus);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, IminusKH, EKF_P_minus, 1.0, KHPminus);
	
	
		//(eye(NumberStates) - K*H)*P_minus*(eye(NumberStates)-K*H)'
		gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, KHPminus, IminusKH, 1.0, EKF_P);
	
		//K*R*K'
		gsl_matrix_set_zero (tempKRKT);
		gsl_matrix_set_zero(KRKT);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, EKF_K, EKF_R, 1.0, tempKRKT);
	
		gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0,tempKRKT, EKF_K, 1.0,KRKT);	
		gsl_matrix_sub(EKF_P, KRKT);

		// apply correction
		Pos_LLH[0] -= gsl_vector_get(EKF_x_hat_out,0);
		Pos_LLH[1] -= gsl_vector_get(EKF_x_hat_out,1);
		Pos_LLH[2] -= gsl_vector_get(EKF_x_hat_out,2);
		Vel_NED[0] -= gsl_vector_get(EKF_x_hat_out,3);
		Vel_NED[1] -= gsl_vector_get(EKF_x_hat_out,4);
		Vel_NED[2] -= gsl_vector_get(EKF_x_hat_out,5);
		
		AccBias[0] -= gsl_vector_get(EKF_x_hat_out,9);
		AccBias[1] -= gsl_vector_get(EKF_x_hat_out,10);
		AccBias[2] -= gsl_vector_get(EKF_x_hat_out,11);
		
		GyroBias[0] -= gsl_vector_get(EKF_x_hat_out,12);
		GyroBias[1] -= gsl_vector_get(EKF_x_hat_out,13);
		GyroBias[2] -= gsl_vector_get(EKF_x_hat_out,14);
		
		// attitude
		del_alpha = -gsl_vector_get(EKF_x_hat_out,6);
		del_beta = -gsl_vector_get(EKF_x_hat_out,7);
		del_gamma = -gsl_vector_get(EKF_x_hat_out,8);
		
		del_att_skew[0][0] = 1.0;
		del_att_skew[0][1] = -del_gamma;
		del_att_skew[0][2] = del_beta;
		del_att_skew[1][0] = del_gamma;
		del_att_skew[1][1] = 1.0;
		del_att_skew[1][2] = -del_alpha;
		del_att_skew[2][0] = -del_beta;
		del_att_skew[2][1] = del_alpha;
		del_att_skew[2][2] = 1.0;
		
		/* C_BN = (eye(3,3) + del_att_skew) * C_BN */
		MatMul333(del_att_skew,C_BN,C_BN_new);
		memcpy(C_BN_new,C_BN,sizeof(C_BN));
		
		// clock
		UserClock[0] -= gsl_vector_get(EKF_x_hat_out,15);
		UserClock[1] -= gsl_vector_get(EKF_x_hat_out,16);
		
		/* re-enable paging */
//		munlockall();

		log("Result: \n");
		fprintf(stdout,"Pos: %f %f %f\n",
			Pos_LLH[0]*RAD2DEG,
			Pos_LLH[1]*RAD2DEG,
			Pos_LLH[2]);
		fprintf(stdout,"Vel: %f %f %f\n",
			Vel_NED[0],
			Vel_NED[1],
			Vel_NED[2]);
		fprintf(stdout,"Clock: %f %f\n",
			UserClock[0],UserClock[1]);	

		fprintf(stdout,"X-out\n");
		gsl_vector_fprintf(stdout,EKF_x_hat_out,"[%f]");
		
//		fprintf(stdout,"delta Z\n");
//		gsl_vector_fprintf(stdout,EKF_deltaz,"[%f]");

		log("Releasing Matrices");

		gsl_vector_free(EKF_z);
		gsl_vector_free(EKF_deltaz);
		gsl_matrix_free(EKF_H);
		gsl_matrix_free(EKF_R);
		gsl_matrix_free(tempV);
		gsl_matrix_free(EKF_V);
		gsl_matrix_free(tempNN);
		gsl_vector_free(SVD_S);
		gsl_vector_free(SVD_Sinv);
		gsl_matrix_free(SVD_SinvMat);
		gsl_matrix_free(SVD_V);
		gsl_vector_free(SVDWork);
		gsl_matrix_free(EKF_K);
		gsl_matrix_free(tempK);
		gsl_matrix_free(KRKT);
		gsl_matrix_free(tempKRKT);


	}

	/* exit */
	fprintf(stdout,"%s GPS-INS Thread Exiting\n",MODULE_NAME);

	gsl_matrix_free(EKF_PHI);
	gsl_matrix_free(EKF_F);
	gsl_matrix_free(EKF_P_minus);
	gsl_matrix_free(tempP);
	gsl_matrix_free(EKF_P);
	gsl_matrix_free(I_17by17);	
	gsl_vector_free(EKF_x_hat);
	gsl_matrix_free(EKF_Q);
	gsl_vector_free(temp1);
	gsl_matrix_free(KH);
	gsl_matrix_free(IminusKH);
	gsl_matrix_free(KHPminus);
	
#ifdef USE_RTAI
	rt_make_soft_real_time();

	if(rttGPSINS)
	{
		rt_task_delete(rttGPSINS);
	}
	stop_rt_timer();
#endif

	pthread_exit(NULL);
}



