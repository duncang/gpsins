/**
 *	\file gpsstore.c - Retrieves measurements from a Novatel GPS receiver and stores in shared memory.
 *
 * 	
 *
 *
 *	Revision Information: $Id: gpsstore.c 340 2007-10-08 11:47:09Z greerd $
 */

#include <signal.h>
#include <stdio.h>
#include <termios.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <strings.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <sys/shm.h>
#include "gpsins.h"
#include "novatel.h"

#define MODULE_NAME	"[GPS Store]"

#define log(string)	fprintf(stdout,"%s %s\n", MODULE_NAME, string)
#define err(string)	fprintf(stderr,"%s ERROR: %s\n",MODULE_NAME, string)

#define MAX_LOG_LENGTH 1024

/** max number of GPS satellites that we can store ephemeris for*/
#define MAX_PRNS	40	

/* example bestpos packet (ascii)
#BESTPOSA,COM2,0,28.0,FINESTEERING,1386,97129.000,00000020,4ca6,1810;SOL_COMPUTED,SINGLE,-27.47735832601,153.02715301736,49.0198,40.7704WGS84,1.6805,1.3796,3.8551,"",0.000,0.000,10,10,0,0,0,0,0,0*aa287bfb
*/

/** Controls whether or not the program should keep running. */
int iShutdown = 0;

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
 *	\fn SendString
 * 	\brief Send string to serial port
 *
 * 	\param iPort The port handler; must be open and initialised
 *	\param string The string to send (null terminated)
 */
int SendString(int iPort,char *string);



/**
 *	
 */
int main(int argc, char *argv[])
{
	char *azcPort = "/dev/ttyS0";
	char *azcBaud = "230400";
	int iBaudRate = 230400;
	int iPort;
	int iBytesRead = 0;
	unsigned char azucBuffer[MAX_LOG_LENGTH];

	unsigned char ucHeaderLength = 0;
	int iTotalBytes = 0;
	unsigned short usMessageLength = 0;
	unsigned short usMessageID = 0;

	unsigned long ulCRC_Calc, ulCRC_Tx;

	/** controls if shared memory is open for use: 0 - No, 1 - Yes */
	int iUseSharedMemory = 0;

	struct termios tPortSettings;

	/** structure to store positions */
	BESTPOS	BestPosData;
	RANGE	RangeData;
	GPSEPHEM	GPSEphemeris[MAX_PRNS];

	int iAbandonPacket = 0;

	long lMeasurementNumber = 0;

	GPSTime GPSTimeStamp;
	unsigned long ulGPSMilliSeconds;
	unsigned long ulPRN;

	/* shared memory objects */
	int shmGPSId;
	key_t shmGPSKey;

	/* GPS Data reference */
	GPSData *ptGPSData;

	/* setup signal handlers */
	signal(SIGHUP,  HandleSignal);
	signal(SIGKILL, HandleSignal);
	signal(SIGTERM, HandleSignal);
	signal(SIGALRM, HandleSignal);
	signal(SIGINT, HandleSignal);

	/* read input parameters */
	switch(argc)
	{
		case 1:
			/* defaults */
			break;
		case 2:
			/* set port */
			azcPort = argv[1];
			break;
		case 3:
			/* set port and baud */
			azcPort = argv[1];
			azcBaud = argv[2];

			switch(atoi(azcBaud))
			{
				case 9600: 	iBaudRate = B9600; break;
				case 19200:	iBaudRate = B19200; break;
				case 38400:	iBaudRate = B38400; break;
				case 57600:	iBaudRate = B57600; break;
				case 115200:	iBaudRate = B115200; break;
				case 230400:	iBaudRate = B230400; break;
				default:	fprintf(stderr,"Unknown Baud - setting default B230400!\n"); iBaudRate = B230400; break;
			}
			break;

		default:
			fprintf(stderr,"%s Unknown Argument List... Exiting\n",MODULE_NAME);
			exit(-1);
			break;
	}

	fprintf(stdout,"%s Using Port: %s at %s\n",MODULE_NAME, azcPort, azcBaud);

	/* open port */
	iPort = open(azcPort, O_RDWR);

	if(iPort < 0)
	{
		fprintf(stderr,"Cannot open port, %s:",azcPort);
		perror("open()");
		exit(-1);
	}

	/* initialise port settings */
	bzero(&tPortSettings,sizeof(tPortSettings));

	/* configure port settings */
	tPortSettings.c_iflag &= ~(IGNPAR|BRKINT|PARMRK|ISTRIP|INLCR|IGNCR|ICRNL|IXON);
	/*tPortSettings.c_iflag |= (IXON|IXOFF);  // enable software flow control */
	tPortSettings.c_cflag &= ~(CSIZE|PARENB);
	/*tPortSettings.c_cflag |= (B115200|CS8|CLOCAL|CREAD); */
	tPortSettings.c_cflag |= (CS8|CLOCAL|CREAD);

	/* Set the Baud Rate (same for input and output) */
	cfsetispeed(&tPortSettings, iBaudRate);
	cfsetospeed(&tPortSettings, iBaudRate);

	tPortSettings.c_oflag &= ~(OPOST);
	tPortSettings.c_lflag &= ~(ECHO|ECHONL|ICANON|ISIG|IEXTEN);
	tPortSettings.c_cc[VTIME] = 10; /* wait for 1 second before returning */
	tPortSettings.c_cc[VMIN] = 0;   /* can return with no data */
	tcflush(iPort,TCIFLUSH);
	tcsetattr(iPort,TCSANOW,&tPortSettings);

	
	/* connect to shared memory for GPS*/
	shmGPSKey = GPSDATA_KEY;
	if((shmGPSId = shmget(shmGPSKey, sizeof(GPSData), 0666)) < 0) 
	{
        	perror("shmget(): cannot find GPS data storage - no shared mem used");
       		iUseSharedMemory = 0;
    	} else
	{
		/* shared memory is intiialised */
		iUseSharedMemory = 1;
	}

	/* attach */
	if(iUseSharedMemory == 1)
	{
		if((ptGPSData = shmat(shmGPSId, NULL, 0)) == (GPSData *) -1) 
		{
			perror("shmat(): attaching GPS Data Storage");
			exit(-1);
		}
	} 

	log("Initialised GPS Store - Hit CTRL-C to exit");
	
	while(iShutdown == 0)
	{
		
		/* clear buffer */
		memset(azucBuffer,0,MAX_LOG_LENGTH);
		iAbandonPacket = 0;

		/* synchronoise packet - look for binary header - 0xAA 0x44 0x12 */
		do
		{
			azucBuffer[0] = azucBuffer[1];
			azucBuffer[1] = azucBuffer[2];
			iBytesRead = read(iPort,&azucBuffer[2],1);
			
			if(iShutdown == 1)
			{	
				iAbandonPacket=1;
				break;
			}
				
		} while(azucBuffer[0] != 0xAA && azucBuffer[1] != 0x44 && azucBuffer[2] != 0x12);
		
		if(iAbandonPacket == 1)
		{
			continue;
		}

		/* read header length */
		iBytesRead = read(iPort,&azucBuffer[3],1);
		if(iBytesRead == 1)
		{
			ucHeaderLength = azucBuffer[3];
			
		} else
		{
			log("Could not read header length");
			break;
		}

		/* sanity check on header length */
		if(ucHeaderLength > 28)
		{
			/* too large - try again */
			err("Header unexpectedly large");
			iAbandonPacket = 1;
			continue;
		}

		/* we have already received the first 4 bytes */
		iTotalBytes = 4;
	
		/* read rest of header */
		while(iTotalBytes<ucHeaderLength)
		{	
			iBytesRead = read(iPort,&azucBuffer[iTotalBytes],MAX_LOG_LENGTH-2);
			iTotalBytes += iBytesRead;
			if(iTotalBytes > MAX_LOG_LENGTH)
			{
				err("Too many bytes received");
				iAbandonPacket = 1;
				break;
			}
		}

		if(iAbandonPacket == 1)
		{
			continue;
		}

		/*  read message length */
		memcpy(&usMessageLength,&azucBuffer[8],2); 

		/* receive rest of message */
		while(iTotalBytes<(ucHeaderLength + usMessageLength + 4)) /* 4 is the length of the CRC */
		{
			iBytesRead = read(iPort,&azucBuffer[iTotalBytes],MAX_LOG_LENGTH-2);
			iTotalBytes += iBytesRead;
			if(iTotalBytes > MAX_LOG_LENGTH)
			{
				err("Too many bytes received");
				iAbandonPacket = 1;
				break;
			}
		}

		if(iAbandonPacket == 1)
		{
			continue;
		}

		/* get the message ID */
		memcpy(&usMessageID,&azucBuffer[4],2);

		/* debug message */
		fprintf(stdout,"Message type: %d received, %d bytes total\n",usMessageID, iTotalBytes);

		/* perform checksum */
		ulCRC_Calc = CalculateBlockCRC32(ucHeaderLength+usMessageLength,azucBuffer);
		
		/* retrieve checksum from message */
		memcpy(&ulCRC_Tx,&azucBuffer[ucHeaderLength+usMessageLength],4);

		if(ulCRC_Tx != ulCRC_Calc)
		{
			/* we have an error */
			log("CRC Check Error");
		} else
		{	
			/* get the time stamp */
			memcpy(&GPSTimeStamp.usGPSWeek,&azucBuffer[14],2);
			memcpy(&ulGPSMilliSeconds,&azucBuffer[16],4);
			GPSTimeStamp.fGPSSecondOfWeek = ((float)ulGPSMilliSeconds) / 1000.0;

			/* extract the message data */
			switch(usMessageID)
			{
				case NOVATEL_BESTPOS:
					/* copy timestamp */
					memcpy(&BestPosData.gtTimeStamp,&GPSTimeStamp,sizeof(GPSTime));	
					
					/*extract data */

					memcpy(&BestPosData.dLatitude,&azucBuffer[ucHeaderLength+8],8);
					memcpy(&BestPosData.dLongitude,&azucBuffer[ucHeaderLength+16],8);
					memcpy(&BestPosData.dHeight,&azucBuffer[ucHeaderLength+24],8);
					memcpy(&BestPosData.fUndulation,&azucBuffer[ucHeaderLength+32],4);
					memcpy(&BestPosData.ulDatumID,&azucBuffer[ucHeaderLength+36],4);
					memcpy(&BestPosData.fLatitudeSigma,&azucBuffer[ucHeaderLength+40],4);
					memcpy(&BestPosData.fLongitudeSigma,&azucBuffer[ucHeaderLength+44],4);
					memcpy(&BestPosData.fHeightSigma,&azucBuffer[ucHeaderLength+48],4);
					memcpy(&BestPosData.azucBaseStationID[0],&azucBuffer[ucHeaderLength+52],4);
					memcpy(&BestPosData.fDifferentialAge,&azucBuffer[ucHeaderLength+56],4);
					memcpy(&BestPosData.fSolutionAge,&azucBuffer[ucHeaderLength+60],4);
					memcpy(&BestPosData.ucNumberObservationsTracked,&azucBuffer[ucHeaderLength+64],1);
					memcpy(&BestPosData.ucNumberL1ObservationsUsed,&azucBuffer[ucHeaderLength+65],1);
					memcpy(&BestPosData.ucNumberL1ObservationsAboveRTKMaskAngle,&azucBuffer[ucHeaderLength+66],1);
					memcpy(&BestPosData.ucNumberL2ObservationsAboveRTKMaskAngle,&azucBuffer[ucHeaderLength+67],1);
					
					fprintf(stdout,"BESTPOS: %lf %lf %lf\n",BestPosData.dLatitude,
										BestPosData.dLongitude,
										BestPosData.dHeight);
					break;

				case NOVATEL_RANGE: /* pseudorange measurements */
					/* copy timestamp */
					memcpy(&RangeData.gtTimeStamp,&GPSTimeStamp,sizeof(GPSTime));	
					
					/* get number of measurements */
					memcpy(&RangeData.lNumberObservations,&azucBuffer[ucHeaderLength],4);
					fprintf(stdout,"RANGE: %ld Measurements\n",RangeData.lNumberObservations);

					/* limit max measurements */
					if(RangeData.lNumberObservations > NOVATEL_MAXCHANNELS)
					{
						RangeData.lNumberObservations = NOVATEL_MAXCHANNELS;
					}

					for(lMeasurementNumber=0;lMeasurementNumber<RangeData.lNumberObservations;lMeasurementNumber++)
					{
						memcpy(&RangeData.usPRN[lMeasurementNumber],
							&azucBuffer[ucHeaderLength + 4 + lMeasurementNumber*44],2);
						memcpy(&RangeData.usGlonassFrequency[lMeasurementNumber],
							&azucBuffer[ucHeaderLength + 6 + lMeasurementNumber*44],2);
						memcpy(&RangeData.dPseudorange[lMeasurementNumber],
							&azucBuffer[ucHeaderLength + 8 + lMeasurementNumber*44],8);
						memcpy(&RangeData.fPseudorangeSigma[lMeasurementNumber],
							&azucBuffer[ucHeaderLength + 16 + lMeasurementNumber*44],4);
						memcpy(&RangeData.dCarrierPhase[lMeasurementNumber],
							&azucBuffer[ucHeaderLength + 20 + lMeasurementNumber*44],8);
						memcpy(&RangeData.fCarrierPhaseSigma[lMeasurementNumber],
							&azucBuffer[ucHeaderLength + 28 + lMeasurementNumber*44],4);
						memcpy(&RangeData.fDoppler[lMeasurementNumber],
							&azucBuffer[ucHeaderLength + 32 + lMeasurementNumber*44],4);
						memcpy(&RangeData.fCNo[lMeasurementNumber],
							&azucBuffer[ucHeaderLength + 36 + lMeasurementNumber*44],4);
						memcpy(&RangeData.fLockTime[lMeasurementNumber],
							&azucBuffer[ucHeaderLength + 40 + lMeasurementNumber*44],4);
						memcpy(&RangeData.ulTrackingStatus[lMeasurementNumber],
							&azucBuffer[ucHeaderLength + 44 + lMeasurementNumber*44],4);
						

					/*	fprintf(stdout,"%d %f: RANGE: PRN:%d %lf\n",
							RangeData.gtTimeStamp.usGPSWeek,
							RangeData.gtTimeStamp.fGPSSecondOfWeek,
							RangeData.usPRN[lMeasurementNumber],
							RangeData.dPseudorange[lMeasurementNumber]); 
					*/
					}
					
					break;

				case NOVATEL_GPSEPHEM:  /* satellite ephemeris */
					
					/* get the PRN */
					memcpy(&ulPRN,&azucBuffer[ucHeaderLength],4);

					/* copy timestamp */
					memcpy(&GPSEphemeris[ulPRN].gtTimeStamp,&GPSTimeStamp,sizeof(GPSTime));	
					
					/* copy ephemeris parameters */
					memcpy(&GPSEphemeris[ulPRN].ulPRN,&azucBuffer[ucHeaderLength],4);	
					memcpy(&GPSEphemeris[ulPRN].dTOW,&azucBuffer[ucHeaderLength+4],8);	
					memcpy(&GPSEphemeris[ulPRN].ulHealth,&azucBuffer[ucHeaderLength+12],4);	
					memcpy(&GPSEphemeris[ulPRN].ulIODE1,&azucBuffer[ucHeaderLength+16],4);	
					memcpy(&GPSEphemeris[ulPRN].ulIODE2,&azucBuffer[ucHeaderLength+20],4);	
					memcpy(&GPSEphemeris[ulPRN].ulGPSWeek,&azucBuffer[ucHeaderLength+24],4);	
					memcpy(&GPSEphemeris[ulPRN].ulZWeek,&azucBuffer[ucHeaderLength+28],4);	
					memcpy(&GPSEphemeris[ulPRN].dTOE,&azucBuffer[ucHeaderLength+32],8);	
					memcpy(&GPSEphemeris[ulPRN].dA,&azucBuffer[ucHeaderLength+40],8);	
					memcpy(&GPSEphemeris[ulPRN].dDeltaN,&azucBuffer[ucHeaderLength+48],8);	
					memcpy(&GPSEphemeris[ulPRN].dM0,&azucBuffer[ucHeaderLength+56],8);	
					memcpy(&GPSEphemeris[ulPRN].dEccentricity,&azucBuffer[ucHeaderLength+64],8);	
					memcpy(&GPSEphemeris[ulPRN].dOmega,&azucBuffer[ucHeaderLength+72],8);	
					memcpy(&GPSEphemeris[ulPRN].dcuc,&azucBuffer[ucHeaderLength+80],8);	
					memcpy(&GPSEphemeris[ulPRN].dcus,&azucBuffer[ucHeaderLength+88],8);	
					memcpy(&GPSEphemeris[ulPRN].dcrc,&azucBuffer[ucHeaderLength+96],8);	
					memcpy(&GPSEphemeris[ulPRN].dcrs,&azucBuffer[ucHeaderLength+104],8);	
					memcpy(&GPSEphemeris[ulPRN].dcic,&azucBuffer[ucHeaderLength+112],8);	
					memcpy(&GPSEphemeris[ulPRN].dcis,&azucBuffer[ucHeaderLength+120],8);	
					memcpy(&GPSEphemeris[ulPRN].dInclination0,&azucBuffer[ucHeaderLength+128],8);	
					memcpy(&GPSEphemeris[ulPRN].dInclination_dot,&azucBuffer[ucHeaderLength+136],8);	
					memcpy(&GPSEphemeris[ulPRN].dOmega0,&azucBuffer[ucHeaderLength+144],8);	
					memcpy(&GPSEphemeris[ulPRN].dOmega_dot,&azucBuffer[ucHeaderLength+152],8);	
					memcpy(&GPSEphemeris[ulPRN].ulIODC,&azucBuffer[ucHeaderLength+160],4);	
					memcpy(&GPSEphemeris[ulPRN].dTOC,&azucBuffer[ucHeaderLength+164],8);	
					memcpy(&GPSEphemeris[ulPRN].dTGD,&azucBuffer[ucHeaderLength+172],8);	
					memcpy(&GPSEphemeris[ulPRN].dA_f0,&azucBuffer[ucHeaderLength+180],8);	
					memcpy(&GPSEphemeris[ulPRN].dA_f1,&azucBuffer[ucHeaderLength+188],8);	
					memcpy(&GPSEphemeris[ulPRN].dA_f2,&azucBuffer[ucHeaderLength+196],8);	
					memcpy(&GPSEphemeris[ulPRN].ulAntiSpoofing,&azucBuffer[ucHeaderLength+204],4);	
					memcpy(&GPSEphemeris[ulPRN].dN,&azucBuffer[ucHeaderLength+208],8);	
					memcpy(&GPSEphemeris[ulPRN].dURA,&azucBuffer[ucHeaderLength+216],8);	
					
					fprintf(stdout,"PRN%ld Ephemeris Stored\n",GPSEphemeris[ulPRN].ulPRN);

					break;
				default:
					log("Unknown Message");
					break;		
			}
			
		}
		
	}

	/* shutdown */
	close(iPort);

	log("Shutting Down");

	/* exit */
	return 0;
}  /* main */



int SendString(int iPort, char *string)
{
	char cWriteData[100];
	int iBytes = 0;

	bzero(&cWriteData,sizeof(cWriteData));
	sprintf(cWriteData,string);
	
	fprintf(stdout,"%s Tx: %s\n",MODULE_NAME, string);

	iBytes = write(iPort,string,strlen(string));

	if(iBytes < 0)
	{
		perror("write()");
	}
	/* send CRLF */
	write(iPort,"\r\n",2);

	return iBytes;
}


