/**
 *	\file gras.h GRAS related data
 *	\author Duncan Greer
 *  $Id: gras.h 362 2007-10-22 03:43:07Z greerd $
 *
 */

#ifndef _GRAS_H
#define _GRAS_H


#define GRAS_MAXMEASUREMENTS 18

typedef struct 
{
	unsigned char ucRangingSourceId;
	unsigned char ucIssueOfData;
	float fPseudorangeCorrection;
	float fRangeRateCorrection;
	float fSigmaPR_GND;
	float fB1;
	float fB2;
	float fB3;
	float fB4;

}GRAS_MEASUREMENT_BLOCK;

typedef struct _GBAS_TYPE101
{
	float fModifiedZCount;
	unsigned char ucNumberMeasurements;
	unsigned char ucMessageType;
	float	fEphemerisDecorrelation;
	unsigned short usEphemerisCRC;
	float fSourceAvailabilityDuration;
	unsigned char ucNumberBParamaters;

	GRAS_MEASUREMENT_BLOCK GRASMeasurementBlock[GRAS_MAXMEASUREMENTS];

} GBAS_TYPE101;


#endif

