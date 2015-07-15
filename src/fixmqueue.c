/**
 * \file fixmqueue.c
 * 
 * Standalone program which unlinks the message queues in use by gpsins if it crashes.
 *
 * $Id: fixmqueue.c 340 2007-10-08 11:47:09Z greerd $
 */

#include <mqueue.h>

#define CONTROL_MQUEUE_NAME "/control-queue"
#define GPSINS_MQUEUE_NAME "/gpsins-queue"
#define	INS_MQUEUE_NAME	"/ins-queue"

void main(void)
{
	mq_unlink(CONTROL_MQUEUE_NAME);
	mq_unlink(GPSINS_MQUEUE_NAME);
	mq_unlink(INS_MQUEUE_NAME);
}
