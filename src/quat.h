/** 
 *    \file quat.h  Implements Quaternion functions
 *
*   $Id: quat.h 365 2007-10-23 23:26:48Z n2523710 $
 */
 
 #ifndef _QUAT_H_
#define _QUAT_H_


#define PI	3.1415926
#define RAD2DEG 180.0/PI
#define DEG2RAD PI/180.0


/*
 * This will construct a direction cosine matrix from 
 * euler angles in the standard rotation sequence 
 * [phi][theta][psi] from NED to body frame
 *
 * body = tBL(3,3)*NED
 */

void eulerDC(double tBL[3][3],double phi,double theta,double psi);

/*
 * This will construct a direction cosine matrix from 
 * quaternions in the standard rotation  sequence
 * [phi][theta][psi] from NED to body frame
 *
 * body = tBL(3,3)*NED
 * q(4,1)
 */

void quatDC(double tBL[3][3],const double q[4]);


/*
 * This will construct the euler omega-cross matrix
 * wx(3,3)
 * p, q, r (rad/sec)
 */

void eulerWx(double Wx[3][3],double p,double q,double r);

/*
 * This will construct the quaternion omega matrix
 * W(4,4)
 * p, q, r (rad/sec)
 */

void quatW(double W[4][4],double p,double q,double r);


/*
 * This will normalize a quaternion vector q
 * q/norm(q)
 * q(4,1)
 */

void normq(double q[4]);


/*
 * This will convert from quaternions to euler angles
 * q(4,1) -> euler[phi;theta;psi] (rad)
 */

void quat2euler(double q[4],double euler[3]);


/*
 * This will convert from euler angles to quaternion vector
 * phi, theta, psi -> q(4,1)
 * euler angles in radians
 */

void euler2quat(double q[4],double phi,double theta,double psi);

/* 
 * quaternion inverse 
 */

void quatinv(double r[4], double q[4]);

/*
 * quaternion multiplication
 */
void quatmul(double r[4], double p[4], double q[4]);

void T_alphaFromQuat(double T_alpha[4][3],double q[4]);


#endif
