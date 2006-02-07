/*                                                 -*- tab-width: 4; -*-
 *  Quaternion.cpp
 *
 *  Copyright (C) 2001-2005 Sarnoff Corporation
 *  All Rights Reserved
 *
 */

#ifndef QUATERNION_H
#include "Quaternion.h"
#endif

#include <cmath>
#include <limits>
#include <cassert>
// #include "MathUtilities.h"

namespace {
char rcsid[] = "$Id$";
char h_rcsid[] = RCSID_QUATERNION_H;
}
using namespace std;

Quaternion::Quaternion( double m[3][3] )
{
	// This implementation is from

	// E. Salamin, Application of quaternions to computation with
	// rotations, Stanford AI Lab, 1979.  This paper incorrectly asserts
	// that if some q_k^2 >= 1/4 then there's no other q_k'^2 larger.
	// This is false.  On the other hand it's true that some q_k^2 >=
	// 1/4 and here we use the first such one as the normalizer.

	double t;
	if ((t = 1 + m[0][0] + m[1][1] + m[2][2]) >= 1) {
		m_w = sqrt(t)/2;
		m_x = (m[2][1] - m[1][2])/(4*m_w);
		m_y = (m[0][2] - m[2][0])/(4*m_w);
		m_z = (m[1][0] - m[0][1])/(4*m_w);
	} else if ((t = 1 + m[0][0] - m[1][1] - m[2][2]) >= 1) {
		m_x = sqrt(t)/2;
		m_w = (m[2][1] - m[1][2])/(4*m_x);
		m_y = (m[0][1] + m[1][0])/(4*m_x);
		m_z = (m[0][2] + m[2][0])/(4*m_x);
	} else if ((t = 1 - m[0][0] + m[1][1] - m[2][2]) >= 1) {
		m_y = sqrt(t)/2;
		m_w = (m[0][2] - m[2][0])/(4*m_y);
		m_x = (m[0][1] + m[1][0])/(4*m_y);
		m_z = (m[1][2] + m[2][1])/(4*m_y);
	} else {
		// Should be >= 1 but may not be because m is not an exact rotation
		t = 1 - m[0][0] - m[1][1] + m[2][2];
		m_z = sqrt(t)/2;
		m_w = (m[1][0] - m[0][1])/(4*m_z);
		m_x = (m[0][2] + m[2][0])/(4*m_z);
		m_y = (m[1][2] + m[2][1])/(4*m_z);
	}
	Normalize();
}

/**
* Create a quaternion that represent a rotation of none
*******************************************************************************/
Quaternion::Quaternion() : m_w(1.0), m_x(0), m_y(0), m_z(0)
{
}

/**
* Create a quaternion based on the specified w, x, y, and z
* components.
*******************************************************************************/
Quaternion::Quaternion(double w, double x, double y, double z) :
	  m_w(w), m_x(x), m_y(y), m_z(z)
{
}

/**
* Create a normalized quaternion from a rotation axis and an angle
*******************************************************************************/
Quaternion::Quaternion(Vector3D<double> axis, double angle)
{
	// normalize the vector
	axis.Normalize();
	axis *= sin(angle/2);
	m_x = axis.x;
	m_y = axis.y;
	m_z = axis.z;
	m_w = cos(angle/2);

	Normalize();
}

// Convert rotate vector into quaternion.  The magnitude of the vector
// is the amount by which to rotate (in radians or "turns"), and the
// direction is the rotation axis.
Quaternion::Quaternion(Vector3D<double> rotate, bool turnp)
{
	double angle = rotate.Length();
	if (angle > 0)
		rotate /= angle;
	if (turnp)
		angle = TurnToAngle(angle);
	rotate *= sin(angle/2);
	m_x = rotate.x;
	m_y = rotate.y;
	m_z = rotate.z;
	m_w = cos(angle/2);

	Normalize();
}

// Return the rotate vector for this quaternion.
Vector3D<double> Quaternion::RotateVector(bool turnp) const
{
	double c = m_w;
	double s = sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
	Vector3D<double> rotate(m_x, m_y, m_z);
	if (c < 0) {
		c = -c;
		rotate *= -1;
	}
	// c and s are non-negative
	double angle = 2 * atan2(s, c); // in [0, pi]
	if (turnp)
		angle = AngleToTurn(angle);
	rotate *= angle/(s == 0 ? 1 : s);
	return rotate;
}

double Quaternion::Magnitude() const
{
	return sqrt(m_w*m_w + m_x*m_x + m_y*m_y + m_z*m_z);
}

/**
* Return the conjugate of this quaternion
*******************************************************************************/
Quaternion Quaternion::Conjugate() const
{
	return Quaternion(m_w, -m_x, -m_y, -m_z);
}

void Quaternion::CircularRotate(size_t i) {
  i %= 4;
  if (i == 1) {
	double t = m_w;
	m_w = m_x;
	m_x = m_y;
	m_y = m_z;
	m_z = t;
  } else if (i == 2) {
	double t = m_w;
	m_w = m_y;
	m_y = t;
	t = m_x;
	m_x = m_z;
	m_z = t;
  } else if (i == 3) {
	double t = m_w;
	m_w = m_z;
	m_z = m_y;
	m_y = m_x;
	m_x = t;
  }
  return;
}

/**
* Normalize this quaternion
*******************************************************************************/
void Quaternion::Normalize()
{
	double t = Magnitude();
	assert(t > 0);
	t = 1/t;

	m_w *= t;
	m_x *= t;
	m_y *= t;
	m_z *= t;
}

void Quaternion::Canonicalize() {
	Normalize();
	if (m_w > 0)
		return;
	if (m_w < 0) {
		m_w *= -1; m_x *= -1; m_y *= -1; m_z *= -1; return;
	}
	if (m_x > 0)
		return;
	if (m_x < 0) {
		m_x *= -1; m_y *= -1; m_z *= -1; return;
	}
	if (m_y > 0)
		return;
	if (m_y < 0) {
		m_y *= -1; m_z *= -1; return;
	}
	if (m_z > 0)
		return;
	if (m_z < 0) {
		m_z *= -1; return;
	}
}

/**
* Returns a rotation matrix that will produce the rotation represented by
* this quaternion.
*******************************************************************************/
/*
RotationMatrix Quaternion::Matrix() const
{
	RotationMatrix m;

	// See http://www.euclideanspace.com/maths/geometry/
	// rotations/conversions/quaternionToMatrix/index.htm
	// Note that the expression below are only valid if the quaternion
	// is normalized
	m.m_matrix[0][0] = 1 - 2*m_y*m_y - 2*m_z*m_z;
	m.m_matrix[0][1] = 	   2*m_x*m_y - 2*m_z*m_w;
	m.m_matrix[0][2] = 	   2*m_x*m_z + 2*m_y*m_w;
	m.m_matrix[1][0] = 	   2*m_x*m_y + 2*m_z*m_w;
	m.m_matrix[1][1] = 1 - 2*m_x*m_x - 2*m_z*m_z;
	m.m_matrix[1][2] = 	   2*m_y*m_z - 2*m_x*m_w;
	m.m_matrix[2][0] = 	   2*m_x*m_z - 2*m_y*m_w;
	m.m_matrix[2][1] = 	   2*m_y*m_z + 2*m_x*m_w;
	m.m_matrix[2][2] = 1 - 2*m_x*m_x - 2*m_y*m_y;

	return m;
}
*/

/**
* Multiply this quaternion by q, which adds the two rotations
* represented by this and q
*******************************************************************************/
Quaternion& Quaternion::operator *=(const Quaternion& q)
{
	double w = m_w;
	double x = m_x;
	double y = m_y;
	double z = m_z;

	m_w = w*q.m_w - x*q.m_x - y*q.m_y - z*q.m_z;
	m_x = w*q.m_x + x*q.m_w + y*q.m_z - z*q.m_y;
	m_y = w*q.m_y + y*q.m_w + z*q.m_x - x*q.m_z;
	m_z = w*q.m_z + z*q.m_w + x*q.m_y - y*q.m_x;
/*
	// Even more optimized
	// From http://www.lboro.ac.uk/departments/ma/gallery/quat/src/quat.cc
	// 8 multiplies + 1 halve + 27 adds instead of 16 multiplies + 12 adds
	// Not sure it's worth it...

	double t0 = (m_z - m_y) * (q.m_y - q.m_z);
	double t1 = (m_w + m_x) * (q.m_w + q.m_x);
	double t2 = (m_w - m_x) * (q.m_y + q.m_z);
	double t3 = (m_z + m_y) * (q.m_w - q.m_x);
	double t4 = (m_z - m_x) * (q.m_x - q.m_y);
	double t5 = (m_z + m_x) * (q.m_x + q.m_y);
	double t6 = (m_w + m_y) * (q.m_w - q.m_z);
	double t7 = (m_w - m_y) * (q.m_w + q.m_z);
	double t8 = t5 + t6 + t7;
	double t9 = (t4 + t8)/2;
	m_w = t0 + t9 - t5;
	m_x = t1 + t9 - t8;
	m_y = t2 + t9 - t7;
	m_z = t3 + t9 - t6;
*/
	Normalize();
	return *this;
}

/**
* Add q to this quaternion (this assumes the equivalence of q and -q)
*******************************************************************************/
Quaternion& Quaternion::operator +=(const Quaternion& q)
{
	double sign = (this->DotProduct(q) >= 0 ? 1 : -1);
	m_w += sign * q.m_w;
	m_x += sign * q.m_x;
	m_y += sign * q.m_y;
	m_z += sign * q.m_z;

	return *this;
}

/**
* Subtract q from this quaternion (this assumes the equivalence of q and -q)
*******************************************************************************/
Quaternion& Quaternion::operator -=(const Quaternion& q)
{
	double sign = (this->DotProduct(q) >= 0 ? 1 : -1);
	m_w -= sign * q.m_w;
	m_x -= sign * q.m_x;
	m_y -= sign * q.m_y;
	m_z -= sign * q.m_z;

	return *this;
}
/**
* Divide (scale) this quaternion by a double
*******************************************************************************/
Quaternion& Quaternion::operator /=(double s)
{
	m_w /= s;
	m_x /= s;
	m_y /= s;
	m_z /= s;

	return *this;
}

/**
* Multiply (scale) this quaternion by a double
*******************************************************************************/
Quaternion& Quaternion::operator *=(double s)
{
	m_w *= s;
	m_x *= s;
	m_y *= s;
	m_z *= s;

	return *this;
}

/**
* Multiply two quaternions, which adds the two rotations they represent
********************************************************************************/
Quaternion Quaternion::operator*(const Quaternion& a) const
{
	Quaternion r = *this;
	return r *= a;
};

/**
* Equality of quaternions
********************************************************************************/
bool Quaternion::operator==(const Quaternion& q) const
{
	return (m_w == q.m_w && m_x == q.m_x && m_y == q.m_y && m_z == q.m_z) ||
		(m_w == -q.m_w && m_x == -q.m_x && m_y == -q.m_y && m_z == -q.m_z);
}

double Quaternion::DotProduct(const Quaternion& q) const
{
	return m_w*q.m_w + m_x*q.m_x + m_y*q.m_y + m_z*q.m_z;
}
Vector3D<double> Quaternion::transformPoint(const Vector3D<double>& pos) const
{
	double A, B, C, D;

	// Form q * (0,p)

	A = - pos.z*m_z - pos.y*m_y - pos.x*m_x;
	B = - pos.y*m_z + pos.z*m_y + pos.x*m_w;
	C =   pos.x*m_z - pos.z*m_x + pos.y*m_w;
	D = - pos.x*m_y + pos.y*m_x + pos.z*m_w;

	// Form ( q * (0,p) ) * conj(q)

	Vector3D<double> newpos;
	newpos.x =  - C*m_z + D*m_y - A*m_x + B*m_w;
	newpos.y =    B*m_z - A*m_y - D*m_x + C*m_w;
	newpos.z =  - A*m_z - B*m_y + C*m_x + D*m_w;
	return newpos;
}

double Quaternion::AngleToTurn(double theta) {
	// Taylor series is t/(6*pi)^(1/3) * (1 - t^2/60 + t^4/8400 + ...)
	if (abs(theta) <=
		// sqrt(sqrt(8400 * numeric_limits<double>::epsilon()))
#if defined(DOUBLE)
		0.00116863766081		// epsilon = 2^-52
#else
		0.1778882843006f		// epsilon = 2^-23
#endif
		)
		// 0.3757..  = 1/(6*pi)^(1/3)
		return 0.375750550595608872207534 * theta * (1 - theta * theta/60);
	double dtheta = theta;
	// If cbrt(x) gets replaced by pow(x, 1/3.0) then need to treat
	// negative x properly: cbrt(-x) = - cbrt(x).
	return cbrt((dtheta-sin(dtheta))/M_PI);
}

// The inverse of the previous function.
double Quaternion::TurnToAngle(double turn) {
	double sign = turn < 0 ? -1 : 1;
	double u = abs(turn);
	double theta;
	if (u <= 1) {
		// Expand about u = 0

		// Taylor series is (6*pi)^(1/3)*u * (1 + (pi^2/6000)^(1/3)*u^2
		// + ((81*pi^4)/171500000)^(1/3)*u^4 + pi^2/700*u^6 + ...)

		double u2 = u * u;
		// 2.661... = (6*pi)^(1/3)
		// 0.118... = (pi^2/6000)^(1/3)
		// 0.035... = ((81*pi^4)/171500000)^(1/3)
		theta =  2.66134007898293758185122 * u *
			(1 + u2 * (0.118045516933348474107871 +
					   0.0358321990321580434428649 * u2));
		// 4.139... = (700/pi^2)^(1/3)
		if (u2 < 4.13935586588440243 * cbrt(numeric_limits<double>::epsilon()))
			return sign * theta;
	} else {
		// Expand about u = 2^(1/3)
		// 1.259... = 2^(1/3)
		double v = 1.25992104989487316476721 - u;
		if (v >= 0) {
			// The case u <= 2^(1/3).
			// 89.765... = 36*pi/2^(1/3)
			//  1.496... = 3*pi/5/2^(1/3)
			theta = (2.0 * M_PI) - cbrt(89.7654146969520912005595 * v) -
				1.49609024494920152000932 * v;
			// 3.495... = 2^(7/12)*3^(1/4)*sqrt(%pi)
			if ( v < 3.495071629824412029 *
				 pow(numeric_limits<double>::epsilon(), 0.75) )
				return sign * theta;
		} else {
			// The general large argument case
			double u3 = u * u * u;
			double n = 2 * floor( 0.5 * (u3 + 1) );
			// 2.661... = (6*pi)^(1/3)
			double delta =  2.66134007898293758185122 * cbrt(u3 - n);
			theta = M_PI * n + delta;
			if (abs(delta) <= sqrt(2 * numeric_limits<double>::epsilon()) ||
				1-cos(theta) == 0)
				return sign * theta;
		}
	}
	double rhs = M_PI * u * u * u;
	double delta;
	bool last = false;
	// If double is float use epsilon, if real is double use 50*epsilon
	// Ask for absolute accuracy for small theta, relative otherwise
	double eps = (50 * numeric_limits<double>::epsilon()) *
		max(1.0, theta/(2.0 * M_PI));
	// Newton's method.  Should only need 1-3 iterations.
	for (size_t i = 0; i < 10; i++) {
		// A standard Newton step on theta - sin(theta) - pi * u^3
		delta = (theta - sin(theta) - rhs)/(1 - cos(theta));
		theta = theta - delta;
		// Do one iteration beyond nominal convergence
		if (last)
			return sign * theta;
		last = abs(delta) < eps;
	}
	assert(delta == 0.0);
	return numeric_limits<double>::signaling_NaN();
}

#if !defined(NDEBUG)
#include <iostream>
void Quaternion::TurnToAngleTest() {
	for (int i = -100; i <= 100; i++) {
		double turn = 0.1 * i;
		double theta = TurnToAngle(turn);
		double turnx = AngleToTurn(theta);
		cout << turn << " " << theta << " "
			 << turnx << " " << turnx - turn
			 << endl;
	}
	for (int i = -100; i <= 10; i++) {
		double turn = pow(2.0, double(i));
		double theta = TurnToAngle(turn);
		double turnx = AngleToTurn(theta);
		cout << turn << " " << theta << " "
			 << turnx << " " << turnx - turn
			 << endl;
	}
	for (int i = -10; i <= 10; i++) {
		double turn0 = cbrt(double(i));
		if (i == 0) continue;
		for (int j = -20; j <= 20; j++) {
			int k = (j == 0) ? 0 :
				(j < 0 ? -1 : 1);
			k = k * (1 << (abs(j) -1));
			double turn = turn0 *
				(1 + double(k) * numeric_limits<double>::epsilon());
			double theta = TurnToAngle(turn);
			double turnx = AngleToTurn(theta);
			cout << turn << " " << theta << " "
				 << turnx << " " << turnx - turn
				 << endl;
		}
	}
}
#endif

// return index to the quaternion in list which is "closest" to this.
// Closest is defined as having the max(abs(this.list[i])).  All these
// quaternions are assumed to have unit magnitude.  The abs() accounts
// for the +/- ambiguity with quaternions.
size_t Quaternion::FindClosest(const vector<Quaternion>& l) const {
	double maxt = -1;
	size_t result = l.size();
	for (size_t i = 0; i < l.size(); i++) {
		double t = abs(DotProduct(l[i]));
		if (t > maxt) {
			maxt = t;
			result = i;
		}
	}
	assert(result < l.size());
	return result;
}
