/*                                                 -*- tab-width: 4; -*-
 * Quaternion.h
 *
 * CopyRight (C) 2001 Sarnoff Corporation
 * All Right Reserved
 *
 */

#ifndef QUATERNION_H
#define QUATERNION_H

#define  RCSID_QUATERNION_H "$Id$"

#include <vector>
using namespace std;

#include "Types.h"

#ifndef ROTATIONMATRIX_H
#include "RotationMatrix.h"
#endif

/**
* A quaternion used for representing rotations. Note that this class
* cannot be used for a general purpose quaternion, since some of the
* methods automatically do normalization and other operations meant
* to restrict the quaternion to valid rotations.
*******************************************************************************/
class Quaternion
{
public:
	Quaternion();
	Quaternion( real m[3][3] );
	Quaternion(real w, real x, real y, real z);
	Quaternion(Vector3D<real> axis, real angle);
	Quaternion(Vector3D<real> rotate, bool turnp = false);
	Vector3D<real> RotateVector(bool turnp = false) const;

	Quaternion& operator*=(const Quaternion& q);
	Quaternion& operator+=(const Quaternion& q);
	Quaternion& operator-=(const Quaternion& q);
	Quaternion& operator*=(real s);
	Quaternion& operator/=(real s);
	Quaternion operator*(const Quaternion& a) const;

	bool operator==(const Quaternion& q) const;
	void Normalize();
	void Canonicalize();
	real DotProduct(const Quaternion& q) const;
	size_t FindClosest(const vector<Quaternion>& l) const;

	real Magnitude() const;
	Quaternion Conjugate() const;
	Vector3D<real> transformPoint(const Vector3D<real>& pos) const;

	RotationMatrix Matrix() const;

	/**
	* get the w (real) component
	************************************/
	real w() const { return m_w; };

	/**
	* get the x (i) component
	************************************/
	real x() const { return m_x; };

	/**
	* get the y (j) component
	************************************/
	real y() const { return m_y; };

	/**
	* get the z (k) component
	************************************/
	real z() const { return m_z; };

#if !defined(NDEBUG)
	static void TurnToAngleTest();
#endif

private:

	real   m_w, m_x, m_y, m_z;
	static real AngleToTurn(real theta);
	static real TurnToAngle(real turn);

};


/**
* Multiply a quaternion by a real
********************************************************************************/
inline Quaternion operator*(const Quaternion& a, real s)
{
	Quaternion r = a;
	return r *= s;
}

/**
* Add two quaternions (this assumes the equivalence of q and -q).
********************************************************************************/
inline Quaternion operator+(const Quaternion& a, const Quaternion& b)
{
	Quaternion r = a;
	return r += b;
};

/**
* Subtract two quaternions (this assumes the equivalence of q and -q).
********************************************************************************/
inline Quaternion operator-(const Quaternion& a, const Quaternion& b)
{
	Quaternion r = a;
	return r -= b;
};

/**
* Divide a quaternion by a real
********************************************************************************/
inline Quaternion operator/(const Quaternion& q, real s)
{
	Quaternion r = q;
	return r /= s;
}

#endif
