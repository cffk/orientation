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

#include "Vector3D.h"
// #include "Types.h"

//#ifndef ROTATIONMATRIX_H
//#include "RotationMatrix.h"
//#endif

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
	Quaternion(double m[3][3] );
	Quaternion(double w, double x, double y, double z);
	Quaternion(Vector3D<double> axis, double angle);
	Quaternion(Vector3D<double> rotate, bool turnp = false);
	Vector3D<double> RotateVector(bool turnp = false) const;

	Quaternion& operator*=(const Quaternion& q);
	Quaternion& operator+=(const Quaternion& q);
	Quaternion& operator-=(const Quaternion& q);
	Quaternion& operator*=(double s);
	Quaternion& operator/=(double s);
	Quaternion operator*(const Quaternion& a) const;

	bool operator==(const Quaternion& q) const;
    void CircularRotate(size_t i);
	void Normalize();
	void Canonicalize();
	double DotProduct(const Quaternion& q) const;
	size_t FindClosest(const vector<Quaternion>& l) const;

	double Magnitude() const;
	Quaternion Conjugate() const;
	Vector3D<double> transformPoint(const Vector3D<double>& pos) const;

  //	RotationMatrix Matrix() const;

	/**
	* get the w (real) component
	************************************/
	double w() const { return m_w; };

	/**
	* get the x (i) component
	************************************/
	double x() const { return m_x; };

	/**
	* get the y (j) component
	************************************/
	double y() const { return m_y; };

	/**
	* get the z (k) component
	************************************/
	double z() const { return m_z; };

#if !defined(NDEBUG)
	static void TurnToAngleTest();
#endif

private:

	double   m_w, m_x, m_y, m_z;
	static double AngleToTurn(double theta);
	static double TurnToAngle(double turn);

};


/**
* Multiply a quaternion by a double
********************************************************************************/
inline Quaternion operator*(const Quaternion& a, double s)
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
* Divide a quaternion by a double
********************************************************************************/
inline Quaternion operator/(const Quaternion& q, double s)
{
	Quaternion r = q;
	return r /= s;
}

#endif
