/*
 * Vector3D.h
 *
 * CopyRight (C) 2001 Sarnoff Corporation
 * All Right Reserved
 *
 */
#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>
#include <cassert>

using namespace std;

#define RCSID_VECTOR3D_H "$Id$"

/**
* vector having x, y, and z components of type T
*******************************************************************************/
template<class T> class Vector3D
{
public:
    Vector3D() : x(0), y(0), z(0) {};
    Vector3D(T newx, T newy, T newz) : x(newx), y(newy), z(newz) {};
    T   x, y, z;

    /**
    * Add a vector to this vector
    *****************************************************/
    Vector3D<T>& operator+=(Vector3D<T> v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    };

    /**
    * Subtract a vector from this vector
    *****************************************************/
    Vector3D<T>& operator-=(Vector3D<T> v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    };

    /**
    * Multiply this vector by a scalar of type T
    *****************************************************/
    Vector3D<T>& operator*=(T n)
    {
        x *= n;
        y *= n;
        z *= n;
        return *this;
    };


    /**
    * Divide this vector by a scalar of type T
    *****************************************************/
    Vector3D<T>& operator/=(T n)
    {
        x /= n;
        y /= n;
        z /= n;
        return *this;
    };

    /**
    * Length of this vector squared.
    * This method exists because its faster to compute the square of the length 
    * than the length, and often comparing lengths squared if just as good as
    * comparing lengths.
    *****************************************************************************/
    T LengthSquared() const
    {
        return x*x + y*y + z*z;
    };

    /**
    * Length of this vector
    *****************************************************************************/
    T Length() const
    {
        return static_cast<T>(sqrt(static_cast<double>(LengthSquared())));
    };


    /**
    * Returns the dot product of this object and b (this dot b )
    ********************************************************************************/
    T DotProduct(Vector3D<T> b) const
    {
        return (x * b.x) + (y * b.y) + (z * b.z);
    };

    /**
    * Returns the angle (in radians) between this vector and v.
    ********************************************************************************/
    T AngleWith(Vector3D<T> v) const
    {
		// a = DotProduct(v) / (Length() * v.Length());
		double a = (x*v.x + y*v.y + z*v.z) /
			sqrt(static_cast<double>((x*x + y*y + z*z) * (v.x*v.x + v.y*v.y + v.z*v.z)));
		if (abs(a) > 1.0) a = 1.0;
		return static_cast<T>(acos(a));
    }

    /**
    * Returns the cross product of this object and b (this cross b)
    ********************************************************************************/
    Vector3D<T> CrossProduct(Vector3D<T> b) const
    {
        return Vector3D<T>( y * b.z - z * b.y,
                            z * b.x - x * b.z,
                            x * b.y - y * b.x);
    };

    /**
    * Normalize this vector
    *******************************************************************************/
    void Normalize()
    {
    T   magnitude;

        magnitude = static_cast<T>(sqrt(static_cast<double>(x*x + y*y + z*z)));
	assert(magnitude > 0);

        x /= magnitude;
        y /= magnitude;
        z /= magnitude;

    }

	/**
	* Subtract two vectors
	********************************************************************************/
	Vector3D<T> operator-(const Vector3D<T>& b) const
	{
		Vector3D<T> r = *this;
		return r -= b;
	}

	/**
	* Add two vectors together.
	********************************************************************************/
	Vector3D<T> operator+(const Vector3D<T>& v2) const
	{
		Vector3D<T> r = *this;
		return r += v2;
	}

	/**
	* Divide a vector by a scalar
	********************************************************************************/
	Vector3D<T> operator/(const T& n) const
	{
		Vector3D<T> r = *this;
		return r /= n;
	}
			
    Vector3D<T> operator-()
    {
        return Vector3D<T>(-x, -y, -z);
    }



	bool operator==(const Vector3D& v) const
	{
		return ((x == v.x) && (y == v.y) && (z == v.z));
	}
};




/**
* Multiply a vector by a scalar
********************************************************************************/
/*
template<class T> __inline Vector3D<T> operator*(Vector3D<T> v, T n)
{
    Vector3D<T> r = v;
    return r *= n;
};
*/



template<class T> Vector3D<T> operator*(const Vector3D<T> &v, double factor  )
{
    return Vector3D<T>((T)(v.x*factor),(T)(v.y*factor),(T)(v.z*factor ));
}
template<class T> Vector3D<T> operator*(const Vector3D<T> &v, float factor  )
{
    return Vector3D<T>((T)(v.x*factor),(T)(v.y*factor),(T)(v.z*factor ));
}
template<class T> Vector3D<T> operator*(double factor, const Vector3D<T> &v )
{
    return Vector3D<T>((T)(v.x*factor),(T)(v.y*factor),(T)(v.z*factor ));
}
template<class T> Vector3D<T> operator*(float factor, const Vector3D<T> &v )
{
    return Vector3D<T>((T)(v.x*factor),(T)(v.y*factor),(T)(v.z*factor ));
}

#endif

