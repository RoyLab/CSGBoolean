#pragma once
#include <typedefs.h>

struct PlaneRepsTriangle
{
	GS::double4 sPlane;
	GS::double4 bPlanes[3];
};


template <class T>
inline T RoundOff(T val)
{
	return round(val * ROUNDOFF_SCALE) / ROUNDOFF_SCALE;
}

template <class T>
inline void RoundOffVec(T& vec)
{
	vec.x = RoundOff(vec.x);
	vec.y = RoundOff(vec.y);
	vec.z = RoundOff(vec.z);
}


template <class Vec1, class Vec2>
void Uniformalize(Vec1& v, Vec2& center, Vec2& scale)
{
	return (v - center) / scale*2.0;
}

template <class Vec1, class Vec2>
void Deuniformalize(Vec1& v, Vec2& center, Vec2& scale)
{
	return v / 2.0*scale + center;
}

void ConvertToPlaneReps(GS::double3& v0, GS::double3& v1, GS::double3& v2, PlaneRepsTriangle& output);
void ConvertToPlaneReps(GS::double3 (&v)[3], PlaneRepsTriangle& output);
