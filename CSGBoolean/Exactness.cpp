#include "precompile.h"
#include "Exactness.h"
#include <adaptive.h>
#include <typedefs.h>
#include <arithmetic.h>

using GS::double3;
using GS::double4;

const double ROUNDOFF_SCALE = std::pow(2.0, 18);
const double EDGE_LENGTH_LIMIT = std::pow(2.0, -3.5);
const double EDGE_LENGTH_OFFSET = std::pow(2.0, -4);


inline void ConvertToPlaneReps(double3& v0, double3& e0, double3& e1, double4& plane)
{
	double3 norm = GS::cross(e0, e1);
	plane = double4(norm, -GS::dot(norm, v0));
}

void ConvertToPlaneReps(GS::double3(&v)[3], PlaneRepsTriangle& output)
{
	ConvertToPlaneReps(v[0], v[1], v[2], output);
}

void ConvertToPlaneReps(double3& v0, double3& v1, double3& v2, PlaneRepsTriangle& output)
{
	double3 e0 = v2 - v1;
	double3 e1 = v0 - v2;
	double3 e2 = v1 - v0;
	double3 norm = GS::cross(e0, e1);
	output.sPlane = double4(norm, -GS::dot(norm, v1));

	double3 normUni = GS::normalize(norm);
	double3 newPoint = normUni*EDGE_LENGTH_OFFSET + v0;
	RoundOffVec(newPoint);

	double3 b0 = newPoint - v0;
	double3 b1 = newPoint - v1;
	double3 b2 = newPoint - v2;

	ConvertToPlaneReps(v0, b1, e0, output.bPlanes[0]);
	ConvertToPlaneReps(v1, b2, e1, output.bPlanes[1]);
	ConvertToPlaneReps(v2, b0, e2, output.bPlanes[2]);
}

void ComputeCoord(double4& p0, double4& p1, double4& p2, double3& v)
{
	GS::double3x3 mat;
	mat[0] = p0.xyz;
	mat[1] = p1.xyz;
	mat[2] = p2.xyz;
	double D = GS::ExactDet3x3(mat);

	mat[0][0] = p0.w;
	mat[1][0] = p1.w;
	mat[2][0] = p2.w;
	double D1 = GS::ExactDet3x3(mat);

	mat[0][1] = p0.w;
	mat[1][1] = p1.w;
	mat[2][1] = p2.w;

	mat[0][0] = p0.x;
	mat[1][0] = p1.x;
	mat[2][0] = p2.x;
	double D2 = GS::ExactDet3x3(mat);

	mat[0][2] = p0.w;
	mat[1][2] = p1.w;
	mat[2][2] = p2.w;

	mat[0][1] = p0.y;
	mat[1][1] = p1.y;
	mat[2][1] = p2.y;
	double D3 = GS::ExactDet3x3(mat);

	v.x = D1 / D;	v.y = D2 / D;	v.z = D3 / D;
}

