#pragma once
#include <Eigen/Core>

using namespace Eigen;

namespace DeformationTransfer
{
	namespace Utility
	{
		struct BoundingRect
		{
			Vector3f min;
			Vector3f max;
		};		
		float SquaredDistance(const Vector3f& a, const Vector3f& b);
		float SquaredDistanceRect(const BoundingRect *rect, const Vector3f& x0);
	}
}