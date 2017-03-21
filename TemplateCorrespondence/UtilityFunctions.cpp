#include "UtilityFunctions.h"

using namespace DeformationTransfer::Utility;

float DeformationTransfer::Utility::SquaredDistance(const Vector3f & a, const Vector3f & b)
{
	float x = a.x() - b.x();
	float y = a.y() - b.y();
	float z = a.z() - b.z();

	return x*x + y* y + z*z;
}
float DeformationTransfer::Utility::SquaredDistanceRect(const BoundingRect *rect, const Vector3f& x0)
{
	int i_dim;
	float dist_sq = 0;

	for (i_dim = 0; i_dim < 3; i_dim++)
	{
		if (x0[i_dim] < rect->min[i_dim]) {
			dist_sq +=
				(rect->min[i_dim] - x0[i_dim]) *
				(rect->min[i_dim] - x0[i_dim]);
		}
		else if (x0[i_dim] > rect->max[i_dim]) {
			dist_sq +=
				(rect->max[i_dim] - x0[i_dim]) *
				(rect->max[i_dim] - x0[i_dim]);
		}
	}

	return dist_sq;
}
