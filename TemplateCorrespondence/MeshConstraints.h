#pragma once
#include <unordered_map>

using namespace std;

namespace DeformationTransfer
{
	namespace Correspondence
	{
		class MeshConstraints
		{
		public:
			MeshConstraints();
			~MeshConstraints();
			int numConstraints;
			unordered_map<int, int> constrainedList;

			bool ReadConstraintFile(string path);

		};
	}
}