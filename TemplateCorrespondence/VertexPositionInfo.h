#pragma once
#include <vector>

using namespace std;

namespace DeformationTransfer
{
	namespace Correspondence
	{
		class CorrespondenceProblem;
		struct VertexInfo
		{
			bool isConstrained;
			int matrixPosition;
		};

		class VertexPositionInfo
		{
		public:
			vector<VertexInfo> info;
			int numberConstrained;
			int numberFree;

			VertexPositionInfo();
			~VertexPositionInfo();

			void PopulateVertexInformationList(CorrespondenceProblem& problem);
		};
	}
}
