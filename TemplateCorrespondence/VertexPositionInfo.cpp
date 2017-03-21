#include "VertexPositionInfo.h"
#include "CorrespondenceProblem.h"

using namespace DeformationTransfer::Correspondence;


VertexPositionInfo::VertexPositionInfo()
{
}


VertexPositionInfo::~VertexPositionInfo()
{
}

void VertexPositionInfo::PopulateVertexInformationList(CorrespondenceProblem& problem)
{
	int freeVertCount = 0;
	int constrainedVertCount = 0;
	MeshData& srcMesh = problem.srcMesh;
	//VertexInformation = new VertexInfo[srcMesh.VertexList.Length];

	for (int i = 0; i < srcMesh.NumberOfVertices(); i++)
	{
		VertexInfo vertInfo;
		unordered_map<int, int>::iterator elem = problem.meshCons.constrainedList.find(i);
		if (elem != problem.meshCons.constrainedList.end())
		{
			//constrained.						
			vertInfo.isConstrained = true;
			vertInfo.matrixPosition = constrainedVertCount;
			constrainedVertCount++;
		}
		else
		{
			//Free.						
			vertInfo.isConstrained = false;
			vertInfo.matrixPosition = freeVertCount;
			freeVertCount++;
		}
		info.push_back(vertInfo);
	}

	numberConstrained = constrainedVertCount;
	numberFree = freeVertCount;
}
