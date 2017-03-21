#pragma once
#include <vector>
#include "MeshData.h"
#include "VertexPositionInfo.h"
#include "MeshConstraints.h"
#include "KDTree.h"

namespace DeformationTransfer
{
	namespace Correspondence
	{

#pragma region DataStructures

		struct ResultContainer
		{			
			int tgtIndex;
			float distance;
		};
#pragma endregion

		class CorrespondenceProblem
		{		
		public:
#pragma region DataContainers
			MeshData srcMesh, tgtMesh;
			VertexPositionInfo vInfo;			
			MeshConstraints meshCons;
#pragma endregion

#pragma region Weights
			float wtSmoothness;
			float wtIdentity;
			float wtClosestStart, wtClosestEnd, wtClosestStep;
#pragma endregion

			//Result Container Declaration 
			unordered_map<int,vector<ResultContainer>> CorrespondenceResultMap;
		private:
			int MatrixPositionForFreeVertex(int vertIdx, int componentIdx);
			int MatrixPositionForTriangle(int triIdx, int vertIdx, int componentIdx);
			void AddToLinearEquationSystemSparse_Unrolled(TripletList & A, int triIdx, int row, float weight);
			void AddToLinearEquationSystemDence_Unrolled(VectorXf& C, int triIdx, int row, float weight);
			void AddToLinearEquationSystemDence_Unrolled(VectorXf& C, Vector91f& c, int row, float weight);
			int MatrixPositionForFourthVertex(int triIdx, int vertIdx, int componentIdx);
		
		public:
			CorrespondenceProblem();
			~CorrespondenceProblem();
			
#pragma region  Methods
			void Solve();
			void SolvePhase1();
			void SolvePhase2();
			void SolveEquations(TripletList & tripA, VectorXf & C, int nRow, int nCol);
			int BuildPhase2Equation(TripletList & tripA, VectorXf & C, vector<int>& correspondenceMap, float wtClosestPoint);
			int BuildPhase1Equation(TripletList& tripA, VectorXf& C);			
			int CreateSmoothnessEquation(TripletList& tripA, VectorXf& C, int row);
			int CreateIdentityEquation(TripletList & tripA, VectorXf & C, int row);
			int CreateClosestPointEquation(TripletList & tripA, VectorXf & C, vector<int>& correspondenceMap, int row, float wtClosest);
			void CreateCorrespondenceMap_Basic(MeshData src, MeshData tgt, vector<int>& correspondenceMap);
			void ResolveTriangleCorrespondence(float threshold);
			void ApplyDeformationToSourceMesh(VectorXf & deformedX);
#pragma endregion

		};
	}
}
