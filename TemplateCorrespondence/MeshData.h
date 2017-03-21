#pragma once

#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>
#include <Eigen\Core>
#include <Eigen\SparseCore>
#include "UtilityFunctions.h"

using namespace Eigen;
using namespace std;

namespace DeformationTransfer
{
	namespace Correspondence
	{
		typedef Eigen::Matrix<float, 9, 4> Matrix94f;
		typedef Eigen::Array<float, 9, 1> Vector91f;

		typedef Eigen::SparseMatrix<float> SparseMat;
		typedef Eigen::Triplet<float> TripletF;
		typedef vector<Eigen::Triplet<float>> TripletList;

		//Forward declaration. 
		class CorrespondenceProblem;
		struct TriangleAdjacency
		{
			int adjacentTriangleCount;
			vector<unordered_set<int>> adjacencyList;
		};
		struct MeshMatrixData
		{
			Matrix94f A;
			Vector91f C;
		};
		struct TriangleMatrixData
		{
			Matrix3f VertexMatrix;
			Matrix3f InvVertexMatrix;
			Vector3f FourthVertex;
		};
		struct TriangleData
		{
			int vertexIdx[3];
			int normalIdx[3];
		};
		class MeshData
		{
		private:
			vector<TriangleData> triangleList;
			vector<Vector3f> vertexList;
			vector<Vector3f> normalList;
			Vector3f CalculateFourthVertex(Vector3f & v2_1, Vector3f & v3_1);
			Vector3f CalculateFourthVertex(Vector3f & v1, Vector3f & v2, Vector3f & v3);
			void CalculateLinearEquationMatricesForTriangle(int index, CorrespondenceProblem& problem, MeshMatrixData& current);
		public:
			vector<MeshMatrixData> LinearEquationForTriangle;
			vector<TriangleMatrixData> MatricesForTriangle;
			vector<Vector3f> centroidList;
			TriangleAdjacency triAdj;

			MeshData();
			~MeshData();

			size_t NumberOfTriangles();
			vector<TriangleData>& TriangleList();
			void TriangleList(vector<TriangleData> triList);

			size_t NumberOfVertices();
			vector<Vector3f>& VertexList();
			void VertexList(vector<Vector3f> vertList);

			size_t NumberofNormals();
			vector<Vector3f>& NormalList();
			void NormalList(vector<Vector3f> normList);

			Utility::BoundingRect ComputeBoundingBox();
			float ComputeTriangleCorrespondenceThreshold();			

			Vector3f CalculateCoefficientV1(Matrix3f invMat);
			void ComputeVertexMatricesForTriangles();
			void ComputeLinearEquationMatricesForTriangles(CorrespondenceProblem& problem);
			void ComputeCentroidsForTriangles();

			void CreateAdjacencyList();

			bool ReadObjFile(string path);
			bool SaveObjFile(string path);
		};
		Vector3f CalculateUnitNorm(MeshData * mesh, int triIndex);
	}
}


