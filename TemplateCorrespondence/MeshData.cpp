#include "MeshData.h"
#include <set>
#include <algorithm>
#include <iterator>
#include <Eigen\LU>
#include <Eigen\Geometry>
#include "CorrespondenceProblem.h"

using namespace DeformationTransfer::Correspondence;
using namespace DeformationTransfer::Utility;
using namespace std;


MeshData::MeshData()
{
}

MeshData::~MeshData()
{
}

size_t MeshData::NumberOfTriangles()
{
	return triangleList.size();
}

inline vector<TriangleData>& MeshData::TriangleList()
{
	return triangleList;
}

inline void MeshData::TriangleList(vector<TriangleData> triList)
{
	triangleList = triList;
}

size_t MeshData::NumberOfVertices()
{
	return vertexList.size();
}

vector<Vector3f>& MeshData::VertexList()
{
	return vertexList;
}

void MeshData::VertexList(vector<Vector3f> vertList)
{
	vertexList = vertList;
}

size_t MeshData::NumberofNormals()
{
	return normalList.size();
}

vector<Vector3f>& MeshData::NormalList()
{
	return normalList;
}

void MeshData::NormalList(vector<Vector3f> normList)
{
	normalList = normalList;
}
BoundingRect MeshData::ComputeBoundingBox()
{
	BoundingRect bRect; 
	Vector3f vMax = vertexList[0];
	Vector3f vMin = vertexList[0];

	/* find the bounding box of the model */
	for (int iv = 0; iv < NumberOfVertices(); iv++)
	{
		Vector3f cV = vertexList[iv];
		vMax.x() = max(cV.x(), vMax.x());  vMin.x() = min(cV.x(), vMin.x());
		vMax.y() = max(cV.y(), vMax.y());  vMin.y() = min(cV.y(), vMin.y());
		vMax.z() = max(cV.z(), vMax.z());  vMin.z() = min(cV.z(), vMin.z());
	}
	bRect.max = vMax;
	bRect.min = vMin;
	return bRect;
}
float MeshData::ComputeTriangleCorrespondenceThreshold()
{
	int iv = 1;
	BoundingRect bBox = ComputeBoundingBox();

	/* width, height and depth of model's bounding box */
	float dx = bBox.max.x() - bBox.min.x();
	float dy = bBox.max.y() - bBox.min.y();
	float dz = bBox.max.z() - bBox.min.z();	

	return sqrt(4 * (dx*dy + dy*dz + dx*dz) / NumberOfTriangles());	
}
Vector3f DeformationTransfer::Correspondence::CalculateUnitNorm(MeshData* mesh , int triIndex)
{
	Vector3f& v0 = mesh->VertexList()[mesh->TriangleList()[triIndex].vertexIdx[0]];
	Vector3f& v1 = mesh->VertexList()[mesh->TriangleList()[triIndex].vertexIdx[1]];
	Vector3f& v2 = mesh->VertexList()[mesh->TriangleList()[triIndex].vertexIdx[2]];

	Vector3f unitNorm = (v1 - v0).cross(v2 - v0);
	return unitNorm;
}
Vector3f MeshData::CalculateCoefficientV1(Matrix3f invMat)
{

	Vector3f CoeffV1;
	
	CoeffV1.x() = invMat(0, 0) + invMat(0, 1) + invMat(0, 2);
	CoeffV1.y() = invMat(1, 0) + invMat(1, 1) + invMat(1, 2);
	CoeffV1.z() = invMat(2, 0) + invMat(2, 1) + invMat(2, 2);

	return -1 * CoeffV1;
}
Vector3f  MeshData::CalculateFourthVertex(Vector3f& v2_1, Vector3f& v3_1)
{
	auto v4_1 = v2_1.cross(v3_1);
	float norm = v4_1.norm();
	v4_1 = v4_1 / sqrt(norm);

	return v4_1;
}

Vector3f  MeshData::CalculateFourthVertex(Vector3f& v1, Vector3f& v2, Vector3f& v3)
{
	return CalculateFourthVertex(Vector3f(v2 - v1), Vector3f(v3 - v1));
}
void MeshData::ComputeVertexMatricesForTriangles()
{
	for (int i = 0; i < NumberOfTriangles(); i++)
	{
		int v1i = triangleList[i].vertexIdx[0];
		int v2i = triangleList[i].vertexIdx[1];
		int v3i = triangleList[i].vertexIdx[2];

		Vector3f v2_1 = vertexList[v2i] - vertexList[v1i];
		Vector3f v3_1 = vertexList[v3i] - vertexList[v1i];
		Vector3f v4_1 = CalculateFourthVertex(v2_1, v3_1);

		TriangleMatrixData triMatData;
		triMatData.FourthVertex = v4_1;

		triMatData.VertexMatrix << v2_1.x(), v2_1.y(), v2_1.z(),
			v3_1.x(), v3_1.y(), v3_1.z(),
			v4_1.x(), v4_1.y(), v4_1.z();

		triMatData.InvVertexMatrix = triMatData.VertexMatrix.inverse();

		MatricesForTriangle.push_back(triMatData);
	}
}

void MeshData::CalculateLinearEquationMatricesForTriangle(int index, CorrespondenceProblem& problem, MeshMatrixData& current)
{

	//var A = srcMesh.matricesForTriangle[index].A;
	//var C = srcMesh.matricesForTriangle[index].C;

	// = LinearEquationForTriangle[index];

	Matrix3f inverseMat = problem.srcMesh.MatricesForTriangle[index].InvVertexMatrix;
	Vector3f coefv1 = problem.srcMesh.CalculateCoefficientV1(inverseMat);
	float coef = 0;

	int v1 =  0;
	int v2 =  1;
	int v3 =  2;
	unordered_map<int, int>::iterator consIter;
	for (int i = 0, row = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++, row++)
		{
			current.C[row] = 0;

			/* coefficient of v1: -(a[0,i] + a[1,i] + a[2,i]) */
			coef = coefv1[j];
			int key = (problem.srcMesh.TriangleList())[index].vertexIdx[v1];
			consIter = problem.meshCons.constrainedList.find(key);
			if ((consIter == problem.meshCons.constrainedList.end()))
			{
				current.A(row, 0) = coef;
			}
			else
			{
				int tgtVertIdx = consIter->second;
				current.C[row] -= coef * problem.tgtMesh.VertexList()[tgtVertIdx][i];
			}


			/* coefficient of v2: a[0,i]*/
			coef = inverseMat(j, 0);
			key = (problem.srcMesh.TriangleList())[index].vertexIdx[v2];
			consIter = problem.meshCons.constrainedList.find(key);
			if ((consIter == problem.meshCons.constrainedList.end()))
			{
				current.A(row, 1) = coef;
			}
			else
			{
				int tgtVertIdx = consIter->second;
				current.C[row] -= coef * problem.tgtMesh.VertexList()[tgtVertIdx][i];
			}

			/* coefficient of v3: a[1,i] */
			coef = inverseMat(j, 1);
			key = (problem.srcMesh.TriangleList())[index].vertexIdx[v3];
			consIter = problem.meshCons.constrainedList.find(key);
			if ((consIter == problem.meshCons.constrainedList.end()))
			{
				current.A(row, 2) = coef;
			}
			else
			{
				int tgtVertIdx = consIter->second;
				current.C[row] -= coef * (problem.tgtMesh.VertexList()[tgtVertIdx][i]);
			}

			/* coefficient of v4: a[2,i] */			
			current.A(row, 3) = inverseMat(j, 2);
			//std::cout << current.A(row, 0) << "\t" << current.A(row, 1) << "\t" << current.A(row, 2) << "\t" << current.A(row, 3) << endl;
		}
	}
	//for (int i = 0; i< 9; i++)
	//	std::cout << "[" << i << "] : " << current.C[i] << endl;
}

void MeshData::ComputeLinearEquationMatricesForTriangles(CorrespondenceProblem& problem)
{
	int triIndex = 0;
	//initialize the list of matrices. 
	LinearEquationForTriangle.reserve(NumberOfTriangles());

	for (triIndex = 0; triIndex < NumberOfTriangles(); triIndex++)
	{		
		MeshMatrixData mData;
		mData.A = Matrix94f::Zero();
		mData.C = Vector91f::Zero();
		CalculateLinearEquationMatricesForTriangle(triIndex, problem, mData);
		LinearEquationForTriangle.push_back(mData);
		//if (display)
		//	DisplayMatrices(triIndex);
	}
}

void DeformationTransfer::Correspondence::MeshData::ComputeCentroidsForTriangles()
{
	Vector3f v0;
	Vector3f v1;
	Vector3f v2;	
	for (int i = 0; i < NumberOfTriangles(); i++)
	{
		v0 = vertexList[triangleList[i].vertexIdx[0]];
		v1 = vertexList[triangleList[i].vertexIdx[1]];
		v2 = vertexList[triangleList[i].vertexIdx[2]];
		Vector3f centroid;

		/* cartesian coordinates of centroid are the means of the coordinates of
		the three vertices */
		centroid.x() = (v0.x() + v1.x() + v2.x()) / 3;
		centroid.y() = (v0.y() + v1.y() + v2.y()) / 3;
		centroid.z() = (v0.z() + v1.z() + v2.z()) / 3;
		centroidList.push_back(centroid);
	}
	
}

void MeshData::CreateAdjacencyList()
{
	cout << "Creating Adjacency List" << endl;
	const int numVertices = NumberOfVertices();		
	const int numTris = NumberOfTriangles();
	
	vector<set<int>> vertMembershipList;
	vertMembershipList.reserve(numVertices);
	triAdj.adjacencyList.reserve(numTris);

	int triangleIndex = 0;
	int i = 0;

	//initialize arrays
	for (i = 0; i < numTris; i++)
	{
		unordered_set<int> temp;
		triAdj.adjacencyList.push_back(temp);
	}
		
	for (i = 0; i < numVertices; i++)
	{
		set<int> temp;
		vertMembershipList.push_back(temp);
	}
		

	// Prepare data for faster triangle adjacency list creation. 
	for (triangleIndex = 0; triangleIndex < numTris; triangleIndex++)
	{
		int v1i = 0;
		int v2i = 1;
		int v3i = 2;

		//Check if index is within bounds.                          
		vertMembershipList[triangleList[triangleIndex].vertexIdx[v1i]].insert(triangleIndex);
		vertMembershipList[triangleList[triangleIndex].vertexIdx[v2i]].insert(triangleIndex);
		vertMembershipList[triangleList[triangleIndex].vertexIdx[v3i]].insert(triangleIndex);
	}
	triangleIndex = 0;
	int totalAdjTriangles = 0;
	// find adjacency
	for (triangleIndex = 0; triangleIndex < numTris; triangleIndex++)
	{
		int v1i = 0;
		int v2i = 1;
		int v3i = 2;

		vector<int> v1v2; 
		vector<int> v1v3;
		vector<int> v2v3;
		
		set<int>& v1 = vertMembershipList[triangleList[triangleIndex].vertexIdx[v1i]];
		set<int>& v2 = vertMembershipList[triangleList[triangleIndex].vertexIdx[v2i]];
		set<int>& v3 = vertMembershipList[triangleList[triangleIndex].vertexIdx[v3i]];

		set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v1v2));
		set_intersection(v1.begin(), v1.end(), v3.begin(), v3.end(), back_inserter(v1v3));
		set_intersection(v2.begin(), v2.end(), v3.begin(), v3.end(), back_inserter(v2v3));
/*
		List<IEnumerable<int>> adjTriangleIndexList = new List<IEnumerable<int>>();
		var v1v2 = vertMembershipList[triangles[v1i]].Intersect(vertMembershipList[triangles[v2i]]);
		var v1v3 = vertMembershipList[triangles[v1i]].Intersect(vertMembershipList[triangles[v3i]]);
		var v2v3 = vertMembershipList[triangles[v2i]].Intersect(vertMembershipList[triangles[v3i]]);*/

		/*adjTriangleIndexList.Add(v1v2);
		adjTriangleIndexList.Add(v1v3);
		adjTriangleIndexList.Add(v2v3);*/

		if (v1v2.size() > 1)
		{
			//max two indices should be available
			//one edge can only be shared between two triangles, hence max two. 
			totalAdjTriangles++;
			triAdj.adjacencyList[v1v2[0]].insert(v1v2[1]);
			triAdj.adjacencyList[v1v2[1]].insert(v1v2[0]);
		}
		if (v2v3.size() > 1)
		{
			//max two indices should be available
			//one edge can only be shared between two triangles, hence max two. 
			totalAdjTriangles++;
			triAdj.adjacencyList[v2v3[0]].insert(v2v3[1]);
			triAdj.adjacencyList[v2v3[1]].insert(v2v3[0]);
		}
		if (v1v3.size() > 1)
		{
			//max two indices should be available
			//one edge can only be shared between two triangles, hence max two. 
			totalAdjTriangles++;
			triAdj.adjacencyList[v1v3[0]].insert(v1v3[1]);
			triAdj.adjacencyList[v1v3[1]].insert(v1v3[0]);
		}

		//for (int j = 0; j < adjTriangleIndexList.Count; j++)
		//{
		//	var triangleForEdge = adjTriangleIndexList[j];
		//	if (triangleForEdge.Count() >  1)
		//	{
		//		//max two indices should be available
		//		//one edge can only be shared between two triangles, hence max two. 
		//		totalAdjTriangles++;
		//		adjList[triangleForEdge.ElementAt(0)].Add(triangleForEdge.ElementAt(1));
		//		adjList[triangleForEdge.ElementAt(1)].Add(triangleForEdge.ElementAt(0));
		//	}
		//}
	}	
	triAdj.adjacentTriangleCount = totalAdjTriangles;
	cout << "Adjacency List Complete! - Total adjacent triangles : "<<totalAdjTriangles << endl;
}

bool MeshData::ReadObjFile(string path)
{
	FILE* file = fopen(path.c_str(), "r");
	char prefix[3];   // prefix reader - 3 chars should be sufficient. 

	if (file == nullptr)
	{
		cout << "Could not open file :" << path << endl;
		fclose(file);
		return false;
	}
	else
	{
		while (fscanf(file, "%2s", prefix) != EOF)
		{
			if (strcmp(prefix, "v") == 0)
			{
				Vector3f vert;
				fscanf(file, "%f %f %f\n", &(vert.x()), &(vert.y()), &(vert.z()));
				vertexList.push_back(vert);
			}
			else if (strcmp(prefix, "vn") == 0)
			{
				Vector3f normal;
				fscanf(file, "%f %f %f\n", &(normal.x()), &(normal.y()), &(normal.z()));
				normalList.push_back(normal);
			}
			else if (strcmp(prefix, "f") == 0)
			{
				TriangleData triData;
				for (int i = 0; i < 3; i++)
				{
					fscanf(file, "%d/", &(triData.vertexIdx[i]));
					if ((prefix[0] = (char)fgetc(file)) != '/')
					{
						//skip texture coordinates if found. 
						ungetc(prefix[0], file);  fscanf(file, "%*d/");
					}
					fscanf(file, "%d/", &(triData.normalIdx[i]));
				}
				//Change obj 1 index to c++ 0 index. 
				triData.vertexIdx[0]--;
				triData.vertexIdx[1]--;
				triData.vertexIdx[2]--;

				triData.normalIdx[0]--;
				triData.normalIdx[1]--;
				triData.normalIdx[2]--;

				triangleList.push_back(triData);
			}
			else
			{
				int readChar;
				while ((readChar = fgetc(file)) != '\n' && readChar != EOF);
			}
		}			
		fclose(file);
		return true;
	}
}

bool MeshData::SaveObjFile(string path)
{
	//const  dtVertex   *vertex = model->vertex;
	//const  dtVector   *normvec = model->normvec;
	//const  dtTriangle *triangle = model->triangle;
	//dt_index_type ind;

	FILE *fd = fopen(path.c_str(), "w");
	if (fd != NULL)
	{
		/* save vertex information */
		for (int  i = 0; i < NumberOfVertices(); i++)
		{
			fprintf(fd, "v   %12.9f   %12.9f   %12.9f\n",
				vertexList[i].x(), vertexList[i].y(), vertexList[i].z());
		}

		/* save normal vector information */
		for (int i  = 0; i < NumberofNormals(); i++)
		{
			fprintf(fd, "vn   %12.9f   %12.9f   %12.9f\n",
				normalList[i].x(), normalList[i].y(), normalList[i].y());
		}

		/* save triangular unit information

		Notice that vertex/normvec indexes in .obj file are one-based, so
		they need a incrementation before being written to filestream. */
		for (int i = 0; i < NumberOfTriangles(); i++)
		{
			fprintf(fd, "f %d//%d %d//%d %d//%d\n",
				triangleList[i].vertexIdx[0] + 1, triangleList[i].normalIdx[0] + 1,
				triangleList[i].vertexIdx[1] + 1, triangleList[i].normalIdx[1] + 1,
				triangleList[i].vertexIdx[2] + 1, triangleList[i].normalIdx[2] + 1 );
		}

		fclose(fd);
		return true;
	}
	else 
	{
		return false;
	}
}

