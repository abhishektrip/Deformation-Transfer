#include <iostream>
#include "CorrespondenceProblem.h"
#include "KDTree.h"
//#include "MeshData.h"

using namespace DeformationTransfer::Correspondence;

int main()
{
	CorrespondenceProblem problem; 
	problem.wtSmoothness = 1.0;
	problem.wtIdentity = 0.01;
	problem.wtClosestStart = 1; 
	problem.wtClosestEnd = 10;
	problem.wtClosestStep = 10;

	//if (!problem.srcMesh.ReadObjFile("C:\\Users\\abhis\\Documents\\3DStyleSwitching_CPP\\TemplateCorrespondence\\TestCases\\Src.obj"))
	//	cout << "Src load Failed" << endl;
	//if (!problem.tgtMesh.ReadObjFile("C:\\Users\\abhis\\Documents\\3DStyleSwitching_CPP\\TemplateCorrespondence\\TestCases\\Tgt.obj"))
	//	cout << "Src load Failed" << endl;

	//problem.meshCons.ReadConstraintFile("C:\\Users\\abhis\\Documents\\3DStyleSwitching_CPP\\TemplateCorrespondence\\TestCases\\Src_Tgt_Cons.txt");

	problem.srcMesh.ReadObjFile("C:\\Users\\abhis\\Documents\\3DStyleSwitching\\TemplateCorrespondence\\Assets\\Models\\horse_ref.obj");
	problem.tgtMesh.ReadObjFile("C:\\Users\\abhis\\Documents\\3DStyleSwitching\\TemplateCorrespondence\\Assets\\Models\\camel_ref.obj");

	problem.meshCons.ReadConstraintFile("C:\\Users\\abhis\\Documents\\3DStyleSwitching\\TemplateCorrespondence\\Assets\\Models\\Horse_Camel_Corres.txt");

	//std::cout << "Vertices : " << problem.srcMesh.NumberOfVertices() << endl;
	//std::cout << "Triangles : " << problem.srcMesh.NumberOfTriangles() << endl;

	//std::cout << "Vertices : " << problem.tgtMesh.NumberOfVertices() << endl;
	//std::cout << "Triangles : " << problem.tgtMesh.NumberOfTriangles() << endl;

	//std::cout << "Constraints : " << problem.meshCons.numConstraints << endl;

	problem.Solve();

	//DeformationTransfer::Utility::KDTree kdTree(3);
	//kdTree.BuildTree(problem.srcMesh.VertexList());

	int n;
	cin >> n;
}
