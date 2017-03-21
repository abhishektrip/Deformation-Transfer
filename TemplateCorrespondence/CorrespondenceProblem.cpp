#include "CorrespondenceProblem.h"
#include <iostream>
#include <Eigen/Sparse>
#include "UtilityFunctions.h"
using namespace DeformationTransfer::Correspondence;


CorrespondenceProblem::CorrespondenceProblem()
{
}


CorrespondenceProblem::~CorrespondenceProblem()
{
}
void CorrespondenceProblem::SolveEquations(TripletList& tripA, VectorXf& C , int nRow, int nCol)
{
	SparseMat A(nRow, nCol);
	/*for(int i = 0 ;i< nRow;i++)
		std::cout << "["<<i<<"] : " << C[i] << endl;*/
	A.setFromTriplets(tripA.begin(), tripA.end());
	/*for(auto iter = tripA.begin(); iter != tripA.end(); iter++)
		cout << (*iter).row() << "\t" << (*iter).col() << "\t" << (*iter).value() <<endl;*/
	//cout<<"Sparse Matrix" <<endl << A << endl;
	//std::cout << "Computing Transpose " << endl;
	SparseMat Atrans = A.transpose();
	VectorXf b = Atrans*C;

	SparseMat AtA = Atrans*A;

	std::cout << "Solving  Equations " << endl;
	AtA.makeCompressed();
	Eigen::SparseLU<SparseMatrix<float, ColMajor>, COLAMDOrdering<int> > solver;
	solver.analyzePattern(AtA);
	solver.factorize(AtA);
	if (solver.info() != Success)
	{
		cout << "Factorize failed - info = " << solver.info()<< endl;
	}
	else
	{
		VectorXf x = solver.solve(b);
		std::cout << "Finished Solving Equations " << endl;
		//std::cout << "Result Mat \n" << x;
		ApplyDeformationToSourceMesh(x);
		srcMesh.SaveObjFile("Test.obj");
	}
	/*Eigen::SimplicialCholesky<SparseMat> chol;
	chol.compute(AtA);
	if (chol.info() != Success)
	{
		cout << "Decompostion failed. " << endl;
	}
	else
	{
		VectorXf x = chol.solve(b);
		std::cout << "Finished Solving Equations " << endl;
		std::cout << "Result Mat \n" << x;
		ApplyDeformationToSourceMesh(x);
		srcMesh.SaveObjFile("Test.obj");
	}*/
}
void CorrespondenceProblem::Solve()
{
	srcMesh.CreateAdjacencyList();
	std::cout << "Computing the vertex information" << endl;
	vInfo.PopulateVertexInformationList(*this);

	std::cout << "Solving Phase 1 " << endl;
	SolvePhase1();
	SolvePhase2();


	float srcThresh = srcMesh.ComputeTriangleCorrespondenceThreshold();
	float tgtThresh = tgtMesh.ComputeTriangleCorrespondenceThreshold();

	float thresh = max(srcThresh, tgtThresh);
	ResolveTriangleCorrespondence(thresh);
}
void CorrespondenceProblem::SolvePhase1()
{
	int nRow = 9 * (srcMesh.triAdj.adjacentTriangleCount + srcMesh.NumberOfTriangles());
	int nCol = 3 * (vInfo.numberFree + srcMesh.NumberOfTriangles());

	TripletList tripA; 
	VectorXf C = VectorXf::Zero(nRow);
	std::cout << "Building Equations " << endl;
	BuildPhase1Equation(tripA, C);

	std::cout << "Building Sparse Matrix " << endl;
	SolveEquations(tripA, C, nRow, nCol);

	//SparseMat A(nRow, nCol);
	////for(int i = 0 ;i< nRow;i++)
	////	std::cout << "["<<i<<"] : " << C[i] << endl;
	//A.setFromTriplets(tripA.begin(), tripA.end());

	////cout<<"Sparse Matrix" <<endl << A << endl;
	////std::cout << "Computing Transpose " << endl;
	//SparseMat Atrans = A.transpose();
	//VectorXf b = Atrans*C;

	//SparseMat AtA = Atrans*A;
	////
	////MatrixXf MtransC = Mtrans*C;

	//std::cout << "Solving  Equations " << endl;	

	////AtA.makeCompressed();
	////Eigen::SparseLU<SparseMatrix<float, ColMajor>, COLAMDOrdering<int> > solver;
	////solver.analyzePattern(AtA);
	////solver.factorize(AtA);
	////if (solver.info() != Success)
	////{
	////	cout << "Factorize failed - info = " << solver.info()<< endl;
	////}
	////else
	////{
	////	VectorXf x = solver.solve(b);
	////	std::cout << "Finished Solving Equations " << endl;
	////	//std::cout << "Result Mat \n" << x;
	////	ApplyDeformationToSourceMesh(x);
	////	srcMesh.SaveObjFile("Test.obj");
	////}
	//
	//

	///*Eigen::LeastSquaresConjugateGradient<SparseMat> lscg;
	//lscg.compute(M);
	//if (lscg.info() != Success)
	//{
	//	cout << "LSCG failed. " << endl;
	//}
	//else
	//{
	//	VectorXf x = lscg.solve(C);
	//	std::cout << "Finished Solving Equations " << endl;
	//	ApplyDeformationToSourceMesh(x);
	//	srcMesh.SaveObjFile("Test.obj");
	//}*/



	//Eigen::SimplicialCholesky<SparseMat> chol;
	//chol.compute(AtA);
	//if (chol.info() != Success)
	//{
	//	cout << "Decompostion failed. " << endl;
	//}
	//else
	//{
	//	VectorXf x = chol.solve(b);
	//	std::cout << "Finished Solving Equations " << endl;
	//	ApplyDeformationToSourceMesh(x);
	//	srcMesh.SaveObjFile("Test.obj");
	//}
}
void CorrespondenceProblem::SolvePhase2()
{
	cout << "Phase 2 Starting..." << endl;

	int nRow = 9 * (srcMesh.triAdj.adjacentTriangleCount + srcMesh.NumberOfTriangles()) + 3 * (vInfo.numberFree);
	int nCol = 3 * (vInfo.numberFree + srcMesh.NumberOfTriangles());

	cout << "Rows : " << nRow << " Cols : " << nCol << endl;
	TripletList tripA;
	VectorXf C;
	for (float wt = wtClosestStart; wt < wtClosestEnd; wt += wtClosestStep)
	{
		//Debug.Log("Iteration Step " + wt);
		//Debug.Log("Creating Correspondence map..");
		vector<int> correspondenceMap;

		CreateCorrespondenceMap_Basic(srcMesh, tgtMesh, correspondenceMap);
		tripA.clear();
		C = VectorXf::Zero(nRow);
		//A = Matrix<float>.Build.Sparse(nRow, nCol, 0);
		//C = Vector<float>.Build.Dense(nRow, 1);

		cout << "Building Equations" << endl;
		BuildPhase2Equation(tripA, C, correspondenceMap, wt);

		cout << "Linear System Solve..." << endl;

		//Debug.Log("A : " + A.ToString());
		//Debug.Log("C : " + C.ToString());
		SolveEquations(tripA, C, nRow, nCol);				
	}

}
int  CorrespondenceProblem::BuildPhase1Equation(TripletList& tripA, VectorXf& C)
{
	int rowIdx = 0;
	srcMesh.MatricesForTriangle.clear();
	srcMesh.LinearEquationForTriangle.clear();
	std::cout << "Computing Vertex Martices" << endl;
	srcMesh.ComputeVertexMatricesForTriangles();
	std::cout << "Computing A & C Matrices " << endl;
	srcMesh.ComputeLinearEquationMatricesForTriangles(*this);

	rowIdx = CreateSmoothnessEquation(tripA, C, rowIdx);
	cout << "RowIdx after smoothness : " << rowIdx << endl;
	rowIdx = CreateIdentityEquation(tripA, C, rowIdx);
	return rowIdx;
}
int  CorrespondenceProblem::BuildPhase2Equation(TripletList& tripA, VectorXf& C, vector<int>& correspondenceMap, float wtClosestPoint)
{
	int rowIdx = 0;

	rowIdx = BuildPhase1Equation(tripA, C);
	cout << "RowIdx after Phase 1 : " << rowIdx << endl;
	rowIdx = CreateClosestPointEquation(tripA, C, correspondenceMap, rowIdx, wtClosestPoint);

	return rowIdx;
}
int  CorrespondenceProblem::MatrixPositionForFreeVertex(int vertIdx, int componentIdx)
{
	int position = -1;
	if (!vInfo.info[vertIdx].isConstrained)
		position = vInfo.info[vertIdx].matrixPosition * 3 + componentIdx;

	return position;
}
int  CorrespondenceProblem::MatrixPositionForFourthVertex(int triIdx, int vertIdx, int componentIdx)
{
	int position = -1;
	if (vertIdx == 3)
	{
		position = (vInfo.numberFree + triIdx) * 3 + componentIdx;
	}
	return position;
}
int  CorrespondenceProblem::MatrixPositionForTriangle(int triIdx, int vertIdx, int componentIdx)
{
	int vIdx = 0;
	int position = -1;
	if (vertIdx < 3)
	{
		vIdx = srcMesh.TriangleList()[triIdx].vertexIdx[vertIdx];
		position = MatrixPositionForFreeVertex(vIdx, componentIdx);
		//Debug.Log(String.Format("Vert index : {0} -> position : {1}", vIdx, position));
	}
	else if (vertIdx == 3)
	{
		position = MatrixPositionForFourthVertex(triIdx, vertIdx, componentIdx);
	}
	return position;
}

void CorrespondenceProblem::AddToLinearEquationSystemSparse_Unrolled(TripletList& A, int triIdx, int row, float weight)
{
	int iRow = row;
	int vertIdx[] = { 0, 1, 2, 3 };
	int compIdx[] = { 0, 0, 0, 1, 1, 1, 2, 2, 2 };
	for (int jRow = 0; jRow < 9; jRow++, iRow++)
	{
		int position = MatrixPositionForTriangle(triIdx, vertIdx[0], compIdx[jRow]);
		if (position != -1)
		{
			float value = weight * srcMesh.LinearEquationForTriangle[triIdx].A(jRow, vertIdx[0]);
			TripletF temp(iRow, position, value);
			A.push_back(temp);
		}

		position = MatrixPositionForTriangle(triIdx, vertIdx[1], compIdx[jRow]);
		if (position != -1)
		{
			//A[iRow, position] = weight * srcMesh.matricesForTriangle[triIdx].A[jRow, vertIdx[1]];
			float value = weight * srcMesh.LinearEquationForTriangle[triIdx].A(jRow, vertIdx[1]);
			TripletF temp(iRow, position, value);
			A.push_back(temp);
		}

		position = MatrixPositionForTriangle(triIdx, vertIdx[2], compIdx[jRow]);
		if (position != -1)
		{
			//A[iRow, position] = weight * problem.srcMesh.matricesForTriangle[triIdx].A[jRow, vertIdx[2]];
			float value = weight * srcMesh.LinearEquationForTriangle[triIdx].A(jRow, vertIdx[2]);
			TripletF temp(iRow, position, value);
			A.push_back(temp);
		}
		position = MatrixPositionForTriangle(triIdx, vertIdx[3], compIdx[jRow]);
		if (position != -1)
		{
			//A[iRow, position] = weight * problem.srcMesh.matricesForTriangle[triIdx].A[jRow, vertIdx[3]];
			float value = weight * srcMesh.LinearEquationForTriangle[triIdx].A(jRow, vertIdx[3]);
			TripletF temp(iRow, position, value);
			A.push_back(temp);
		}
	}
}
void CorrespondenceProblem::AddToLinearEquationSystemDence_Unrolled(VectorXf& C, int triIdx, int row, float weight)
{
	int iRow = row;
	int jRow = 0;
	for (int i = 0; i < 9; i++, iRow++, jRow++)
	{
		C[iRow] = C[iRow] + weight * srcMesh.LinearEquationForTriangle[triIdx].C[jRow];
	}
}
void CorrespondenceProblem::AddToLinearEquationSystemDence_Unrolled(VectorXf& C, Vector91f& c, int row, float weight)
{
	int iRow = row;
	int jRow = 0;
	for (int i = 0; i < 9; i++, iRow++, jRow++)
	{
		C[iRow] = C[iRow] + weight * c[jRow];
	}
}

int  CorrespondenceProblem::CreateSmoothnessEquation(TripletList& tripA, VectorXf& C, int row)
{
	int matRow = row;
	//srcMesh.TriangleList();
	int numTris = srcMesh.NumberOfTriangles();
	int triangleIndex = 0;
	float sqrtWt = sqrtf(wtSmoothness);
	for (triangleIndex = 0; triangleIndex < numTris; triangleIndex++)
	{
		unordered_set<int>& adjTriangles = srcMesh.triAdj.adjacencyList[triangleIndex];
		int adjTrisCt = adjTriangles.size();
		if (adjTrisCt > 0)
		{
			
			for (auto tri = adjTriangles.begin(); tri != adjTriangles.end() ; tri++)
			{
				//cout<< "Triangle Pair #1 : "<<triangleIndex << endl;
				//cout << "Triangle Pair #2 : " << *tri << endl;
				AddToLinearEquationSystemDence_Unrolled(C, triangleIndex, matRow, sqrtWt);
				AddToLinearEquationSystemDence_Unrolled(C, *tri, matRow, -sqrtWt);
				
				AddToLinearEquationSystemSparse_Unrolled(tripA, triangleIndex, matRow, sqrtWt);
				AddToLinearEquationSystemSparse_Unrolled(tripA, *tri, matRow, -sqrtWt);
				matRow += 9;
			}
		}
	}
	/*for (int i = 0; i< C.rows(); i++)
		std::cout << "[" << i << "] : " << C[i] << endl;*/
	return matRow;	
}
int  CorrespondenceProblem::CreateIdentityEquation(TripletList& tripA, VectorXf& C, int row)
{
	int matRow = row;
	// srcTriangles = problem.srcMesh.TriangleList;
	int numTris = srcMesh.NumberOfTriangles();
	int triangleIndex = 0;
	float sqrtWt = sqrtf(wtIdentity);

	for (triangleIndex = 0; triangleIndex < numTris; triangleIndex++)
	{
		//Create Identity vector. 
		Vector91f cIdentity;
		cIdentity <<1, 0, 0, 0, 1, 0, 0, 0, 1 ;
		
		// T - I
		for (int iden = 0; iden<9; iden++)
		{
			cIdentity[iden] += srcMesh.LinearEquationForTriangle[triangleIndex].C[iden];
		}
		// get A matrix for triangle. 
		//var a = srcMesh.LinearEquationForTriangle[triangleIndex].A;
		// Add to linear system. 
		AddToLinearEquationSystemSparse_Unrolled(tripA, triangleIndex, matRow, sqrtWt);
		AddToLinearEquationSystemDence_Unrolled(C, cIdentity, matRow, sqrtWt);
		//Utils.UtilityFunctions.AddToLinearEquationSystem(A, C, a, cIdentity, problem, triangleIndex, matRow, sqrtWt);
		matRow += 9;
	}
	return matRow;
}

int CorrespondenceProblem::CreateClosestPointEquation(TripletList& tripA, VectorXf& C, vector<int>& correspondenceMap, int row, float wtClosest)
{

	for (int i = 0; i < srcMesh.NumberOfVertices(); i++)
	{
		if (vInfo.info[i].isConstrained == false)
		{
			int iVx = MatrixPositionForFreeVertex(i, 0);
			int iVy = MatrixPositionForFreeVertex(i, 1);
			int iVz = MatrixPositionForFreeVertex(i, 2);

			int tgtIdx = correspondenceMap[i];

			tripA.push_back(TripletF(row, iVx, wtClosest));
			C[row] = wtClosest * tgtMesh.VertexList()[tgtIdx].x();
			row++;

			tripA.push_back(TripletF(row, iVy, wtClosest));
			C[row] = wtClosest * tgtMesh.VertexList()[tgtIdx].y();
			row++;

			tripA.push_back(TripletF(row, iVz, wtClosest));
			C[row] = wtClosest * tgtMesh.VertexList()[tgtIdx].z();
			row++;
		}
	}
	return row;
}
void CorrespondenceProblem::CreateCorrespondenceMap_Basic(MeshData src, MeshData tgt, vector<int>& correspondenceMap)
{
	correspondenceMap.reserve(src.NumberOfVertices());
	int closestVertIndex = -1;
	float distance = 0.0f, minDistance = LONG_MAX;
	int numSrcVerts = src.NumberOfVertices();
	int numTgtVerts = tgt.NumberOfVertices();
	for (int i = 0; i < numSrcVerts; i++)
	{
		Vector3f srcVertex = src.VertexList()[i];
		Vector3f srcNormal = src.NormalList()[i];
		closestVertIndex = -1;
		for (int j = 0; j < numTgtVerts; j++)
		{
			Vector3f tgtVertex = tgt.VertexList()[j];
			Vector3f tgtNormal = tgt.NormalList()[j];

			if (srcNormal.dot(tgtNormal) > 0)
			{
				distance = DeformationTransfer::Utility::SquaredDistance(srcVertex, tgtVertex);
				if (closestVertIndex == -1 || distance < minDistance)
				{
					minDistance = distance;
					closestVertIndex = j;
				}
			}
		}
		correspondenceMap.push_back(closestVertIndex);
	}	
}
void CorrespondenceProblem::ResolveTriangleCorrespondence(float threshold)
{
	Utility::KDTree centroidTree(3);
	
	//Compute centroids for each triangle.
	srcMesh.ComputeCentroidsForTriangles();
	tgtMesh.ComputeCentroidsForTriangles();

	// build a centroid tree for src mesh. 
	centroidTree.BuildTree(srcMesh.centroidList);
	vector<Utility::KDNode*> resultNodes;
	vector<float> resultDist;
	vector<Vector3f>& srcCentroidList = srcMesh.centroidList;
	MeshData& srcMeshPtr = srcMesh;
	for (int i = 0 ; i < tgtMesh.NumberOfTriangles(); i++)
	{
		Vector3f x0 = tgtMesh.centroidList[i];
		Vector3f tgtNorm = CalculateUnitNorm(&tgtMesh, i);

		int nResult = centroidTree.RangeSearch(centroidTree.root, x0, threshold, resultNodes, resultDist, [i, tgtNorm , &srcMeshPtr ](Utility::KDNode* node)->bool
		{
			Vector3f unitNormSrc = CalculateUnitNorm(&srcMeshPtr, node->data.elementIdx);
			float dotprod = tgtNorm.dot(unitNormSrc);
			return (dotprod> 0);
		});

		for (int j = 0; j < resultNodes.size(); j++)
		{
			ResultContainer res; 
			res.tgtIndex = i; 
			res.distance = resultDist[j];
			
			auto iter = CorrespondenceResultMap.find(resultNodes[j]->data.elementIdx);

			if (iter == CorrespondenceResultMap.end())
			{
				int idx = resultNodes[j]->data.elementIdx;
				vector<ResultContainer> resList{ res };
				CorrespondenceResultMap.insert({ idx, resList});
			}
			else
			{
				iter->second.push_back(res);
			}
		}
		resultNodes.clear();
		resultDist.clear();
	}

	srcMesh.centroidList.clear();
	tgtMesh.centroidList.clear();
}



void CorrespondenceProblem::ApplyDeformationToSourceMesh(VectorXf& deformedX )
{
	for (int i = 0; i < srcMesh.NumberOfVertices(); i++)
	{
		if (vInfo.info[i].isConstrained)
		{
			int tgtIdx = meshCons.constrainedList[i];
			srcMesh.VertexList()[i] = tgtMesh.VertexList()[tgtIdx];
		}
		else
		{
			//Vector3f& deformedVert = srcMesh.VertexList()[i];

			srcMesh.VertexList()[i].x() = deformedX[MatrixPositionForFreeVertex(i, 0)];
			srcMesh.VertexList()[i].y() = deformedX[MatrixPositionForFreeVertex(i, 1)];
			srcMesh.VertexList()[i].z() = deformedX[MatrixPositionForFreeVertex(i, 2)];

			//.x() = deformedVert;
		}
	}
}


