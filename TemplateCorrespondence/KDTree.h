#pragma once
#include <vector>
#include <Eigen/Core>
#include "UtilityFunctions.h"
using namespace std;
using namespace Eigen;

namespace DeformationTransfer
{
	namespace Utility
	{
		//enum SplitDimension
		//{
		//	dimX = 0,
		//	dimY, 
		//	dimZ
		//};

		struct NodeInfo
		{
			int splitDim;
			int elementIdx;
			Vector3f element;

			NodeInfo(int idx, Vector3f elem, int currentDim)
				: elementIdx(idx), element(elem), splitDim(currentDim){}
		};
		class KDNode
		{
		public:
			NodeInfo data;
			KDNode* leftTree;
			KDNode* rightTree;
		
			KDNode(int idx, Vector3f element, int cd);
		};
		class KDTree
		{
			const int dimK; 
			 
			
			KDNode* insert(const Vector3f& element ,const int elemIdx, KDNode* t, int cd);
			bool SearchElement(const Vector3f & elem, KDNode * topNode, NodeInfo & result, int cd);
			void NearestNeighbor(KDNode * topNode, const Vector3f & x, BoundingRect & boundRect, KDNode** nearest, float & distance, std::function<bool(KDNode*)> conditionOperator);
		public:			
			KDNode* root;
			
			KDTree(int k);
			~KDTree();
			void BuildTree(vector<Vector3f> listOfElements);
			int RangeSearch(KDNode * topNode, const Vector3f & x, float & range, vector<KDNode*>& result, vector<float>& distance, std::function<bool(KDNode*)> conditionOperator);
		};
	}
}


