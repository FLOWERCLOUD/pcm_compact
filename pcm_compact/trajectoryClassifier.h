#ifndef _TRAJECTORY_CLASSIFIER_H
#define _TRAJECTORY_CLASSIFIER_H

#include <set>
#include "sample_set.h"
//#include "DeformaleRegistration.h"

#define SAVE_CORRESPONDENCE
#define SAVE_LABELS

class TrajectoryClassifier
{
public:
	TrajectoryClassifier(IndexType cFrame);
	void run() ;
public:
	void derive_rotation_by_svd(VecX& rot,const MatrixX3 &X,  MatrixX3& Y,MatrixXXi& vtx_map);
	void bubleSort(std::vector<IndexType>& oriData,std::vector<IndexType>& labels,IndexType size);
	void bubleSort(std::vector<IndexType>& oriData,std::vector<ScalarType>& labels,IndexType size);
	void visDistor();
	int orderLabels(std::vector<IndexType>& labels);
	void  setParamter(IndexType _trajLen,IndexType _octreeReso,ScalarType _perC,
		ScalarType _thresHold,IndexType _modelT, IndexType _smallL,bool _isEqual,bool _isRigid);
	void setNeigNum(IndexType _neigbNum);
	void setCenterFrame(IndexType _centerNum);
private:
	void derive_rotation_by_svd( VecX& rot,const MatrixX3 &X, const MatrixX3& Y);
	void diff_using_bfs( std::vector<IndexType>& labels,std::vector<IndexType>& centerVtxId,IndexType centerFrame );

	IndexType centerFrame;
	IndexType trajLen;
	IndexType octreeRes;
	ScalarType perC;
	ScalarType threshold;
	IndexType modelT;
	IndexType lifeT;

	bool isEqual;
	bool isRigid;
	IndexType neigborNum;

};
#endif