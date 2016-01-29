#ifndef _GCOP_H
#define _GCOP_H
//------STL
#include <vector>
#include <unordered_map>
#include <map>
#include <set>

#include <basic_types.h>
#include "vertex.h"
#include "color_table.h"


#define  ONLY_THREE
//#define FIT_PLAN
//#define  USING_CUR

#define EXPSWAP 
class GraphNodeCtr;
struct GraphCutNode;
class DualwayPropagation;
class GCoptimizationGeneralGraph ;
class Sample;

class GCop 
{


public:
	void run();

public:
;

public:
	GCop();
	~GCop();

public:
	void emAlgoritm(std::vector<GraphCutNode*>& oriData);

	void refineSegm();

	void setNeighbor();

	void setNeighborSmpTree();

	void setDataItem();

	void setDataItem(MatrixXX&totError,std::unordered_map<IndexType,std::set<IndexType> >& edgePoints );
	void setSmoothItem();

	void visulizationLabels();

	void setParamter(IndexType _nlabel,IndexType _nExpansion,IndexType _nSwap,IndexType _nGNeig,
		IndexType _nCNeig,ScalarType _nTraOf,ScalarType _CSigma,ScalarType _DSigma,
		ScalarType _nDiffu,IndexType centerF);
public:
	//co-segmentation 20150528
	void coSegmentation();
	void visualCosegmentation(char *labels_file);
	void visualCosegOriPointCloud(char *labels_file);

public:
	void ransacAlgorithm(std::vector<GraphCutNode*>& oriData);
	void ransacAllAlgorithm(std::vector<GraphCutNode*>& oriData);
	void initPlan(Matrix3X & oriCoor,PointType & center,NormalType & planNorm,std::set<IndexType>& inliers);
	void calculateCenterNorm(Matrix3X& smpCoor,PointType& center,NormalType & planNorm);
	void realRansac(Matrix3X& oriCoor,PointType& center,NormalType& planNorm,std::set<IndexType>& inliers,ScalarType noiseR,IndexType iterNum);
	void dist2Plan(Matrix3X& oriCoor,PointType& center,NormalType& planNorm,std::set<IndexType>& inliers,ScalarType noiseR);
	void getSmpCoor(Matrix3X&smpCoor,Matrix3X & oriCoor,std::set<IndexType>& inliers);
	void setInliers(std::map<IndexType,IndexType>& id2Id,std::set<IndexType>& inliers);
	void drawPlan(PointType&center,NormalType&planNorm);
	//fits a  gauss distribution
	void initDistribution(MatrixXX& trans,MatrixXX& miu,MatrixXX& sigma,std::set<IndexType>&inliers);
	void calculateMeanVariance(MatrixXX& trans,MatrixXX& miu,MatrixXX& sigma);
	void gaussRansac(MatrixXX& trans,MatrixXX& miu,MatrixXX& sigma,std::set<IndexType>&inliers,ScalarType noiseR,IndexType itNum);
	void dist2Distribution(MatrixXX& trans,MatrixXX& miu,MatrixXX& sigma,std::set<IndexType>&inliers,ScalarType noiseR);
	void getSmpTrans(MatrixXX & smpTrans,MatrixXX&oriTrans,std::set<IndexType>& inliers);
	//
	void findEdgePoint(std::unordered_map<IndexType, std::set<IndexType>> & edgePoints);
	void visualEdgePoints(std::unordered_map<IndexType, std::set<IndexType>> & edgePoints);
	void calculateEdgeMean(std::unordered_map<IndexType, std::set<IndexType>> & edgePoints,std::unordered_map<IndexType,PointType>&edgeMeanCoor);
	// ransacRT---5-20
	void ransacRotTan(std::vector<GraphCutNode*>& oriData);
	void initRotTrans(Matrix3X& trans,MatrixXX& miu,Matrix33& sigma,std::set<IndexType>&inliers,Matrix3X& tCoor);
	void getCorCoor(PointType& pCoor,IndexType gId);
	void rotRansac(Matrix3X& trans,MatrixXX& miu,Matrix33& sigma,std::set<IndexType>&inliers,ScalarType noiseR,
		IndexType itNum,std::set<IndexType>& graphId,Matrix3X& tCoor);
	void dist2Transformation(Matrix3X& srCoor,MatrixXX& miu,Matrix33& sigma,std::set<IndexType>&inliers,ScalarType noiseR,std::set<IndexType>& graphId,Matrix3X& tCoor);
	void getSmpCorCoor(Matrix3X& smpCoor,std::set<IndexType>&inliers,Matrix3X& tCoor);

	void ransacMultiRotTan(std::vector<GraphCutNode*>& oriData);//累计的旋转平移
	void ransacMultiRotTan(std::vector<GraphCutNode*>& oriData,MatrixXX & totDist);//重载
	void smoothError(MatrixXX& totDist,std::map<IndexType,IndexType>& id2Id);//滤波误差并且设置数据项

public:
	void productionSamples(std::vector<GraphCutNode*>& oriData,MatrixXX & samples);

	void calculateTrans(GraphCutNode* node,MatrixXX & trans);
	void point2point(Matrix3X & srCloud,Matrix3X & tgCloud,Matrix33 & rotMat,MatrixXX & transVec);
	void point2point(Matrix3X & srCloud,Matrix3X & tgCloud,Matrix33 & rotMat,MatrixXX & transVec,Vec3&sMean,Vec3&tMean);
	//
	void calculateMean(MatrixXX& sample_,MatrixXX & mean_,MatrixXX& cov_0,MatrixXX& cov_1,MatrixXX& pro,MatrixXX& labels,std::vector<IndexType>& nLabels);
	//void calculateMean(MatrixXX& sample_,MatrixXX & mean_,CvMat* cov_1,CvMat* cov_2,MatrixXX& pro,MatrixXX& labels,std::vector<IndexType>& nLabels);
	void calculateCov(MatrixXX& sample_,MatrixXX& cov_);
	void calculatePro(MatrixXX& sample_,MatrixXX& mean_,MatrixXX& cov_,MatrixXX pro_);
	void getInitLabels(std::vector<GraphCutNode*>& oriData_,std::vector<IndexType> & labels_,std::vector<IndexType>& nLabels);
	ScalarType calculateProb(IndexType nDim,MatrixXX& mean_,MatrixXX& cov,MatrixXX& sValue);
	void rotMat2EAngle(Matrix33& rotMat,Vec3& eulerAg);
	ScalarType rationCurNode(IndexType frame,IndexType index);
	void smoothSampleLabel(Sample& oriPC,std::vector<IndexType>& sampleVtxId,std::vector<IndexType>& label_smp,std::vector<IndexType>& label_smooth);
	void propagateLabel2Orignal(Sample& oriPC,std::vector<IndexType>& sampleVtxId,std::vector<IndexType>& label_smp,std::vector<IndexType>& label_ori);
	void diff_using_bfs( std::vector<IndexType>& labels,std::vector<IndexType>& centerVtxId,IndexType centerFrame );
	IndexType orderLabels(std::vector<IndexType>& labels);
	void set2DataItem(std::vector<GraphCutNode*>& oriData);
	ScalarType dist2Mean(MatrixXX& _sample,MatrixXX& _mean,IndexType Itid,IndexType LabelId);

	void orderLabelsOnly();
	void bubleSort(std::vector<IndexType>& oriData,std::vector<IndexType>& labels,IndexType lSize);
	//void diff_using_bfs( std:: std::vector<IndexType>& labels,std::std::vector<IndexType>& centerVtxId,IndexType centerFrame );

	void mergeFile(char* CORR_FILESNAME ,char* CORR_FILEOUT_NAME , char* LABEL_FILESNAME , char* LABEL_FILEOUT_NAME);
	void mergeFile();
	void appendcorr( const char* file_in ,char* file_out);
	void appendlabel( const char* file_in ,char* file_out);

	void splitProcess(DualwayPropagation& dp_solver);
	void cosegProcessing(DualwayPropagation& dp_solver);
	void mergeProcess(DualwayPropagation& dp_solver);

	void setCurSmpId(IndexType curId){m_centerF = curId;}

private:
	GraphNodeCtr* gcNode;
	GCoptimizationGeneralGraph *m_gc;
	float * m_data;
	float * m_smooth;
	std::vector<bool> m_isSelect;
	std::vector<IndexType> m_smpId;
	//IndexType m_nLabels;
	IndexType m_nDim;
	std::vector<IndexType> m_optLabel;

	std::vector<bool> m_inliers;
	IndexType m_curFId;//当前选中帧
	std::vector<bool> m_isEdgePoint;//用来可视化边界点

	//parameter from Ui
private:
	IndexType m_nLabels;
	IndexType m_nExpansion;
	IndexType m_nSwap;
	IndexType m_nGraphNeig;
	IndexType m_nCurNeig;
	ScalarType m_nTradeOff;
	ScalarType m_cSigma;
	ScalarType m_dSigma;
	ScalarType m_nDiffu;
	IndexType  m_centerF;
};

#endif // !_GCOPTIMIZATION_H
