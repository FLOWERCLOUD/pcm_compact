#ifndef _MULTIWAY_PROPAGATE_H
#define _MULTIWAY_PROPAGATE_H

//STL
#include <queue>
#include <vector>

#include "basic_types.h"
#include "BoostGraph.h"
#include "pool_allocator.h"

#define  COUT_DEBUG 0

class GCoptimizationGeneralGraph;
class Sample;



class DualwayPropagation
{
public:

	DualwayPropagation();
	~DualwayPropagation();

	void propagate_back_with_link(IndexType f);

	void propagate_back(IndexType f);

	void propagate_front(IndexType f);

	void read_data(char *label_name,char *corr_name);


	void init_labeles_graph();

	void getEdgeVertexs( IndexType _CFrameId ,IndexType lLabelId ,IndexType rLabelId , std::map<IndexType, std::map<IndexType ,HVertex*> >& _edgepoints );
	void getEdgeVertexs2( IndexType _CFrameId , distanPriQueue& _PriQuemap, std::map<IndexType, std::map<IndexType ,HVertex*> >& _edgepoints );
	void init_labeles_graph_hier(ScalarType distThr);


	void read_label_file(char *filename);

	void read_corres_file(char *filename);

	void  read_label_file_hier(char *filename);



	void smooth_label(IndexType frame_idx);

	void write_label_file(std::string filename);

	void wirteSplitGraphLables(std::string filename);

	void compute();

	void buildSmapleKDtree(CFrame* smpCloud, Sample* smp);

	void point2point(Matrix3X & srCloud,Matrix3X & tgCloud,Matrix33 & rotMat,MatrixXX & transVec);

	void buildKdTree(IndexType _cframeid);

	//0821
	void split_twoAjacent_graph_next(IndexType srFrame, IndexType tgFrame);

	void split_twoAjacent_graph_prev(IndexType srFrame, IndexType tgFrame);

	void split_nest_graph_prev(IndexType startFrame,IndexType srFrame, IndexType tgFrame);

	void splitAllSquenceGraph(IndexType iterN);

	void getNextCorVtx(std::map<IndexType,HVertex*>& edgePs, std::map<IndexType,HVertex*>& edgePsCor);

	IndexType checkNextLabelBucket(std::map<IndexType,HVertex*>& edgePs, std::map<IndexType,HVertex*>& edgePsCor);

	bool checkPrevLabelBucket(std::map<IndexType,HVertex*>& edgePs, std::map<IndexType,HVertex*>& edgePsCor, IndexType& isSplit);

	void minDistBeTwoParts(IndexType cFrame,std::map<IndexType,HVertex*>& fPart,std::map<IndexType,HVertex*>& sPart, ScalarType& misDis);

	void determinateUnmarkPoints(IndexType cFrame,std::map<IndexType,HVertex*>& unMarkPs,std::vector<HLabel*> oriLabelBucket,IndexType nodeId,IndexType newLabe,IndexType tgSize);

	ScalarType p2PatchAvgDis(IndexType cFrame, PointType& pCoor,std::map<IndexType,HVertex*>& parthPs);//随机点的距离

	ScalarType p2PatchMinDis(IndexType cFrame, PointType& pCoor,std::map<IndexType,HVertex*>& parthPs);//所有点的最小距离

	void copyLabelBucket(std::vector<HLabel*>& leftLabels, const std::vector<HLabel*>& oriLabels);

	void checkPsNewLabelParentPtr(std::vector<HLabel*> oriLabelBucket,IndexType labParSize);

	void smoothAfterSplit();

	void smoothSingleFrame(IndexType frameId);

	void buildMainPatchMatching();

	void buildNextPatchCorr(IndexType srFrame,IndexType tgFrame);

	void buildPrevPatchCorr(IndexType srFrame,IndexType tgFrame);

	void buildSquenceUniqLabel(std::string filename);

	void writeSquenceLabel(std::map<IndexType,IndexType>& patchTrajId,std::string filename);

	std::map<IndexType,HFrame>* getCompents(){return &hier_componets_;}

	void read_label_file_coseg(char *label_file_name);

	void read_corr_file_coseg(char* corr_file_name);

	void buildPatchCorrespondenceByLabel();

	void mergePatchesAfterCoSeg();

	void findBestPatches(IndexType srFrame, IndexType tgFrame, set<IndexType>& srBestSet, set<IndexType>& tgBestSet); 

	void mergeSingleTinyPatches(IndexType vSize);

	//split by ordered edges 0903
	void split_twoAjacent_graph_next_order(IndexType srFrame, IndexType tgFrame);

	void split_twoAjacent_graph_prev_order(IndexType srFrame, IndexType tgFrame);

	void generateOrderededges(IndexType srFrame, IndexType tgFrame);

	void constructPCloudGraph();

	void addGraphVertex(PCloudGraph& pcGraph, IndexType frameId);

	void addGraphEdge(PCloudGraph& pcGraph, IndexType frameId);

	ScalarType p2PatchGeoDis(IndexType cFrame, HVertex& oriP, std::map<IndexType,HVertex*>& parthPs);

	//
	void generateOrderPrevEdges(IndexType srFrame, IndexType tgFrame);

	void show_corresponding(IndexType f);

	void mergePatchTraj();

	void generTrajNodes(std::vector<PatchTraj>& pNodes);// 生成轨迹节点

	void graphCuts(std::vector<PatchTraj>& pNodes, std::vector<IndexType>& labels);

	void mergeSquences(std::vector<IndexType>& labes);

	void setSegNeihbor(std::vector<PatchTraj>& pNodes, GCoptimizationGeneralGraph& segGraphC);

	void setSegDataItem(GCoptimizationGeneralGraph& segGraphC);

	void setSegSmoothItem(GCoptimizationGeneralGraph& segGraphC);

	bool isAdjInSeq(PatchTraj& nodeA, PatchTraj& nodeB);

	void getSegLabels(GCoptimizationGeneralGraph& segGraphC, std::vector<IndexType>& labels);

	ScalarType motionSimilarityBetw2Nodes(IndexType i, IndexType j, std::vector<PatchTraj>& oriTraj);// i connect j

	void calculateTrans(IndexType lab,IndexType sFrame, IndexType tFrame, Matrix34& toTrans, Matrix34& backTrans);

	void getCoorByVtxBucket(IndexType lab, IndexType frameId, Matrix3X& vtxCoor);

	void initGraphAfterCoseg();

	void calculateSimilar2Componet();

	void propagateLabel2Orignal(Sample& oriPC,std::vector<IndexType>& sampleVtxId,std::vector<IndexType>& label_smp,std::vector<ScalarType>& label_ori);

public:

	Sample* downSample;
	std::vector<Eigen::MatrixXd> patchSimilarMea;

private:
	std::map<IndexType, CFrame> components_;
	std::map<IndexType,IndexType> frame_index_map_label_; 
	static  PoolAllocator allocator_;

	std::map<IndexType, HFrame> hier_componets_;

	std::map<IndexType,IndexType> o2d_map; // origin to down
	std::map<IndexType,IndexType> d2o_map;

public:
	OrderEdgeQueue orderedEdgeQ;
};


#endif