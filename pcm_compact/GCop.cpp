#include "GCop.h"
#include "multiway_propagation.h"
#include "graph_cut_node.h"//default current folder
#include <gco-v3.0/GCoptimization.h>
#include "co_segmentation.h"
#include "merge_write_utilities.h"
#include <stdlib.h>
#include <windows.h>
#include <tchar.h>
extern char* filename;
GCop::GCop()
{
	gcNode = NULL;
	m_gc = NULL;
	m_data = NULL;
	m_smooth = NULL;
	m_nLabels = 2;

	m_nGraphNeig = 10;
	m_nCurNeig = 10;
	m_nExpansion = 1;
	m_nSwap = 1;
	m_nTradeOff = 1;
	m_cSigma = 1.;
	m_dSigma = 1.;
	m_nDiffu = 1.;
	m_centerF = 10;

	if (!m_optLabel.empty())
	{
		m_optLabel.clear();
	}
}

GCop::~GCop()
{
	if (gcNode)
	{
		delete gcNode;
		gcNode = NULL;
	}

	if (m_gc)
	{
		delete m_gc;
		m_gc = NULL;
	}

	if (m_data)
	{
		delete [] m_data;
		m_data = NULL;
	}

	if (m_smooth)
	{
		delete [] m_smooth;
		m_smooth = NULL;
	}

	m_optLabel.clear();
	m_smpId.clear();
	m_isSelect.clear();
	m_isEdgePoint.clear();
}

void GCop::setParamter(IndexType _nlabel,IndexType _nExpansion,IndexType _nSwap,IndexType _nGNeig,IndexType _nCNeig,
	ScalarType _nTraOf,ScalarType _CSigma,ScalarType _DSigma,ScalarType _nDiffu,IndexType centerF)
{
	m_nLabels = _nlabel;
	m_nExpansion = _nExpansion;
	m_nSwap = _nSwap;
	m_nGraphNeig = _nGNeig;
	m_nCurNeig = _nCNeig;
	m_nTradeOff = _nTraOf;
	m_cSigma = _CSigma;
	m_dSigma = _DSigma;
	m_nDiffu = _nDiffu;
	m_centerF = centerF;
}

void GCop::run()
{
	bool isrefine = false;
	if (isrefine)
	{
		Logger<<"Start refine single frame.\n";

		refineSegm();

		Logger<<"End refine single frame.\n";
	}else
	{
		Logger<<"Start Split & Coseg & Merge Process!\n";

		//在初始分割之后, 分别进行split,coseg以及merge过程. 20151009
	
		char* in_label_file = ".\\tmp\\results\\cosegOOrder.txt";//smoothing again--12-13
		char* in_corr_file  = ".\\tmp\\results\\hksingle_corr.txt";


		DualwayPropagation dp_solver;

		dp_solver.read_data(in_label_file,in_corr_file);

		dp_solver.init_labeles_graph_hier(m_nTradeOff); //构建分块对应的Graph;

		dp_solver.constructPCloudGraph();    //构建点云对应的Graph;

		splitProcess(dp_solver);

		Logger<<"End Split Process!\n";

		cosegProcessing(dp_solver);

		//Logger<<"End  Coseg Process!\n";

		mergeProcess(dp_solver);

		//Logger<<"End Merge Process!\n";
	}


}

void GCop::refineSegm()
{

#ifdef ONLY_THREE
	m_nDim = 3;
#else
	m_nDim = 6;
#endif // ONLY_THREE

	gcNode = new GraphNodeCtr();//diff_labels(0) 

	char input_label_file[2048];
	char input_cor_file[2048];

	sprintf_s(input_label_file,".\\tmp\\results\\cosegOOrder%.2d.txt",m_centerF); //in order to smoothig after coseg,2015-12-05
	sprintf_s(input_cor_file,".\\tmp\\results\\hksingle_corr%.2d_0.50.txt",m_centerF);

	m_nLabels = gcNode->readnLabelFile(input_label_file);
	gcNode->read_corres_file(input_cor_file); 
	IndexType nodeNum = (IndexType)gcNode->node_vec.size();
	m_gc = new GCoptimizationGeneralGraph(nodeNum,m_nLabels,gcNode);

	setNeighborSmpTree(); 


#ifdef EXPSWAP

#else
	ransacMultiRotTan(gcNode->node_vec);//利用多帧间的平均误差---用来可视化
#endif // Expansion


	visulizationLabels();
}

void GCop::splitProcess(DualwayPropagation& dp_solver)
{
	char* out_label_file = ".\\tmp\\results\\1223SigsplitResultsSmth.txt";
	if(_access(".\\tmp",0) == -1){
		CreateDirectory(_T(".\\tmp"), NULL);
	}
	if(_access(".\\tmp\\results",0) == -1){
		CreateDirectory(_T(".\\tmp\\results"), NULL);
	}

	dp_solver.splitAllSquenceGraph(0);//读取j-linkagelabel文件之后进行前后的分裂操作,参数表示序列分裂的帧数;

	//dp_solver.smoothAfterSplit(); //k =30,分裂之后进行smooth操作

	//平滑处理之后,有些块点个数变为零,或者个数很小,需要进行合并操作!
	//dp_solver.mergeSingleTinyPatches(30); //remove empty segments 加入了循环操作,

	dp_solver.wirteSplitGraphLables(out_label_file);//可视化合并后的结果
	if(_access(".\\tmp",0) == -1){
		CreateDirectory(_T(".\\tmp"), NULL);
	}
	if(_access(".\\tmp\\ply",0) == -1){
		CreateDirectory(_T(".\\tmp\\ply"), NULL);
	}
	if(_access(".\\tmp\\ply\\split",0) == -1){
		CreateDirectory(_T(".\\tmp\\ply\\split"), NULL);
	}
	generateSamplePly( out_label_file ,".\\tmp\\ply\\split\\",filename);


}

void GCop::cosegProcessing(DualwayPropagation& dp_solver)
{
	char* out_label_file = ".\\tmp\\results\\1205cosegResults.txt";
	if(_access(".\\tmp",0) == -1){
		CreateDirectory(_T(".\\tmp"), NULL);
	}
	if(_access(".\\tmp\\results",0) == -1){
		CreateDirectory(_T(".\\tmp\\results"), NULL);
	}

	CoSegmentation coseg_solver(SampleSet::get_instance(),dp_solver.getCompents());

	coseg_solver.hierComponets2Components();

	coseg_solver.compute();

	coseg_solver.components2HierComponets();

	//dp_solver.init_labeles_graph_hier(0.015); //重新构造一个图结构

	////dp_solver.initGraphAfterCoseg();//不需要构造一个图结构, 因为共分割并没有改变原来的图结构.只需要调整vector<HLabel*>到原来的就可以

	//dp_solver.mergeSingleTinyPatches(50);

	dp_solver.wirteSplitGraphLables(out_label_file);
	if(_access(".\\tmp",0) == -1){
		CreateDirectory(_T(".\\tmp"), NULL);
	}
	if(_access(".\\tmp\\ply",0) == -1){
		CreateDirectory(_T(".\\tmp\\ply"), NULL);
	}
	if(_access(".\\tmp\\ply\\cog",0) == -1){
		CreateDirectory(_T(".\\tmp\\ply\\cog"), NULL);
	}
	generateSamplePly( out_label_file ,".\\tmp\\ply\\cog\\",filename);


}

void GCop::mergeProcess(DualwayPropagation& dp_solver)
{
	char* out_label_file = ".\\tmp\\results\\1205mergeResults.txt";
	if(_access(".\\tmp",0) == -1){
		CreateDirectory(_T(".\\tmp"), NULL);
	}
	if(_access(".\\tmp\\results",0) == -1){
		CreateDirectory(_T(".\\tmp\\results"), NULL);
	}

	dp_solver.mergePatchesAfterCoSeg(); // 读取共分割的label文件之后,建立图结构, 然后进行图的merge操作 0831

	////merge squences by Graph cuts

	//dp_solver.mergePatchTraj();

	dp_solver.wirteSplitGraphLables(out_label_file);//可视化合并后的结果
	if(_access(".\\tmp",0) == -1){
		CreateDirectory(_T(".\\tmp"), NULL);
	}
	if(_access(".\\tmp\\ply",0) == -1){
		CreateDirectory(_T(".\\tmp\\ply"), NULL);
	}
	if(_access(".\\tmp\\ply\\merge",0) == -1){
		CreateDirectory(_T(".\\tmp\\ply\\merge"), NULL);
	}
	generateSamplePly( out_label_file,".\\tmp\\ply\\merge\\",filename);

}

void GCop::setNeighborSmpTree()
{
	IndexType * neighbours = new IndexType[m_nGraphNeig];

	ScalarType * dist = new ScalarType [m_nGraphNeig];

	IndexType fId = (gcNode->node_vec)[0]->frame;

	Sample & curF = gcNode->m_smpSet[fId];

	m_isSelect.resize(curF.num_vertices(),false);

	m_smpId.resize(m_gc->numSites());

	unordered_map<IndexType,bool> neigh_pair;

	IndexType nodeId = 0;

	Sample* downSmp = new Sample;

	for (auto iter = gcNode->node_vec.begin(); iter!= gcNode->node_vec.end(); iter ++,nodeId ++)
	{
		Vertex& vtx = curF[(*iter)->index];

		PointType v( vtx.x(), vtx.y(), vtx.z() );
		ColorType cv(vtx.r(), vtx.g(), vtx.b(), vtx.alpha());
		NormalType nv(vtx.nx(), vtx.ny(), vtx.nz());

		downSmp->add_vertex(v,nv,cv);
	}

	downSmp->build_kdtree();

	nodeId = 0;

	for (auto iter = gcNode->node_vec.begin(); iter!= gcNode->node_vec.end(); iter ++,nodeId ++)
	{
		m_isSelect[(*iter)->index] = true;

		m_smpId[nodeId] = (*iter)->index;

		m_gc->setLabel((*iter)->graph_index,(*iter)->label);

		downSmp->neighbours((*iter)->graph_index,m_nGraphNeig,neighbours,dist);

		for (IndexType neig_id = 1; neig_id < m_nGraphNeig; neig_id ++)
		{
			bool temp = neigh_pair[neighbours[neig_id],(*iter)->graph_index];

			if (!temp)
			{
				m_gc->setNeighbors((*iter)->graph_index,neighbours[neig_id] );

				neigh_pair[(*iter)->graph_index,neighbours[neig_id]] = true;
			}

		}

	}

	delete [] neighbours;

	delete [] dist;

	delete downSmp;

	neigh_pair.clear();

}

void GCop::setNeighbor()
{
	IndexType * neighbours = new IndexType[m_nGraphNeig];
	ScalarType * dist = new ScalarType [m_nGraphNeig];

	IndexType fId = (gcNode->node_vec)[0]->frame;

	Sample & curF = gcNode->m_smpSet[fId];

	m_isSelect.resize(curF.num_vertices(),false);
	m_smpId.resize(m_gc->numSites());

	unordered_map<IndexType,bool> neigh_pair;

	IndexType nodeId = 0;

	vector<ScalarType> all_dist;
	vector<ScalarType> all_curvature;

	for (auto iter = gcNode->node_vec.begin(); iter!= gcNode->node_vec.end(); iter ++,nodeId ++)
	{
		m_isSelect[(*iter)->index] = true;

		m_smpId[nodeId] = (*iter)->index;

		curF.neighbours((*iter)->index,m_nGraphNeig,neighbours,dist);


		m_gc->setLabel((*iter)->graph_index,(*iter)->label);



		for (IndexType neig_id = 1; neig_id < m_nGraphNeig; neig_id ++)
		{
			IndexType pId = neighbours[neig_id];
			GraphCutNode *fNode = gcNode->node_map[frame_index_to_key(fId,pId)];

			if(fNode)
			{
				bool temp = neigh_pair[frame_index_to_key(fNode->graph_index,(*iter)->graph_index)];

				if (!temp)//重复添加了邻域边?
				{

					m_gc->setNeighbors((*iter)->graph_index,fNode->graph_index);//因为在优化能量项中 数据从0下标开始

					neigh_pair[frame_index_to_key((*iter)->graph_index,fNode->graph_index)] = true;


				}

			}
		}
	}

	


	delete [] neighbours;
	delete [] dist;

	neigh_pair.clear();
}

void GCop::visulizationLabels()
{
	IndexType centFrame = (gcNode->node_vec)[0]->frame;
	Sample & curF = gcNode->m_smpSet[centFrame];


	IndexType nNode = (IndexType)gcNode->node_vec.size();
	// init prolabel

	IndexType k = 0;
	for (IndexType smpId = 0; smpId < nNode; smpId ++,k++)
	{
		m_optLabel.push_back(m_gc->whatLabel(k));
	}


	// 	IndexType nIter = 1;
	// 	while (nIter-->0)
	// 	{`
	// 		smoothSampleLabel(curF,m_smpId,m_optLabel,m_optLabel);
	// 	}
	//// 
	//diff_using_bfs(m_optLabel,m_smpId,centFrame);

	IndexType finalLabels = orderLabels(m_optLabel);
	//IndexType finalLabels = 2;

#ifdef  SHOW_SAMPLE

#ifdef OUTPUT_LABELS
	char outputLabelName[1024]; 
	sprintf(outputLabelName,"F:\\EG2015\\compar\\diffusionOrder\\labelAfterCoseg\\boundaryLabels%.2d.txt",centFrame);
	FILE *in_label = fopen(outputLabelName,"w");
	fprintf(in_label,"%d\n",finalLabels);
#endif // OUTPUT_LABELS

	IndexType i = 0;
	IndexType ik = 0;
	for (Sample::vtx_iterator v_iter = curF.begin();v_iter != curF.end();v_iter++,i++)
	{
		if (m_isSelect[i])
		{
#ifdef EXPSWAP
			(*v_iter)->set_visble(true);
			// 					if (m_isEdgePoint[i])
			// 					{
			// 						(*v_iter)->set_label(m_nLabels +1);
			// 					}else
			// 					{
			(*v_iter)->set_label( m_optLabel[ik]);
			//					}

#else
			(*v_iter)->set_visble(true);

			if (m_inliers[i] == false)//可视化拟合平面--剔除噪声
			{
				(*v_iter)->set_label( m_optLabel[ik]);

			}else
			{
				(*v_iter)->set_label(m_nLabels + 1);	
			}
#endif // EXPSWAP

#ifdef OUTPUT_LABELS
			fprintf(in_label,"%d %d %d\n",centFrame,m_optLabel[ik],i);
#endif // OUTPUT_LABELS

			ik++;
		}else
		{
			(*v_iter)->set_visble(false);
		}
	}

#ifdef OUTPUT_LABELS
	fclose(in_label);
#endif // OUTPUT_LABELS

#else

	vector<IndexType> result_labels(curF.num_vertices(),0);
	propagateLabel2Orignal(curF,m_smpId,m_optLabel,result_labels);

#ifdef OUTPUT_LABELS
	char outputLabelName[256];
	sprintf(outputLabelName,"G:\\Data\\horse\\quaEva1215\\J-linkage threshold test\\res-0.50\\J-0.50-splitMerge\\cosegfinal_labels%.2d.txt",centFrame);
	FILE *in_label = fopen(outputLabelName,"w");
	fprintf(in_label,"%d\n",finalLabels);
#endif // OUTPUT_LABELS

	IndexType i = 0;
	for (Sample::vtx_iterator v_iter = curF.begin();v_iter != curF.end();v_iter++,i++)
	{
		(*v_iter)->set_label( result_labels[i]);

#ifdef OUTPUT_LABELS
		//fprintf(in_label,"%d %d %d\n",centFrame,m_gc->whatLabel(k),i);
		fprintf(in_label,"%d %d %d\n",centFrame,result_labels[i],i);
#endif // OUTPUT_LABELS

	}

#ifdef OUTPUT_LABELS
	fclose(in_label);
#endif // OUTPUT_LABELS

#endif //  SHOW_SAMPLE


}

IndexType GCop::orderLabels(vector<IndexType>& labels)
{
	if (labels.size() == 0)
	{
		Logger<<" label's vector is empty!.\n";
		return 0;
	}

	map<IndexType,IndexType> recordLabel;
	map<IndexType,IndexType>::iterator fIt;


	IndexType temp = 0;
	auto iter = labels.begin();

	recordLabel[(*iter)] = temp;

	for (;iter != labels.end(); iter ++)
	{
		fIt = recordLabel.find(*iter);
		if (fIt != recordLabel.end())
		{
			(*iter) = fIt->second;
		}else
		{
			temp++;
			recordLabel[(*iter)] = temp;
			(*iter) = temp;
		}
	}

	recordLabel.clear();

	return temp + 1;
}

void GCop::propagateLabel2Orignal(Sample& oriPC,vector<IndexType>& sampleVtxId,vector<IndexType>& label_smp,vector<IndexType>& label_ori)
{
	Logger<<"start.\n";

	IndexType nCluster = 35;

	map<IndexType,IndexType> smpLabel;
	map<IndexType,IndexType>::iterator IsValidIter;
	for (int i = 0; i < label_smp.size(); i++)
	{
		smpLabel.insert(make_pair(sampleVtxId[i],label_smp[i]));
	}

	const IndexType k = 250;
	IndexType neighbours[k];
	ScalarType dist[k];

	IndexType vtx_num = (IndexType)oriPC.num_vertices();
	IndexType result_label;


	for(IndexType vtx_id = 0; vtx_id < vtx_num; vtx_id ++)
	{
		vector<IndexType> recordLabelTime(nCluster,0);
		result_label = -1;
		oriPC.neighbours(vtx_id, k, neighbours, dist);
		for(IndexType neig_id = 0; neig_id < k; neig_id ++)
		{
			IsValidIter = smpLabel.find(neighbours[neig_id]);
			if(IsValidIter != smpLabel.end())
			{
				recordLabelTime[IsValidIter->second] += 1;
			}
		}
		for (int i = 0; i < nCluster; i++)
		{
			if(result_label < recordLabelTime[i])
			{
				label_ori[vtx_id] = i;
				result_label = recordLabelTime[i];
			}
		}

	}

	smpLabel.clear();
	Logger<<"end.\n";
}