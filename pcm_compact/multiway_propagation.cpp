#include "multiway_propagation.h"
#include "GraphMatching.h"

PoolAllocator DualwayPropagation::allocator_;

void DualwayPropagation::constructPCloudGraph()
{
	Logger<<"  Begin initialize point cloud graphs.\n";

	for (auto citer = hier_componets_.begin(); citer!=hier_componets_.end(); citer++)
	{

		IndexType nodeSize = (IndexType)citer->second.label_of_vtx.size();//.size();

		PCloudGraph* new_pcGraph_space = allocator_.allocate<PCloudGraph>();

		PCloudGraph* new_pcGraph = new (new_pcGraph_space)PCloudGraph;

		addGraphVertex(*new_pcGraph,citer->first);

		addGraphEdge(*new_pcGraph,citer->first);

		citer->second.pcGraph = new_pcGraph;

	}

	Logger<<"  End initialize point cloud graphs.\n";
}

void DualwayPropagation::addGraphVertex(PCloudGraph& pcGraph, IndexType frameId)
{
	IndexType nSize = (IndexType)hier_componets_[frameId].label_of_vtx.size();

	auto vIter = hier_componets_[frameId].label_of_vtx.begin();
	auto vEnd = hier_componets_[frameId].label_of_vtx.end();

	IndexType gIndex = 0;
	map<IndexType,IndexType> gNodeIdx;

	for (; vIter != vEnd; ++ vIter,++ gIndex)
	{
		IndexType labelId = (*vIter).second;
		IndexType vtxId = (*vIter).first;

		IndexType idxm = hier_componets_[frameId].hier_label_vtxBucket_index[0][labelId];
		HVertex* vtx = hier_componets_[frameId].hier_label_bucket[0][idxm]->vertex_bucket[vtxId];

		gNodeIdx[vtxId] = gIndex;

		PCVertexProperty vp;
		vp.index = gIndex;
		vp.vtxSite = vtx;
		boost::add_vertex(vp,pcGraph);	
	}

	hier_componets_[frameId].gId_of_vtx = gNodeIdx;

}

void DualwayPropagation::init_labeles_graph_hier(ScalarType distThr)
{

	Logger<<"  Begin initialize graphs.\n";

	if (distThr >= 1.0 || distThr < 0.0)
	{
		distThr = 0.05;
	}

	for (auto citer = hier_componets_.begin(); citer!=hier_componets_.end(); citer++)
	{

		IndexType lbsize = (IndexType)citer->second.hier_label_bucket.size();
		assert( lbsize > 0);
		IndexType nodeSize = (IndexType)citer->second.hier_label_bucket[lbsize -1].size();//.size();

		LabelsGraph* new_labelGraph_space = allocator_.allocate<LabelsGraph>();

		LabelsGraph* new_labelGraph = new (new_labelGraph_space)LabelsGraph;

		GraphVertexProperty gvp;

		map<IndexType,IndexType> labelLevel = citer->second.hier_label_vtxBucket_index[lbsize - 1];

		auto labIndex = labelLevel.begin();

		IndexType gep_count  = 0;

		for( IndexType i = 0 ; i< nodeSize && labIndex != labelLevel.end(); ++i,++ labIndex)
		{
			gvp.index = i ;
			gvp.prev =  -1;
			gvp.next = -1;
			gvp.label_id = labIndex->first;
			add_vertex( gvp , *new_labelGraph);
		}

		for (IndexType i = 0; i < nodeSize-1; i ++)
		{
			for (IndexType j = i + 1; j < nodeSize; j++)
			{
				GraphEdgeProperty   gep;

				distanPriQueue PriQue;   //
				while ( ! PriQue.empty() ) PriQue.pop(); 

				HLabel* fir = citer->second.hier_label_bucket[lbsize -1][i];

				HLabel* sec = citer->second.hier_label_bucket[lbsize -1][j];

				map<IndexType,HVertex*> fir_vtx = fir->vertex_bucket;

				map<IndexType,HVertex*> sec_vtx = sec->vertex_bucket;

				//����������֮�����̾���
				ScalarType minDis = 1e5;
				SampleSet& sample_set = SampleSet::get_instance();
				map<IndexType,HVertex*>::iterator biter1 , eiter1 ,biter2,eiter2;
				eiter1 = fir_vtx.end();
				eiter2 = sec_vtx.end();
				ScalarType dia ;

				for( biter1 = fir_vtx.begin() ;biter1 != eiter1 ;++biter1)
				{
					for( biter2 = sec_vtx.begin() ;biter2 != eiter2 ;++biter2)
					{
						IndexType index1 = biter1->first;
						Sample& s = sample_set[citer->first];
						dia = s.getBox().diag();

						PointType point1 =  s.vertices_matrix().col(index1);

						IndexType index2 = biter2->first;

						PointType point2 =  s.vertices_matrix().col(index2);

						ScalarType distance = (point1 - point2).norm();

						//��dijstra��̾�������������֮��ľ���

						if(distance < minDis)minDis = distance;
						IndexType label1 = hier_componets_[citer->first].label_of_vtx[ index1];
						IndexType label2 = hier_componets_[citer->first].label_of_vtx[ index2];
						if( distance < dia*distThr)
						{
							pointdistance tpd( index1 ,  label1 , index2 ,label2 ,distance );
							PriQue.push( tpd);
						}

					}

				}

				if (minDis <  dia * distThr )//0.0588
				{
					getEdgeVertexs2( citer->first ,PriQue ,gep.edgePoints);

					//�߽�㼯����̫��Ҳ����Ҫ���,������������.

					gep.index = gep_count;
					if( i < j)
					{
						gep.start_ = i;
						gep.end_ = j;
						boost::add_edge(i,j,gep ,*new_labelGraph);
						++gep_count;
					}else if( i> j)
					{
						gep.start_ = j;
						gep.end_ = i;
						boost::add_edge( j, i ,gep ,*new_labelGraph);
						++gep_count;
					}else
					{
						// i == j ʱ������
					}

				}

			}

		}

		//citer->second.hier_graph.clear();

		citer->second.hier_graph.push_back( new_labelGraph);

	}//֡��������

	Logger<<"  End initialize graphs.\n";
}

void DualwayPropagation::getEdgeVertexs2( IndexType _CFrameId , distanPriQueue& _PriQuemap, map<IndexType, map<IndexType ,HVertex*> >& _edgepoints )
{
	map<IndexType,HVertex*> edgedps;
	pointdistance pd;


	IndexType sr_labelLevel,tg_labelLevel;
	for( IndexType i  = 10 ; i >0 ; --i)
	{
		if (_PriQuemap.empty())
		{
			Logger<<"Queue empty.\n";
			break;
		}
		pd = _PriQuemap.top();
		_PriQuemap.pop();

		sr_labelLevel = hier_componets_[_CFrameId].hier_label_vtxBucket_index[0][pd.labelofvtx1_];
		tg_labelLevel = hier_componets_[_CFrameId].hier_label_vtxBucket_index[0][pd.labelofvtx2_];

		HVertex* tmpVtx;


		tmpVtx = hier_componets_[_CFrameId].hier_label_bucket[0][sr_labelLevel]->vertex_bucket[pd.vtx1Id_];
		edgedps[ pd.vtx1Id_] = tmpVtx;

		// 		tmpVtx = hier_componets_[_CFrameId].hier_label_bucket[0][pd.labelofvtx2_]->vertex_bucket[pd.vtx2Id_];
		// 		edgedps[ pd.vtx2Id_] = tmpVtx;



	}

	IndexType seKey = frame_label_to_key(sr_labelLevel, tg_labelLevel);
	_edgepoints[ seKey] = edgedps;

}

void DualwayPropagation::read_data(char *label_name,char *corr_name)
{
	//read_label_file(label_name); //��ԭʼ�ĺ���

	//read_label_file_hier(label_name); //��ȡ�����Ϣ

	//read_label_file_coseg(label_name); //��ȡ���ָ��label�ļ�

	//201501010

	read_label_file_coseg(label_name); //��ȡ���ָ��label�ļ�,����֮֡�䲻����label����

	read_corr_file_coseg(corr_name);

}

void DualwayPropagation::read_label_file_coseg(char *label_file_name)
{
	FILE* in_file = fopen(label_file_name, "r");
	if (in_file==NULL)
	{
		return;
	}
	IndexType frame, label, vtx_idx;
	while ( true )
	{
		int stat =  fscanf(in_file, "%d %d %d\n",&frame, &label, &vtx_idx);
		if (stat == EOF)
			break;


		if ( hier_componets_.find(frame)== hier_componets_.end() )
		{
			hier_componets_.insert(make_pair(frame, HFrame()));
			hier_componets_[frame].frame_id = frame;	
			hier_componets_[frame].hier_label_bucket.resize(1);
		}
		if ( label >= hier_componets_[frame].hier_label_bucket[0].size() )
		{
			hier_componets_[frame].hier_label_bucket[0].resize( label+1 );
		}
		if (   nullptr==hier_componets_[frame].hier_label_bucket[0][label] )
		{
			HLabel* new_label_space = allocator_.allocate<HLabel>();
			HLabel* new_label = new (new_label_space)HLabel;
			hier_componets_[frame].hier_label_bucket[0][label] = new_label;
			hier_componets_[frame].hier_label_bucket[0][label]->frame_parent = &hier_componets_[frame];
			hier_componets_[frame].hier_label_bucket[0][label]->label_id = label;
		}
		HVertex* new_space = allocator_.allocate<HVertex>();

		HVertex* new_vtx = new (new_space)HVertex(vtx_idx, hier_componets_[frame].hier_label_bucket[0][label]);
		hier_componets_[frame].hier_label_bucket[0][label]->vertex_bucket.insert( make_pair(vtx_idx,new_vtx) );
		hier_componets_[frame].label_of_vtx.insert( make_pair(vtx_idx, label) );

	}

}

void DualwayPropagation::read_corr_file_coseg(char* corr_file_name)
{
	FILE* in_file = fopen(corr_file_name,"r");
	if(in_file==NULL)
	{
		return;
	}

	IndexType cur_frame, cur_vtx_idx, next_frame, next_vtx_idx;
	while (true)
	{
		int stat = fscanf(in_file,"%d %d %d %d\n",&cur_frame, &cur_vtx_idx, &next_frame, &next_vtx_idx);
		if(stat==EOF)break;

		if( hier_componets_.find(cur_frame)==hier_componets_.end() || hier_componets_.find(next_frame)==hier_componets_.end())
			continue;

		IndexType label = hier_componets_[cur_frame].label_of_vtx[cur_vtx_idx];

		IndexType next_label = hier_componets_[next_frame].label_of_vtx[next_vtx_idx];

		HVertex& cur_vtx = *hier_componets_[cur_frame].hier_label_bucket[0][label]->vertex_bucket[cur_vtx_idx];

		if ( cur_frame+1 == next_frame  )
		{
			cur_vtx.next_corr = hier_componets_[next_frame].hier_label_bucket[0][ next_label]->vertex_bucket[next_vtx_idx];
		}
		else if (cur_frame-1 == next_frame)
		{
			cur_vtx.prev_corr = hier_componets_[next_frame].hier_label_bucket[0][ next_label]->vertex_bucket[next_vtx_idx];
		}
	}


	//ɾ��label_buckert��Ϊ�յ�Ԫ��
	for (auto fiter = hier_componets_.begin(); fiter != hier_componets_.end(); ++ fiter)
	{
		IndexType frame_id = fiter->first;
		vector<HLabel*> fHierBucket = hier_componets_[frame_id].hier_label_bucket[0];

		vector<HLabel*>::iterator  hlabelIter = fHierBucket.begin();

		map<IndexType,IndexType> vtxIndex;
		IndexType vtxLevel = 0;
		for (; hlabelIter != fHierBucket.end(); )
		{
			if (*hlabelIter == NULL)
			{
				hlabelIter = fHierBucket.erase(hlabelIter);
			}else
			{
				IndexType labelId = (*hlabelIter)->label_id;
				vtxIndex[labelId] = vtxLevel;

				++ vtxLevel;
				++ hlabelIter;
			}
		}

		hier_componets_[frame_id].hier_label_bucket[0] = fHierBucket;

		hier_componets_[frame_id].hier_label_vtxBucket_index.resize(1);

		hier_componets_[frame_id].hier_label_vtxBucket_index[0] = vtxIndex;

	}

}
void DualwayPropagation::wirteSplitGraphLables(std::string filename)
{

	FILE* outfile = fopen(filename.c_str(),"w");

	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); fIter ++)
	{
		IndexType gLevel = (IndexType)fIter->second.hier_label_bucket.size();//������߲��label_bucket
		IndexType fId = fIter->first;
		vector<HLabel*>& label_buctet = fIter->second.hier_label_bucket[gLevel - 1];

		for (auto lIter = label_buctet.begin(); lIter != label_buctet.end(); ++ lIter)
		{

			IndexType pId = (*lIter)->label_id;
			
			map<IndexType,HVertex*>& vtx_bucket = (*lIter)->vertex_bucket;

			for (auto vIt = vtx_bucket.begin(); vIt != vtx_bucket.end(); ++vIt)
			{
				IndexType vtx_id = vIt->first;

				fprintf(outfile, "%d %d %d\n", fId, pId, vtx_id);
			}
		}
	}

	fclose(outfile);
}

void DualwayPropagation::splitAllSquenceGraph(IndexType iterN)
{
	IndexType fSize = (IndexType)hier_componets_.size();

	map<IndexType,HFrame>::iterator  cIter= hier_componets_.begin();

	map<IndexType,HFrame>::iterator  cEnd = hier_componets_.end();


	// 	IndexType startF = 4;
	// 	while (startF -- > 0)
	// 	{
	// 		++cIter;
	// 	}

	--cEnd;

	IndexType iterNum = 0;

	IndexType startFrameId = cIter->first;



	for (; cIter != cEnd && iterNum < iterN; ++ cIter++, ++iterNum)
	{
		//���һ���ڵ㲻��������
		IndexType srFrame = cIter->first;
		IndexType tgFrame = srFrame + 1;

		//split_twoAjacent_graph_next(srFrame,tgFrame );

		//show_corresponding(srFrame);

		split_twoAjacent_graph_next_order(srFrame,tgFrame );

		split_nest_graph_prev(startFrameId,srFrame,tgFrame);

	}

	//show_corresponding(1);

	//���ÿ֡�ķָ����

	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); ++ fIter)
	{
		IndexType gLevel = (IndexType)fIter->second.hier_label_bucket.size();
		IndexType labeSize =(IndexType) fIter->second.hier_label_bucket[gLevel - 1].size();
		Logger<<"  ��"<<fIter->first<<"֡���ָ��"<<labeSize<<"����.\n";
	}

}

//�Ա�����֮����з���

void DualwayPropagation::split_twoAjacent_graph_next_order(IndexType srFrame, IndexType tgFrame)
{
	//��ǰ����

	//get the new graph of tgGrame--��Ҫ���

	Logger<<" .......\n";
	Logger<<"  Start next split.\n";
	IndexType tgGraSize = (IndexType)hier_componets_[tgFrame].hier_graph.size();
	LabelsGraph* oriGra = hier_componets_[tgFrame].hier_graph[tgGraSize - 1];
	LabelsGraph* new_graph = new LabelsGraph(*oriGra);

	//new_graph = oriGra;

	//vector<HLabel* > new_label_bucket =  hier_componets_[tgFrame].hier_label_bucket[tgGraSize - 1];
	vector<HLabel* > new_label_bucket;
	copyLabelBucket(new_label_bucket,hier_componets_[tgFrame].hier_label_bucket[tgGraSize - 1] );

	//
	IndexType gLevel = 0;

	IndexType srGraSize = (IndexType)hier_componets_[srFrame].hier_graph.size();

	gLevel  = srGraSize - 1;//��ȡ���µĲ�
	LabelsGraph* srGraLat = hier_componets_[srFrame].hier_graph[gLevel];

	IndexType labParentsize = tgGraSize + 1; //���ɵĲ���

	Logger<<srFrame<<"֡�ĵ�"<<gLevel<<"��߽�ָ�"<<tgFrame<<"��"<<tgGraSize - 1 <<"��"<<endl;

	pair<EdgeIterator,EdgeIterator> ei = boost::edges(*srGraLat);

	generateOrderededges(srFrame,tgFrame);

	while (orderedEdgeQ.size() > 0)
	{
		EdgeSplitOrder oEdge = orderedEdgeQ.top();
		orderedEdgeQ.pop();

		EdgeDescriptor ed = oEdge.EdgeDec;

		GraphEdgeProperty& ep = (*srGraLat)[ed];

		map<IndexType,HVertex*> edgePoints;

		auto ePsIt = ep.edgePoints.begin();

		edgePoints.insert(ePsIt->second.begin(),ePsIt->second.end() );  

		if (edgePoints.size() < 1)
		{
			Logger<<"���ϵĶ�����̫��,�޷�����.\n";
			continue;
		}

		map<IndexType,HVertex*> edgeCorrNextVtx;

		IndexType newGraphEdgeSize = (IndexType)new_graph->m_edges.size();

		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//��ñ߽������һ֡��Ӧ�Ŀ�Ͷ�Ӧ��
		//IndexType nodeId = edgeCorrNextVtx.size();


		//��Ӧ��ȥ,��ǩ����ʼ��,�򱣳ֲ���
		HLabel* splitedLabel = new_label_bucket[nodeId];

		IndexType eS = ep.start_;
		IndexType eE = ep.end_;

		Logger<<"  �ߵ����Ϊ"<<eS<<"�յ�Ϊ"<<eE<<endl;


		//�����ѳ����ĵ������һ�����ݺ���,��ñ߲����ѱ�

		if ( oEdge.srCorNum < 10 || oEdge.tgCorNum < 10 )
		{
			Logger<<"�߽��̫����,����Ҫ����.\n";
			continue;
		}

		if ( oEdge.unMarkedRation < 0.2)
		{
			Logger<<"Unmark���ֵ̫��,��ʱ������.\n";
			continue;
		}

		//����nodeId �������ӵı�

		pair<VertexIterator, VertexIterator> vi = boost::vertices(*new_graph);

		VertexIterator nodeIter = (vi.first + nodeId);

		VertexDescriptor nodeDesc = *nodeIter;

		//�ڵ��Ӧ�����г���
		pair<OutEdgeIterator,OutEdgeIterator> nodeEiter = boost::out_edges(nodeDesc,*new_graph);

		map<IndexType,map<IndexType,HVertex*> > recordColapseEdges;

		set<GraphEdgeProperty> collapseEdges;

		OutEdgeIterator oit,nextIt;

		oit = nodeEiter.first;

		for (nextIt = oit; oit != nodeEiter.second; oit = nextIt )
		{
			++nextIt;

			EdgeDescriptor nextEdgeD = *oit;

			GraphEdgeProperty& nextEP = (*new_graph)[nextEdgeD];

			collapseEdges.insert(nextEP);

			boost::remove_edge(*oit,*new_graph);//ɾ��������
		}

		// 		//����һ���ڵ�// ���Ӷ����������ǰ��.
		// 
		IndexType nSize = (IndexType)boost::num_vertices(*new_graph);
		// 
		// 		GraphVertexProperty vp(nSize,-1,-1);
		// 
		// 		boost::add_vertex(vp,*new_graph);

		//���·ָ����Ϣ,�����ӵ�Label���ΪnSize. �����ѵĵ�ΪnodeId
		IndexType new_label = nSize;

		new_label_bucket.push_back((HLabel*)0 );
		HLabel* new_label_space = allocator_.allocate<HLabel>();
		HLabel* new_label_obj = new (new_label_space)HLabel;
		new_label_bucket[new_label] = new_label_obj;
		new_label_bucket[new_label]->label_id = new_label;
		new_label_bucket[new_label]->frame_parent = &hier_componets_[tgFrame];

		// 		//��Ӧ��ȥ,��ǩ����ʼ��,�򱣳ֲ���

		map<IndexType,HVertex*> unMarkPs;

		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); )
		{
			IndexType vtx_id = iter->first;

			IndexType pId = iter->second->prev_corr->vtx_id;

			IndexType prev_id = iter->second->prev_corr->label_parent[gLevel]->label_id; //�ò������µ�label_parent��ַ.

			if (prev_id == eS)
			{

				iter->second->label_parent.resize(tgGraSize + 1);
				iter->second->label_parent[tgGraSize] =  new_label_bucket[nodeId];
				++iter;

			}else if(prev_id == eE)
			{

				iter->second->label_parent.resize(tgGraSize + 1);
				iter->second->label_parent[tgGraSize] = new_label_bucket[new_label] ;
				new_label_bucket[new_label]->vertex_bucket.insert(*iter);
				hier_componets_[tgFrame].label_of_vtx[vtx_id ] = new_label;
				iter = splitedLabel->vertex_bucket.erase(iter);


			}else//һЩ��ȷ��label�ĵ�
			{
				unMarkPs.insert(*iter);
				iter = splitedLabel->vertex_bucket.erase(iter);
			}

		}

		if (!unMarkPs.empty())
		{
			if ( (!new_label_bucket[new_label]->vertex_bucket.empty() ) && (!new_label_bucket[nodeId]->vertex_bucket.empty() ) )
			{
				//�����ȡ���������С�������жϲ�ȷ���������ĸ���.unmark Ҫô����nodeid Ҫô����new_label
				determinateUnmarkPoints(tgFrame,unMarkPs,new_label_bucket,nodeId,new_label,tgGraSize);

			}else if (new_label_bucket[new_label]->vertex_bucket.empty())
			{
				Logger<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ѳ����ĵ�̫��.\n";

				new_label_bucket.pop_back();

				for (auto iter = unMarkPs.begin(); iter != unMarkPs.end();)//�Ż�ԭ���Ŀ���--��û���ӽڵ�
				{
					new_label_bucket[nodeId]->vertex_bucket.insert(*iter);

					hier_componets_[tgFrame].label_of_vtx[iter->first] = nodeId;

					iter->second->label_parent.resize(tgGraSize + 1);
					iter->second->label_parent[tgGraSize] = new_label_bucket[nodeId];

					iter = unMarkPs.erase(iter); 

				}

				continue;

			}else if(new_label_bucket[nodeId]->vertex_bucket.empty() )
			{
				Logger<<"ori parch empty!.\n";
				//unmarked��new_label���Żص�ԭ���Ŀ���.

				for (auto iter = unMarkPs.begin(); iter != unMarkPs.end();)//�Ż�ԭ���Ŀ���--��û���ӽڵ�
				{
					new_label_bucket[nodeId]->vertex_bucket.insert(*iter);

					hier_componets_[tgFrame].label_of_vtx[iter->first] = nodeId;

					iter->second->label_parent.resize(tgGraSize + 1);
					iter->second->label_parent[tgGraSize] = new_label_bucket[nodeId];

					iter = unMarkPs.erase(iter); 

				}

				//�������ɵĿ�Ҳ�Ż�ԭ����

				auto bMapVtx = new_label_bucket[new_label]->vertex_bucket.begin();
				auto eMapVtx = new_label_bucket[new_label]->vertex_bucket.end();

				for (; bMapVtx != eMapVtx; ++ bMapVtx)
				{
					new_label_bucket[nodeId]->vertex_bucket.insert(*bMapVtx);

					hier_componets_[tgFrame].label_of_vtx[bMapVtx->first] = nodeId;

					bMapVtx->second->label_parent.resize(tgGraSize + 1);

					bMapVtx->second->label_parent[tgGraSize] = new_label_bucket[nodeId];

				}

				new_label_bucket.pop_back();

				break;
			}

		}


		//����һ���ڵ� --�����Ӷ����������ȷ������һ��������.

		//IndexType nSize = boost::num_vertices(*new_graph);

		GraphVertexProperty vp(nSize,-1,-1);

		boost::add_vertex(vp,*new_graph);


		//������������мӱ߲���,
		//node<--->new_node
		GraphEdgeProperty newEP;
		newEP.start_ = nodeId;
		newEP.end_ = nSize;
		newEP.index = newGraphEdgeSize;

		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
		//��û��index��ֵ

		newEP.edgePoints[edgeKey].insert(edgeCorrNextVtx.begin(),edgeCorrNextVtx.end());

		boost::add_edge(nodeId,nSize,newEP,*new_graph);

		for (auto iter = collapseEdges.begin(); iter != collapseEdges.end(); iter ++)
		{
			GraphEdgeProperty glueEdge;

			glueEdge = (*iter);
			IndexType  sEdgeId = glueEdge.start_;
			IndexType  eEdgeId = glueEdge.end_;

			map<IndexType,HVertex*>  edgePoints;

			auto bIter = glueEdge.edgePoints.begin();
			edgePoints.insert(bIter->second.begin(),bIter->second.end() );

			map<IndexType,HVertex*> startVtx = new_label_bucket[sEdgeId]->vertex_bucket;
			map<IndexType,HVertex*> endVtx = new_label_bucket[eEdgeId]->vertex_bucket;
			map<IndexType,HVertex*> nodeVtx = new_label_bucket[nodeId]->vertex_bucket;
			map<IndexType,HVertex*> nSizeVtx = new_label_bucket[nSize]->vertex_bucket;

			ScalarType minNode = 0.0;
			ScalarType minSize = 0.0;


			if (sEdgeId < nodeId) //nodeidΪ�յ�
			{
				minDistBeTwoParts(tgFrame,startVtx,nodeVtx,minNode);
				minDistBeTwoParts(tgFrame,startVtx,nSizeVtx,minSize);

				if (minNode < minSize)
				{
					glueEdge.end_ = nodeId;
				}else
				{
					glueEdge.end_ = nSize;
				}

			}else//nodeidΪ���
			{
				minDistBeTwoParts(tgFrame,endVtx,nodeVtx,minNode);
				minDistBeTwoParts(tgFrame,endVtx,nSizeVtx,minSize);

				if (minNode < minSize)
				{
					glueEdge.start_ = nodeId;
					glueEdge.end_ = eEdgeId;

				}else
				{
					glueEdge.start_ = eEdgeId;
					glueEdge.end_ = nSize;
				}

			}

			IndexType eKey = frame_index_to_key(glueEdge.start_,glueEdge.end_);

			glueEdge.edgePoints[eKey] = edgePoints;

			boost::add_edge(glueEdge.start_,glueEdge.end_,glueEdge,*new_graph);

		}//����collapse�ı�


	} //end  while

	checkPsNewLabelParentPtr(new_label_bucket,labParentsize);//next dirction

	map<IndexType,IndexType> labelIndex;
	IndexType kk=0;
	for (auto iter = new_label_bucket.begin(); iter != new_label_bucket.end(); ++ iter,++kk)
	{
		IndexType label = (*iter)->label_id;
		labelIndex[label] = kk;
	}

	hier_componets_[tgFrame].hier_label_bucket.push_back(new_label_bucket);

	hier_componets_[tgFrame].hier_graph.push_back(new_graph);//�������µ�graph

	hier_componets_[tgFrame].hier_label_vtxBucket_index.push_back(labelIndex);




	Logger<<"  End next split.\n";

	Logger<<" .......\n";
}
void DualwayPropagation::split_nest_graph_prev(IndexType startFrame,IndexType srFrame, IndexType tgFrame)
{
	if(srFrame == startFrame)
	{
		split_twoAjacent_graph_prev(srFrame,tgFrame);
		return;
	}

	while (srFrame > startFrame -1 && (startFrame > -1) )
	{
		split_twoAjacent_graph_prev(srFrame,tgFrame);

		srFrame --;

		tgFrame --;
	}

	return;
}

void DualwayPropagation::split_twoAjacent_graph_prev(IndexType srFrame, IndexType tgFrame)
{

	//��ȡ��Ҫ���µ�graph
	Logger<<" .......\n";
	Logger<<"  Start prev split.\n";

	IndexType srGraphSize = (IndexType)hier_componets_[srFrame].hier_graph.size();
	IndexType tgGraphSize = (IndexType)hier_componets_[tgFrame].hier_graph.size();

	assert(srGraphSize > 0 && tgGraphSize > 0);

	LabelsGraph* oriSpGra = hier_componets_[srFrame].hier_graph[srGraphSize - 1];
	LabelsGraph* shouldSplitGraph = new LabelsGraph(*oriSpGra);

	vector<HLabel* > new_label_bucket;
	copyLabelBucket(new_label_bucket,hier_componets_[srFrame].hier_label_bucket[srGraphSize - 1] );

	//��ȡָ���ָ��graph
	LabelsGraph* guideSplitGraph;
	IndexType graphLevel = 0;

	assert(tgGraphSize > 1);

	if (tgGraphSize > 2)
	{
		graphLevel = tgGraphSize - 1;
	}else
	{
		graphLevel = tgGraphSize - 2;
	}
	guideSplitGraph = hier_componets_[tgFrame].hier_graph[graphLevel];

	//�Ա߽������ȼ�����




	// guideSplitGraph��ÿ��������һ�ηָ�
	pair<EdgeIterator,EdgeIterator> ei = boost::edges(*guideSplitGraph);

	for (EdgeIterator eit = ei.first; eit != ei.second; ++eit)
	{
		EdgeDescriptor ed = *eit;

		GraphEdgeProperty& ep = (*guideSplitGraph)[ed];

		map<IndexType,HVertex*> edgePoints;

		auto ePsIt = ep.edgePoints.begin();

		edgePoints.insert(ePsIt->second.begin(),ePsIt->second.end() );  

		if (edgePoints.size() < 1)
		{
			Logger<<"���ϵĶ�����̫��,�޷�����.\n";
			continue;
		}

		Logger<<tgFrame<<"֡�ı����"<<ep.start_<<"���յ�"<<ep.end_<<"��ʼ����"<<endl;

		map<IndexType,HVertex*> edgeCorrPrevVtx;

		IndexType nSize = (IndexType)boost::num_vertices(*shouldSplitGraph);

		IndexType newGraphEdgeSize = (IndexType)shouldSplitGraph->m_edges.size(); //Ϊ�˸������ӵı�������

		bool isSplit = true;

		IndexType edgePsCorNode = 0;

		isSplit = checkPrevLabelBucket(edgePoints,edgeCorrPrevVtx,edgePsCorNode);

		if ( edgePsCorNode < 0 || edgePsCorNode > nSize - 1)
		{
			Logger<<"�߽��ҵ��Ŀ�Խ��!.\n";
			break;
		}

		if (!isSplit)
		{
			Logger<<"�߽��ָ����������,���ָ�.\n";
			continue;
		}

		HLabel* splitedLabel = new_label_bucket[edgePsCorNode]; //������Ҫ���ѵĿ�

		IndexType curEdgeStart = ep.start_;
		IndexType curEdgeEnd   = ep.end_;
		IndexType strCorPsSzie = 0;
		IndexType endCorPsSize = 0;

		//����ǰ��һ���򵥵�Ԥ�ж�
		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); iter ++)
		{
			if (iter->second->next_corr->label_parent[graphLevel] != NULL)
			{

				IndexType nextVtx_label = iter->second->next_corr->label_parent[graphLevel]->label_id; 

				if (nextVtx_label == curEdgeStart)
				{
					strCorPsSzie ++;

				}else if(nextVtx_label == curEdgeEnd)
				{
					endCorPsSize ++;
				}
			}

		}

		IndexType  vtxBSize = (IndexType)splitedLabel->vertex_bucket.size();

		ScalarType ration = (ScalarType)(strCorPsSzie + endCorPsSize)/vtxBSize;

		if ( strCorPsSzie < 10 || endCorPsSize < 10)
		{
			Logger<<"�߽�㿿���߽�(��Ӧ�����),����Ҫ����.\n";
			continue;
		}

		if ( ration < 0.2)
		{
			Logger<<"������ǧ��?���˰�!.\n";
			continue;
		}

		//��eit �߿�ʼ����
		pair<VertexIterator, VertexIterator> vi = boost::vertices(*shouldSplitGraph);

		VertexIterator curNodeIter = (vi.first + edgePsCorNode);

		VertexDescriptor curNodeDesc = *curNodeIter;

		//�ڵ��Ӧ�����г���
		pair<OutEdgeIterator,OutEdgeIterator> nodeEiter = boost::out_edges(curNodeDesc,*shouldSplitGraph);

		map<IndexType,map<IndexType,HVertex*> > recordColapseEdges;

		set<GraphEdgeProperty> collapseEdges;

		OutEdgeIterator oit,nextIt;

		oit = nodeEiter.first;

		for (nextIt = oit; oit != nodeEiter.second; oit = nextIt )
		{
			++nextIt;

			EdgeDescriptor nextEdgeD = *oit;

			GraphEdgeProperty& nextEP = (*shouldSplitGraph)[nextEdgeD];

			collapseEdges.insert(nextEP);

			boost::remove_edge(*oit,*shouldSplitGraph);//ɾ��������
		}

		//����һ���ڵ�

		GraphVertexProperty vp(nSize,-1,-1);
		boost::add_vertex(vp,*shouldSplitGraph);

		//���·ָ����Ϣ,�����ӵ�Label���ΪnSize. �����ѵĵ�ΪedgePsCorNode
		IndexType new_label = nSize;
		new_label_bucket.push_back((HLabel*)0 );
		HLabel* new_label_space = allocator_.allocate<HLabel>();
		HLabel* new_label_obj = new (new_label_space)HLabel;
		new_label_bucket[new_label] = new_label_obj;
		new_label_bucket[new_label]->label_id = new_label;
		new_label_bucket[new_label]->frame_parent = &hier_componets_[tgFrame];

		//��ʼ����
		map<IndexType,HVertex*> unMarkPs;

		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); )
		{
			if (iter->second->next_corr->label_parent[graphLevel])
			{
				IndexType nextVtx_label = iter->second->next_corr->label_parent[graphLevel]->label_id; //�ò������µ�label_parent��ַ.

				IndexType vtx_id = iter->first;

				if (nextVtx_label == curEdgeStart)
				{
					iter->second->label_parent.resize(srGraphSize + 1);
					iter->second->label_parent[srGraphSize] = new_label_bucket[edgePsCorNode];
					++iter;

				}else if(nextVtx_label == curEdgeEnd)
				{
					iter->second->label_parent.resize(srGraphSize + 1);

					iter->second->label_parent[srGraphSize] = new_label_bucket[new_label];

					new_label_bucket[new_label]->vertex_bucket.insert(*iter);

					hier_componets_[srFrame].label_of_vtx[vtx_id ] = new_label;

					iter = splitedLabel->vertex_bucket.erase(iter);

				}else//һЩ��ȷ��label�ĵ�
				{
					unMarkPs.insert(*iter);
					iter = splitedLabel->vertex_bucket.erase(iter);
				}
			}

		} //������Ҫ���ѿ��ÿ����

		if (!unMarkPs.empty())
		{
			if ( (!new_label_bucket[new_label]->vertex_bucket.empty() ) && (!new_label_bucket[edgePsCorNode]->vertex_bucket.empty() ) )
			{
				//�����ȡ���������С�������жϲ�ȷ���������ĸ���.unmark Ҫô����nodeid Ҫô����new_label
				determinateUnmarkPoints(srFrame,unMarkPs,new_label_bucket,edgePsCorNode,new_label,srGraphSize);

			}else if (new_label_bucket[new_label]->vertex_bucket.empty())
			{
				Logger<<"���ѳ����ĵ�̫��.\n";

				new_label_bucket.pop_back();

				for (auto iter = unMarkPs.begin(); iter != unMarkPs.end();)//�Ż�ԭ���Ŀ���
				{
					new_label_bucket[edgePsCorNode]->vertex_bucket.insert(*iter);

					hier_componets_[tgFrame].label_of_vtx[iter->first] = edgePsCorNode;

					iter->second->label_parent.resize(srGraphSize + 1);
					iter->second->label_parent[srGraphSize] = new_label_bucket[edgePsCorNode];

					iter = unMarkPs.erase(iter); 

				}

				continue;

			}else if(new_label_bucket[edgePsCorNode]->vertex_bucket.empty() )
			{
				Logger<<"ori parch empty!.\n";

				break;
			}

		}

		//�����ȡ���������С�������жϲ�ȷ���������ĸ���.unmark Ҫô����nodeid Ҫô����new_label
		//determinateUnmarkPoints(srFrame,unMakePs,new_label_bucket,edgePsCorNode,new_label,srGraphSize);


		//������������мӱ߲���,
		//node<--->new_node
		GraphEdgeProperty newEP;
		newEP.start_ = edgePsCorNode;
		newEP.end_ = nSize;
		newEP.index = newGraphEdgeSize;
		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
		newEP.edgePoints[edgeKey].insert(edgeCorrPrevVtx.begin(),edgeCorrPrevVtx.end());
		boost::add_edge(edgePsCorNode,nSize,newEP,*shouldSplitGraph);


		//�϶����������ڵ��������ڵ�������߲���ֻ����� recordColapseEdges.size()����.
		for (auto iter = collapseEdges.begin(); iter != collapseEdges.end(); iter ++)
		{
			GraphEdgeProperty glueEdge;
			glueEdge = (*iter);
			IndexType  sEdgeId = glueEdge.start_;
			IndexType  eEdgeId = glueEdge.end_;

			map<IndexType,HVertex*>  edgePoints;

			auto bIter = glueEdge.edgePoints.begin();
			edgePoints.insert(bIter->second.begin(),bIter->second.end() );

			map<IndexType,HVertex*> startVtx = new_label_bucket[sEdgeId]->vertex_bucket;
			map<IndexType,HVertex*> endVtx = new_label_bucket[eEdgeId]->vertex_bucket;
			map<IndexType,HVertex*> nodeVtx = new_label_bucket[edgePsCorNode]->vertex_bucket;
			map<IndexType,HVertex*> nSizeVtx = new_label_bucket[nSize]->vertex_bucket;

			ScalarType minNode = 0.0;
			ScalarType minSize = 0.0;


			if (sEdgeId < edgePsCorNode) //nodeidΪ�յ�
			{
				minDistBeTwoParts(srFrame,startVtx,nodeVtx,minNode);
				minDistBeTwoParts(srFrame,startVtx,nSizeVtx,minSize);

				if (minNode < minSize)
				{
					glueEdge.end_ = edgePsCorNode;
				}else
				{
					glueEdge.end_ = nSize;
				}

			}else//nodeidΪ���
			{
				minDistBeTwoParts(srFrame,endVtx,nodeVtx,minNode);
				minDistBeTwoParts(srFrame,endVtx,nSizeVtx,minSize);

				if (minNode < minSize)
				{
					glueEdge.start_ = edgePsCorNode;
					glueEdge.end_ = eEdgeId;

				}else
				{
					glueEdge.start_ = eEdgeId;
					glueEdge.end_ = nSize;
				}

			}

			IndexType eKey = frame_index_to_key(glueEdge.start_,glueEdge.end_);

			glueEdge.edgePoints[eKey] = edgePoints;

			boost::add_edge(glueEdge.start_,glueEdge.end_,glueEdge,*shouldSplitGraph);

		}//����collapse�ı�

	} //���������ָ�ͼ��ÿ����

	checkPsNewLabelParentPtr(new_label_bucket,srGraphSize + 1);

	map<IndexType,IndexType> labelIndex;
	IndexType kk=0;
	for (auto iter = new_label_bucket.begin(); iter != new_label_bucket.end(); ++ iter,++kk)
	{
		IndexType label = (*iter)->label_id;
		labelIndex[label] = kk;
	}

	hier_componets_[srFrame].hier_label_bucket.push_back(new_label_bucket);

	hier_componets_[srFrame].hier_graph.push_back(shouldSplitGraph);//�������µ�graph

	hier_componets_[srFrame].hier_label_vtxBucket_index.push_back(labelIndex);

	Logger<<"  End prev split.\n";
	Logger<<" .......\n";
}
void DualwayPropagation::addGraphEdge(PCloudGraph& pcGraph, IndexType frameId)
{
	IndexType nSize = (IndexType)hier_componets_[frameId].label_of_vtx.size();

	auto vIter = hier_componets_[frameId].label_of_vtx.begin();
	auto vEnd = hier_componets_[frameId].label_of_vtx.end();

	unordered_map<IndexType,bool> recordEdges;

	buildKdTree(frameId);

	IndexType gIndex = 0;
	const IndexType k = 6;
	IndexType neighbours[k];
	ScalarType dist[k];

	IndexType edNum = 0;
	for (; vIter != vEnd; ++ vIter,++ gIndex)
	{
		IndexType vtxId = (*vIter).first;
		downSample->neighbours(gIndex,k,neighbours,dist);

		for (IndexType i = 1; i < k; ++i)
		{
			//�������㷨��н��ж�,��������180��,�����һ����;������ӱ߲���.

			bool temp = recordEdges[frame_index_to_key(gIndex,neighbours[i]) ];

			if (!temp)
			{
				PCEdgeProperty ep;
				ep.index = edNum;
				ep.start_ = gIndex;
				ep.end_ = neighbours[i];
				ep.dist = dist[i];

				boost::add_edge(gIndex,neighbours[i],ep,pcGraph);

				recordEdges[frame_index_to_key(neighbours[i],gIndex)] = true;

				++ edNum;
			}
		}
	}

}

// ȷ��ÿ�������labparentָ����Ŀ������ͬ
void DualwayPropagation::checkPsNewLabelParentPtr(vector<HLabel*> oriLabelBucket,IndexType labParSize)
{
	IndexType vId = 0;
	for (auto vIter = oriLabelBucket.begin(); vIter != oriLabelBucket.end(); vIter ++,vId ++)
	{
		map<IndexType,HVertex*>& vtxBucket = (*vIter)->vertex_bucket;

		if (vtxBucket.empty())
		{
			continue;
		}

		HVertex* fVtx = (*vtxBucket.begin()).second;

		IndexType vlfSize = (IndexType)fVtx->label_parent.size();

		if ( vlfSize < labParSize)
		{
			for (auto lIter =vtxBucket.begin(); lIter != vtxBucket.end(); lIter ++ )
			{
				(*lIter).second->label_parent.push_back(oriLabelBucket[vId] );
			}
		}

	}
}

void DualwayPropagation::minDistBeTwoParts(IndexType cFrame, map<IndexType,HVertex*>& fPart,map<IndexType,HVertex*>& sPart, ScalarType& misDis)
{

	Sample& smp = SampleSet::get_instance()[cFrame];

	ScalarType diag = smp.getBox().diag();

	vector<ScalarType> totDis;

	for (auto fIter = fPart.begin(); fIter != fPart.end(); ++fIter)
	{
		for (auto sIter = sPart.begin(); sIter != sPart.end(); ++sIter)
		{
			PointType fVtx = smp.vertices_matrix().col(fIter->first);
			PointType sVtx = smp.vertices_matrix().col(sIter->first);

			ScalarType dis = (fVtx - sVtx).norm();

			if (dis < 0.5 * diag)
			{
				totDis.push_back(dis);
			}
		}
	}

	if (totDis.size() == 0)
	{
		misDis =  0.5*diag;
	}else
	{
		misDis = * min_element(totDis.begin(),totDis.end() );
	}


}

void DualwayPropagation::determinateUnmarkPoints(IndexType cFrame,map<IndexType,HVertex*>& unMarkPs,vector<HLabel*> oriLabelBucket,IndexType nodeId,IndexType newLabe,IndexType tgSize)
{
	if (unMarkPs.size() < 1)
	{
		Logger<<"û��Unmark ��.\n";
		return;
	}

	Sample& smp = SampleSet::get_instance()[cFrame];

	for (auto iter = unMarkPs.begin(); iter != unMarkPs.end();)
	{
		PointType pCoor = smp.vertices_matrix().col(iter->first);

		// 		ScalarType d2node = p2PatchAvgDis(cFrame,pCoor,oriLabelBucket[nodeId]->vertex_bucket);//���ȡֵ����Сֵ
		// 		ScalarType d2size = p2PatchAvgDis(cFrame,pCoor,oriLabelBucket[newLabe]->vertex_bucket );


		// 		ScalarType d2node = p2PatchMinDis(cFrame,pCoor,oriLabelBucket[nodeId]->vertex_bucket);//�����ϸ���С
		// 		ScalarType d2size = p2PatchMinDis(cFrame,pCoor,oriLabelBucket[newLabe]->vertex_bucket );


		// 		//����߾���
		// 		ScalarType d2node = 0.;
		// 		ScalarType d2size = 0.;
		// 
		// 		if (oriLabelBucket[nodeId]->vertex_bucket.empty() )
		// 		{
		// 			d2node = 1e5;
		// 		}else
		// 		{
		// 			d2node = p2PatchGeoDis(cFrame,*(*iter).second,oriLabelBucket[nodeId]->vertex_bucket);//�����ϸ���С
		// 		}
		// 
		// 		if (oriLabelBucket[newLabe]->vertex_bucket.empty() )
		// 		{
		// 			d2size = 1e5;
		// 		}else
		// 		{
		// 			d2size = p2PatchGeoDis(cFrame,*(*iter).second,oriLabelBucket[newLabe]->vertex_bucket );
		// 		}

		ScalarType d2node = p2PatchGeoDis(cFrame,*(*iter).second,oriLabelBucket[nodeId]->vertex_bucket);//�����ϸ���С
		ScalarType d2size = p2PatchGeoDis(cFrame,*(*iter).second,oriLabelBucket[newLabe]->vertex_bucket );

		if (d2node <= d2size)
		{

			oriLabelBucket[nodeId]->vertex_bucket.insert(*iter);

			hier_componets_[cFrame].label_of_vtx[iter->first] = nodeId;

			iter->second->label_parent.resize(tgSize + 1);
			iter->second->label_parent[tgSize] = oriLabelBucket[nodeId];

		}else
		{

			oriLabelBucket[newLabe]->vertex_bucket.insert(*iter);

			hier_componets_[cFrame].label_of_vtx[iter->first] = newLabe;

			//iter->second->label_parent.push_back(oriLabelBucket[newLabe] );
			iter->second->label_parent.resize(tgSize + 1);
			iter->second->label_parent[tgSize] = oriLabelBucket[newLabe] ;
		}

		iter = unMarkPs.erase(iter); 
	}

}

IndexType DualwayPropagation::checkNextLabelBucket(map<IndexType,HVertex*>& edgePs, map<IndexType,HVertex*>& edgePsCor)
{
	assert(edgePs.size() > 0);

	set<IndexType > labelS;
	labelS.clear();

	for (auto iter = edgePs.begin(); iter != edgePs.end(); iter ++)
	{
		IndexType corPId = iter->second->next_corr->vtx_id;
		//

		//IndexType TestnextVtxLabid = hier_components_[2].label_of_vtx[corPId];
		IndexType labParSize = (IndexType)iter->second->next_corr->label_parent.size();
		IndexType updateLevel = labParSize - 1;

		IndexType nextVtxLabid = iter->second->next_corr->label_parent[updateLevel]->label_id;

		labelS.insert(nextVtxLabid);

		HVertex* nextP = iter->second->next_corr;

		edgePsCor[corPId] = nextP;

	}

	if (labelS.size() > 1)
	{
		Logger<<"�߽��ָ�����¸��������Label!.\n";

	}

	return (*labelS.begin() );
}

void DualwayPropagation::generateOrderededges(IndexType srFrame, IndexType tgFrame)
{
	Logger<<"�Ա߽�������.\n";

	while (!orderedEdgeQ.empty() )
	{
		orderedEdgeQ.pop();
	}

	//��ǰ����

	//get the new graph of tgGrame--��Ҫ���


	Logger<<" .......\n";
	Logger<<"  Start next split.\n";
	IndexType tgGraSize = (IndexType)hier_componets_[tgFrame].hier_graph.size();
	LabelsGraph* oriGra = hier_componets_[tgFrame].hier_graph[tgGraSize - 1];
	LabelsGraph* new_graph = new LabelsGraph(*oriGra);

	//new_graph = oriGra;

	//vector<HLabel* > new_label_bucket =  hier_componets_[tgFrame].hier_label_bucket[tgGraSize - 1];
	vector<HLabel* > new_label_bucket;
	copyLabelBucket(new_label_bucket,hier_componets_[tgFrame].hier_label_bucket[tgGraSize - 1] );


	//
	IndexType gLevel = 0;

	IndexType srGraSize = (IndexType)hier_componets_[srFrame].hier_graph.size();

	gLevel  = srGraSize - 1;//��ȡ���µĲ�
	LabelsGraph* srGraLat = hier_componets_[srFrame].hier_graph[gLevel];

	IndexType labParentsize = tgGraSize + 1; //���ɵĲ���

	//Logger<<srFrame<<"֡�ĵ�"<<gLevel<<"��߽�ָ�"<<tgFrame<<"��"<<tgGraSize - 1 <<"��"<<endl;

	pair<EdgeIterator,EdgeIterator> ei = boost::edges(*srGraLat);


	EdgeSplitOrder tempEdge;

	for (EdgeIterator eit = ei.first; eit != ei.second; ++eit)
	{
		EdgeDescriptor ed = *eit;

		GraphEdgeProperty& ep = (*srGraLat)[ed];

		map<IndexType,HVertex*> edgePoints;

		auto ePsIt = ep.edgePoints.begin();

		edgePoints.insert(ePsIt->second.begin(),ePsIt->second.end() );  

		if (edgePoints.size() < 3)
		{
			Logger<<"���ϵĶ�����̫��,�޷�����.\n";
			continue;
		}

		map<IndexType,HVertex*> edgeCorrNextVtx;

		IndexType newGraphEdgeSize = (IndexType)new_graph->m_edges.size();

		//��û�и���ͼ�ṹ

		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//��ñ߽������һ֡��Ӧ�Ŀ�Ͷ�Ӧ��
		//IndexType nodeId = edgeCorrNextVtx.size();


		//��Ӧ��ȥ,��ǩ����ʼ��,�򱣳ֲ���
		HLabel* splitedLabel = new_label_bucket[nodeId];

		IndexType eS = ep.start_;
		IndexType eE = ep.end_;

		//Logger<<"  �ߵ����Ϊ"<<eS<<"�յ�Ϊ"<<eE<<endl;

		IndexType recordS = 0;       
		IndexType recordE = 0;

		IndexType vtxBSzie = (IndexType)splitedLabel->vertex_bucket.size();
		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); iter ++)
		{
			IndexType prev_id = iter->second->prev_corr->label_parent[gLevel]->label_id; //

			if (prev_id == eS)
			{
				recordS ++;

			}else if(prev_id == eE)
			{
				recordE ++;
			}
		}

		ScalarType ration = (ScalarType)(recordS + recordE)/vtxBSzie;


		tempEdge.EdgeDec = ed;
		tempEdge.unMarkedRation = ration;
		tempEdge.srCorNum = recordS;
		tempEdge.tgCorNum = recordE;


		orderedEdgeQ.push(tempEdge);
	}

	Logger<<"����������.\n";
	//��Ҫɾ��������ڴ�!!!
}

void DualwayPropagation::copyLabelBucket(vector<HLabel*>& leftLabels, const vector<HLabel*>& oriLabels)
{
	assert(oriLabels.size() > 0);

	for (auto iter = oriLabels.begin(); iter != oriLabels.end(); iter ++)
	{
		IndexType LId = (*iter)->label_id;

		HLabel* new_label_space = allocator_.allocate<HLabel>();
		HLabel* new_label = new (new_label_space)HLabel((*iter)->label_id,(*iter)->frame_parent,(*iter)->vertex_bucket,(*iter)->prev_corr,(*iter)->next_corr);

		// 		for (auto cvIter = new_label->vertex_bucket.begin(); cvIter!= new_label->vertex_bucket.end(); cvIter++)
		// 		{
		// 			cvIter->second->label_parent = new_label;
		// 		}

		leftLabels.push_back(new_label);
	}
}

bool DualwayPropagation::checkPrevLabelBucket(map<IndexType,HVertex*>& edgePs, map<IndexType,HVertex*>& edgePsCor, IndexType& isSplit)
{
	assert(edgePs.size() > 0);

	set<IndexType > labelS;
	labelS.clear();

	for (auto iter = edgePs.begin(); iter != edgePs.end(); iter ++)
	{
		IndexType labParSize =  (IndexType)iter->second->prev_corr->label_parent.size();

		IndexType updateLevel = labParSize - 1;

		IndexType corPId = iter->second->prev_corr->vtx_id;

		IndexType nextVtxLabid = iter->second->prev_corr->label_parent[updateLevel]->label_id;

		labelS.insert(nextVtxLabid);

		HVertex* nextP = iter->second->prev_corr;

		edgePsCor[corPId] = nextP;

	}

	if (labelS.size() > 1)
	{
		Logger<<"�߽��ָ�����¸��������Label!.\n";

	}

	isSplit = (* labelS.begin() );

	return (labelS.size() == 1);
}

void DualwayPropagation::buildKdTree( IndexType _cframeId)
{
	if (downSample != NULL)
	{
		delete downSample;
	}

	downSample =new Sample();
	Sample& smp = SampleSet::get_instance()[_cframeId];
	auto iter = hier_componets_[_cframeId].label_of_vtx.begin();
	IndexType i = 0;
	for(; iter != hier_componets_[_cframeId].label_of_vtx.end(); ++iter,i++ )
	{
		IndexType vtx_id = (*iter).first;
		Vertex& vtx = smp[vtx_id];

		PointType v( vtx.x(), vtx.y(), vtx.z() );
		ColorType cv(vtx.r(), vtx.g(), vtx.b(), vtx.alpha());
		NormalType nv(vtx.nx(), vtx.ny(), vtx.nz());

		downSample->add_vertex(v,nv,cv);
		o2d_map[vtx_id] = i;
		d2o_map[i] = vtx_id;

	}

	downSample->build_kdtree();

}

ScalarType DualwayPropagation::p2PatchGeoDis(IndexType cFrame, HVertex& oriP,map<IndexType,HVertex*>& parthPs)
{

	ScalarType avgDis = 1e5;

	ScalarType randDis = 0.0;

	IndexType pSize = (IndexType)parthPs.size();

	assert(pSize > 0);

	IndexType srPId = oriP.vtx_id;
	IndexType nodeId = hier_componets_[cFrame].gId_of_vtx[srPId];

	ScalarType geoDis = 0.0 ;

	PCloudGraph* pg = hier_componets_[cFrame].pcGraph;

	IndexType nSize = (IndexType)boost::num_vertices(*pg);

	//����������
	pair<pcVertexIterator, pcVertexIterator> vi = boost::vertices(*pg);

	pcVertexIterator nodeIter = (vi.first + nodeId);

	pcVertexDescriptor nodeDesc = *nodeIter;

	vector<ScalarType> djDis(boost::num_vertices(*pg), 0.);
	vector<pcVertexDescriptor> parents(boost::num_vertices(*pg) );
	//djDis.resize();

	auto p_map = boost::make_iterator_property_map(&parents[0], boost::get(boost::vertex_index, *pg));
	auto w_map = boost::get(&PCEdgeProperty::dist,*pg);
	auto d_map = boost::make_iterator_property_map(&djDis[0],boost::get(boost::vertex_index,*pg));

	boost::dijkstra_shortest_paths(*pg,nodeDesc,boost::weight_map(w_map).predecessor_map(p_map).distance_map(d_map) );



	auto bIter = parthPs.begin();

	auto eIter = parthPs.end();


	for (; bIter != eIter; ++bIter)
	{
		IndexType tgPId = bIter->first;
		IndexType tgNodeId = hier_componets_[cFrame].gId_of_vtx[tgPId];

		randDis = djDis[tgNodeId];

		if(randDis < avgDis)
		{
			avgDis = randDis;
		}

	}

	return  avgDis;
}

void DualwayPropagation::mergePatchesAfterCoSeg() //0831
{

	//����merge��ֹ���� ,��ȡCo-segmentation��label�ļ�֮������ú���

	buildPatchCorrespondenceByLabel();

	IndexType itNum = 9;
	IndexType i = 0;

	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end() &&  i < itNum ; ++ fIter,++i)
	{
		IndexType frameId = fIter->first;
		if (hier_componets_.find(frameId + 1) != hier_componets_.end() )
		{
			IndexType combineNum = 1;//�ϲ�����

			while (combineNum -- > 0) 
			{
				set<IndexType> srBestSet,tgBestSet;
				srBestSet.clear();
				tgBestSet.clear();

				IndexType srLevel = (IndexType)hier_componets_[frameId].hier_graph.size();
				IndexType tgLevel = (IndexType)hier_componets_[frameId + 1].hier_graph.size();

				GraphMatch pairFrameSimilar(SampleSet::get_instance(), hier_componets_, frameId, frameId + 1);

				//pairFrameSimilar.findBestPatches(srLevel - 1, tgLevel - 1, srBestSet,tgBestSet);

				pairFrameSimilar.mergePatchesAfterCoSeg(srLevel - 1, tgLevel - 1,srBestSet,tgBestSet);

			}
		}
	}

}

void DualwayPropagation::buildPatchCorrespondenceByLabel()
{
	for (auto fiter = hier_componets_.begin(); fiter != hier_componets_.end(); ++fiter )
	{
		IndexType frame_id = fiter->first;
		IndexType gLevel = (IndexType)fiter->second.hier_label_bucket.size();
		--gLevel;

		vector<HLabel*>& labelBucket = fiter->second.hier_label_bucket[gLevel];

		if (hier_componets_.find(frame_id + 1) != hier_componets_.end() )
		{

			map<IndexType,IndexType> sr_labelIndexMap = fiter->second.hier_label_vtxBucket_index[gLevel];

			IndexType tLevel = (IndexType)hier_componets_[frame_id + 1].hier_label_vtxBucket_index.size();
			--tLevel;

			map<IndexType,IndexType> tg_labelIndexMap = hier_componets_[frame_id + 1].hier_label_vtxBucket_index[tLevel];

			vector<HLabel*>& tg_labelBucket = hier_componets_[frame_id + 1].hier_label_bucket[tLevel];


			for ( auto mIter = sr_labelIndexMap.begin(); mIter != sr_labelIndexMap.end(); ++mIter )
			{
				IndexType labelId = mIter->first;
				IndexType realSrLabIdx = mIter->second;

				if (tg_labelIndexMap.find(labelId) != tg_labelIndexMap.end() )//������ͬ��label
				{
					IndexType corLabelId = tg_labelIndexMap[labelId];

					labelBucket[realSrLabIdx]->next_corr = tg_labelBucket[corLabelId]; //next correspondence

					tg_labelBucket[corLabelId]->prev_corr = labelBucket[realSrLabIdx]; //prev correspondence
				}
			}

		}
	}
}

DualwayPropagation::DualwayPropagation()
{
	downSample = NULL;
}

DualwayPropagation::~DualwayPropagation()
{
	if( downSample){
		delete downSample;
	}
}
