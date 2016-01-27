#include "TrajectoryClassifier.h"
#include "DeformableRegistration.h"
#include "J_linkage_Matlab.hpp"
#include "traj2AffineModel_distance.h"
#include "traj2model_distance.h"
#include "bfs_classifier.hpp"
using namespace  std;
void TrajectoryClassifier::run()
{
	SampleSet& set = SampleSet::get_instance();
	if (set.empty())
	{
		Logger<< " End Clustering(Warning: Empty Set).\n";
		return;
	}
	//Step 1: compute rotate feature
	DeformableRegistration nonrigid;

	vector<PCloudTraj> totalTraj;

	vector<IndexType> sampleCenterVtxId;

	totalTraj.clear();

	sampleCenterVtxId.clear();
	if (isEqual)
	{

		Logger<<"Using Equal Length Traj"<<endl;

		nonrigid.calculateFixedLengthTrajWithTracingAlong(totalTraj,centerFrame,sampleCenterVtxId,trajLen,octreeRes);

	}else
	{
		Logger<<"Using Unequal Length Traj"<<endl;
		nonrigid.calculateDownSmpLifeSpanTrajCenter(totalTraj,centerFrame,sampleCenterVtxId,trajLen,octreeRes,threshold,lifeT); 
	}


	//calculate motion models

	IndexType modeNum = modelT * sampleCenterVtxId.size();
	vector<IndexType> labels;
	labels.resize(totalTraj.size(),0);

	if (totalTraj.size() < 50)
	{
		Logger<<"Traj size is small!.\n";
		Logger<< " End Clustering.\n";
		return;
	}

	//rigid or affine motion model

	if (isRigid)
	{
		Logger<<"Using Rigid motion  model"<<endl;

		vector<PCloudModel> totalModel;

		totalModel.clear();

		nonrigid.sampleModel(totalTraj,totalModel,modeNum);

		J_LinkageAdapter_Matlab<PCloudTraj, PCloudModel, Traj2ModDistanceFunc>	algo_adaptor(totalTraj, totalModel, labels, Traj2ModDistanceFunc(threshold),perC);

		algo_adaptor.compute();

	}else
	{
		Logger<<"Using Affine motion  model"<<endl;

		vector<PCloudAffModel> totalModel;

		totalModel.clear();

		nonrigid.sampleAffineModel(totalTraj,totalModel,modeNum);

		J_LinkageAdapter_Matlab<PCloudTraj, PCloudAffModel, Traj2AffModelDistFunc>	algo_adaptor(totalTraj, totalModel, labels, Traj2AffModelDistFunc(threshold),perC);

		algo_adaptor.compute();

	}

	Logger<<"centerFrame = "<<centerFrame<<endl;


	Sample& sample0 = set[centerFrame];

	vector<bool> isSelector(sample0.num_vertices(),false);

	//in order to visual sample vertexes 

	for (int i = 0; i < sampleCenterVtxId.size(); i++)
	{
		isSelector[sampleCenterVtxId[i]] = true;
	}

	if (isEqual)
	{
#ifdef SAVE_CORRESPONDENCE

		IndexType smp_range[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,
			51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,
			101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,
			141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173};

		char corr_file_name[1024];



		sprintf(corr_file_name,"G:\\Data\\horse\\quaEva1215\\J-linkage threshold test\\hksingle_corr%.2d_%.2f.txt",centerFrame,perC);

		FILE *in_correspond = fopen(corr_file_name,"w");

		for ( int i=0; i<sampleCenterVtxId.size();i++ )
		{
			for (int j=0; j<(sizeof(smp_range)/sizeof(IndexType));j++)
			{
				if(centerFrame==smp_range[j] || totalTraj[i].trajLifeSpan.start > smp_range[j] )continue;
				if(totalTraj[i].trajLifeSpan.end < smp_range[j])break;
				IndexType v = sampleCenterVtxId[i];

				fprintf(in_correspond,"%d %d %d %d\n", centerFrame, v, smp_range[j], totalTraj[i].trajNode[smp_range[j] - totalTraj[i].trajLifeSpan.start]);

			}
		}


		fclose(in_correspond);

#endif

		bubleSort(sampleCenterVtxId,labels,labels.size());

		vector<IndexType> label_smooth(labels.size(),0);

		diff_using_bfs(labels,sampleCenterVtxId,centerFrame);

		nonrigid.smoothSmapleLabel_KDTree(sample0,sampleCenterVtxId,labels,label_smooth);

		IndexType nLabels = 1;

		// in order to set the label system continuous and the 0 as the first label
		nLabels = orderLabels(label_smooth);

		Logger<<"seg size = "<<nLabels; 

		vector<IndexType> result_label(sample0.num_vertices(),0);

		// propagate the sample vertexes to original point cloud.
		nonrigid.propagateLabel2Orignal(sample0,sampleCenterVtxId,label_smooth,result_label); 

#ifdef SAVE_LABELS

		char label_labsmooth[1024];


		sprintf(label_labsmooth,"G:\\Data\\horse\\quaEva1215\\J-linkage threshold test\\orilabels%.2d_%.2f.txt",centerFrame,perC);

		FILE *in_label_smooth = fopen(label_labsmooth, "w");

		IndexType tpd = 0;

		for ( int i=0; i<set[centerFrame].num_vertices(); i++ )
		{
			fprintf( in_label_smooth, "%d %d %d\n", centerFrame, result_label[i], i );
		}

		fclose(in_label_smooth);

#endif

#ifdef SAVE_CORRESPONDENCE

		//visualize the segmentation results in point cloud.
		IndexType i = 0;
		IndexType k = 0;
		for (Sample::vtx_iterator v_iter = sample0.begin();
			v_iter != sample0.end();
			v_iter++,i++ )
		{
			if (isSelector[i])
			{
				(*v_iter)->set_visble(true);

				(*v_iter)->set_label( label_smooth[k] );//smooth label
				k++;
			}else
			{
				(*v_iter)->set_visble(false);
			}
		}
#endif 

	}




	//IndexType frames[] = {1,2,3,4,5,6,7,8,9,10/*,11,12,13,14,15,16,17,18,19,20*/};
	if ( !isEqual)
	{
		for ( IndexType i= (centerFrame - trajLen/2 ); i<= (centerFrame + trajLen/2 ); i++ )
		{
			for (IndexType j=0; j<set[i].num_vertices(); j++)
			{
				set[i][j].set_visble(false);
			}
		}

		map<IndexType, IndexType> inter_labels;
		IndexType trajSize = totalTraj.size();
		for (IndexType tId = 0; tId < trajSize; tId ++) 
		{
			IndexType nodeId = 0;

			for (IndexType stF = totalTraj[tId].trajLifeSpan.start; stF <= totalTraj[tId].trajLifeSpan.end; stF ++, nodeId ++)
			{
				if ( centerFrame == stF)
				{

					inter_labels[totalTraj[tId].trajNode[nodeId] ] = labels[tId];

				}

			}
		}

		vector<IndexType> order_labels;
		order_labels.resize(sampleCenterVtxId.size(),0);

		for (IndexType i= 0; i < sampleCenterVtxId.size(); i ++)
		{
			if (inter_labels.find(sampleCenterVtxId[i]) != inter_labels.end() )
			{
				order_labels[i] = inter_labels[sampleCenterVtxId[i] ];

			}else
			{
				Logger<<"Unmark point.\n";
			}

		}

		bubleSort(sampleCenterVtxId,order_labels,order_labels.size());
		vector<IndexType> label_smooth(order_labels.size(),0);

		diff_using_bfs(order_labels,sampleCenterVtxId,centerFrame); 
		nonrigid.smoothSmapleLabel_KDTree(sample0,sampleCenterVtxId,order_labels,label_smooth);

		IndexType nLabels = 0;

		nLabels = orderLabels(label_smooth);

		Logger<<"seg size = "<<nLabels; 

		IndexType i = 0;
		IndexType k = 0;
		for (Sample::vtx_iterator v_iter = sample0.begin();
			v_iter != sample0.end();
			v_iter++,i++ )
		{
			if (isSelector[i])
			{
				(*v_iter)->set_visble(true);
				(*v_iter)->set_label( label_smooth[k] );//smooth label
				k++;
			}else
			{
				(*v_iter)->set_visble(false);
			}
		}

	}
	Logger<< " End Clustering.\n";

}
void TrajectoryClassifier::setParamter(IndexType _trajLen,IndexType _octreeReso,ScalarType _perC,
	ScalarType _thresHold,IndexType _modelT, IndexType _smallL,bool _isEqual, bool _isRigid)
{
	trajLen = _trajLen;
	octreeRes = _octreeReso;
	perC = _perC;
	threshold = _thresHold;
	modelT = _modelT;
	lifeT = _smallL;
	isEqual = _isEqual;
	isRigid = _isRigid;

}
void TrajectoryClassifier::setNeigNum(IndexType _neigbNum)
{
	neigborNum = _neigbNum;
}

void TrajectoryClassifier::bubleSort(vector<IndexType>& oriData,vector<IndexType>& labels,IndexType lSize)
{
	int temp;
	IndexType labelTemp = 0;
	bool flag=false;
	for (int i=0;i<lSize;i++)
	{
		flag=true;
		for (int j=0;j<lSize-i-1;j++)
		{
			if(oriData[j]>oriData[j+1])
			{
				temp=oriData[j];
				labelTemp = labels[j];

				oriData[j]=oriData[j+1];
				labels[j] = labels[j+1];

				oriData[j+1]=temp;
				labels[j+1] = labelTemp;

				flag = false;
			}
		}
		if(flag) break;
	}
}

void TrajectoryClassifier::diff_using_bfs(std:: vector<IndexType>& labels,std::vector<IndexType>& centerVtxId,IndexType centerFrame )
{
	SampleSet& set = SampleSet::get_instance();
	IndexType max_label= *max_element(labels.begin(),labels.end());
	vector< std::set<IndexType> > label_bucket(max_label+1);

	//IndexType centerFrame = 10;
	for ( int i=0; i<labels.size(); i++ )
	{
		label_bucket[labels[i]].insert( i );
	}
	IndexType new_label = max_label;
	Logger<<"max label before post:"<<new_label<<endl;
	for ( IndexType l=0; l<label_bucket.size(); l++ )
	{
		std::set<IndexType>& idx_set = label_bucket[l];
		if (idx_set.size()==0)
		{
			continue;
		}
		IndexType idx_set_size = idx_set.size();
		Matrix3X vtx_data;
		vtx_data.setZero( 3, idx_set_size );
		size_t i=0;
		for (std::set<IndexType>::iterator iter = idx_set.begin();
			iter != idx_set.end();
			iter++,i++)
		{
			IndexType a = *iter;
			vtx_data(0,i)  = set[centerFrame][centerVtxId[*iter]].x();
			vtx_data(1,i)  = set[centerFrame][centerVtxId[*iter]].y();
			vtx_data(2,i)  = set[centerFrame][centerVtxId[*iter]].z();
		}


#ifdef USE_RADIUS_NEAREST
		ScalarType rad = set[centerFrame].box_diag();
		BFSClassifier<ScalarType> classifier(vtx_data,rad);
#else
		ScalarType rad = set[centerFrame].box_diag()/10;
		BFSClassifier<ScalarType> classifier(vtx_data,rad,10);
#endif

		classifier.run();
		int *sep_label = classifier.get_class_label();
		i=0;
		for (std::set<IndexType>::iterator iter = idx_set.begin();
			iter != idx_set.end();
			iter++,i++)
		{
			if(sep_label[i]==0)continue;
			labels[*iter] = new_label + sep_label[i];
		}
		new_label += classifier.get_num_of_class()-1;


	}
	Logger<<"max label after post:"<<new_label<<endl;
}

int TrajectoryClassifier::orderLabels(vector<IndexType>& labels)
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
