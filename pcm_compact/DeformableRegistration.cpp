#include "DeformableRegistration.h"
#include"sample_properity.h"
#include "Eigen/Dense"
using namespace std;



//sampleVtxId were oridered!
void DeformableRegistration::propagateLabel2Orignal(Sample& oriPC,vector<IndexType>& sampleVtxId,vector<IndexType>& label_smp,vector<IndexType>& label_ori)
{
	IndexType nCluster = -1;
	for (int i = 0; i < label_smp.size(); i++)
	{
		if(nCluster < label_smp[i])
		{
			nCluster = label_smp[i];
		}
	}
	Logger<<"The number of segments are "<<nCluster + 1<<endl;

	map<IndexType,IndexType> smpLabel;
	map<IndexType,IndexType>::iterator IsValidIter;
	for (int i = 0; i < label_smp.size(); i++)
	{
		smpLabel.insert(make_pair(sampleVtxId[i],label_smp[i]));
	}

	const IndexType k = 300;
	IndexType neighbours[k];
	ScalarType dist[k];

	IndexType vtx_num = (IndexType)oriPC.num_vertices();
	IndexType result_label;


	for(IndexType vtx_id = 0; vtx_id < vtx_num; vtx_id ++)
	{
		vector<IndexType> recordLabelTime(nCluster + 1,0);
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
		for (int i = 0; i <= nCluster; i++)
		{
			if(result_label < recordLabelTime[i])
			{
				label_ori[vtx_id] = i;
				result_label = recordLabelTime[i];
			}
		}

	}
	Logger<<"k-neig's k = "<<k<<endl;
	Logger<<"End Propagate!"<<endl;
}

//------------------------
void DeformableRegistration::smoothSmapleLabel_KDTree(Sample& oriPC,vector<IndexType>& sampleVtxId,vector<IndexType>& label_smp,vector<IndexType>& label_smooth)
{

	Sample* downSmp = new Sample;

	for (auto iter = sampleVtxId.begin(); iter != sampleVtxId.end(); iter++)
	{
		Vertex& vtx = oriPC[*iter];

		PointType v( vtx.x(), vtx.y(), vtx.z() );
		ColorType cv(vtx.r(), vtx.g(), vtx.b(), vtx.alpha());
		NormalType nv(vtx.nx(), vtx.ny(), vtx.nz());

		downSmp->add_vertex(v,nv,cv);
	}

	downSmp->build_kdtree();


	IndexType nCluster = -1;

	for (int i = 0; i < label_smp.size(); i++)
	{
		if(nCluster < label_smp[i])
		{
			nCluster = label_smp[i];
		}
	}

	Logger<<"The number of segments are "<<nCluster + 1<<endl;


	const IndexType k = 15;

	IndexType neighbours[k];

	ScalarType dist[k];

	IndexType result_label;

	IndexType s_id = 0;


	for (IndexType i = 0; i < sampleVtxId.size(); i ++)
	{
		vector<IndexType> recordLabelTime(nCluster + 1,0);

		result_label = -1;

		downSmp->neighbours(i, k, neighbours, dist);

		for(IndexType neig_id = 0; neig_id < k; neig_id ++)
		{
			recordLabelTime[ label_smp [neighbours[neig_id] ] ] += 1;
		}

		for (int j = 0; j <= nCluster; j++)
		{
			if(result_label < recordLabelTime[j])
			{
				label_smooth[i] = j;

				result_label = recordLabelTime[j];
			}
		}
	}

	Logger<<"k-neig's k = "<<k<<endl;

	Logger<<"End  Sample Propagate!"<<endl;


	delete downSmp;

}


//---------------------------Affine transformation
bool DeformableRegistration::sampleAffineModel(vector<PCloudTraj>& oriTraj,vector<PCloudAffModel>& smpOriModel,IndexType modelNum)
{
	LifeSpan commonLifeSp;
	vector<PCloudTraj> resTraj;
	PCloudAffModel resModel;

	while (modelNum-- > 0) 
	{
		while (true)
		{
			if (sampleValidTraj(oriTraj,resTraj,commonLifeSp,4) )
			{
				//Logger<<"sample success!.\n";
				break;
			}
			resTraj.clear();
			commonLifeSp.start = 0;
			commonLifeSp.end   = 0;
		}
		/*resModel = constrcutRigidMotionMode(resTraj,commonLifeSp);*/
		resModel = constrcutAffineMotionMode(resTraj,commonLifeSp);
		smpOriModel.push_back(resModel);
		resTraj.clear();
		commonLifeSp.start = 0;
		commonLifeSp.end   = 0;
	}
	// 
	Logger<<"sample mode size = "<<smpOriModel.size()<<endl;
	return true;
}
//---------------------------

bool DeformableRegistration::sampleModel(vector<PCloudTraj>& oriTraj,vector<PCloudModel>& smpOriModel,IndexType modelNum)
{
	LifeSpan commonLifeSp;
	vector<PCloudTraj> resTraj;
	PCloudModel resModel;

	while (modelNum-- > 0) 
	{
		while (true)
		{
			if (sample3validTraj(oriTraj,resTraj,commonLifeSp))
			{
				//Logger<<"sample success!.\n";
				break;
			}
			resTraj.clear();
			commonLifeSp.start = 0;
			commonLifeSp.end   = 0;
		}
		resModel = constrcutRigidMotionMode(resTraj,commonLifeSp);
		smpOriModel.push_back(resModel);
		resTraj.clear();
		commonLifeSp.start = 0;
		commonLifeSp.end   = 0;
	}
	// 
	Logger<<"sample mode size = "<<smpOriModel.size()<<endl;
	return true;
}

//------------------------------------------从中间开始,向前后帧传播, 生成不等长轨迹
void DeformableRegistration::calculateDownSmpLifeSpanTrajCenter(vector<PCloudTraj>& trajs, IndexType smp_idx,vector<IndexType>& centerVtxMap, IndexType length,IndexType octreeRes ,
	ScalarType dfThreshold, IndexType LifeMin)
{
	Logger<<"start Calculate  down sample life-span Traj from  Center!.\n";

	SampleSet& smp_set = SampleSet::get_instance();
	IndexType start_idx = smp_idx - length/2 ;
	IndexType end_idx = smp_idx + length/2;

	if( length%2!=0 )
		++end_idx;
	if( start_idx<0 )
	{
		start_idx = 0;
		end_idx = length;
	}
	if( end_idx > smp_set.size() - 1  )
	{
		end_idx = (IndexType)smp_set.size() - 1;
		start_idx = end_idx - length;
	}

	vector<Sample*> p_subsmp_set( end_idx - start_idx + 1);
	vector< vector<IndexType> > vtx_map(end_idx - start_idx + 1);

	for (int i = 0; i < end_idx - start_idx+1; i++)
	{
		p_subsmp_set[i] = new Sample;
		SampleManipulation::SubSampleSolution solv( smp_set[start_idx+i], *p_subsmp_set[i],octreeRes);
		vtx_map[i] = solv.compute();
	}

	centerVtxMap = vtx_map[smp_idx - start_idx];

	set<IndexType> validEnd;//valid end index
	set<IndexType> nextValidEnd;//product new valid end index 
	vector<PCloudTraj> willGoOnTraj;
	vector<PCloudTraj> nextwillGoOnTraj;
	set<IndexType>::iterator isInValidIter;

	IndexType vtx_cor = 0;
	IndexType willGoSetSize = 0;
	PCloudTraj tempTraj;
	validEnd.clear();
	Matrix3X srVtxMatrix;
	vector<bool> onceLabel;

	IndexType invalidTraj = 0;


	// forward trace! life spans trajectories  ------>>

	for ( int i = (smp_idx - start_idx); i < (end_idx - start_idx); i++ )
	{

		Sample* srGraph = p_subsmp_set[i]; //对应原始的 ( start_idx + i) 帧

		Sample* tgGraph = p_subsmp_set[i + 1];//对应原始的 ( start_idx + i + 1) 帧

		Matrix3X srVtxMatrix = srGraph->vertices_matrix();

		MatrixXXi   global_vtx_map;

		MatrixXX    isDeformable;

		std::set<IndexType> distoredLargeSet;


		non_rigidRegister(*srGraph,*tgGraph,global_vtx_map,isDeformable,distoredLargeSet);


		srGraph->vertices_matrix() = srVtxMatrix; 

		srGraph->update();//back to original position after registeration

		IndexType srGraphVtxN = (IndexType)srGraph->num_vertices();

		for ( int j = 0; j < srGraphVtxN; j++)
		{
			//j号点为i帧中的实际顶点号码为: vtx_map[i][j]
			//j 点在下一帧的对应点为 vtx_map[i + 1][global_vtx_map(0,j)]
			vtx_cor = vtx_map[i + 1][global_vtx_map(0,j)];
			isInValidIter = validEnd.find(vtx_map[i][j]);

			ScalarType defValue = isDeformable(0,(__int64)j);

			if (isInValidIter != validEnd.end())
			{

				vector<PCloudTraj>::iterator v_iter =willGoOnTraj.begin();
				for (; v_iter != willGoOnTraj.end();)
				{
					IndexType trajLife = (IndexType)(*v_iter).trajNode.size();	
					IndexType trajEndIndex = (*v_iter).trajNode[trajLife - 1];
					IndexType endF = (*v_iter).trajLifeSpan.end;

					if (endF == (start_idx + i) && trajEndIndex == vtx_map[i][j])
					{

						if (defValue <= dfThreshold)
						{

							(*v_iter).trajLifeSpan.end += 1;//Ti.end = Ti.end + 1;
							(*v_iter).trajNode.push_back(vtx_cor);//轨迹的最终顶点号码追加为对应号码
							nextwillGoOnTraj.push_back(*v_iter);// 没有断,继续走下去
							nextValidEnd.insert(vtx_cor);
							v_iter = willGoOnTraj.erase(v_iter); //discard

						}else
						{
							if ( (*v_iter).trajLifeSpan.start == smp_idx )//若起点是中心帧,则保留,为后向传播做准备
							{

								v_iter ++;

							}else
							{

								if (((*v_iter).trajLifeSpan.end - (*v_iter).trajLifeSpan.start) >= LifeMin)//加入轨迹长度判断
								{
									trajs.push_back(*v_iter);//保存轨迹        
								}

								v_iter = willGoOnTraj.erase(v_iter); //前向不能延续,直接抛弃该轨迹!

							}

						}

					}else
					{

						v_iter ++;

					}
				}

			}else
			{
				if (defValue <= dfThreshold)
				{
					tempTraj.trajLifeSpan.start = (start_idx + i);
					tempTraj.trajLifeSpan.end = (start_idx + i + 1);
					tempTraj.trajNode.clear();

					tempTraj.trajNode.push_back(vtx_map[i][j]);

					tempTraj.trajNode.push_back(vtx_cor);

					nextwillGoOnTraj.push_back(tempTraj);

					nextValidEnd.insert(vtx_cor);

				}else
				{
					invalidTraj++;
				}
			}
		}//end for very vtx in current frame

		Logger<<i <<"Forward! frame begin is end, never stay this world num size = " <<invalidTraj<<endl;
		validEnd.clear();
		validEnd = nextValidEnd;
		nextValidEnd.clear();

		willGoOnTraj.clear();
		willGoOnTraj  = nextwillGoOnTraj;
		nextwillGoOnTraj.clear();

		invalidTraj = 0;
		// save willGOTraj for backward trace!
		// 		if ((end_idx - 1) == (i + start_idx) )
		// 		{		
		// 			vector<PCloudTraj>::iterator willStart = willGoOnTraj.begin();
		// 			vector<PCloudTraj>::iterator willEnd = willGoOnTraj.end();
		// 			for (; willStart != willEnd; willStart ++)
		// 			{
		// 				if (((*willStart).trajLifeSpan.end - (*willStart).trajLifeSpan.start) >= LifeMin)
		// 				{
		// 					trajs.push_back(*willStart);
		// 				}
		// 			}           
		// 		}

	}


	// Backward !!!!!!!      <<-----------------------------------

	validEnd.clear();

	vector<PCloudTraj>::iterator vw_iter =willGoOnTraj.begin();

	for (; vw_iter != willGoOnTraj.end(); vw_iter ++) //保存仍然需要延续轨迹起始的顶点索引,使得在遍历中心帧时需不要重新生成新的轨迹
	{
		if ((*vw_iter).trajLifeSpan.start == smp_idx)
		{
			validEnd.insert((*vw_iter).trajNode[0]);
		}
	}

	for ( int i = (smp_idx - start_idx); i > 0; i --)
	{
		Sample* srGraph = p_subsmp_set[i]; //对应原始的 ( start_idx + i) 帧

		Sample* tgGraph = p_subsmp_set[i - 1];//对应原始的 ( start_idx + i + 1) 帧

		Matrix3X srVtxMatrix = srGraph->vertices_matrix();

		MatrixXXi   global_vtx_map;

		MatrixXX    isDeformable;

		std::set<IndexType> distoredLargeSet;


		non_rigidRegister(*srGraph,*tgGraph,global_vtx_map,isDeformable,distoredLargeSet);


		srGraph->vertices_matrix() = srVtxMatrix; 

		srGraph->update();//back to original position after registeration

		IndexType srGraphVtxN = (IndexType)srGraph->num_vertices();

		for ( int j = 0; j < srGraphVtxN; j++)
		{
			//j号点为i帧中的实际顶点号码为: vtx_map[i][j]
			//j 点在下一帧的对应点为 vtx_map[i + 1][global_vtx_map(0,j)]
			vtx_cor = vtx_map[i - 1][global_vtx_map(0,j)];
			isInValidIter = validEnd.find(vtx_map[i][j]);

			ScalarType defValue = isDeformable(0,(__int64)j);


			if (isInValidIter != validEnd.end())
			{
				vector<PCloudTraj>::iterator v_iter =willGoOnTraj.begin();
				for (; v_iter != willGoOnTraj.end();)
				{
					// 					IndexType trajLife = (*v_iter).trajNode.size();	
					// 					IndexType trajEndIndex = (*v_iter).trajNode[trajLife - 1];

					IndexType trajStartIndex = (*v_iter).trajNode[0]; //获得轨迹的起点顶点索引
					IndexType strF = (*v_iter).trajLifeSpan.start;

					if ( strF == (start_idx + i)  && trajStartIndex == vtx_map[i][j])
					{

						if (defValue <= dfThreshold)
						{

							(*v_iter).trajLifeSpan.start -= 1;//向前推进一帧

							vector<IndexType>::iterator itB = (*v_iter).trajNode.begin();

							(*v_iter).trajNode.insert(itB,vtx_cor);//在起始位置插入对应点索引

							nextwillGoOnTraj.push_back(*v_iter);// 没有断,继续走下去

							nextValidEnd.insert(vtx_cor);

							v_iter = willGoOnTraj.erase(v_iter); //discard

						}else
						{

							if (((*v_iter).trajLifeSpan.end - (*v_iter).trajLifeSpan.start) >= LifeMin)//加入轨迹长度判断
							{
								trajs.push_back(*v_iter);//保存轨迹        
							}

							v_iter = willGoOnTraj.erase(v_iter); //前向不能延续,直接抛弃该轨迹!
						}

					}else //该轨迹不以该点为结束且不以(start+ i)为起始帧,"那么它永远不会续约下去,进行判断,要么放入最终的轨迹中,要么抛弃." 上面已经抛弃了一次!
					{
						v_iter ++;
					}

				}// end for

			}else //生成一条新的轨迹
			{
				if (defValue <= dfThreshold)
				{
					tempTraj.trajLifeSpan.start = (start_idx + i - 1);
					tempTraj.trajLifeSpan.end = (start_idx + i);
					tempTraj.trajNode.clear();

					tempTraj.trajNode.push_back(vtx_cor);

					tempTraj.trajNode.push_back(vtx_map[i][j]);

					nextwillGoOnTraj.push_back(tempTraj);

					nextValidEnd.insert(vtx_cor);

				}else
				{
					invalidTraj++;
				}
			}
		}//end for very vtx in current frame

		Logger<<i <<"Backward! frame begin is end, never stay this world num size = " <<invalidTraj<<endl;
		validEnd.clear();
		validEnd = nextValidEnd;
		nextValidEnd.clear();

		willGoOnTraj.clear();
		willGoOnTraj  = nextwillGoOnTraj;
		nextwillGoOnTraj.clear();


		if ( 1 == i ) //???
		{		
			vector<PCloudTraj>::iterator willStart = willGoOnTraj.begin();
			vector<PCloudTraj>::iterator willEnd = willGoOnTraj.end();
			for (; willStart != willEnd; willStart ++)
			{
				if (((*willStart).trajLifeSpan.end - (*willStart).trajLifeSpan.start) >= LifeMin)
				{
					trajs.push_back(*willStart);
				}
			}           
		}

		invalidTraj = 0;

	}

	//draw life-spans  trajectories
//	drawLifespansTraj(trajs);

	Logger<<"End Calculate  down sample life-span Traj.\n";

}
//------------------------------------------

void DeformableRegistration::calculateFixedLengthTrajWithTracingAlong(vector<PCloudTraj>& trajs, IndexType smp_idx,vector<IndexType>& centerVtxMap,IndexType length, IndexType octreeRes )
{

	Logger<<"start Calculate Traj,using concenstive rigid.\n";



	SampleSet& smp_set = SampleSet::get_instance();

	IndexType start_idx = smp_idx - length/2 ;

	IndexType end_idx = smp_idx + length/2;



	if( length%2!=0 )

		++end_idx;

	if( start_idx<0 )

	{

		start_idx = 0;

		end_idx = length;

	}

	if( end_idx > smp_set.size() - 1  )

	{

		end_idx = (IndexType)smp_set.size() - 1;

		start_idx = end_idx - length;

	}



	vector<Sample*> p_subsmp_set( end_idx - start_idx + 1);

	vector< vector<IndexType> > vtx_map(end_idx - start_idx + 1);



	for (int i = 0; i < end_idx - start_idx+1; i++)

	{

		p_subsmp_set[i] = new Sample;

		SampleManipulation::SubSampleSolution solv( smp_set[start_idx+i], *p_subsmp_set[i],octreeRes);

		vtx_map[i] = solv.compute();

	}



	centerVtxMap = vtx_map[smp_idx - start_idx];



	trajs.clear();

	trajs.resize( p_subsmp_set[smp_idx - start_idx]->num_vertices() );

	for ( int i=0; i<trajs.size();i++ )

	{

		trajs[i].trajLifeSpan.start = start_idx;

		trajs[i].trajLifeSpan.end = end_idx;

		trajs[i].trajNode.resize( end_idx-start_idx+1 );

	}



	for ( int i=0; i<trajs.size();i++ )

	{

		trajs[i].trajNode[smp_idx-start_idx] = vtx_map[smp_idx-start_idx][i];

	}



	vector<int> pre_idx_map(trajs.size());

	for ( int i=0; i<trajs.size(); ++i )

	{

		pre_idx_map[i] = i;

	}



	//backward trace

	for ( int j = smp_idx; j>start_idx; --j )

	{

		Sample* srGraph = p_subsmp_set[j- start_idx];

		Sample* tgGraph = p_subsmp_set[j - start_idx - 1];

		Matrix3X srVtxMatrix = srGraph->vertices_matrix();

		MatrixXXi   global_vtx_map;

		MatrixXX    isDeformable;

		std::set<IndexType> distoredLargeSet;



		non_rigidRegister(*srGraph,*tgGraph,global_vtx_map,isDeformable,distoredLargeSet);



		srGraph->vertices_matrix() = srVtxMatrix; 

		srGraph->update();//back to original position after registeration

		for ( int i=0; i<trajs.size();i++ )

		{

			trajs[i].trajNode[j-start_idx-1] = vtx_map[j-start_idx-1][global_vtx_map(0,pre_idx_map[i])];

			pre_idx_map[i] = global_vtx_map(0,pre_idx_map[i]);

		}



	}


	for ( int i=0; i<trajs.size(); ++i )

	{
		pre_idx_map[i] = i;
	}

	//forward trace

	for ( int j = smp_idx+1; j<=end_idx; ++j )

	{

		Sample* srGraph = p_subsmp_set[j- start_idx-1];

		Sample* tgGraph = p_subsmp_set[j - start_idx];

		Matrix3X srVtxMatrix = srGraph->vertices_matrix();

		MatrixXXi   global_vtx_map;

		MatrixXX    isDeformable;

		std::set<IndexType> distoredLargeSet;



		non_rigidRegister(*srGraph,*tgGraph,global_vtx_map,isDeformable,distoredLargeSet);



		srGraph->vertices_matrix() = srVtxMatrix; 

		srGraph->update();//back to original position after registeration

		for ( int i=0; i<trajs.size();i++ )

		{

			trajs[i].trajNode[j-start_idx] = vtx_map[j-start_idx][global_vtx_map(0,pre_idx_map[i])];

			pre_idx_map[i] = global_vtx_map(0,pre_idx_map[i]);

		}



	}

	//draw traj

//	drawLifespansTraj(trajs);

}

//affine model 7-30
PCloudAffModel DeformableRegistration::constrcutAffineMotionMode(vector<PCloudTraj> & oriTraj,LifeSpan & comTime)
{
	//oriTraj size is 3,and  begin and end  of those trajectories even known
	assert(comTime.end > comTime.start);

	IndexType totFrame = comTime.end - comTime.start + 1;
	PCloudAffModel pModel;
	pModel.modelLifeSpan = comTime;


	Vector3Type mean;
	Matrix34 f_transf;
	Matrix34 b_transf;
	Matrix3X srPoints;
	Matrix3X tgPoints;    

	//calculate rotmate matrix  between frame (_start and _start + 1)
	// rotmate * srPos = tgPos  ----- [3×3] * [3×3] = [3×3];column-major storage vertex position
	getMultiTrajNodeCoordinate(oriTraj,comTime.start,srPoints);


	for (IndexType frameIndex = comTime.start; frameIndex < comTime.end; frameIndex ++)//nodes  number = (end - start)
	{
		getMultiTrajNodeCoordinate(oriTraj,frameIndex + 1,tgPoints);

		f_transf.setZero();
		b_transf.setZero();

		calculateAffineTrans(srPoints,tgPoints,f_transf);

		pModel.fNode.push_back(f_transf);

		calculateAffineTrans(tgPoints,srPoints,b_transf);

		pModel.bNode.push_back(b_transf);

		//affine-2015-7-30--ldlt

		srPoints = tgPoints;
	}
	return pModel; 
}
PCloudModel DeformableRegistration::constrcutRigidMotionMode(PCloudTraj& traj)
{
	/*NormalType sNorm,tNorm;*/
	ScalarType Nx,Ny,Nz;
	IndexType i = 0;
	PCloudModel tempModel;
	tempModel.modeLifeSpan.start = 18;
	tempModel.modeLifeSpan.end = 22;

	for (IndexType startF = traj.trajLifeSpan.start; startF < traj.trajLifeSpan.end; startF ++,i ++)
	{
		Nx = m_smpSet[startF][traj.trajNode[i]].nx();
		Ny = m_smpSet[startF][traj.trajNode[i]].ny();
		Nz = m_smpSet[startF][traj.trajNode[i]].nz();
		NormalType sNorm(Nx,Ny,Nz);

		Nx = m_smpSet[startF + 1][traj.trajNode[i + 1]].nx();
		Ny = m_smpSet[startF + 1][traj.trajNode[i + 1]].ny();
		Nz = m_smpSet[startF + 1][traj.trajNode[i + 1]].nz();
		NormalType tNorm(Nx,Ny,Nz);

		// 		if (sNorm(0,0) != tNorm(0,0) ||sNorm(1,0) != tNorm(1,0) ||sNorm(2,0) != tNorm(2,0) )
		// 		{
		//             Logger<<"debug.\n";
		// 		}

		Matrix33 rotMat;
		calculateNormalRotMat(sNorm,tNorm,rotMat);
		tempModel.fNode.push_back(rotMat);

		calculateNormalRotMat(tNorm,sNorm,rotMat);
		tempModel.bNode.push_back(rotMat);
	}
	return tempModel;
}
//9-30 the follow codes are not test!
bool DeformableRegistration::sample3validTraj(vector<PCloudTraj> & totTraj,vector<PCloudTraj>& resTraj,LifeSpan & lifesp)
{
	//rand()

	IndexType firstTraj = rand()%totTraj.size();
	IndexType secondTraj = rand()%totTraj.size();
	IndexType thirdTraj = rand()%totTraj.size();

	resTraj.push_back(totTraj[firstTraj]);
	resTraj.push_back(totTraj[secondTraj]);
	resTraj.push_back(totTraj[thirdTraj]);
	//
	//Logger<<"first sample traj is "<<firstTraj<<endl;
	//Logger<<"second sample traj is "<<secondTraj<<endl;
	//Logger<<"third sample traj is "<<thirdTraj<<endl;
	//
	return getCommonTime(resTraj,lifesp);

}



void DeformableRegistration::non_rigidRegister(Sample& srGraph,Sample& tgGraph,MatrixXXi & vtxMap,
	MatrixXX &isDeformable,std::set<IndexType> & distLrgIndex)
{
	IndexType num_vtx = (IndexType)srGraph.num_vertices();
	Matrix3X resPosition;
	isDeformable.setZero(1,(__int64)num_vtx);
	resPosition.setZero(3,num_vtx);
	vtxMap.setZero(1,num_vtx);

	//do  global ICP
	SICP::Parameters pa(false,2.0f,10,1.2f,1e5,100,100,1,0.00001f); //release
	//SICP::Parameters pa(false,2.0,10,1.2,1e5,20,20,1,1e-5);     //debug  
	Eigen::Matrix3Xf srgMatrix,tgMatrix ;
	srgMatrix = srGraph.vertices_matrix().cast<float>();
	tgMatrix = tgGraph.vertices_matrix().cast<float>();
	SICP::point_to_point(srgMatrix ,tgMatrix,vtxMap,pa);
	srGraph.update();//为了查找点的邻域，需要更新

	//local ICP and mix coordinate with rest  and correspondence position
	ScalarType rigidRate = 0.95;
	IndexType iterMax = 1;

	while (iterMax-- > 0)
	{
		//localICP(srGraph,tgGraph,vtxMap,resPosition,isDeformable);
		localICPDisto(srGraph,tgGraph,vtxMap,resPosition,isDeformable);//用对应点的距离比值计算扭曲量
		optimizeCoor(srGraph,tgGraph,vtxMap);
		combinationCoor(rigidRate,srGraph,tgGraph,resPosition,vtxMap,isDeformable,distLrgIndex);
	}
}


bool DeformableRegistration::sampleValidTraj(vector<PCloudTraj> & totTraj,vector<PCloudTraj>& resTraj,LifeSpan & lifesp,IndexType num)
{
	assert(num > 3);

	for (IndexType i = 0; i < num; i ++)
	{
		IndexType indexTraj = rand()%totTraj.size();
		resTraj.push_back(totTraj[indexTraj]);
	}

	return getCommonTime(resTraj,lifesp);

}
PCloudModel DeformableRegistration::constrcutRigidMotionMode(vector<PCloudTraj> & oriTraj,LifeSpan & comTime)
{
	//oriTraj size is 3,and  begin and end  of those trajectories even known
	assert(comTime.end > comTime.start);

	IndexType totFrame = comTime.end - comTime.start + 1;
	PCloudModel pModel;
	pModel.modeLifeSpan = comTime;
	pModel.Centers.setZero(3,totFrame);

	Vector3Type mean;
	Matrix33 f_rotmate;
	Matrix33 b_rotmate;
	Matrix3X srPoints;
	Matrix3X tgPoints;    

	//calculate rotmate matrix  between frame (_start and _start + 1)
	// rotmate * srPos = tgPos  ----- [3×3] * [3×3] = [3×3];column-major storage vertex position
	getMultiTrajNodeCoordinate(oriTraj,comTime.start,srPoints);

	IndexType centerId = 0;
	calculateCenter(srPoints,mean);
	pModel.Centers.col(centerId) = mean;

	for (IndexType frameIndex = comTime.start; frameIndex < comTime.end; frameIndex ++)//nodes  number = (end - start)
	{
		getMultiTrajNodeCoordinate(oriTraj,frameIndex + 1,tgPoints);	
		// calculate center 2014-10-24
		centerId ++;
		mean.setZero();
		calculateCenter(tgPoints,mean);
		pModel.Centers.col(centerId) = mean;
		//
		f_rotmate.setIdentity();
		b_rotmate.setIdentity();

		//affine-9-30
		//calculateTransformation2PClouds(srPoints,tgPoints,f_rotmate);//calculate  front rotmate --least squares
		//calculateTransformation2PClouds(tgPoints,srPoints,b_rotmate);//calculate  back  rotmate --least squares

		//rigid-10-8--svd-factorization
		point2point(srPoints,tgPoints,f_rotmate);
		point2point(tgPoints,srPoints,b_rotmate);

		pModel.fNode.push_back(f_rotmate);
		pModel.bNode.push_back(b_rotmate);

		srPoints = tgPoints;
	}
	return pModel; 
}

void DeformableRegistration::calculateAffineTrans(Matrix3X& srCoor, Matrix3X& tgCoor,Matrix34& affineMat)
{
	MatrixXX srOriCoor = srCoor;

	MatrixXX tgOriCoor = tgCoor;

	IndexType pSzie = (IndexType)srCoor.cols();
	MatrixXX rA;

	rA.setOnes(4,pSzie);

	rA.block(0,0,3,pSzie) = srOriCoor.block(0,0,3,pSzie);

	MatrixXX rAT = rA.transpose();


	//MatrixXX res = rAT.ldlt().solve(tgOriCoor.transpose() );
	MatrixXX tmp = tgOriCoor.transpose().cast<ScalarType>(); 
	MatrixXX res = rAT.fullPivHouseholderQr().solve(tmp);

	affineMat = res.transpose();
}
////
void DeformableRegistration::getMultiTrajNodeCoordinate(vector<PCloudTraj>& multiTraj,IndexType frameId,Matrix3X &resCoorMat)
{//获取不同轨迹在同一帧的顶点坐标<假设它们有共同的life-time>

	//assert(each trajectory have the common instance of frameId)

	IndexType trajNum = (IndexType)multiTraj.size();
	resCoorMat.setZero(3,trajNum);
	PointType tempPoint;
	for (IndexType trajIndex = 0; trajIndex < trajNum; trajIndex ++)
	{
		getTrajNodeCoordinate(multiTraj[trajIndex],frameId,tempPoint);
		resCoorMat.col(trajIndex) = tempPoint;
	}
}
//
void DeformableRegistration::calculateNormalRotMat(NormalType& sN,PointType& tN,Matrix33& rotMat)
{
	ScalarType sL,sA,tL,tA;
	//ScalarType PI = 3.141592653;

	sL = acos(sN.z());
	tL = acos(tN.z());

	//assert(sN.x() != 0 && tN.x() != 0);
	if(sN.x() > -1e-5 && sN.x() < 1e-5 ) 
	{
		if (sN.y() > 1e-5)
		{
			sA = 0.5 *PI;
		}
		sA = 1.5*PI;
	}
	else
	{
		sA = atan(sN.y()/sN.x());
	}


	if(tN.x() > -1e-5 && tN.x() < 1e-5 ) 
	{
		if (tN.y() > 1e-5)
		{
			tA = 0.5 *PI;
		}
		tA = 1.5*PI;
	}
	else
	{
		tA = atan(tN.y()/tN.x());
	}


	tA -= sA;
	tL -= sL;
	//rotate with A axis
	Matrix33 rotZ;

	rotZ(0,0) = cos(2*PI + tA);rotZ(0,1) =-sin(2*PI + tA);rotZ(0,2) =0; 
	rotZ(1,0) = sin(2*PI + tA);rotZ(1,1) =cos(2*PI + tA) ;rotZ(1,2) =0; 
	rotZ(2,0) = 0;             rotZ(2,1) =0 ;            rotZ(2,2) =1; 

	//Logger<<"Z axis = "<<rotZ<<endl;
	Matrix33 rotX;

	rotX(0,0) = 1;rotX(0,1) = 0;rotX(0,2) =0; 
	rotX(1,0) = 0;rotX(1,1) =cos(2*PI + tL) ;rotX(1,2) =-sin(2*PI +tL); 
	rotX(2,0) = 0;rotX(2,1) =sin(2*PI +tL) ;rotX(2,2) =cos(2*PI + tL); 

	//Logger<<"X axis = "<<rotX<<endl;
	rotMat = rotZ * rotX;
}

bool DeformableRegistration::getCommonTime(vector<PCloudTraj> & eleTraj,LifeSpan & lifesp)
{
	IndexType lifeStart = (IndexType)0;
	IndexType lifeEnd   = (IndexType)1e5;
	IndexType trajSize = (IndexType)eleTraj.size();
	PCloudTraj temp;
	for (IndexType trajId = 0; trajId < trajSize; trajId ++)
	{
		temp = eleTraj[trajId];
		if (temp.trajLifeSpan.start > lifeStart)
		{
			lifeStart = temp.trajLifeSpan.start;
		}

		if (temp.trajLifeSpan.end < lifeEnd)
		{
			lifeEnd = temp.trajLifeSpan.end;
		}
	}

	if ( lifeEnd > lifeStart)
	{
		lifesp.start = lifeStart;
		lifesp.end   = lifeEnd;
		return true;
	}else
	{
		return false;
	}
}
//
void DeformableRegistration::combinationCoor(ScalarType stiffCoff,Sample & srGraph,Sample & trGraph,Matrix3X & restPos,MatrixXXi& vtx_cor)
{
	IndexType	num_vtx	= (IndexType)srGraph.num_vertices();
	PointType temp ;
	temp.setZero();
	for (IndexType v_it = 0; v_it < num_vtx; v_it++)
	{
		temp = stiffCoff * restPos.col(v_it) + (1- stiffCoff) * trGraph.vertices_matrix().col(vtx_cor(0,v_it));
		srGraph.vertices_matrix().col(v_it) = temp;
	}
}
//
void DeformableRegistration::optimizeCoor(Sample&srGraph,Sample&tgGraph,MatrixXXi & res_cor)
{
	IndexType	num_vtx	= (IndexType)srGraph.num_vertices();
	Matrix3X    distance(3,num_vtx);
	ScalarType  total_dis = 0.0;

	for (IndexType v_idx = 0; v_idx < num_vtx;v_idx ++)
	{
		unsigned int mp = res_cor(0,v_idx);
		distance.col(v_idx) = tgGraph.vertices_matrix().col(mp) - srGraph.vertices_matrix().col(v_idx);
	}

	//calculate total distance
	total_dis = totalDistance(srGraph,m_neigNum,distance);

	//update correspondence

	IndexType *tg_neig = new IndexType[m_neigNum];
	PointType srPoint;
	IndexType vtx_cor = 0;
	Vec3 srMean;
	srPoint.setZero();
	srMean.setZero();
	IndexType newCor = 0;
	ScalarType criterion;
	IndexType iterMax = 2;
	Matrix3X temp;
	temp.setZero(3,num_vtx);
	while (iterMax-- > 0)
	{
		criterion = total_dis;
		for (IndexType v_idx = 0; v_idx < num_vtx; v_idx++)
		{
			srPoint = srGraph.vertices_matrix().col(v_idx);
			vtx_cor = res_cor(0,v_idx);
			tgGraph.neighbours(vtx_cor,m_neigNum,tg_neig);             //find neighbor of vtx_cor in tgGraph
			avgDistancePoint(v_idx,srGraph,m_neigNum,distance,srMean); 
			res_cor(0,v_idx) = findMinCorr(srPoint,srMean,tg_neig,tgGraph);//update correspondence		
			temp.col(v_idx) = tgGraph.vertices_matrix().col(newCor) - srPoint;
		}
		distance = temp;
		total_dis = totalDistance(srGraph,m_neigNum,distance);
		if (total_dis > criterion) break;	//???	
	}
}
//
void DeformableRegistration::localICPDisto(Sample & srGraph,Sample& tgGraph,MatrixXXi & vtx_cor, Matrix3X& restPos,MatrixXX& isDeformable)
{
	IndexType num_vtx = (IndexType)srGraph.num_vertices();
	IndexType * neigSrGraph = new IndexType[m_neigNum];
	IndexType * neigTgGraph = new IndexType[m_neigNum];
	Matrix3X neigSrCoor(3,m_neigNum);
	Matrix3X neigTgCoor(3,m_neigNum);
	Matrix3X OrineigSrCoor;

	VecX neigDis;     //record ori distance
	VecX resNeigDis;  // record target distances 
	neigDis.setZero(m_neigNum,1);
	resNeigDis.setZero(m_neigNum,1);

	//loacl ICP 
	SICP::Parameters pa(false,2,10,1.2f,1e5,20,20,1,1e-5f);
	MatrixXXi localCor;
	localCor.setZero(1,m_neigNum);


	for (IndexType v_idx = 0; v_idx < num_vtx; v_idx++)
	{
		neigSrCoor.setZero();
		neigTgCoor.setZero();

		srGraph.neighbours(v_idx,m_neigNum,neigSrGraph);
		tgGraph.neighbours(vtx_cor(0,v_idx),m_neigNum,neigTgGraph);// in order to local ICP

		getConstantNeigCoorByIndex(srGraph,neigSrGraph,neigSrCoor);
		getConstantNeigCoorByIndex(tgGraph,neigTgGraph,neigTgCoor);// in order to local ICP
		calculateNeigDis(neigSrCoor,neigDis);//calculate the ori distance around the v_idx 为了计算扭曲量用到的距离

		OrineigSrCoor = neigSrCoor;
		Eigen::Matrix3Xf neigSrMatrix ,neigTgMatrix;
		neigSrMatrix = neigSrCoor.cast<float>();
		neigTgMatrix = neigTgCoor.cast<float>();
		SICP::point_w_point(neigSrMatrix,neigTgMatrix,localCor,pa);

		restPos.col(v_idx) = neigSrCoor.col(0);//record the res position(rigid ICP)
		vtx_cor(0,v_idx) = neigTgGraph[localCor(0,0)];//update the correspondence of v_idx in tgGraph



		//计算中心点到邻域边
		alignTargetCoor(neigTgCoor,localCor);
		//calculateNeigDis(neigTgCoor,resNeigDis);
		//isDeformable(0,v_idx) = deformableValue(neigDis,resNeigDis);//

		isDeformable(0,(__int64)v_idx) = deformableTotal(OrineigSrCoor,neigTgCoor);
	}

	delete []neigSrGraph;
	delete []neigTgGraph;	
}
//
void DeformableRegistration::combinationCoor(ScalarType stiffCoff,Sample & srGraph,Sample & trGraph,
	Matrix3X & restPos,MatrixXXi& vtx_cor,MatrixXX& isDeformable,std::set<IndexType> & distLrgIndex)
{
	ScalarType threshold = 0.2;
	IndexType	num_vtx	= (IndexType)srGraph.num_vertices();
	//IndexType recordDistortionNum05 = 0;
	PointType temp ;
	temp.setZero();
	// 	 	std::fstream op("E:\\distor.txt",std::ios::app);
	// 		if (op.is_open())
	// 		{
	// 		}
	// 		else
	// 		{
	// 			cout<<"open file error！"<<endl;
	// 			return ;
	// 		}
	for (IndexType v_it = 0; v_it < num_vtx; v_it++)
	{
		if (isDeformable(0,(__int64)v_it) < threshold)
		{
			temp = stiffCoff * restPos.col(v_it) + (1- stiffCoff) * trGraph.vertices_matrix().col(vtx_cor(0,v_it));
			srGraph.vertices_matrix().col(v_it) = temp;
		}else
		{
			//recordDistortionNum05++;
			// op<<v_it<<endl;
			srGraph.vertices_matrix().col(v_it) = restPos.col(v_it);
			distLrgIndex.insert(v_it);
		}

	}
	// 	op<<"this iter is  over  wait the next time !!!!!!.\n";
	// 	op.close();
	//	Logger<<"out of distortion number = "<<recordDistortionNum05<<endl;
}
//
void DeformableRegistration::point2point(Matrix3X & srCloud,Matrix3X & tgCloud,Matrix33 & rotMat)
{
	Vector3Type X_mean, Y_mean;
	for(int i=0; i<3; ++i) //计算两点云的均值
	{
		X_mean(i) = (ScalarType)srCloud.row(i).sum()/srCloud.cols();
		Y_mean(i) = (ScalarType)tgCloud.row(i).sum()/tgCloud.cols();
	}

	srCloud.colwise() -= X_mean;
	tgCloud.colwise() -= Y_mean;

	 //Compute transformation
	Eigen::Affine3f transformation;
	Eigen::Matrix3f sigma = (srCloud * tgCloud.transpose()).cast<float>();
	Eigen::JacobiSVD<Eigen::Matrix3f> svd(sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);
	if(svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0)//contains reflection
	{
		Eigen::Vector3f S = Eigen::Vector3f::Ones(); S(2) = -1.0;
		transformation.linear().noalias() = svd.matrixV()*S.asDiagonal()*svd.matrixU().transpose();
	} else 
	{
		transformation.linear().noalias() = svd.matrixV()*svd.matrixU().transpose();//计算旋转矩阵
	}
	//transformation.translation().noalias() = Y_mean - transformation.linear()*X_mean;//计算平移向量
	rotMat = transformation.linear().cast<ScalarType>() ;

	srCloud.colwise() += X_mean;
	tgCloud.colwise() += Y_mean;
}
//
void DeformableRegistration::calculateCenter(Matrix3X& oriCoordinate,Vector3Type& mean)
{
	assert(oriCoordinate.cols() > 0);
	for(int i=0; i<3; ++i) //calculate center
	{
		mean(i) = oriCoordinate.row(i).sum()/oriCoordinate.cols();
	}
}
void DeformableRegistration::getTrajNodeCoordinate(PCloudTraj & traj,IndexType frameIndex,PointType & resCoor)
{
	IndexType vtxId = 0;
	assert(traj.trajLifeSpan.end > traj.trajLifeSpan.start);
	if (frameIndex <= traj.trajLifeSpan.end && frameIndex >= traj.trajLifeSpan.start)
	{
		vtxId = traj.trajNode[frameIndex - traj.trajLifeSpan.start];
		resCoor = m_smpSet[frameIndex].vertices_matrix().col(vtxId);
	}else
	{
		Logger<<" index out of range.\n";
	}
}
//
IndexType DeformableRegistration::findMinCorr(PointType & srPoint,Vec3& mean, int * neigIndex,Sample & smp)
{
	PointType tgPoint;
	IndexType minIndex = 0;
	ScalarType diff = 1e7;
	srPoint += mean;
	for (IndexType neig_it = 0; neig_it < m_neigNum;neig_it++)
	{
		tgPoint = smp.vertices_matrix().col(neigIndex[neig_it]);
		tgPoint -= srPoint;
		if (tgPoint.squaredNorm() < diff)
		{
			diff = tgPoint.squaredNorm();
			minIndex = neigIndex[neig_it];
		}
	}

	return minIndex;
}
//
ScalarType DeformableRegistration::avgDistancePoint(IndexType pointIndex,Sample & smp,IndexType neigNum,Matrix3X & dst,Vec3 & mean)
{
	assert(neigNum > 0);
	ScalarType dis = 0.0;
	IndexType *sr_neig = new IndexType[neigNum];
	Vec3 ori_point;
	ori_point.setZero();
	mean.setZero();

	ori_point = smp.vertices_matrix().col(pointIndex);
	smp.neighbours(pointIndex,neigNum,sr_neig);

	for (IndexType neigIndex = 0; neigIndex < neigNum; neigIndex++)
	{
		mean += dst.col(sr_neig[neigIndex]);
	}
	mean /= neigNum;

	return (ori_point - mean).squaredNorm();
}
//
ScalarType DeformableRegistration::totalDistance(Sample &smp,IndexType neigNum,Matrix3X & dst)
{
	ScalarType tot_dis = 0.0;
	IndexType vtx_num = (IndexType)dst.cols();
	Vec3 mean_;
	for (IndexType v_idx = 0; v_idx < vtx_num; v_idx++)
	{
		tot_dis += avgDistancePoint(v_idx,smp,neigNum,dst,mean_);
	}
	return tot_dis;
}
//
ScalarType DeformableRegistration::deformableTotal(Matrix3X& srCoor, Matrix3X& tgCoor)
{
	ScalarType distor = 0.0;
	IndexType pSize = (IndexType)srCoor.cols();
	IndexType itTime = 0;
	for (IndexType i = 0; i < pSize - 1; i++ )
	{
		for (IndexType j = i + 1; j < pSize; j ++)
		{
			ScalarType srDis = ( srCoor.col(i) - srCoor.col(j) ).norm();
			ScalarType tgDis = ( tgCoor.col(i) - tgCoor.col(j) ).norm();

			ScalarType denominator = srDis;
			ScalarType molecule = abs(srDis - tgDis);
			assert(denominator > 1e-7);
			distor += molecule/denominator;
			itTime ++;
		}
	}

	assert(itTime > 0);

	distor /= itTime;

	return distor;
}
//
void DeformableRegistration::alignTargetCoor(Matrix3X & tgCloud,MatrixXXi & vertexMap)
{
	MatrixXX temp = tgCloud;
	int verN = (int)temp.cols();
	for (int i = 0; i < verN; i++)
	{
		tgCloud.col(i) = temp.col(vertexMap(0,i));
	}
}
//
void DeformableRegistration::calculateNeigDis(Matrix3X & neigCoor,VecX& resDis)
{
	IndexType vtx_num = (IndexType)neigCoor.cols();
	PointType ori = neigCoor.col(0);
	PointType diff;
	diff.setZero();
	for (IndexType vtx_iter = 1; vtx_iter < vtx_num;vtx_iter ++)
	{
		diff = neigCoor.col(vtx_iter) - ori;
		resDis(vtx_iter,0) = diff.norm();
	}
}
//
void DeformableRegistration::getConstantNeigCoorByIndex(Sample & smp,IndexType * neigIndex,Matrix3X & res_coor)
{
	for (IndexType v_it = 0; v_it < m_neigNum;v_it++)
	{
		res_coor.col(v_it) = smp.vertices_matrix().col(neigIndex[v_it]);
	}
}
DeformableRegistration::DeformableRegistration():
	m_smpSet(SampleSet::get_instance())
{
	m_neigNum = 60;
	//Logger<<"Local ICP using neighbor K ="<<m_neigNum<<endl;
}
