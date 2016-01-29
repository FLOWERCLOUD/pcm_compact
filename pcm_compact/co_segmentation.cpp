#include "co_segmentation.h"
#define  INF_LOCAL 1000000
PoolAllocator CoSegmentation::allocator_;

void CoSegmentation::components2HierComponets()
{
	map<IndexType, map<IndexType,IndexType> > v_Label_of_vtx;
	map<IndexType, vector<HLabel*> > v_label_bucket;

	IndexType cluster_count=0;
	for ( map<IndexType,set<IndexType> >::iterator ii=mint_.union_find_set.begin();	ii!=mint_.union_find_set.end(); ii++,cluster_count++)
	{
		set<IndexType>& ss = ii->second;
		for ( set<IndexType>::iterator jj=ss.begin(); jj!=ss.end(); jj++ )
		{
			IndexType frame = components_[*jj].frame;

			IndexType old_label = components_[*jj].label;

			HLabel* new_label_space = allocator_.allocate<HLabel>();

			HLabel* new_label = new (new_label_space)HLabel;

			new_label->label_id = cluster_count;

			new_label->frame_parent = &((*hier_componets)[frame]);

			if (!components_[*jj].vtx_corr_next_frame.empty())
			{
				for ( auto viter=components_[*jj].vtx_corr_next_frame.begin(); viter!=components_[*jj].vtx_corr_next_frame.end();viter++)
				{
					IndexType vtx_idx = viter->first;
					if (v_Label_of_vtx.find(frame) == v_Label_of_vtx.end() )
					{
						v_Label_of_vtx.insert( make_pair(frame,map<IndexType,IndexType>() ) );

						v_Label_of_vtx[frame].insert(make_pair(vtx_idx,cluster_count) );
					}else
					{
						v_Label_of_vtx[frame].insert(make_pair(vtx_idx,cluster_count) );
					}

					IndexType hg = (IndexType)(*hier_componets)[frame].hier_label_bucket.size();
					--hg;

					IndexType idMap = (*hier_componets)[frame].hier_label_vtxBucket_index[hg][old_label];

					HVertex* getVtx = (*hier_componets)[frame].hier_label_bucket[hg][idMap]->vertex_bucket[vtx_idx];

					getVtx->label_parent.push_back(new_label);

					new_label->vertex_bucket.insert(make_pair(vtx_idx, getVtx) );
				}
			}else
			{
				for ( auto viter=components_[*jj].vtx_corr_prev_frame.begin(); viter!=components_[*jj].vtx_corr_prev_frame.end();viter++)
				{
					IndexType vtx_idx = viter->first;
					if (v_Label_of_vtx.find(frame) == v_Label_of_vtx.end() )
					{
						v_Label_of_vtx.insert( make_pair(frame,map<IndexType,IndexType>() ) );

						v_Label_of_vtx[frame].insert(make_pair(vtx_idx,cluster_count) );
					}else
					{
						v_Label_of_vtx[frame].insert(make_pair(vtx_idx,cluster_count) );
					}

					IndexType hg = (IndexType)(*hier_componets)[frame].hier_label_bucket.size();
					--hg;

					HVertex* getVtx = (*hier_componets)[frame].hier_label_bucket[hg][old_label]->vertex_bucket[vtx_idx];

					getVtx->label_parent.push_back(new_label);

					new_label->vertex_bucket.insert(make_pair(vtx_idx, getVtx) );
				}
			}



			if (v_label_bucket.find(frame) == v_label_bucket.end() ) 
			{
				v_label_bucket.insert( make_pair(frame, vector<HLabel*>() ));

				v_label_bucket[frame].push_back(new_label);

			}else
			{
				v_label_bucket[frame].push_back(new_label);
			}

		}
	}


	for (auto fIter = v_label_bucket.begin(); fIter != v_label_bucket.end(); ++ fIter)
	{
		IndexType frame = fIter->first;

		vector<HLabel*>& label_buctet = fIter->second;

		map<IndexType,IndexType> bucket_index;

		IndexType idx = 0;

		for (auto vIter = label_buctet.begin(); vIter != label_buctet.end(); ++ vIter,++ idx)
		{
			IndexType lab_id = (*vIter)->label_id;

			bucket_index[lab_id] = idx;
		}

		(*hier_componets)[frame].hier_label_bucket.push_back(fIter->second);

		(*hier_componets)[frame].label_of_vtx = v_Label_of_vtx[frame];

		(*hier_componets)[frame].hier_label_vtxBucket_index.push_back(bucket_index);
	}

}

void CoSegmentation::compute()
{
	calc_similarity();
	cluster_component();
	co_segment();
}
void CoSegmentation::calc_similarity()
{
	for(IndexType i=0;  i<components_.size()-1; ++i)
	{
		// 			if (components_[i].frame == 8 && components_[i].label == 3)
		// 			{
		// 				Logger<<" Start to Debug!.\n";
		// 			}

		for (IndexType j=i+1;j<components_.size();++j)
		{
			if(components_[j].frame != components_[i].frame+1)
				continue;

			// 				Logger<<"Frame = "<<components_[j].frame<<endl;
			// 				Logger<<"Label = "<<components_[j].label<<endl;

			//ScalarType d = similarity_between_component(i,j);

			ScalarType d = similarity_between_component_back(i,j);

			if(d!=INF_LOCAL)
			{
				cq_.push( compo_dist_node(i,j,d) );
			}
		}

	}
}
void CoSegmentation::cluster_component()
{
	mint_.init( (IndexType)components_.size() );
	while (!cq_.empty())
	{
		compo_dist_node nod = cq_.top();
		cq_.pop();

		/*check whether to unite*/
		IndexType x = mint_.find(nod.i);
		IndexType y = mint_.find(nod.j);
		set<IndexType>& sx  =mint_.union_find_set[x];
		set<IndexType>& sy = mint_.union_find_set[y];
		bool same_frame = false;
		for ( set<IndexType>::iterator iter1 = sx.begin(); iter1!=sx.end(); iter1++ )
		{
			for (set<IndexType>::iterator iter2=sy.begin(); iter2!=sy.end(); iter2++)
			{
				Component& c1 = components_[*iter1];
				Component& c2 = components_[*iter2];
				if(c1.frame == c2.frame)
				{
					same_frame = true;
					break;
				}
			}
			if(same_frame)break;
		}

		if(!same_frame)
			mint_.unite(x, y);
	}
}
void CoSegmentation::co_segment()
{
	IndexType cluster_count=0;
	for ( map<IndexType,set<IndexType> >::iterator ii=mint_.union_find_set.begin();
		ii!=mint_.union_find_set.end(); ii++)
	{
		set<IndexType>& ss = ii->second;
		Logger<<"cluster "<<cluster_count++<<endl;
		for ( set<IndexType>::iterator jj=ss.begin(); jj!=ss.end(); jj++ )
		{
			IndexType frame = components_[*jj].frame;
			Logger<<frame<<"  "<<components_[*jj].label<<endl;
			for ( auto viter=components_[*jj].vtx_corr_next_frame.begin();
				viter!=components_[*jj].vtx_corr_next_frame.end();viter++)
			{
				IndexType vtx_idx = viter->first;
				SampleSet::get_instance()[frame][vtx_idx].set_label( cluster_count );
			}
		}
	}

	//build_mini_seg();

}

ScalarType CoSegmentation::similarity_between_component_back(IndexType i, IndexType j)
{
	assert( components_[i].frame==components_[j].frame-1 );
	ScalarType toDis = 0.;
	ScalarType backDis = 0.;
	map<IndexType,IndexType> com_part;
	map<IndexType,IndexType> com_back_part;

	IndexType com_num = compute_common(components_[i], components_[j],com_part);

	IndexType com_back_num = compute_common_back(components_[i], components_[j], com_back_part);

	if (com_num==0 && com_back_num == 0) // should bigger than ..vertexes 
	{
		return INF_LOCAL;
	}

	const IndexType k = 300;
	IndexType neighbours[k];
	ScalarType dist[k];
	map<IndexType, IndexType> valid_neigs;

	//to distances

	if ( com_num > 0)
	{
		for ( map<IndexType,IndexType>::iterator iter = com_part.begin();
			iter != com_part.end(); iter++)
		{
			valid_neigs.clear();
			Sample& cur_frame = SampleSet::get_instance()[components_[i].frame];
			cur_frame.neighbours( iter->first, k, neighbours, dist );
			for (IndexType kk=0; kk<k; kk++)
			{
				auto fiter = com_part.find( neighbours[kk] );
				if( fiter!=com_part.end() )
				{
					assert(fiter->second!=-1);
					valid_neigs.insert( *fiter );
				}
			}

			IndexType neig_siz = (IndexType)valid_neigs.size();
			if(neig_siz<3)continue; 
			Matrix3X s_coord, t_coord;
			s_coord.setZero(3, neig_siz);
			t_coord.setZero(3, neig_siz);

			Sample& next_frame = SampleSet::get_instance()[ components_[j].frame ];
			IndexType vi = 0;
			IndexType kk = 0;
			for ( map<IndexType,IndexType>::iterator neig_iter = valid_neigs.begin();
				neig_iter != valid_neigs.end(); ++neig_iter )
			{
				s_coord.col(kk) = cur_frame.vertices_matrix().col(neig_iter->first );
				t_coord.col(kk) = next_frame.vertices_matrix().col(neig_iter->second);
			}
			Matrix33 rot_mat;
			MatrixXX tran_vec;

			point2point(s_coord, t_coord, rot_mat, tran_vec);

			auto bias = rot_mat*cur_frame.vertices_matrix().col(iter->first) + tran_vec 
				- next_frame.vertices_matrix().col(iter->second);

			toDis += bias.norm(); //local ICP 

		}
	}


	

	return  toDis + backDis;

}
//
void CoSegmentation::point2point(Matrix3X & srCloud,Matrix3X & tgCloud,Matrix33 & rotMat,MatrixXX & transVec)
{
	Eigen::Vector3d X_mean, Y_mean;

	for(int i=0; i<3; ++i) //计算两点云的均值
	{
		X_mean(i) = srCloud.row(i).sum()/srCloud.cols();
		Y_mean(i) = tgCloud.row(i).sum()/tgCloud.cols();
	}

	srCloud.colwise() -= X_mean;
	tgCloud.colwise() -= Y_mean;

	/// Compute transformation
	Eigen::Affine3d transformation;
	Eigen::Matrix3d sigma = srCloud * tgCloud.transpose();
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);
	if(svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0)//contains reflection
	{
		Eigen::Vector3d S = Eigen::Vector3d::Ones(); S(2) = -1.0;
		transformation.linear().noalias() = svd.matrixV()*S.asDiagonal()*svd.matrixU().transpose();
	} else 
	{
		transformation.linear().noalias() = svd.matrixV()*svd.matrixU().transpose();//计算旋转矩阵
	}

	transVec = Y_mean - transformation.linear()*X_mean;
	rotMat = transformation.linear() ;

	srCloud.colwise() += X_mean;
	tgCloud.colwise() += Y_mean;

}
//
IndexType CoSegmentation::compute_common_back(Component & c1, Component & c2, map<IndexType,IndexType> & common_part )
{
	IndexType count = 0;
	for ( map<IndexType,IndexType>::iterator iter = c2.vtx_corr_prev_frame.begin();
		iter != c2.vtx_corr_prev_frame.end(); iter++)
	{
		if( c1.vtx_corr_next_frame.find( iter->second )!=c1.vtx_corr_next_frame.end() )
		{
			count++;
			common_part.insert( *iter );
		}
	}
	return count;
}
//
IndexType CoSegmentation::compute_common( Component & c1, Component & c2, map<IndexType,IndexType> & common_part )
{
	IndexType count = 0;
	for ( map<IndexType,IndexType>::iterator iter = c1.vtx_corr_next_frame.begin();
		iter != c1.vtx_corr_next_frame.end(); iter++)
	{
		if( c2.vtx_corr_next_frame.find( iter->second )!=c2.vtx_corr_next_frame.end() )
		{
			count++;
			common_part.insert( *iter );
		}
	}
	return count;
}
//
void CoSegmentation::hierComponets2Components()
{
	//把hier_component最上层的数据赋值给Components

	assert(!hier_componets->empty() );

	if (!components_.empty() )
	{
		components_.clear();
	}

	for (auto fIter = hier_componets->begin(); fIter != hier_componets->end(); ++ fIter)
	{
		IndexType gLevel = (IndexType)fIter->second.hier_label_bucket.size();//访问最高层的label_bucket

		IndexType fId = fIter->first;

		vector<HLabel*>& label_buctet = fIter->second.hier_label_bucket[gLevel - 1];

		IndexType compo_idx;

		for (auto lIter = label_buctet.begin(); lIter != label_buctet.end(); ++ lIter)
		{

			IndexType pId = (*lIter)->label_id;

			map<IndexType,HVertex*>& vtx_bucket = (*lIter)->vertex_bucket;

			compo_idx = (IndexType)components_.size();

			components_.push_back( Component(fId, pId) );

			for (auto vIt = vtx_bucket.begin(); vIt != vtx_bucket.end(); ++vIt)
			{
				IndexType vtx_id = vIt->first;

				HVertex* vtx = vIt->second;

				if (vtx->prev_corr != NULL)
				{
					IndexType prevId = vtx->prev_corr->vtx_id;

					components_[compo_idx].vtx_corr_prev_frame.insert(make_pair(vtx_id,prevId) );
				}else
				{
					components_[compo_idx].vtx_corr_prev_frame.insert(make_pair(vtx_id, -1) );
				}

				if (vtx->next_corr != NULL)
				{
					IndexType nextId = vtx->next_corr->vtx_id;

					components_[compo_idx].vtx_corr_next_frame.insert(make_pair(vtx_id, nextId) );

				}else
				{
					components_[compo_idx].vtx_corr_next_frame.insert(make_pair(vtx_id, -1) );
				}

			}
		}
	}

}