#include "graph_cut_node.h"
using namespace std;
void GraphNodeCtr::read_corres_file(char *filename)
{
	FILE *in_file = fopen(filename,"r");
	if (in_file==NULL)
	{
		return;
	}

	while (true)
	{
		int frame,index,cor_frame,cor_index;
		int status = fscanf(in_file,"%d %d %d %d\n",&frame,&index,&cor_frame,&cor_index);
		if(status==EOF)break;
		add_corresponding_relation(frame,index,cor_frame,cor_index);
	}
}

IndexType GraphNodeCtr::readnLabelFile(char *filename)
{
	FILE *in_file = fopen(filename,"r");
	if (in_file==NULL)
	{
		return 0;
	}

	IndexType nLabels = 0;
	fscanf(in_file,"%d\n",&nLabels);

	/// original sample labels
	IndexType num = 0;
	while (true)
	{
		int frame,label,index;
		int status = fscanf(in_file,"%d %d %d\n",&frame,&label,&index);
		if(status==EOF)break;
		add_node(frame, label, index);
		label_bucket[frame_label_to_key(frame,label)].insert(index);

		num++;
	}

	return nLabels;
}

void GraphNodeCtr::add_node(IndexType frame, IndexType label, IndexType index)
{

	GraphCutNode *new_space = allocator_.allocate<GraphCutNode>();
	GraphCutNode *new_node = new(new_space) GraphCutNode(frame,label,index,cur_graph_index_++);
	node_vec.push_back(new_node);
	node_map[frame_index_to_key(frame,index)] = new_node;
}


void GraphNodeCtr::add_corresponding_relation( IndexType frame, IndexType index, IndexType cor_frame, IndexType cor_idx )
{
	node_map[frame_index_to_key(frame,index)]->cor_frame_index.insert(make_pair(cor_frame,cor_idx));
}