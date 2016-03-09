#include "basic_types.h"
#include "file_io.h"
#include "TrajectoryClassifier.h"
#include "GCop.h"
#include "merge_write_utilities.h"
#include <stdlib.h>
#include <windows.h>
#include <tchar.h>
#define step1
#define step2
#define  MERGE
#define generatePly
#define step3
char* filename = "panther";
int main(int argc, char *argv[])
{
#ifdef step1
	//1. read ply files
	FileIO::load_point_cloud_files( "D:\\point_data\\plystandard\\panda20150918\\",FileIO::PLY );
	IndexType sampleNum = SampleSet::get_instance().size();
#endif
#ifdef step2
	//2. set parameter
	IndexType neigNum, trajLen,resolution,modelT ,smallLife;
	ScalarType _perC,thresHold;
	bool isEqual,isRigid;

	neigNum = 11;
	trajLen =2;
	resolution =32;
	_perC =0.55;  // lamda
	thresHold = 0.7;
	modelT =1;
	smallLife =2;
	isEqual = true;
	isRigid = false;	
	IndexType selectedSmpIdx = 2;  //set the central frame
	//step1
	TrajectoryClassifier* cluster = new TrajectoryClassifier(selectedSmpIdx);
	cluster->setNeigNum(neigNum);
	cluster->setParamter(trajLen,resolution,_perC,thresHold,modelT,smallLife,isEqual,isRigid);
	cluster->run();

	for(IndexType frameId = selectedSmpIdx ; frameId >0 ; --frameId)
	{
		cluster->setCenterFrame(frameId-1);
		cluster->setParamter(trajLen,resolution,_perC,thresHold,modelT,smallLife,isEqual,isRigid);
		cluster->run();
	}
	for(IndexType frameId = selectedSmpIdx ; frameId < sampleNum-1 ; ++frameId)
	{
		cluster->setCenterFrame(frameId+1);
		cluster->setParamter(trajLen,resolution,_perC,thresHold,modelT,smallLife,isEqual,isRigid);
		cluster->run();
	}
#endif
#ifdef MERGE
	mergeFile( 
		".\\tmp\\corr\\" ,
		".\\tmp\\results\\hksingle_corr.txt",
		".\\tmp\\label\\",
		".\\tmp\\results\\cosegOOrder.txt");

	//return 0;
#endif
#ifdef generatePly
	if(_access(".\\tmp",0) == -1){
		CreateDirectory(_T(".\\tmp"), NULL);
	}
	if(_access(".\\tmp\\ply",0) == -1){
		CreateDirectory(_T(".\\tmp\\ply"), NULL);
	}
	if(_access(".\\tmp\\ply\\init",0) == -1){
		CreateDirectory(_T(".\\tmp\\ply\\init"), NULL);
	}
	generateSamplePly( ".\\tmp\\results\\cosegOOrder.txt" ,".\\tmp\\ply\\init\\",filename);
//	return 0;
#endif
#ifdef step3
	//step2
	GCop * graphCut = new GCop;
	IndexType m_nLabels,m_nExpansion,m_nSwap,m_nGraphNeig,m_nCurNeig,m_centerF;
	ScalarType m_nTradeOff,m_cSigma,m_dSigma,m_nDiffu;
	m_nLabels = 1; //read label file to determine
	m_nExpansion =1;
	m_nSwap =1 ;
	m_nGraphNeig = 10;
	m_nCurNeig = 10;
	m_nTradeOff = 1;
	m_cSigma = 1.;
    m_dSigma = 1.;
	m_nDiffu=  1.;
	m_centerF = 2;
	graphCut->setParamter(m_nLabels,m_nExpansion,m_nSwap,m_nGraphNeig,m_nCurNeig,m_nTradeOff,m_cSigma,m_dSigma,m_nDiffu,m_centerF);
	graphCut->run();
#endif	

}