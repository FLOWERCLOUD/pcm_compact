#include "basic_types.h"
#include "file_io.h"
#include "TrajectoryClassifier.h"
#include "GCop.h"

int main(int argc, char *argv[])
{
	//1. read ply files
	//2. set parameter
	IndexType neigNum, trajLen,resolution,modelT ,smallLife;
	ScalarType _perC,thresHold;
	bool isEqual,isRigid;

	neigNum = 11;
	trajLen =2;
	resolution =32;
	_perC =0.45;
	thresHold = 0.7;
	modelT =1;
	smallLife =2;
	isEqual = true;
	isRigid = false;	
	IndexType selectedSmpIdx = 6;  //set the central frame
	//step1
	TrajectoryClassifier* cluster = new TrajectoryClassifier(selectedSmpIdx);
	cluster->setNeigNum(neigNum);
	cluster->setParamter(trajLen,resolution,_perC,thresHold,modelT,smallLife,isEqual,isRigid);
	cluster->run();
	return 0;
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
	m_centerF = 10;
	graphCut->setParamter(m_nLabels,m_nExpansion,m_nSwap,m_nGraphNeig,m_nCurNeig,m_nTradeOff,m_cSigma,m_dSigma,m_nDiffu,m_centerF);
	graphCut->run();
	

}