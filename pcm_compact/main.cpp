#include "basic_types.h"
#include "TrajectoryClassifier.h"
#include "GCop.h"
int main(int argc, char *argv[])
{
	IndexType neigNum, trajLen,resolution,_octreeReso,modelT,_smallL ,smallLife;
	IndexType selectedSmpIdx = 0;
	ScalarType _perC,thresHold,lambda;
	bool isEqual,isRigid;
	neigNum =trajLen =resolution =_octreeReso=modelT=_smallL =smallLife =0;
	_perC=thresHold=lambda=0.0;
	isEqual = isRigid =false;
	//step1
	TrajectoryClassifier* cluster = new TrajectoryClassifier(selectedSmpIdx);
	cluster->setNeigNum(neigNum);
	cluster->setParamter(trajLen,resolution,lambda,thresHold,modelT,smallLife,isEqual,isRigid);
	cluster->run();
	//step2
	GCop * graphCut = new GCop;
	IndexType m_nLabels,m_nExpansion,m_nSwap,m_nGraphNeig,m_nCurNeig,m_centerF;
	m_nLabels=m_nExpansion=m_nSwap=m_nGraphNeig=m_nCurNeig=m_centerF =0;
	ScalarType m_nTradeOff,m_cSigma,m_dSigma,m_nDiffu;
	m_nTradeOff=m_cSigma=m_dSigma=m_nDiffu =0.0;
	graphCut->setParamter(m_nLabels,m_nExpansion,m_nSwap,m_nGraphNeig,m_nCurNeig,m_nTradeOff,m_cSigma,m_dSigma,m_nDiffu,m_centerF);
	graphCut->run();
	//step3

}