#include "merge_write_utilities.h"
#include "sample_set.h"
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include<io.h>
#include <stdlib.h>
#include <windows.h>
#include <tchar.h>
using namespace std;
struct mvertex
{
	ScalarType x,y,z,nx,ny,nz;
	IndexType r,g,b;
};


//************************************
// Method:    mergeFile
// Parameter: char * CORR_FILESNAME     input files path of directory of corr files 
// Parameter: char * CORR_FILEOUT_NAME   output path path of merged corr files
// Parameter: char * LABEL_FILESNAME      input files path of directory of label files
// Parameter: char * LABEL_FILEOUT_NAME   output path of merged label files
//************************************
void mergeFile(char* CORR_FILESNAME ,char* CORR_FILEOUT_NAME , char* LABEL_FILESNAME , char* LABEL_FILEOUT_NAME)
{
	std::vector<std::string> labelfilenames,corrfilenames;
	struct _finddata_t fa;
	intptr_t fHandle;
	char labeldirname[256] ,corrdirname[256];
	strcpy(labeldirname ,LABEL_FILESNAME);
	strcpy(corrdirname ,CORR_FILESNAME);
	strcat( labeldirname ,"*.txt");  //in order to eliminate paths: "." , ".."
	strcat(corrdirname ,"*.txt");
	if(_access(LABEL_FILESNAME,0) == -1){
		printf( "input label files path not exits\n");
		return;
	}
	if(_access(CORR_FILESNAME,0) == -1){
		printf("input corr files path not exits\n");
		return;
	}
	if( (fHandle=_findfirst( labeldirname, &fa ))==-1L )//这里可以改成需要的目录 
	{
		printf( "there are no ply file\n");
	}
	else{
		do
		{
			char cdirname[256];
			strcpy(cdirname , LABEL_FILESNAME);
			strcat(cdirname, fa.name);
			printf( "find:%s\n",cdirname );
			labelfilenames.push_back(cdirname);
		}while( _findnext(fHandle,&fa)==0 );
		_findclose( fHandle );
	}
	if( (fHandle=_findfirst( corrdirname, &fa ))==-1L )//这里可以改成需要的目录 
	{
		printf( "there are no ply file\n");
	}
	else{
		do
		{
			char cdirname[256];
			strcpy(cdirname , CORR_FILESNAME);
			strcat(cdirname, fa.name);
			printf( "find:%s\n",cdirname );
			corrfilenames.push_back( cdirname);
		}while( _findnext(fHandle,&fa)==0 );
		_findclose( fHandle );
	}

	char* corrout_ = CORR_FILEOUT_NAME;
	FILE* corrTMPfile = fopen( corrout_ ,"w");
	if(!corrTMPfile)
	{ 
		if(_access(".\\tmp",0) == -1){
			CreateDirectory(_T(".\\tmp"), NULL);
		}
		if(_access(".\\tmp\\results",0) == -1){
			CreateDirectory(_T(".\\tmp\\results"), NULL);
		}	
		corrout_=".\\tmp\\results\\\hksingle_corr.txt";
		corrTMPfile = fopen( corrout_ ,"w");
	}
	if(corrTMPfile)fclose( corrTMPfile);
	for (IndexType file_idx = 0; file_idx < corrfilenames.size(); file_idx++)
	{		
		appendcorr(  corrfilenames[file_idx].c_str() , corrout_);
	}
	std::cout<<"merge corr_files success"<<std::endl;

	char* labelout_ = LABEL_FILEOUT_NAME;
	FILE* labelTMPfile = fopen( labelout_ ,"w");

	if(!labelTMPfile)
	{
		if(_access(".\\tmp",0) == -1){
			CreateDirectory(_T(".\\tmp"), NULL);
		}
		if(_access(".\\tmp\\results",0) == -1){
			CreateDirectory(_T(".\\tmp\\results"), NULL);
		}	
		labelout_=".\\tmp\\results\\cosegOOrder.txt";
		labelTMPfile = fopen( labelout_ ,"w");
	}
	if(labelTMPfile)fclose( labelTMPfile);
	for (IndexType file_idx = 0; file_idx < labelfilenames.size(); file_idx++)
	{		
		appendlabel( labelfilenames[file_idx].c_str() , labelout_);
	}
	std::cout<<"merge label_files success"<<std::endl;
}

void mergeFile()
{
	char* CORR_FILESNAME = "H:\\povay_pointcloud\\point_data\\shiqun\\Propagation_Info\\Propagation_Info\\corrs";
	char* CORR_FILEOUT_NAME= "H:\\povay_pointcloud\\point_data\\shiqun\\Propagation_Info\\Propagation_Info\\corroutput0_29center5.txt";
	char* LABEL_FILESNAME = "H:\\povay_pointcloud\\point_data\\shiqun\\Propagation_Info\\Propagation_Info\\labels";
	char* LABEL_FILEOUT_NAME = "H:\\povay_pointcloud\\point_data\\shiqun\\Propagation_Info\\Propagation_Info\\labeloutput0_29center5.txt";
	mergeFile(CORR_FILESNAME ,CORR_FILEOUT_NAME , LABEL_FILESNAME , LABEL_FILEOUT_NAME);
	//QString dir_in = QFileDialog::getExistingDirectory( &(main_window::getInstance()) ,QString("Import point cloud files"),".");
}

void appendcorr(const char* file_in ,char* file_out)
{
	FILE* in_file = fopen( file_in,"r");
	FILE* outfile = fopen( file_out,"a");
	//	ofstream fout(file_out);
	if(!in_file){
		printf("input file not exist\n");
		return;
	}
	if(!outfile){
		printf("output file not exist\n");
		return;
	}
	std::string str;
	int buf[4];
	while(true){

		int stat = fscanf(in_file,"%d %d %d %d\n",  &buf[0] ,&buf[1] , &buf[2] ,&buf[3]);	 
		if(stat==EOF)break;
		fprintf(outfile, "%d %d %d %d\n", buf[0] ,buf[1] , buf[2] ,buf[3]);
	}

	fclose(outfile);	
	fclose(in_file);
}

void appendlabel(const char* file_in ,char* file_out)
{
	FILE* in_file = fopen( file_in,"r");
	FILE* outfile = fopen( file_out,"a");
	//	ofstream fout(file_out);
	if(!in_file){
		printf("input file not exist\n");
		return;
	}
	if(!outfile){
		printf("output file not exist\n");
		return;
	}
	std::string str;

	int buf[3];
	//fscanf(in_file,"%d\n",  &buf[0] );	 //略过第一行
	while(true){

		int stat = fscanf(in_file,"%d %d %d\n",  &buf[0] ,&buf[1] , &buf[2]);	 
		if(stat==EOF)break;
		fprintf(outfile, "%d %d %d\n", buf[0] ,buf[1] , buf[2] );
	}

	fclose(outfile);	
	fclose(in_file);
}

ColorType getLabelColor(void* pvtx)
{

	ColorType pClr = Color_Utility::span_color_from_hy_table( ((Vertex*)pvtx)->label() );
	return pClr;
}



void generateSamplePly( char* _label_filename ,char* _output_file_path,char* _prefix )
{
	char* file_sample = _label_filename;

	FILE* in_smpfile = fopen( file_sample ,"r");
	if(!in_smpfile){
		printf( "input label files path not exits\n");
		return ;
	}
	IndexType frameId, labelId, vtxId;

	IndexType cframe = -1;
	IndexType endindex = 0;

	Sample* smp;
	std::ofstream* coutStream =NULL;
	std::map<IndexType , vector<mvertex> > verticeMap;
	vector<mvertex> vtxvec;
	while(true){



		int stat = fscanf( in_smpfile,"%d %d %d\n",  &frameId , &labelId ,&vtxId);
		if(stat==EOF){ //fprintf(outfile ,"#end\n");
			if( NULL!=coutStream)coutStream->close();
			if( cframe != -1){
				verticeMap.insert(  make_pair( cframe ,vtxvec ));
				vtxvec.clear();
			}
			break;
		}

		if( cframe !=frameId){


			smp = &(SampleSet::get_instance()[frameId]);
			if( cframe != -1){
				verticeMap.insert(  make_pair( cframe ,vtxvec ));
				vtxvec.clear();
			}
			cframe = frameId;			
		}

		Vertex& cvtx = SampleSet::get_instance()[frameId][vtxId];

		cvtx.set_label( labelId);
		ColorType pClr = getLabelColor(&cvtx );
		mvertex mvtx;
		mvtx.x = cvtx.x();
		mvtx.y = cvtx.y();
		mvtx.z = cvtx.z();
		mvtx.nx = cvtx.nx();
		mvtx.ny = cvtx.ny();
		mvtx.nz = cvtx.nz();
		mvtx.r = pClr(0 , 0);
		mvtx.g = pClr(1 ,0);
		mvtx.b = pClr(2 ,0);
		vtxvec.push_back( mvtx);

	}
	auto framebitr = verticeMap.begin();
	auto frameeitr = verticeMap.end();
	for ( ; framebitr != frameeitr ; ++framebitr){
		cframe  = (*framebitr).first;
		char path[100];
		strcpy(path ,_output_file_path);
		char fullPath[250];
		sprintf( fullPath ,"%s%s%.3d%s",path ,_prefix, cframe ,".ply");     //必须加入.3d ，使得文件排序正常

		if( NULL!=coutStream)coutStream->close();
		coutStream  =  new std::ofstream( fullPath , std::ofstream::out);
		std::ofstream& outfile = *coutStream;

		outfile<<"ply"<<std::endl;
		outfile<<"format ascii 1.0"<<std::endl;
		outfile<<"element vertex "<<  framebitr->second.size()<<std::endl;
		outfile<<"property   float   x"<<std::endl;
		outfile<<"property   float   y "<<std::endl;
		outfile<<"property   float   z "<<std::endl;
		outfile<<"property   float   nx"<<std::endl;
		outfile<<"property   float   ny "<<std::endl;
		outfile<<"property   float   nz "<<std::endl;
		outfile<<"property   uchar red "<<std::endl;
		outfile<<"property   uchar   green"<<std::endl;
		outfile<<"property   uchar  blue"<<std::endl;
		outfile<<"end_header"<<std::endl;


		for (int i=0; i<framebitr->second.size(); i++)
		{				
			mvertex& v =(framebitr->second)[i];
			std::ofstream& outfile = *coutStream;
			outfile<<v.x<<" "<<v.y<<" "<<v.z<<" "<<v.nx<<" "<<v.ny<<" "<<v.nz<<" "<<v.r<<" "<<v.g<<" "<<v.b<<std::endl;
		}
		coutStream->close();

	} 
	fclose(in_smpfile);
}