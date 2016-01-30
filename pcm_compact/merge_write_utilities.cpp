#include "merge_write_utilities.h"




void mergeFile(char* CORR_FILESNAME ,char* CORR_FILEOUT_NAME , char* LABEL_FILESNAME , char* LABEL_FILEOUT_NAME)
{
	std::vector<std::string> labelfilenames,corrfilenames;
	struct _finddata_t fa;
	intptr_t fHandle;
	char labeldirname[256] ,corrdirname[256];
	strcpy(labeldirname ,LABEL_FILESNAME);
	strcpy(corrdirname ,CORR_FILESNAME);
	strcat( labeldirname ,"*");
	strcat(corrdirname ,"*");
	if( (fHandle=_findfirst( labeldirname, &fa ))==-1L )//这里可以改成需要的目录 
	{
		printf( "there are no ply file\n");
	}
	else{
		do
		{
			printf( "find:%s\n", fa.name );
			labelfilenames.push_back( fa.name);
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
			printf( "find:%s\n", fa.name );
			corrfilenames.push_back( fa.name);
		}while( _findnext(fHandle,&fa)==0 );
		_findclose( fHandle );
	}

	char* corrout_ = CORR_FILEOUT_NAME;
	FILE* corrTMPfile = fopen( corrout_ ,"w");
	fclose( corrTMPfile);
	for (IndexType file_idx = 0; file_idx < corrfilenames.size(); file_idx++)
	{		
		appendcorr(  corrfilenames[file_idx].c_str() , corrout_);
	}
	std::cout<<"merge corr_files success"<<std::endl;

	char* labelout_ = LABEL_FILEOUT_NAME;
	FILE* labelTMPfile = fopen( labelout_ ,"w");
	fclose( labelTMPfile);
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
