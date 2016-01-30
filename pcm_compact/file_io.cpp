#include "file_io.h"
#include "basic_types.h"
#include "sample.h"
#include "sample_set.h"

#include <stdio.h>
#include<io.h>
#include<vector>
#include<string>

namespace FileIO
{
	Sample* load_point_cloud_file( std::string filename, FILE_TYPE type, IndexType sample_idx )
	{
		//read .xyz or .ply files
		FILE* in_file = fopen(filename.c_str(), "r");

		if(in_file==NULL)
			return nullptr;

		Sample* new_sample = new Sample;
		if( type == FileIO::XYZ)
		{
			while ( true )
			{

				ScalarType vx,vy,vz,nx,ny,nz,cx,cy,cz;
#ifdef USE_FLOAT64

				int status = fscanf(in_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&vx,&vy,&vz,
#else
				int status = fscanf(in_file, "%f %f %f %f %f %f %f %f %f\n",&vx,&vy,&vz,
#endif // 
					&nx,&ny,&nz,
					&cx,&cy,&cz);
				if(status==EOF)break;
				PointType v(vx,vy,vz);
				ColorType cv(cx/255.,cy/255.,cz/255.,1.);
				NormalType nv(nx,ny,nz);

				new_sample->add_vertex(v, nv, cv);


			}

		}else if( type == FileIO::PLY){
			int vcount=0;
			int fcount=0;
			fscanf(in_file,"ply\nformat ascii 1.0\nelement vertex %d\n",&vcount);
			fscanf(in_file,"property float x\nproperty float y\nproperty float z\nproperty float nx\nproperty float ny\nproperty float nz\nproperty uchar red\nproperty uchar green\nproperty uchar blue\n");
			fscanf(in_file,"end_header\n");
			while(true){

				ScalarType vx,vy,vz,nx,ny,nz,cx,cy,cz;
#ifdef USE_FLOAT64

				int status = fscanf(in_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&vx,&vy,&vz,
#else
				int status = fscanf(in_file, "%f %f %f %f %f %f %f %f %f\n",&vx,&vy,&vz,
#endif // 
					&nx,&ny,&nz,
					&cx,&cy,&cz);
				if(status==EOF)break;

				PointType v(vx,vy,vz);
				ColorType cv(cx/255.,cy/255.,cz/255.,1.);
				NormalType nv(nx,ny,nz);

				new_sample->add_vertex(v, nv, cv);
			}

		}else{
			return NULL;
		}
		fclose(in_file);
		//give the new sample a color
		new_sample->set_visble(false);
		new_sample->set_color( Color_Utility::span_color_from_table( sample_idx ) );
		new_sample->build_kdtree();
		//add by huayun
		Box box = new_sample->getBox();
		ScalarType dia_distance = box.diag();
		return new_sample;


	}
	void load_point_cloud_files(std::string filesDir ,FILE_TYPE type)
	{
		std::vector<std::string> filenames;
		struct _finddata_t fa;
		intptr_t fHandle;
		char dirname[256];
		strcpy(dirname ,filesDir.c_str());
		strcat( dirname ,"*.ply");
		if( (fHandle=_findfirst( dirname, &fa ))==-1L )//这里可以改成需要的目录 
		{
			printf( "there are no ply file\n");
		}
		else{
			do
			{
				printf( "find:%s\n", fa.name );
				char cdirname[256];
				strcpy(cdirname ,filesDir.c_str());
				strcat(cdirname, fa.name);
				filenames.push_back( cdirname);
			}while( _findnext(fHandle,&fa)==0 );
			_findclose( fHandle );
		}

		IndexType fileNum = (IndexType)filenames.size();
		SampleSet::get_instance().clear();
		for(IndexType fileid = 0; fileid< fileNum ;++fileid)
		{
			Sample* cur_sample = load_point_cloud_file( filenames[fileid] , FileIO::PLY ,fileid);
			if(cur_sample) SampleSet::get_instance().push_back(cur_sample);
		}
	}
}