#ifndef _FILE_IO_H
#define _FILE_IO_H
//stl
#include <string>
#include "basic_types.h"
class Sample;
class SampleSet;
namespace FileIO
{
	enum FILE_TYPE{ XYZ,PLY,OBJ,NONE };

	Sample* load_point_cloud_file( std::string filename, FILE_TYPE type, IndexType );
	void load_point_cloud_files(std::string filesDir ,FILE_TYPE type);

}

#endif