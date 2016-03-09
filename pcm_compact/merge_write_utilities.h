#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include<io.h>
#include "basic_types.h"
typedef int IndexType;
void appendcorr( const char* file_in ,char* file_out);
void appendlabel( const char* file_in ,char* file_out);

void mergeFile(char* CORR_FILESNAME ,char* CORR_FILEOUT_NAME , char* LABEL_FILESNAME , char* LABEL_FILEOUT_NAME);
void mergeFile( );
void generateSamplePly( char* _label_filename ,char* _output_file_path_,char* _prefix);
ColorType getLabelColor(void* pvtx);