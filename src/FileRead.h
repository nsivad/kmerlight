/*
 * FileRead.h
 *
 *  Part of software kmerlight distributed under GNU GPL 3 license.
 *  Created on: 15-May-2016
 *  Author: Naveen Sivadasan
 *
 */
#ifndef SRC_FILEREAD_H_
#define SRC_FILEREAD_H_
#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include "../libs/zlib.h"
#include "types.h"
#include "constants.h"

#define PLAIN   0
#define GZIP    1
#define BZIP2	3

#define FASTA		0
#define FASTQ		1
#define FUNKNOWN	-1



using namespace std;

struct FileRead{

	FILE * fp;
    gzFile_s * fpgz;
    bool zip;
    size_t offset;
    char ** infiles;
    int nfiles;
    int currfile;

	static int ftype;
	static ulong_t fsize;

	FileRead(){
		fp = NULL;
		fpgz = NULL;
		zip = false;
		offset = 0;
		infiles = NULL;
		nfiles = 0;
		currfile=0;
	}

	int openFiles(char ** infiles, int nfiles){
		this->infiles = infiles;
		this->nfiles = nfiles;
		openFile(infiles[currfile]);
		return 1;
	}

	bool hasMoreFile(){
		if(currfile == (nfiles-1)) return false;
		openFile(infiles[++currfile]);
		return true;
	}

	char * getfname(){
		return infiles[currfile];
	}

	int openFile(const char * fname){ //check different file formats and compression...

		if(fpgz != NULL) gzclose(fpgz);
		if(fp != NULL) fclose(fp);

		std::string fs(fname);
		std::string fq(".fq");
		std::string fq1(".fastq");
		std::string fa(".fa");
		std::string fa1(".fasta");
		std::string faz(".fa.gz");
		std::string faz1(".fasta.gz");
		std::string fqz(".fq.gz");
		std::string fqz1(".fastq.gz");

		zip = false;
		offset = 0;

		if(has_suffix(fs, fa) || has_suffix(fs, fa1)){
			ftype = FASTA;
		}else if (has_suffix(fs, fq) || has_suffix(fs, fq1)){
			ftype = FASTQ;
		}else if (has_suffix(fs,faz) ||has_suffix(fs,faz1) ){
			ftype = FASTA;
			zip = true;
		}else if(has_suffix(fs, fqz) || has_suffix(fs, fqz1)){
			ftype = FASTQ;
			zip = true;
		}else{
			ftype = FUNKNOWN;
		}

		if(ftype == FUNKNOWN){
			cout << "Unknown file type for : " << fname << "\n";
			exit(1);
		}

		struct stat fstat;
		stat(fname,&fstat);
		fsize = fstat.st_size;


		if(!zip){
			fp = fopen(fname, "r");
			if(fp == NULL){
				cout << "Cannot open file " << fname << "\n";
				exit(1);
			}
		}else{
			fpgz = gzopen(fname,"rb");
			if(fpgz == NULL){
				cout << "Cannot open file " << fname << "\n";
				exit(1);
			}
			gzbuffer(fpgz, Constants::seq_buf_size);
		}

		return 1;
	}

	size_t fileRead(byte_t * ptr){
		if(!zip){
			offset = fread(ptr, 1, Constants::seq_buf_size,fp);
			return offset;
		}
		return  gzread(fpgz, ptr, Constants::seq_buf_size);
	}

	INLINE size_t getReadsize(){
		if(!zip) {
			return offset;
		}

		size_t tmp =  offset; //prev one..
		offset = gzoffset(fpgz);
		return (offset - tmp);
	}

	int getType(){
		return ftype;
	}
	static bool has_suffix(const std::string &str, const std::string &suffix)
	{
	    return str.size() >= suffix.size() &&
	           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
	}
};


int FileRead::ftype =  FASTA;
ulong_t FileRead::fsize = 0;


#endif /* SRC_FILEREAD_H_ */
