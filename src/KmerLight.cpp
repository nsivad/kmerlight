/*
 * KmerLight.cpp
 *
 *  Part of software kmerlight distributed under GNU GPL 3 license.
 *  Created on: 15-May-2016
 *  Author: Naveen Sivadasan
 *
 */
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include "KmerLight.h"
using namespace std;

int main(int argc, char * argv[]){

	int k = 21;
	char ** infiles=NULL;
	int nfiles=0;
	//char * fname = (char*)"test.fq";
	char * ofile = (char*) "out";
	int num_freq = 8;


	string  arg;
	if(argc <= 1){
		cout << "\nUsage: kmerlight  -k kmerlength -f freqs  -o outputfile -i <.fa/.fa.gz/.fasta/fasta.gz OR .fq/.fq.gz/.fastq/.fastq.gz file(s)>\n";
		return 0;
	}


	int i=1;
	while(i < argc){
		arg.assign(argv[i]);
		if(arg.compare("-k") == 0){
			++i;
			k = atoi(argv[i]);
		}else if(arg.compare("-f") == 0){
			++i;
			num_freq = atoi(argv[i]);
		}else if(arg.compare("-o") == 0){
			++i;
			ofile = argv[i];
		}else if(arg.compare("-i")== 0){
			for(int j=i+1; j<argc; ++j){
				if(argv[j][0] != '-') {
					++nfiles;
				}else{
					break;
				}
			}
			infiles = new char *[nfiles];
			for(int j=0; j<nfiles; ++j) infiles[j] = argv[i+j+1];
		}
		++i;
	}

	cout << "Parms :  k = " << k << "; nfreq = " << num_freq << "; infile(s) = ";
	for(int j=0; j< nfiles; ++j) cout << infiles[j] << ", ";
	cout << "; outfile = " << ofile << "; copies = " << Constants::copies << endl;


	ulong_t * F;
	KmerLight kl(k, infiles, nfiles);

	time_t timer1, timer2;
	time(&timer1);
	struct tm * timeinfo;
	timeinfo = localtime(&timer1);
	cout << asctime(timeinfo) <<  flush;

	F = kl.processAndEstimate(num_freq);

	time(&timer2);
	timeinfo = localtime(&timer2);
	cout << asctime(timeinfo) <<  flush;

	ofstream ofs;
	ofs.open (ofile);
	ofs << "#Parms :  k = " << k << "; nfreq = " << num_freq << "; infile(s) = ";
	for(int j=0; j< nfiles; ++j) ofs << infiles[j] << ", ";
	ofs << "; outfile = " << ofile << "; copies = " << Constants::copies << endl;


	ofs << "#F0 = " << F[0] << endl;
	ofs << "f\tk="<< k<<endl;
	for(int i=1; i<= num_freq; ++i){
		ofs << "f" << i << "\t" << F[i] << endl;
	}
	ofs.close();
	free (F);

	cout << "Time taken ~ " << difftime(timer2, timer1) << " Sec\n";
	cout << "RAM Use : ~" << (int)(KmerLight::memsize/1000000) + 1 << " MB (total); ~" << (int)(((1+Constants::buckets)*4*Constants::R*Constants::copies)/1000000) << " MB (for Sketch)." << endl;
	cout << "Output in \"" << ofile << "\"\n";

	return 0;
}


