/*
 * CountSketch.h
 *
 *  Part of software kmerlight distributed under GNU GPL 3 license.
 *  Created on: 15-May-2016
 *  Author: Naveen Sivadasan
 *
 */

#ifndef SRC_COUNTSKETCH_H_
#define SRC_COUNTSKETCH_H_

#include <algorithm>
#include <pthread.h>
#include "CountSketchInstance.h"

struct CSthrdData{
	int id;
	static int hashbufidx;
};

struct CountSketch{

	static CountSketchInstance *sketch_arr;
	static CSthrdData* thrd_data_arr;


	static pthread_t CSthreads[Constants::nthreads];


	CountSketch(){
		sketch_arr = new CountSketchInstance[Constants::nthreads];

	}
	~CountSketch(){
	}

	static void *storeValsInSketchMTWrap(void *myinp){
		CSthrdData * inp = (CSthrdData*) (myinp);
		sketch_arr[inp->id].storeValsInSketch(inp->hashbufidx);
		pthread_exit(NULL);
	}

	void storeValsInSketchMT(int hashbufidx){

		if(thrd_data_arr == NULL){
			thrd_data_arr = new CSthrdData[Constants::nthreads];
			for(int i=0; i< Constants::nthreads; ++i){
				thrd_data_arr[i].id = i;
			}
		}
		CSthrdData::hashbufidx = hashbufidx;
		for(int i=0; i< Constants::nthreads; ++i){
			pthread_create(&CSthreads[i], NULL, storeValsInSketchMTWrap, &thrd_data_arr[i]);
		}

		void * status;
		for(int i=0; i < Constants::nthreads; ++i){
			pthread_join(CSthreads[i], &status);
		}

	}

	void storeValsInSketch(int hashbufidx){


		for(int i=0; i< Constants::copies; ++i){
			sketch_arr[i].storeValsInSketch(hashbufidx);
		}
	}

	static void * analyzeSketchMTWrap(void *in){
			int * idp = (int *)in;
			sketch_arr[(*idp)].analyzeSketch();
			pthread_exit(NULL);
	}


	static void analyzeSketchMT(){
		int idx[Constants::nthreads];
		for(int i=0; i< Constants::nthreads; ++i){
			idx[i] = i;
			pthread_create(&CSthreads[i], NULL, analyzeSketchMTWrap, &idx[i]);
		}

		void * status;
		for(int i=0; i < Constants::nthreads; ++i){
			pthread_join(CSthreads[i], &status);
		}
	}


	static void analyzeSketch(){
		for(int i=0; i < Constants::copies; ++i){
			sketch_arr[i].analyzeSketch();
		}
	}


	static ulong_t computeF0(){
		ulong_t vals[Constants::copies];
		for (int i=0; i< Constants::copies; ++i){
			vals[i] = sketch_arr[i].computeF0();
		}
		std::sort(vals, vals + Constants::copies);
		return vals[Constants::copies >> 1]; //middle value
	}

	static ulong_t computeFJ(int j, ulong_t F0){
		ulong_t vals[Constants::copies];
		for (int i=0; i< Constants::copies; ++i){
			vals[i] = sketch_arr[i].computeFJ(j, F0);

		}
		std::sort(vals, vals + Constants::copies);
		return vals[Constants::copies >> 1]; //middle value
	}

	static ulong_t * computeAllF(int num_freq){
		ulong_t F0 = computeF0();
		ulong_t *vals = new ulong_t[num_freq + 1];
		vals[0] = F0;
		for(int j=1; j<= num_freq; ++j){
			vals[j] = computeFJ(j, F0);
		}
		return vals;
	}

	static void setnfreq(int nfreq){
		CountSketchInstance::num_freq = nfreq;
	}

};

CountSketchInstance *CountSketch::sketch_arr;
CSthrdData* CountSketch::thrd_data_arr;
pthread_t CountSketch::CSthreads[Constants::nthreads];
int CSthrdData::hashbufidx = 0;

#endif /* SRC_COUNTSKETCH_H_ */
