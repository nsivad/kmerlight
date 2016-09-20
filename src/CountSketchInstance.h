/*
 * CountSketchInstance.h
 *
 *  Part of software kmerlight distributed under GNU GPL 3 license.
 *  Created on: 15-May-2016
 *  Author: Naveen Sivadasan
 *
 */

#ifndef SRC_COUNTSKETCHINSTANCE_H_
#define SRC_COUNTSKETCHINSTANCE_H_

#include <iostream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "types.h"
#include "constants.h"
#include "Rands.h"
#include "KmerHash.h"
#include "murmur3.h"

#define POSMASK		0xffff0000
#define VALUEMASK	0x0000ffff
#define VALUEINVALID 0xffff

struct Finfo{
	ulong_t t;
	int w_opt;
	ulong_t t_opt;
	ulong_t F;
	Finfo(){
		t=0;
		w_opt=-1;
		t_opt = 0;
		F = 0;
	}
};

struct CountSketchInstance{

	static const int R;
	static const int U1;
	static const int buckets;


	uint32 * sketch[Constants::buckets+1];
	Finfo* finfo_arr;
	Rands rand;

	static int num_freq;//default is only F0


	CountSketchInstance(){
		for(int i=0; i<=Constants::buckets; ++i){
			if( i <= Constants::initial_buckets){
				sketch[i] = new uint32[R];
				memset(sketch[i], 0, (R<<2));
			}else{
				sketch[i] = NULL;
			}
		}

		finfo_arr=NULL;

	}

	~CountSketchInstance(){
		for(int i=0; i<= buckets; ++i) {
			if(sketch[i]  != NULL) free(sketch[i]);
		}
		free (finfo_arr);
	}


	/**
	 * @author Martin Läuter (1997)
	 *         Charles E. Leiserson
	 *         Harald Prokop
	 *         Keith H. Randall
	 * "Using de Bruijn Sequences to Index a 1 in a Computer Word"
	 */

	const uint64_t index64[64] = {
	   63,  0, 58,  1, 59, 47, 53,  2,
	   60, 39, 48, 27, 54, 33, 42,  3,
	   61, 51, 37, 40, 49, 18, 28, 20,
	   55, 30, 34, 11, 43, 14, 22,  4,
	   62, 57, 46, 52, 38, 26, 32, 41,
	   50, 36, 17, 19, 29, 10, 13, 21,
	   56, 45, 25, 31, 35, 16,  9, 12,
	   44, 24, 15,  8, 23,  7,  6,  5
	};
	const uint64_t debruijn64 = 0x07EDD5E59A4E28C2ULL;




	void storeValsInSketch(int hashbufidx) {

			ulong_t val;
			unsigned int W, col, pos, position;
			uint32 value;
			ulong_t  out[2];

			const int incr = KmerHash::km_word_len;
			const int upper = KmerHash::km_hash_idx_arr[hashbufidx];
			const int murmsize = KmerHash::km_word_len << 3;
			const unsigned int seed = rand.murmseed;
				for (int i=0; i<= upper; i+=incr) {

					MurmurHash3_x64_128 ( &(KmerHash::km_hash_vals_arr[hashbufidx][i]), murmsize, seed,  &out );
					val = out[0];

					W = index64[((val & -val) * debruijn64) >> 58];




					val = ((ulong_t)(rand.rand_odd_RU * val + rand.rand_small_RU))  >> Constants::rshift_RU;

					col = (unsigned int) (val >> Constants::Uexp);
					pos = (unsigned int) (val & (U1));



					if(sketch[W] != NULL) {
						value = sketch[W][col] & VALUEMASK;
						position = (sketch[W][col] & POSMASK) >> 16;

						if(value != VALUEINVALID){
							if(value == 0){
								sketch[W][col] = (pos << 16) | 1; 
							}else if(position != pos){
								sketch[W][col] = VALUEINVALID;
							}else {
								++sketch[W][col];
							}
						}

					}else{
						sketch[W] = new uint32[R];
						memset(sketch[W], 0, (R<<2));
						sketch[W][col] = (pos << 16) | 1; 
					}
				}
			}


	void analyzeSketch(){ 

		finfo_arr = new Finfo[num_freq+1];

		const int nf = num_freq;

		for (int i=0; i<=buckets; i=i+1) {

			//Scan the array
			if(sketch[i] == NULL){ 
				finfo_arr[0].t = R;
			}else{
				for (int j=0; j<R; j=j+1) {
					int b = sketch[i][j] & VALUEMASK;
					if ((b>= 0) &&(b <= num_freq)) { 
						++finfo_arr[b].t;
					}
				}
			}
			int R2 = R>>1;


			for (int b=0; b<=nf; b=b+1) {
				if(b == 0){ 
					if(labs(finfo_arr[b].t_opt - R2) > labs(finfo_arr[b].t - R2)){//better one..
						finfo_arr[b].t_opt = finfo_arr[b].t;
						finfo_arr[b].w_opt = i+1;
					}
				}else{
					if(finfo_arr[b].t_opt <  finfo_arr[b].t){//better one..
						finfo_arr[b].t_opt = finfo_arr[b].t;
						finfo_arr[b].w_opt = i+1;
					}
				}
				finfo_arr[b].t = 0;
			}
		}
	}

	ulong_t computeF0(){

		if(finfo_arr[0].t_opt == 0) {
			finfo_arr[0].F = 0;
			return finfo_arr[0].F;
		}

		double p0 = (double) ((double) finfo_arr[0].t_opt / (double) R) ;
		double x = (double) R;
		x = 1/x;
		x = 1 - x;

		double num = (double) log(p0);
		double den = (double) log(x);

		finfo_arr[0].F = (ulong_t) ceil(((double) (1<< finfo_arr[0].w_opt)) * num/den) ;
		return finfo_arr[0].F;
	}

	ulong_t computeFJ(int j, long F0){ 

		int w_opt = finfo_arr[j].w_opt;
		long t_opt = finfo_arr[j].t_opt;

		if(w_opt <= 0 || (t_opt == 0)) { //No estimate
			finfo_arr[j].F = 0;
			return finfo_arr[j].F;
		}

		double x = (double) R;
		x = 1 - (1/x);
		double exp = (double) ((double)F0/ (double)(1 << w_opt));
		double p0 = pow(x, exp);
		double pj = (double) ((double) t_opt/ (double) R);
		double val = (double)(1L<< w_opt);
		val = val * ((double) (R -1)) * (pj / p0);

		finfo_arr[j].F = (ulong_t) ceil(val);

		return finfo_arr[j].F;
	}

	ulong_t * computeAllF(long F0){ 

		ulong_t* F = new ulong_t[num_freq];
		F[0] = F0;
		for(int j=0; j<= num_freq; ++j) F[j] = computeFJ(j, F0);
		return F;
	}



};
const int CountSketchInstance::R = Constants::R;
const int CountSketchInstance::U1 = Constants::U-1;
const int CountSketchInstance::buckets = Constants::buckets;
int CountSketchInstance::num_freq = 0;


#endif /* SRC_COUNTSKETCHINSTANCE_H_ */
