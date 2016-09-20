/*
 * KmerHash.h
 *
 *  Part of software kmerlight distributed under GNU GPL 3 license.
 *  Created on: 15-May-2016
 *  Author: Naveen Sivadasan
 *
 */

#ifndef SRC_KMERHASH_H_
#define SRC_KMERHASH_H_

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <bitset>
#include <string.h>
#include "types.h"
#include "Rands.h"
#include "constants.h"
#include "Kmer.h"

#define SEG_VALID			3
#define ON_COMMENT_LINE		1
#define ON_COMMENT_LINE_NL	4
#define FQ_SEEN_PLUS		5
#define FQ_IN_QUAL			6
#define FQ_SEEN_PLUS_NL		7
#define FQ_SEEN_QUAL_NL		8
#define SEEN_VALID_NL		9




struct KmerHash{

		static ulong_t * km_hash_vals_arr[2]; 
		static int km_hash_idx_arr[2];
		static ulong_t km_hash_live_count;

		static ulong_t * km_buf_arr; 
		static ulong_t * km_rev_buf_arr; 

		static int last_segment_status;

		static unsigned int km_len; //kmer length
		static unsigned int km_word_len;

		static ulong_t lmask;
		static ulong_t rev_mask;
		static int rev_shift;

		static int (*lShiftSeqHashPart) (byte_t*  seqbuf, int & start, int & end, int & idx, int hashbufidx);

		static void (*restartKMbufNewSeg)();

		static unsigned __int128 x_1_128;



		static void init(unsigned int k){
			km_len = k;
			km_word_len = (km_len >> 5) + ((km_len & 31) > 0 ? 1 : 0);

			last_segment_status = ON_COMMENT_LINE;


			km_hash_vals_arr[0] = new ulong_t[Constants::seq_buf_size * km_word_len];
			km_hash_vals_arr[1] = new ulong_t[Constants::seq_buf_size * km_word_len];


			km_buf_arr = new ulong_t[2*km_word_len](); //initialize to zero
			km_rev_buf_arr = km_buf_arr + km_word_len;
			km_hash_idx_arr[0] = km_hash_idx_arr[1]  =  -km_word_len;
			km_hash_live_count = 0;

			computeLmask();

			if(km_len <= 32){
				lShiftSeqHashPart = lShiftSeqHashPartSingle;
				restartKMbufNewSeg = restartKMbufNewSegSingle;

			}else if(km_len <= 64){
				lShiftSeqHashPart = lShiftSeqHashPartTwo;
				restartKMbufNewSeg = restartKMbufNewSegTwo;
			}else{
				lShiftSeqHashPart = lShiftSeqHashPartLong;
				restartKMbufNewSeg = restartKMbufNewSegLong;
			}

		}


		~KmerHash(){
			free(km_hash_vals_arr[0]);
			free(km_hash_vals_arr[1]);
		}

		static void computeLmask(){
			unsigned int k1 = km_len & 31; 
			rev_shift = (km_len<<1) -2;
			lmask = (k1 == 0)? (0xffffffffffffffffUL) : (ulong_t) ( (((ulong_t)1) << (2*k1)) -1 );
		}




		INLINE static int lShiftSeqHashPartSingle(byte_t*  seqbuf, int & start, int & end, int & idx, int hashbufidx){
			register ulong_t c=0;
			ulong_t * km_hash_vals = km_hash_vals_arr[hashbufidx];
			register ulong_t & km_buf = km_buf_arr[0];
			register ulong_t & km_rev_buf = km_rev_buf_arr[0];
			int & km_hash_idx =  km_hash_idx_arr[hashbufidx];


			while((idx <= end)  && ((c = Kmer::bpval(seqbuf[idx])) <= 3)){

				km_buf = ((km_buf << 2) & lmask) | c;
				km_rev_buf = (km_rev_buf >> 2)  | ((~c & 3) << rev_shift);

				if(km_hash_live_count >= (km_len -1)) {

					km_hash_vals[++km_hash_idx] = (km_buf < km_rev_buf) ? km_buf : km_rev_buf;
				}
				++km_hash_live_count;
				++idx;
			}

			return ((idx > end)? _valEOFBUF: c);
		}


//####################  TWO  #####################################


		INLINE static int lShiftSeqHashPartTwo(byte_t*  seqbuf, int & start, int & end, int & idx, int hashbufidx){
			register ulong_t c=0;
			ulong_t * km_hash_vals = km_hash_vals_arr[hashbufidx];
			register ulong_t & km_buf0 = km_buf_arr[0];
			register ulong_t & km_buf1 = km_buf_arr[1];

			register ulong_t & km_rev_buf0 = km_rev_buf_arr[0];
			register ulong_t & km_rev_buf1 = km_rev_buf_arr[1];

			int & km_hash_idx =  km_hash_idx_arr[hashbufidx];
			int bufidx = 0;

			while((idx <= end)  && ((c = Kmer::bpval(seqbuf[idx])) <= 3)){


				km_buf1 = ((km_buf1 << 2) & lmask) | (km_buf0 >> 62); 
				km_buf0 = (km_buf0 << 2) | c;

				km_rev_buf0 = (km_rev_buf0 >> 2) | (km_rev_buf1 << 62);
				km_rev_buf1 = (km_rev_buf1 >> 2)  | ((~c & 3) << rev_shift);

				if(km_hash_live_count >= (km_len -1)) {
					km_hash_idx += 2;

					bufidx = (km_buf1 != km_rev_buf1) ? ((km_buf1 < km_rev_buf1) ? 0: 2) : ((km_buf0 < km_rev_buf0)? 0: 2);
					km_hash_vals[km_hash_idx] = km_buf_arr[bufidx];
					km_hash_vals[km_hash_idx+1] = km_buf_arr[bufidx+1];
				}
				++km_hash_live_count;
				++idx;
			}

			return ((idx > end)? _valEOFBUF: c); //go back for skip
		}





		//####################  LONG  #####################################


				INLINE static int lShiftSeqHashPartLong(byte_t*  seqbuf, int & start, int & end, int & idx, int hashbufidx){ //for km_word_len > 2....
					register ulong_t c=0;
					//ulong_t tmp;
					ulong_t * km_hash_vals = km_hash_vals_arr[hashbufidx];

					const  int nwords = km_word_len;
					int & km_hash_idx =  km_hash_idx_arr[hashbufidx];
					int bufidx = 0;

					while((idx <= end)  && ((c = Kmer::bpval(seqbuf[idx])) <= 3)){

						km_buf_arr[nwords-1] = ((km_buf_arr[nwords-1] << 2) & lmask) | (km_buf_arr[nwords-2] >> 62);
						for( int q=nwords-2; q > 0; --q)  { km_buf_arr[q] = (km_buf_arr[q] << 2)  | (km_buf_arr[q-1] >> 62);} 
						km_buf_arr[0] = (km_buf_arr[0] << 2) | c;

						for( int q=0; q < nwords-1; ++q) {km_rev_buf_arr[q] = (km_rev_buf_arr[q] >> 2) | (km_rev_buf_arr[q+1] << 62);}
						km_rev_buf_arr[nwords-1] = (km_rev_buf_arr[nwords-1] >> 2)  | ((~c & 3) << rev_shift);


						if(km_hash_live_count >= (km_len -1)) {
							km_hash_idx += nwords;
							bufidx = 0;
							for( int q=nwords-1; q >=0; --q){ //compute the lesser one...
								if(km_buf_arr[q] != km_rev_buf_arr[q]){
									if(km_rev_buf_arr[q] < km_buf_arr[q]) {
										bufidx = nwords;
									}
									break;
								}
							}
							for( int q=0; q<nwords; ++q){km_hash_vals[km_hash_idx+q] = km_buf_arr[bufidx+q];}
						}
						++km_hash_live_count;
						++idx;
					}

					return ((idx > end)? _valEOFBUF: c); //go back for skip
				}







		INLINE static void resetKMBuf(int hashbufidx){ 


			if(km_hash_live_count > km_len - 1)  km_hash_live_count = km_len - 1;
			km_hash_idx_arr[hashbufidx] = -KmerHash::km_word_len;

		}




		INLINE static void restartKMbufNewSegSingle(){ 

			if(km_hash_live_count > 0){
				km_buf_arr[0] = 0;
				km_rev_buf_arr[0] = 0;
				km_hash_live_count = 0;
			}
		}

		INLINE static void restartKMbufNewSegTwo(){ 

			if(km_hash_live_count > 0){
				km_buf_arr[0] = km_buf_arr[1] = 0;
				km_rev_buf_arr[0] = km_rev_buf_arr[1] = 0;
				km_hash_live_count = 0;
			}
		}

		INLINE static void restartKMbufNewSegLong(){ 

			if(km_hash_live_count > 0){
				memset(km_buf_arr, 0, km_word_len << 2);
				km_hash_live_count = 0;
			}
		}


};

ulong_t * KmerHash::km_buf_arr; 
ulong_t * KmerHash::km_rev_buf_arr;
ulong_t * KmerHash::km_hash_vals_arr[2]; 
int KmerHash::km_hash_idx_arr[2];
ulong_t KmerHash::km_hash_live_count;

int KmerHash::last_segment_status;

unsigned int KmerHash::km_len; 
unsigned int KmerHash::km_word_len; 

ulong_t KmerHash::lmask;
ulong_t KmerHash::rev_mask;
int KmerHash::rev_shift;


int (*KmerHash::lShiftSeqHashPart) (byte_t*  seqbuf, int & start, int & end, int & idx, int hashbufidx);
void (*KmerHash::restartKMbufNewSeg)();




#endif /* SRC_KMERHASH_H_ */
