/*
 * constants.h
 *
 *  Part of software kmerlight distributed under GNU GPL 3 license.
 *  Created on: 15-May-2016
 *  Author: Naveen Sivadasan
 *
 */

#ifndef SRC_CONSTANTS_H_
#define SRC_CONSTANTS_H_

#include "types.h"

struct Constants{

	static const ulong_t prime;
	static  const int primeexp;

	static const int A;
	static const int T;
	static const int G;
	static const int C;

	static const int Rexp;
	static const int Uexp;
	static const int RUexp;
	static const int R;  //pow of 2 is for faster hashing...
	static const int U;  //pow of 2 is for faster hashing...
	static const int ULSIZE;
	static const unsigned int rshift_RU;


	static const int buckets;
	static const int copies; //should be an odd number.. its median is the final estimate.
	static const int nthreads;

	static const int seqlist_size;
	static const int seqsize;
	static const int initial_buckets; //for lazy start.. <= buckets

	static const bool rev_compl;

	static const int seq_buf_size; // ~1MB; need to experiment,

	static const int seq_seg_size_max;

	static const bool streamProcessMT;



};

const int Constants::primeexp = 31;
const ulong_t Constants::prime = (1L << Constants::primeexp) - 1;//mersanne prime. This is only for kmer string to 64bit long hashing.

const int Constants::A = 1;
const int Constants::T = 2;
const int Constants::G = 3;
const int Constants::C = 4;

const int Constants::Rexp = 19; //18 ~ 500 MB RAM, 19 ~ 1GB RAM
const int Constants::Uexp = 13;
const int Constants::RUexp = Rexp + Uexp;
const int Constants::R = 1<< Rexp;  //pow of 2 is for faster hashing...
const int Constants::U = 1<< Uexp;  //pow of 2 is for faster hashing...
const int Constants::ULSIZE = (8 * sizeof(ulong_t));
const unsigned int Constants::rshift_RU = Constants::ULSIZE - Constants::RUexp;


const int Constants::buckets = 64;
const int Constants::copies = 7; //odd number.. its median is the final estimate.
const int Constants::nthreads = Constants::copies;

const int Constants::seqlist_size = 1000;
const int Constants::seqsize = 2000;
const int Constants::initial_buckets = 5; //for lazy start.. <= buckets

const bool Constants::rev_compl = true;

const int Constants::seq_buf_size = (1<<22); // ~1MB; need to experiment,

const int Constants::seq_seg_size_max = (int)(Constants::seq_buf_size / Constants::nthreads) + 1;

const bool Constants::streamProcessMT = true;


#endif /* SRC_CONSTANTS_H_ */
