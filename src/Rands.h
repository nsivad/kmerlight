/*
 * Rands.h
 *
 *  Part of software kmerlight distributed under GNU GPL 3 license.
 *  Created on: 15-May-2016
 *  Author: Naveen Sivadasan
 *
 */

#ifndef SRC_RANDS_H_
#define SRC_RANDS_H_
#include <iostream>
#include <random>
#include <cmath>
#include "types.h"

#define RANDMAX		0xffffffffffffffffUL
struct Rands{

	ulong_t rand_bucket;
	ulong_t rand_bucket2;
	unsigned __int128 a_bucket_rand_128;
	unsigned __int128 b_bucket_rand_128;


	ulong_t rand_odd_RU;
	ulong_t rand_small_RU; 

	static unsigned __int128 a_kmer_rand_128;
	static unsigned __int128 b_kmer_rand_128;
	static bool set_kmer_rands;

	unsigned int murmseed;



	Rands(){

		std::random_device rd;    
		std::mt19937_64 rng(rd());   
		std::uniform_int_distribution<ulong_t> uni(0, RANDMAX); // uniform

		if(!set_kmer_rands){
			a_kmer_rand_128 = uni(rng);
			a_kmer_rand_128 = a_kmer_rand_128 << 64;
			a_kmer_rand_128 |= uni(rng);

			b_kmer_rand_128 = uni(rng);

			b_kmer_rand_128 = b_kmer_rand_128 << 64;
			b_kmer_rand_128 |= uni(rng);
			set_kmer_rands = true;
		}

		a_bucket_rand_128 = uni(rng);
		a_bucket_rand_128 = a_bucket_rand_128 << 64;
		a_bucket_rand_128 |= uni(rng);

		b_bucket_rand_128 = uni(rng);

		b_bucket_rand_128 = b_bucket_rand_128 << 64;
		b_bucket_rand_128 |= uni(rng);




		rand_bucket = uni(rng);

		rand_bucket2 = uni(rng);
		if(rand_bucket2 == 0){
			++rand_bucket2;
		}else if((rand_bucket2 & 1) != 1){
			--rand_bucket2;
		}


		rand_odd_RU = uni(rng);



		if(rand_odd_RU == 0){
			++rand_odd_RU;
		}else if((rand_odd_RU & 1) != 1){
			--rand_odd_RU;
		}

		rand_small_RU = uni(rng);

		rand_small_RU = rand_small_RU >> Constants::RUexp;

		murmseed = (unsigned int) uni(rng);

	}
};

unsigned __int128 Rands::a_kmer_rand_128;
unsigned __int128 Rands::b_kmer_rand_128;
bool Rands::set_kmer_rands = false;
#endif /* SRC_RANDS_H_ */
