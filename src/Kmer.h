/*
 * Kmer.h
 *
 *  Part of software kmerlight distributed under GNU GPL 3 license.
 *  Created on: 15-May-2016
 *  Author: Naveen Sivadasan
 *
 */

#ifndef SRC_KMER_H_
#define SRC_KMER_H_
#include "types.h"

#define _valX  15  //others
#define _valGT  14 //>
#define _valAT 13  //@
#define _valPL  12 //+
#define _valNL  11 //\n
#define _valCR  10 //\r
#define _valEOFBUF  9 //End of stream part.. to be used inside.
#define _valEOSEG  8 //End of stream part.. to be used inside.

#define _valA	0
#define _valT	3
#define _valC	1
#define _valG	2
#define _valN   4 //N/n


struct Kmer{

	static ulong_t bpval_tab[16]; 
	Kmer(){
			for(int i =255; i>=0; --i){
				unsigned char t = (unsigned char) i;

				int idx = t>>4; 
				bpval_tab[idx] <<= 4;

				if((t == 'A')||(t == 'a')){
					bpval_tab[idx] |= _valA; //A

				}else if ((t == 't')||(t == 'T')){
					bpval_tab[idx] |= _valT; //T

				}else if((t == 'g')||(t == 'G')){
					bpval_tab[idx] |= _valG; //G

				}else if((t=='c')||(t=='C')){
					bpval_tab[idx] |= _valC; //C

				}else if((t=='n')||(t=='N')){
					bpval_tab[idx] |= _valN; //N

				}else if(t=='>'){
					bpval_tab[idx] |= _valGT; //>

				}else if(t=='@'){
					bpval_tab[idx] |= _valAT; //@

				}else if(t=='+'){
					bpval_tab[idx] |= _valPL; //+

				}else if(t=='\n'){
					bpval_tab[idx] |= _valNL; //\n

				}else if(t=='\r'){
					bpval_tab[idx] |= _valCR; //\r

				}else{
					bpval_tab[idx] |= _valX; //X
				}
			}
		}

	static inline ulong_t bpval(byte_t b){
		return (((ulong_t)(bpval_tab[b >> 4] & ( ((ulong_t)15) << ((b & 15)<<2) ))) >> ((b & 15)<<2));
	}

}kmertab;

ulong_t Kmer::bpval_tab[16];

#endif /* SRC_KMER_H_ */
