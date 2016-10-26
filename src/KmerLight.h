/*
 * KmerLight.h
 *
 *  Part of software kmerlight distributed under GNU GPL 3 license.
 *  Created on: 15-May-2016
 *  Author: Naveen Sivadasan
 *
 */

#ifndef SRC_KMERLIGHT_H_
#define SRC_KMERLIGHT_H_
#include <pthread.h>
#include <time.h>
#include "types.h"
#include "Kmer.h"
#include "constants.h"
#include "FileRead.h"
#include "KmerHash.h"
#include "CountSketch.h"

#define nbufs 2 

struct KMLThrdData{
	byte_t * buf;
	int hashbufidx;
	int start;
	int end;
	unsigned int size;
	int id;
	CountSketch * cs = NULL;
};


class KmerLight{
	static byte_t * seqbufs[nbufs];

	static Rands rands[Constants::copies];


public:


	CountSketch CS;

	static void (*partAdjust) (byte_t * & seqbuf, int & start, int & end, int & idx );
	static void (*bufAdjust) (byte_t * & seqbuf,  int & start,  int end);

	static int (*streamSkip) (byte_t * & seqbuf, int & start, int & end, int & idx);


	static FileRead filep;

	static pthread_t thrds[3];
	static KMLThrdData thrd_data[3];

	static ulong_t memsize;
	static int readBeg;




	KmerLight(int k,  char ** infiles, int nfiles){
		KmerHash::init(k);
		for(int i=0; i< nbufs; ++i) {
			seqbufs[i] = (byte_t *) new byte_t[Constants::seq_buf_size];
			memsize += Constants::seq_buf_size;
		}
		KmerLight::memsize += (1+Constants::buckets)*4*Constants::R*Constants::copies;

		KmerLight::memsize += (16*Constants::seq_buf_size * KmerHash::km_word_len);



		filep.openFiles(infiles, nfiles);

		switch(filep.getType()){
			case FASTA:
				partAdjust = partAdjustFasta;
				bufAdjust = bufAdjustFasta;
				streamSkip = streamSkipFasta;
				break;
			case FASTQ:
				partAdjust = partAdjustFastQ;
				streamSkip = streamSkipFastaQ;
				bufAdjust = bufAdjustFastQ;
				break;
		}
	}
	~KmerLight(){
		for(int i=0; i< nbufs; ++i) {
			free (seqbufs[i]);
		}
	}

	static size_t loadStream(byte_t* seqbuf){
		return filep.fileRead(seqbuf);
	}

	static void *loadStreamWrap(void *inp){
		KMLThrdData * data = (KMLThrdData *)inp;
		unsigned int size;
		size = (unsigned int) loadStream(data->buf);
		data->size = size;
		data->start = 0;
		pthread_exit(NULL);
	}

	static void *storeValsWrap(void *inp){
		KMLThrdData * data = (KMLThrdData *)inp;
		(*data->cs).storeValsInSketchMT(data->hashbufidx);
		pthread_exit(NULL);
	}

	static void *streamHashPartWrap(void *inp){
			KMLThrdData * data = (KMLThrdData *)inp;
			 int start = data -> start;
			int size = data -> end;
			bufAdjust(data->buf, start, size);
			data-> start = start;
			data-> end = size;
			streamHashPart(data->buf,data->start, data->end, data->hashbufidx);
			pthread_exit(NULL);
	}

	void processStreamMT(){ // this should be MT
		int size;
		int sizeprev;
		int bufidx =0;
		int hashbufidx = 0;
		int start=0;
		bool flag = true;
		ulong_t progress = 0;
		const ulong_t prog_delta = (ulong_t) (FileRead::fsize / 20);
		ulong_t nextupd = prog_delta;
		int progcntr = 0;



		size = (unsigned int) loadStream(seqbufs[bufidx]);

		bufAdjust(seqbufs[bufidx], start, size-1);
		streamHashPart(seqbufs[bufidx],start, size - 1, hashbufidx);
		sizeprev = filep.getReadsize();
		start = 0;
		bufidx = (bufidx + 1) & 1;//toggle
		size = (unsigned int) loadStream(seqbufs[bufidx]);
		while(flag){
			if(size == 0) {
				CS.storeValsInSketchMT(hashbufidx);
				flag = false;
				cout << "\r100%\n" << flush;
				break;
			}
			thrd_data[0].id = 0;
			thrd_data[0].hashbufidx = hashbufidx;
			thrd_data[0].cs = &CS;
			pthread_create(&thrds[0], NULL, storeValsWrap, &thrd_data[0]);


			hashbufidx  = (hashbufidx + 1) & 1;//toggle
			KmerHash::resetKMBuf(hashbufidx);
			thrd_data[1].id = 1;
			thrd_data[1].buf = seqbufs[bufidx];
			thrd_data[1].start = start;
			thrd_data[1].end = size-1;
			thrd_data[1].hashbufidx = hashbufidx;
			pthread_create(&thrds[1], NULL, streamHashPartWrap, &thrd_data[1]);

			bufidx = (bufidx + 1) & 1;//toggle
			thrd_data[2].buf = seqbufs[bufidx];
			thrd_data[2].id = 2;

			pthread_create(&thrds[2], NULL, loadStreamWrap, &thrd_data[2]); 

			void * status;
			pthread_join(thrds[0], &status);

			pthread_join(thrds[1], &status);

			pthread_join(thrds[2], &status);
			progress += sizeprev;
			if(progress >= nextupd){
				progcntr += 5;
				cout << "\r" << progcntr << "%" << flush;
				nextupd += prog_delta;
			}


			sizeprev = filep.getReadsize();
			size = thrd_data[2].size;

		}
	}



	void processStream(){ // this should be MT
		 int size;
		 int start=0;


		while((size = (  int) loadStream(seqbufs[0])) > 0){
			start = 0;
			bufAdjust(seqbufs[0], start, size-1);
			streamHashPart(seqbufs[0],start, size - 1, 0);


			CS.storeValsInSketch(0);

			KmerHash::resetKMBuf(0);
		}
	}

	static inline void streamHashPart(byte_t * & seqbuf, int  start, int  end, int hashbufidx){//hash a large part of the stream - MT
		int ret= _valX;
		int idx = start;

		partAdjust(seqbuf, start, end, idx);

		while(ret != _valEOFBUF){
			ret = KmerHash::lShiftSeqHashPart(seqbuf, start, end, idx, hashbufidx );

			if(ret == _valEOFBUF) break;
			ret = streamSkip(seqbuf, start, end, idx);
			if(idx > end) break;
		}
	}

	static INLINE void skipCommentFasta(byte_t * & seqbuf,  int &idx,  int &end){
		int val;
		bool sawNL = (KmerHash::last_segment_status == ON_COMMENT_LINE_NL) ? true : false;
		while(idx <= end){
			val = Kmer::bpval(seqbuf[idx]);
			if((val == _valNL) || (val == _valCR)){
				sawNL = true;
				KmerHash::last_segment_status = ON_COMMENT_LINE_NL;
			}
			if(sawNL && ((val != _valNL) && (val != _valCR))) {
				KmerHash::last_segment_status = SEG_VALID;
				break;
			}
			++idx;
		}
	}

	static INLINE void bufAdjustFasta (byte_t * & seqbuf,  int & start,  int  end){

		if(KmerHash::last_segment_status == SEG_VALID) return;
		skipCommentFasta(seqbuf, start, end);

	}

	static inline void partAdjustFasta (byte_t * & seqbuf, int & start, int & end, int & idx ){

	}

	static inline void bufAdjustFastQ (byte_t * & seqbuf,  int & start,  int  end ){
		if((KmerHash::last_segment_status == SEG_VALID) || (KmerHash::last_segment_status == SEEN_VALID_NL)) {
			readBeg -= Constants::seq_buf_size; 
			if( KmerHash::last_segment_status == SEG_VALID) return;
		}
		skipFastQLines(seqbuf, start, end);
	}

	static inline void partAdjustFastQ (byte_t * & seqbuf, int & start, int & end, int & idx ){

	}

	static INLINE int streamSkipFasta (byte_t * & seqbuf, int & start, int & end, int & idx){ 
		int val = Kmer::bpval(seqbuf[idx]);

		if((val == _valNL)|| (val == _valCR)) {
			++idx;
			if(idx <= end){
				val = Kmer::bpval(seqbuf[idx]);

				if((idx <= end) && ((val == _valNL)|| (val == _valCR))) {
					++idx;
				}
			}

			return 0;
		}
		if((val == _valN) || (val== _valX)){//reset on symbols such as N, U, M, S etc..

			KmerHash::restartKMbufNewSeg();
			++idx;
		}else if(val == _valGT){
			KmerHash::restartKMbufNewSeg();
			KmerHash::last_segment_status = ON_COMMENT_LINE;
			++idx;
			skipCommentFasta(seqbuf, idx, end);
		}

		return 1;
	}

	static inline void skipFastQLines(byte_t * & seqbuf, int & idx, int & end){
		int val;

		if((KmerHash::last_segment_status == SEEN_VALID_NL) && ((idx + idx - readBeg + 2) <= end )){ 
			idx += (idx - readBeg + 2);
			KmerHash::last_segment_status = ON_COMMENT_LINE;
		}else{
			if(idx > end) return;

			int val = Kmer::bpval(seqbuf[idx]);


			if(KmerHash::last_segment_status == SEEN_VALID_NL){
				KmerHash::last_segment_status = FQ_SEEN_PLUS;
				++idx;
				if(idx > end) return;
			}


			if(KmerHash::last_segment_status == FQ_SEEN_PLUS){

				KmerHash::last_segment_status = FQ_SEEN_PLUS_NL;
				++idx;
				if(idx > end) return;
			}


			if(KmerHash::last_segment_status == FQ_SEEN_PLUS_NL){

				KmerHash::last_segment_status = FQ_IN_QUAL;
			}


			if(KmerHash::last_segment_status == FQ_IN_QUAL){
				while(idx <= end){
					val = Kmer::bpval(seqbuf[idx]);
					if(val == _valNL){
						KmerHash::last_segment_status = FQ_SEEN_QUAL_NL;
						++idx;
						break;
					}
					++idx;
				}
				if(idx > end) return;
			}


			if(KmerHash::last_segment_status == FQ_SEEN_QUAL_NL){//get to next valid
				KmerHash::last_segment_status = ON_COMMENT_LINE;
			}



		}

		if(KmerHash::last_segment_status == ON_COMMENT_LINE){//get to next valid
			while(idx <= end){
				val = Kmer::bpval(seqbuf[idx]);
				if(val == _valNL){
					KmerHash::last_segment_status = SEG_VALID;
					++idx;
					readBeg = idx;
					break;
				}
				++idx;
			}
			if(idx > end) return;
		}

	}


	static inline int streamSkipFastaQ (byte_t * & seqbuf, int & start, int & end, int & idx){ 
		int val = Kmer::bpval(seqbuf[idx]);
		if((val == _valNL)) {
			KmerHash::last_segment_status = SEEN_VALID_NL;
			KmerHash::restartKMbufNewSeg();
			++idx;
			skipFastQLines(seqbuf, idx, end);
			return 0;
		}


		KmerHash::restartKMbufNewSeg();
		++idx;

		return 1;
	}

	ulong_t * processAndEstimate(int num_freq){

		do{
			cout << "Processing " << filep.getfname() << " ...\n";
			if(Constants::streamProcessMT){
				processStreamMT();
			}else{
				processStream();
			}

			KmerHash::last_segment_status = ON_COMMENT_LINE;
			KmerHash::resetKMBuf(0);
			KmerHash::resetKMBuf(1);
			KmerHash::restartKMbufNewSeg();

		}while(filep.hasMoreFile());

		CS.setnfreq(num_freq);



		if(Constants::streamProcessMT){
			CS.analyzeSketchMT();
		}else{
			CS.analyzeSketch();
		}

		return CS.computeAllF(num_freq);

	}
};


byte_t * KmerLight::seqbufs[nbufs];
Rands KmerLight::rands[Constants::copies];

void (*KmerLight::partAdjust) (byte_t * & seqbuf, int & start, int & end, int & idx );
void (*KmerLight::bufAdjust) (byte_t * & seqbuf,  int & start,  int  end);

int (*KmerLight::streamSkip) (byte_t * & seqbuf, int & start, int & end, int & idx);

pthread_t KmerLight::thrds[3];
KMLThrdData KmerLight::thrd_data[3];
ulong_t KmerLight::memsize=0;
int KmerLight::readBeg = -1;


FileRead KmerLight::filep;



#endif /* SRC_KMERLIGHT_H_ */
