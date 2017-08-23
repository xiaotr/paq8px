#include "../filter.h"

// Detect WAVE, AIFF, S3M formats
class DetectAUDIO {
  Filter *f;
  int wavi, wavsize, wavch, wavbps, wavm, wavtype, wavlen, wavlist;  // WAVE
  int aiff, aiffm, aiffs;  // AIFF
  int s3mi, s3mno, s3mni;  // S3M
  U32 buf1, buf0;
public:
  DetectAUDIO(Filter *filter) : f(filter),
    wavi(0), wavsize(0), wavch(0), wavbps(0), wavm(0), wavtype(0), wavlen(0), wavlist(0),
    aiff(0), aiffm(0), aiffs(0),
    s3mi(0), s3mno(0), s3mni(0) {};
  bool detect(int c, int i, long start, FILE *in);
};
