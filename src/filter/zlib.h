#include "../filter.h"

// Detect ZLIB stream
class DetectZLIB {
  Filter *f;
  unsigned char zbuf[32], zin[1<<16], zout[1<<16];
  int zbufpos,zzippos;
  int pdfim,pdfimw,pdfimh,pdfimb,pdfimp;
  U32 buf3, buf2, buf1, buf0;
public:
  DetectZLIB(Filter *filter) : f(filter), zbufpos(0), zzippos(-1), pdfim(0), pdfimw(0), pdfimh(0), pdfimb(0), pdfimp(0) {};
  bool detect(int c, int i, long start, int n, FILE *in);
};

int encode_zlib(FILE* in, FILE* out, int len, int info, int &hdrsize);
int decode_zlib(Encoder& en, FILE* in, int size, int info, FILE *out, FMode mode, int &diffFound);
