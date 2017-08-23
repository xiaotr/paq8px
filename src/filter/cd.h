#include "../filter.h"

// CD sectors detection (mode 1 and mode 2 form 1+2 - 2352 bytes)
class DetectCD {
  Filter *f;
  int cdi, cda, cdm;
  U32 cdf;
  U32 buf1, buf0;
public:
  DetectCD(Filter *filter) : f(filter), cdi(0), cda(0), cdm(0), cdf(0) {};
  bool detect(int c, int i, long start, int n, FILE *in, Filetype type);
};

int encode_cd(FILE* in, FILE* out, int len, int info, int &hdrsize);
int decode_cd(Encoder& en, FILE *in, int size, int info, FILE *out, FMode mode, int &diffFound);
