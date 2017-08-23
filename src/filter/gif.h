#include "../filter.h"

// Detect GIF format
class DetectGIF {
  Filter *f;
  int gif, gifa, gifi, gifw, gifc, gifb;
  U32 buf1, buf0;
public:
  DetectGIF(Filter *filter) : f(filter), gif(0), gifa(0), gifi(0), gifw(0), gifc(0), gifb(0) {};
  bool detect(int c, int i, Filetype type, Filetype &nextBlockType);
};

int encode_gif(FILE* in, FILE* out, int len, int info, int &hdrsize);
int decode_gif(Encoder &en, FILE* in, int size, int info, FILE *out, FMode mode, int &diffFound);
