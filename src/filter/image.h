#include "../filter.h"

// Detect BMP, RGB, TGA, PBM, PGM, PPM, PAM formats
class DetectIMAGE {
  Filter *f;
  int bmp, imgbpp, bmpx, bmpy, bmpof, hdrless;  // For BMP detection
  int rgbi, rgbx, rgby;  // For RGB detection
  int tga, tgax, tgay, tgaz, tgat;  // For TGA detection
  int pgm, pgmcomment, pgmw, pgmh, pgm_ptr, pgmc, pgmn, pamatr, pamd;  // For PBM, PGM, PPM, PAM detection
  char pgm_buf[32];
  U32 buf2, buf1, buf0;
public:
  DetectIMAGE(Filter *filter) : f(filter),
  bmp(0), imgbpp(0), bmpx(0), bmpy(0), bmpof(0), hdrless(0), 
  rgbi(0), rgbx(0), rgby(0), 
  tga(0), tgax(0), tgay(0), tgaz(0), tgat(0), 
  pgm(0), pgmcomment(0), pgmw(0), pgmh(0), pgm_ptr(0), pgmc(0), pgmn(0), pamatr(0), pamd(0) {};
  bool detect(int c, int i, long start, int n, FILE *in);
};

int encode_bmp(FILE* in, FILE* out, int len, int width, int &hdrsize);
int decode_bmp(Encoder& en, FILE *in, int size, int width, FILE *out, FMode mode, int &diffFound);
int encode_im32(FILE* in, FILE* out, int len, int width, int &hdrsize);
int decode_im32(Encoder& en, FILE *in, int size, int width, FILE *out, FMode mode, int &diffFound);
