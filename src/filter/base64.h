#include "../misc.h"
#include "../filter.h"

// Detect base64 data
class DetectBASE64 {
  Filter *f;
  int b64s, b64i, b64line, b64nl;
  U32 buf1, buf0;
public:
  DetectBASE64(Filter *filter) : f(filter), b64s(0), b64i(0), b64line(0), b64nl(0) {};
  bool detect(int c, int i);
};

int encode_base64(FILE* in, FILE* out, int len, int info, int &hdrsize);
int decode_base64(Encoder& en, FILE* in, int size, int info, FILE *out, FMode mode, int &diffFound);