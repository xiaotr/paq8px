#include "../filter.h"

// Detect JPEG format
class DetectJPEG {
  Filter *f;
  int soi, sof, sos, app;
  U32 buf1, buf0;
public:
  DetectJPEG(Filter *filter) : f(filter), soi(0), sof(0), sos(0), app(0) {};
  bool detect(int c, int i, Filetype type);
};
