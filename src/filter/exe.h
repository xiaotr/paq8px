#include "../misc.h"
#include "../filter.h"

// Detect executable (EXE, DLL) format
class DetectEXE {
  Filter *f;
  Array<int> abspos,  // CALL/JMP abs. addr. low byte -> last offset
             relpos;  // CALL/JMP relative addr. low byte -> last offset
  int e8e9count;  // number of consecutive CALL/JMPs
  int e8e9pos;    // offset of first CALL or JMP instruction
  int e8e9last;   // offset of most recent CALL or JMP
  U32 buf1, buf0;
public:
  DetectEXE(Filter *filter) : f(filter), abspos(256), relpos(256), e8e9count(0), e8e9pos(0), e8e9last(0) {};
  bool detect(int c, int i, int start, Filetype type);
};


int encode_exe(FILE* in, FILE* out, int len, int begin, int &hdrsize);
int decode_exe(Encoder& en, FILE* in, int size, int info, FILE *out, FMode mode, int &diffFound);