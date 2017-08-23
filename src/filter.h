/////////////////////////// Filters /////////////////////////////////
//
// Before compression, data is encoded in blocks with the following format:
//
//   <type> <size> <encoded-data>
//
// Type is 1 byte (type Filetype): DEFAULT=0, JPEG, EXE
// Size is 4 bytes in big-endian format.
// Encoded-data decodes to <size> bytes.  The encoded size might be
// different.  Encoded data is designed to be more compressible.
//
//   void encode(FILE* in, FILE* out, int n);
//
// Reads n bytes of in (open in "rb" mode) and encodes one or
// more blocks to temporary file out (open in "wb+" mode).
// The file pointer of in is advanced n bytes.  The file pointer of
// out is positioned after the last byte written.
//
//   en.setFile(FILE* out);
//   int decode(Encoder& en);
//
// Decodes and returns one byte.  Input is from en.decompress(), which
// reads from out if in COMPRESS mode.  During compression, n calls
// to decode() must exactly match n bytes of in, or else it is compressed
// as type 0 without encoding.
//
//   Filetype detect(FILE* in, int n, Filetype type);
//
// Reads n bytes of in, and detects when the type changes to
// something else.  If it does, then the file pointer is repositioned
// to the start of the change and the new type is returned.  If the type
// does not change, then it repositions the file pointer n bytes ahead
// and returns the old type.
//
// For each type X there are the following 2 functions:
//
//   void encode_X(FILE* in, FILE* out, int n, ...);
//
// encodes n bytes from in to out.
//
//   int decode_X(Encoder& en);
//
// decodes one byte from en and returns it.  decode() and decode_X()
// maintain state information using static variables.

#ifndef FILTER_H_INCLUDED
#define FILTER_H_INCLUDED

#include "misc.h"
#include "encoder.h"

typedef enum {FDECOMPRESS, FCOMPARE, FDISCARD} FMode;

class Filter {
  Filetype type, blockType, nextBlockType;
  int info, nextBlockOffset, nextBlockLength, blockLength;

public:
  Filter(): nextBlockType(DEFAULT), nextBlockOffset(0) {};
  void compressRecursive(FILE *in, long n, Encoder &en, char *blstr, int it=0, float p1=0.0, float p2=1.0);
  int decompressRecursive(FILE *out, long n, Encoder& en, FMode mode, int it=0);
  bool detect(FILE* in, Filetype type, int n);
  void direct_encode_block(Filetype type, FILE *in, int len, Encoder &en, int info=-1);
  void transform_encode_block(Filetype type, FILE *in, int len, Encoder &en, int info, char *blstr, int it, float p1, float p2, long begin);
  
  bool detectBlock(Filetype foundType, int foundPosition, int foundLength = 0, int blockInfo = 0, int hdrLength = 0) {
    blockType = hdrLength > 0 ? HDR : foundType;
    blockLength = foundPosition;
    nextBlockLength=foundLength;
    info=blockInfo;
    if (hdrLength > 0) {
      nextBlockType = foundType;
      nextBlockOffset = hdrLength;
    }
    return true;
  }
};

#endif // #ifndef FILTER_H_INCLUDED