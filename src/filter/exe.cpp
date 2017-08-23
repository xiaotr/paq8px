#include "exe.h"

// EXE transform: <encoded-size> <begin> <block>...
// Encoded-size is 4 bytes, MSB first.
// begin is the offset of the start of the input file, 4 bytes, MSB first.
// Each block applies the e8e9 transform to strings falling entirely
// within the block starting from the end and working backwards.
// The 5 byte pattern is E8/E9 xx xx xx 00/FF (x86 CALL/JMP xxxxxxxx)
// where xxxxxxxx is a relative address LSB first.  The address is
// converted to an absolute address by adding the offset mod 2^25
// (in range +-2^24).

int encode_exe(FILE* in, FILE* out, int len, int begin, int &hdrsize) {
  (void)(hdrsize);
  const int BLOCK=0x10000;
  Array<U8> blk(BLOCK);
  fprintf(out, "%c%c%c%c", begin>>24, begin>>16, begin>>8, begin);

  // Transform
  for (int offset=0; offset<len; offset+=BLOCK) {
    int size=min(len-offset, BLOCK);
    int bytesRead=fread(&blk[0], 1, size, in);
    if (bytesRead!=size) quit("encode_exe read error");
    for (int i=bytesRead-1; i>=5; --i) {
      if ((blk[i-4]==0xe8 || blk[i-4]==0xe9 || (blk[i-5]==0x0f && (blk[i-4]&0xf0)==0x80))
         && (blk[i]==0||blk[i]==0xff)) {
        int a=(blk[i-3]|blk[i-2]<<8|blk[i-1]<<16|blk[i]<<24)+offset+begin+i+1;
        a<<=7;
        a>>=7;
        blk[i]=a>>24;
        blk[i-1]=a^176;
        blk[i-2]=(a>>8)^176;
        blk[i-3]=(a>>16)^176;
      }
    }
    fwrite(&blk[0], 1, bytesRead, out);
  }
  return 1;
}

int decode_exe(Encoder& en, FILE* in, int size, int info, FILE *out, FMode mode, int &diffFound) {
  (void)(in);
  (void)(info);
  const int BLOCK=0x10000;  // block size
  int begin, offset=6, a;
  U8 c[6];
  begin=en.decompress()<<24;
  begin|=en.decompress()<<16;
  begin|=en.decompress()<<8;
  begin|=en.decompress();
  size-=4;
  for (int i=4; i>=0; i--) c[i]=en.decompress();  // Fill queue

  while (offset<size+6) {
    memmove(c+1, c, 5);
    if (offset<=size) c[0]=en.decompress();
    // E8E9 transform: E8/E9 xx xx xx 00/FF -> subtract location from x
    if ((c[0]==0x00 || c[0]==0xFF) && (c[4]==0xE8 || c[4]==0xE9 || (c[5]==0x0F && (c[4]&0xF0)==0x80))
     && (((offset-1)^(offset-6))&-BLOCK)==0 && offset<=size) { // not crossing block boundary
      a=((c[1]^176)|(c[2]^176)<<8|(c[3]^176)<<16|c[0]<<24)-offset-begin;
      a<<=7;
      a>>=7;
      c[3]=a;
      c[2]=a>>8;
      c[1]=a>>16;
      c[0]=a>>24;
    }
    if (mode==FDECOMPRESS) putc(c[5], out);
    else if (mode==FCOMPARE && c[5]!=getc(out) && !diffFound) diffFound=offset-6+1;
    if (mode==FDECOMPRESS && !(offset&0xfff)) en.print_status();
    offset++;
  }
  return size;
}


bool DetectEXE::detect(int c, int i, int start, Filetype type) {
  buf1=buf1<<8|buf0>>24;
  buf0=buf0<<8|c;

  // Detect EXE if the low order byte (little-endian) XX is more
  // recently seen (and within 4K) if a relative to absolute address
  // conversion is done in the context CALL/JMP (E8/E9) XX xx xx 00/FF
  // 4 times in a row.  Detect end of EXE at the last
  // place this happens when it does not happen for 64KB.

  if (((buf1&0xfe)==0xe8 || (buf1&0xfff0)==0x0f80) && ((buf0+1)&0xfe)==0) {
    int r=buf0>>24;  // relative address low 8 bits
    int a=((buf0>>24)+i)&0xff;  // absolute address low 8 bits
    int rdist=i-relpos[r];
    int adist=i-abspos[a];
    if (adist<rdist && adist<0x800 && abspos[a]>5) {
      e8e9last=i;
      ++e8e9count;
      if (e8e9pos==0 || e8e9pos>abspos[a]) e8e9pos=abspos[a];
    }
    else e8e9count=0;
    if (type==DEFAULT && e8e9count>=4 && e8e9pos>5)
      return f->detectBlock(EXE, e8e9pos-5);
    abspos[a]=i;
    relpos[r]=i;
  }
  if (i-e8e9last>0x4000) {
    if (type==EXE) return f->detectBlock(DEFAULT, e8e9last, 0, start);
    e8e9count=e8e9pos=0;
  }
  return false;
}