#include "misc.h"
#include "encoder.h"
#include "filter/zlib.h"
#include "filter/cd.h"
#include "filter/jpeg.h"
#include "filter/audio.h"
#include "filter/image.h"
#include "filter/gif.h"
#include "filter/exe.h"
#include "filter/base64.h"
#include "filter.h"

typedef struct
{
    Filetype type;
    const char *name;
    int (*encode_func)(FILE*,FILE*,int,int,int&);
    int (*decode_func)(Encoder&,FILE*,int,int,FILE*,FMode,int&);
} BlockDescription;

const BlockDescription blockDescription[] = {
  { DEFAULT   , "default"         , NULL         , NULL },
  { JPEG      , "jpeg"            , NULL         , NULL },
  { HDR       , "hdr"             , NULL         , NULL },
  { IMAGE1    , "1b-image"        , NULL         , NULL },
  { IMAGE8    , "8b-image"        , NULL         , NULL },
  { IMAGE8GRAY, "8b-img-grayscale", NULL         , NULL },
  { IMAGE24   , "24b-image"       , encode_bmp   , decode_bmp },
  { IMAGE32   , "32b-image"       , encode_im32  , decode_im32 },
  { AUDIO     , "audio"           , NULL         , NULL },
  { EXE       , "exe"             , encode_exe   , decode_exe },
  { CD        , "cd"              , encode_cd    , decode_cd },
  { ZLIB      , "zlib"            , encode_zlib  , decode_zlib },
  { BASE64    , "base64"          , encode_base64, decode_base64 },
  { GIF       , "gif"             , encode_gif   , decode_gif },
};

// Detect blocks
bool Filter::detect(FILE* in, Filetype type, int n) {
  DetectZLIB detectZLIB(this);
  DetectCD detectCD(this);
  DetectJPEG detectJPEG(this);
  DetectAUDIO detectAUDIO(this);
  DetectIMAGE detectIMAGE(this);
  DetectGIF detectGIF(this);
  DetectEXE detectEXE(this);
  DetectBASE64 detectBASE64(this);

  long start=ftell(in);
  for (int i=0; i<n; ++i) {
    int c = getc(in);
    if (c == EOF) break;
    if (detectZLIB.detect(c, i, start, n, in)) return true;
    if (detectCD.detect(c, i, start, n, in, type)) return true;
    if (type == CD) continue;
    if (detectJPEG.detect(c, i, type)) return true;
    if (detectAUDIO.detect(c, i, start, in)) return true;
    if (detectIMAGE.detect(c, i, start, n, in)) return true;
    if (detectGIF.detect(c, i, type, nextBlockType)) return true;
    if (detectEXE.detect(c, i, start, type)) return true;
    if (detectBASE64.detect(c, i)) return true;
  }
  return detectBlock(type, n);
}

void Filter::direct_encode_block(Filetype type, FILE *in, int len, Encoder &en, int info) {
  en.compress(type);
  en.compress(len>>24);
  en.compress(len>>16);
  en.compress(len>>8);
  en.compress(len);
  if (info!=-1) {
    en.compress(info>>24);
    en.compress(info>>16);
    en.compress(info>>8);
    en.compress(info);
  }
  printf("Compressing... ");
  for (int j=0; j<len; ++j) {
    if (!(j&0xfff)) en.print_status(j, len);
    en.compress(getc(in));
  }
  printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
}


void Filter::transform_encode_block(Filetype type, FILE *in, int len, Encoder &en, int info, char *blstr, int it, float p1, float p2, long begin) {
  if (blockDescription[type].encode_func!=NULL) {
    FILE* tmp=tmpfile();  // temporary encoded file
    if (!tmp) perror("tmpfile"), quit();
    int diffFound=0, hdrsize=0;
    diffFound=blockDescription[type].encode_func(in, tmp, len, info, hdrsize)?0:1;
    const long tmpsize=ftell(tmp);
    fseek(tmp, tmpsize, SEEK_SET);
    if (!diffFound) {
      rewind(tmp);
      en.setFile(tmp);
      fseek(in, begin, SEEK_SET);
      blockDescription[type].decode_func(en, tmp, tmpsize, info, in, FCOMPARE, diffFound);
    }
    // Test fails, compress without transform
    if (diffFound || fgetc(tmp)!=EOF) {
      printf("Transform fails at %d, skipping...\n", diffFound-1);
      fseek(in, begin, SEEK_SET);
      direct_encode_block(DEFAULT, in, len, en);
    } else {
      rewind(tmp);
      if (hasRecursion(type)) {
        en.compress(type), en.compress(tmpsize>>24), en.compress(tmpsize>>16);
        en.compress(tmpsize>>8), en.compress(tmpsize);
        Filetype type2=(Filetype)(info>>24);
        if (type2!=DEFAULT) {
          direct_encode_block(HDR, tmp, hdrsize, en);
          transform_encode_block(type2, tmp, tmpsize-hdrsize, en, info&0xffffff, blstr, it, p1, p2, hdrsize);
        } else {
          compressRecursive(tmp, tmpsize, en, blstr, it+1, p1, p2);
        }
      } else {
        direct_encode_block(type, tmp, tmpsize, en, hasInfo(type)?info:-1);
      }
    }
    fclose(tmp);  // deletes
  } else {
    direct_encode_block(type, in, len, en, hasInfo(type)?info:-1);
  }
}

void Filter::compressRecursive(FILE *in, long n, Encoder &en, char *blstr, int it, float p1, float p2) {
  static const char* audiotypes[4]={"8b mono", "8b stereo", "16b mono", "16b stereo"};
  blockType=DEFAULT;
  int blnum=0;  // image width or audio type
  long begin=ftell(in), end0=begin+n;
  char b2[32];
  strcpy(b2, blstr);
  if (b2[0]) strcat(b2, "-");
  if (it==5) {
    direct_encode_block(DEFAULT, in, n, en);
    return;
  }
  float pscale=n>0?(p2-p1)/n:0;

  // Transform and test in blocks
  while (n>0) {
    Filetype type = blockType;
    if (nextBlockOffset) { blockLength=nextBlockOffset; blockType=nextBlockType; nextBlockOffset=0; }
    else if (nextBlockLength) { blockLength=nextBlockLength; blockType=DEFAULT; nextBlockLength=0; }
    else detect(in, type, n);
    long end=begin+blockLength;
    if (end>end0) {  // if some detection reports longer then actual size file is
      end=begin+1;
      type=DEFAULT;
    }
    int len=int(end-begin);
    if (len>0) {
      en.set_status_range(p1,p2=p1+pscale*len);
      sprintf(blstr,"%s%d",b2,blnum++);
      assert(blockDescription[type].type==type);
      printf(" %-11s | %-16s |%10d bytes [%ld - %ld]",blstr,blockDescription[type].name,len,begin,begin+len-1);
      if (type==AUDIO) printf(" (%s)", audiotypes[info%4]);
      else if (type==IMAGE1 || type==IMAGE8 || type==IMAGE8GRAY || type==IMAGE24 || type==IMAGE32) printf(" (width: %d)", info);
      else if (type==CD) printf(" (m%d/f%d)", info==1?1:2, info!=3?1:2);
      else if (type==ZLIB && info>0) printf(" (image)");
      printf("\n");
      fseek(in, begin, SEEK_SET);
      transform_encode_block(type, in, len, en, info, blstr, it, p1, p2, begin);
      p1=p2;
      n-=len;
    }
    begin=end;
  }
}

int Filter::decompressRecursive(FILE *out, long n, Encoder& en, FMode mode, int it) {
  Filetype type;
  long len, i=0;
  int diffFound=0, info=0;
  while (i<n) {
    type=(Filetype)en.decompress();
    len=en.decompress()<<24;
    len|=en.decompress()<<16;
    len|=en.decompress()<<8;
    len|=en.decompress();

    if (hasInfo(type)) {
      info=0; for (int i=0; i<4; ++i) { info<<=8; info+=en.decompress(); }
    }
    if (hasRecursion(type)) {
      FILE *tmp=tmpfile();
      if (!tmp) perror("tmpfile"), quit();
      decompressRecursive(tmp, len, en, FDECOMPRESS, it+1);
      if (mode!=FDISCARD) {
        rewind(tmp);
        if (blockDescription[type].decode_func != NULL) len=blockDescription[type].decode_func(en, tmp, len, info, out, mode, diffFound);
      }
      fclose(tmp);
    } else if (blockDescription[type].decode_func != NULL) {
      len=blockDescription[type].decode_func(en, NULL, len, info, out, mode, diffFound);
    } else
    {
      for (int j=0; j<len; ++j) {
        if (!(j&0xfff)) en.print_status();
        if (mode==FDECOMPRESS) putc(en.decompress(), out);
        else if (mode==FCOMPARE) {
          if (en.decompress()!=fgetc(out) && !diffFound) {
            mode=FDISCARD;
            diffFound=i+j+1;
          }
        } else en.decompress();
      }
    }
    i+=len;
  }
  return diffFound;
}


