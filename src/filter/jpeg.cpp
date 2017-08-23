#include "jpeg.h"

bool DetectJPEG::detect(int c, int i, Filetype type) {
  buf1=buf1<<8|buf0>>24;
  buf0=buf0<<8|c;

  // Detect JPEG by code SOI APPx (FF D8 FF Ex) followed by
  // SOF0 (FF C0 xx xx 08) and SOS (FF DA) within a reasonable distance.
  // Detect end by any code other than RST0-RST7 (FF D9-D7) or
  // a byte stuff (FF 00).
  if (!soi && i>=3 && (buf0&0xffffff00)==0xffd8ff00 && ((U8)buf0==0xC0 || (U8)buf0==0xC4 || ((U8)buf0>=0xDB && (U8)buf0<=0xFE) )) soi=i, app=i+2, sos=sof=0;
  if (soi) {
    if (app==i && (buf0>>24)==0xff &&
       ((buf0>>16)&0xff)>0xc0 && ((buf0>>16)&0xff)<0xff) app=i+(buf0&0xffff)+2;
    if (app<i && (buf1&0xff)==0xff && (buf0&0xff0000ff)==0xc0000008) sof=i;
    if (sof && sof>soi && i-sof<0x1000 && (buf0&0xffff)==0xffda) {
      sos=i;
      if (type!=JPEG) return f->detectBlock(JPEG, soi-3);
    }
    if (i-soi>0x40000 && !sos) soi=0;
  }
  if (type==JPEG && sos && i>sos && (buf0&0xff00)==0xff00
      && (buf0&0xff)!=0 && (buf0&0xf8)!=0xd0) return f->detectBlock(DEFAULT, i);

  return false;
}