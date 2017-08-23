#include "image.h"

bool IsGrayscalePalette(FILE* in, int n = 256, int isRGBA = 0){
  long offset = ftell(in);
  int stride = 3+isRGBA, res = (n>0)<<8, order=1;
  for (int i = 0; (i < n*stride) && (res>>8); i++) {
    int b = getc(in);
    if (b==EOF){
      res = 0;
      break;
    }
    if (!i) {
      res = 0x100|b;
      order = 1-2*(b>0);
      continue;
    }

    //"j" is the index of the current byte in this color entry
    int j = i%stride;
    if (!j)
      res = (res&((b-(res&0xFF)==order)<<8))|b; // load first component of this entry
    else if (j==3)
      res&=((!b || (b==0xFF))*0x1FF); // alpha/attribute component must be zero or 0xFF
    else
      res&=((b==(res&0xFF))<<9)-1;
  }
  fseek(in, offset, SEEK_SET);
  return res>>8;
}



// 24-bit image data transform:
// simple color transform (b, g, r) -> (g, g-r, g-b)

int encode_bmp(FILE* in, FILE* out, int len, int width, int &hdrsize) {
  (void)(hdrsize);
  int r,g,b;
  for (int i=0; i<len/width; i++) {
    for (int j=0; j<width/3; j++) {
      b=fgetc(in), g=fgetc(in), r=fgetc(in);
      fputc(g, out);
      fputc(g-r, out);
      fputc(g-b, out);
    }
    for (int j=0; j<width%3; j++) fputc(fgetc(in), out);
  }
  return 1;
}

int decode_bmp(Encoder& en, FILE *in, int size, int width, FILE *out, FMode mode, int &diffFound) {
  (void)(in);
  int r,g,b,p;
  for (int i=0; i<size/width; i++) {
    p=i*width;
    for (int j=0; j<width/3; j++) {
      b=en.decompress(), g=en.decompress(), r=en.decompress();
      if (mode==FDECOMPRESS) {
        fputc(b-r, out);
        fputc(b, out);
        fputc(b-g, out);
        if (!j && !(i&0xf)) en.print_status();
      }
      else if (mode==FCOMPARE) {
        if (((b-r)&255)!=getc(out) && !diffFound) diffFound=p+1;
        if (b!=getc(out) && !diffFound) diffFound=p+2;
        if (((b-g)&255)!=getc(out) && !diffFound) diffFound=p+3;
        p+=3;
      }
    }
    for (int j=0; j<width%3; j++) {
      if (mode==FDECOMPRESS) {
        fputc(en.decompress(), out);
      }
      else if (mode==FCOMPARE) {
        if (en.decompress()!=getc(out) && !diffFound) diffFound=p+j+1;
      }
    }
  }
  return size;
}

// 32-bit image
int encode_im32(FILE* in, FILE* out, int len, int width, int &hdrsize) {
  (void)(hdrsize);
  int r,g,b,a;
  for (int i=0; i<len/width; i++) {
    for (int j=0; j<width/4; j++) {
      b=fgetc(in), g=fgetc(in), r=fgetc(in); a=fgetc(in);
      fputc(g, out);
      fputc(g-r, out);
      fputc(g-b, out);
      fputc(a, out);
    }
    for (int j=0; j<width%4; j++) fputc(fgetc(in), out);
  }
  return 1;
}

int decode_im32(Encoder& en, FILE *in, int size, int width, FILE *out, FMode mode, int &diffFound) {
  (void)(in);
  int r,g,b,a,p;
  bool rgb = (width&(1<<31))>0;
  if (rgb) width^=(1<<31);
  for (int i=0; i<size/width; i++) {
    p=i*width;
    for (int j=0; j<width/4; j++) {
      b=en.decompress(), g=en.decompress(), r=en.decompress(), a=en.decompress();
      if (mode==FDECOMPRESS) {
        fputc(b-r, out); fputc(b, out); fputc(b-g, out); fputc(a, out);
        if (!j && !(i&0xf)) en.print_status();
      }
      else if (mode==FCOMPARE) {
        if (((b-r)&255)!=getc(out) && !diffFound) diffFound=p+1;
        if (b!=getc(out) && !diffFound) diffFound=p+2;
        if (((b-g)&255)!=getc(out) && !diffFound) diffFound=p+3;
        if (((a)&255)!=getc(out) && !diffFound) diffFound=p+4;
        p+=4;
      }
    }
    for (int j=0; j<width%4; j++) {
      if (mode==FDECOMPRESS) {
        fputc(en.decompress(), out);
      }
      else if (mode==FCOMPARE) {
        if (en.decompress()!=getc(out) && !diffFound) diffFound=p+j+1;
      }
    }
  }
  return size;
}

bool DetectIMAGE::detect(int c, int i, long start, int n, FILE *in) {
  buf2=buf2<<8|buf1>>24;
  buf1=buf1<<8|buf0>>24;
  buf0=buf0<<8|c;

  // Detect .bmp image
  if ( !bmp && (((buf0&0xffff)==16973) || (!(buf0&0xFFFFFF) && ((buf0>>24)==0x28))) ) //possible 'BM' or headerless DIB
    imgbpp=bmpx=bmpy=0,hdrless=!(U8)buf0,bmpof=hdrless*68,bmp=i-hdrless*16;
  if (bmp) {
    const int p=i-bmp;
    if (p==12) bmpof=bswap(buf0);
    else if (p==16 && buf0!=0x28000000) bmp=0; //BITMAPINFOHEADER (0x28)
    else if (p==20) bmpx=bswap(buf0),bmp=((bmpx==0||bmpx>0x30000)?0:bmp); //width
    else if (p==24) bmpy=abs((int)bswap(buf0)),bmp=((bmpy==0||bmpy>0x10000)?0:bmp); //height
    else if (p==27) imgbpp=c,bmp=((imgbpp!=1 && imgbpp!=8 && imgbpp!=24 && imgbpp!=32)?0:bmp);
    else if ((p==31) && buf0) bmp=0;
    // check number of colors in palette (4 bytes), must be 0 (default) or <= 1<<bpp.
    // also check if image is too small, since it might not be worth it to use the image models
    else if (p==48){
      if ( (!buf0 || ((bswap(buf0)<=(U32)(1<<imgbpp)) && (imgbpp<=8))) && (((bmpx*bmpy*imgbpp)>>3)>512) ) {
        // possible icon/cursor?
        if (hdrless && (bmpx*2==bmpy) && (
          (bmpx==8)   || (bmpx==10) || (bmpx==14) || (bmpx==16) || (bmpx==20) ||
          (bmpx==22)  || (bmpx==24) || (bmpx==32) || (bmpx==40) || (bmpx==48) ||
          (bmpx==60)  || (bmpx==64) || (bmpx==72) || (bmpx==80) || (bmpx==96) ||
          (bmpx==128) || (bmpx==256)
        ))
          bmpy=bmpx;

        // if DIB and not 24bpp, we must calculate the data offset based on BPP or num. of entries in color palette
        if (hdrless && (imgbpp<=24))
          bmpof+=(buf0)?bswap(buf0)*4:4<<imgbpp;

        if (imgbpp==1) f->detectBlock(IMAGE1, bmp-1, (((bmpx-1)>>5)+1)*4, (((bmpx-1)>>5)+1)*4*bmpy, bmpof);
        else if (imgbpp==8){
          fseek(in, start+bmp+53, SEEK_SET);
          f->detectBlock((IsGrayscalePalette(in, (buf0)?bswap(buf0):1<<imgbpp, 1))?IMAGE8GRAY:IMAGE8, bmp-1, (bmpx+3)&-4, ((bmpx+3)&-4)*bmpy, bmpof);
        }
        else if (imgbpp==24) f->detectBlock(IMAGE24, bmp-1, ((bmpx*3)+3)&-4, (((bmpx*3)+3)&-4)*bmpy, bmpof);
        else if (imgbpp==32) f->detectBlock(IMAGE32, bmp-1, bmpx*4, bmpx*4*bmpy, bmpof);
      }
      bmp=0;
    }
  }

  // Detect .pbm .pgm .ppm .pam image
  if ((buf0&0xfff0ff)==0x50300a) {
    pgmn=(buf0&0xf00)>>8;
    if ((pgmn>=4 && pgmn<=6) || pgmn==7) pgm=i,pgm_ptr=pgmw=pgmh=pgmc=pgmcomment=pamatr=pamd=0;
  }
  if (pgm) {
    if (i-pgm==1 && c==0x23) pgmcomment=1; //pgm comment
    if (!pgmcomment && pgm_ptr) {
      int s=0;
      if (pgmn==7) {
         if ((buf1&0xdfdf)==0x5749 && (buf0&0xdfdfdfff)==0x44544820) pgm_ptr=0, pamatr=1; // WIDTH
         if ((buf1&0xdfdfdf)==0x484549 && (buf0&0xdfdfdfff)==0x47485420) pgm_ptr=0, pamatr=2; // HEIGHT
         if ((buf1&0xdfdfdf)==0x4d4158 && (buf0&0xdfdfdfff)==0x56414c20) pgm_ptr=0, pamatr=3; // MAXVAL
         if ((buf1&0xdfdf)==0x4445 && (buf0&0xdfdfdfff)==0x50544820) pgm_ptr=0, pamatr=4; // DEPTH
         if ((buf2&0xdf)==0x54 && (buf1&0xdfdfdfdf)==0x55504c54 && (buf0&0xdfdfdfff)==0x59504520) pgm_ptr=0, pamatr=5; // TUPLTYPE
         if ((buf1&0xdfdfdf)==0x454e44 && (buf0&0xdfdfdfff)==0x4844520a) pgm_ptr=0, pamatr=6; // ENDHDR
         if (c==0x0a) {
           if (pamatr==0) pgm=0;
           else if (pamatr<5) s=pamatr;
           if (pamatr!=6) pamatr=0;
         }
      } else if (c==0x20 && !pgmw) s=1;
      else if (c==0x0a && !pgmh) s=2;
      else if (c==0x0a && !pgmc && pgmn!=4) s=3;
      if (s) {
        pgm_buf[pgm_ptr++]=0;
        int v=atoi(pgm_buf);
        if (s==1) pgmw=v; else if (s==2) pgmh=v; else if (s==3) pgmc=v; else if (s==4) pamd=v;
        if (v==0 || (s==3 && v>255)) pgm=0; else pgm_ptr=0;
      }
    }
    if (!pgmcomment) pgm_buf[pgm_ptr++]=c;
    if (pgm_ptr>=32) pgm=0;
    if (pgmcomment && c==0x0a) pgmcomment=0;
    if (pgmw && pgmh && !pgmc && pgmn==4) f->detectBlock(IMAGE1, pgm-2, (pgmw+7)/8, (pgmw+7)/8*pgmh, i-pgm+3);
    if (pgmw && pgmh && pgmc && (pgmn==5 || (pgmn==7 && pamd==1 && pamatr==6))) f->detectBlock(IMAGE8GRAY, pgm-2, pgmw, pgmw*pgmh, i-pgm+3);
    if (pgmw && pgmh && pgmc && (pgmn==6 || (pgmn==7 && pamd==3 && pamatr==6))) f->detectBlock(IMAGE24, pgm-2, pgmw*3, pgmw*3*pgmh, i-pgm+3);
    if (pgmw && pgmh && pgmc && (pgmn==7 && pamd==4 && pamatr==6)) f->detectBlock(IMAGE32, pgm-2, pgmw*4, pgmw*4*pgmh, i-pgm+3);
  }

  // Detect .rgb image
  if ((buf0&0xffff)==0x01da) rgbi=i,rgbx=rgby=0;
  if (rgbi) {
    const int p=i-rgbi;
    if (p==1 && c!=0) rgbi=0;
    else if (p==2 && c!=1) rgbi=0;
    else if (p==4 && (buf0&0xffff)!=1 && (buf0&0xffff)!=2 && (buf0&0xffff)!=3) rgbi=0;
    else if (p==6) rgbx=buf0&0xffff,rgbi=(rgbx==0?0:rgbi);
    else if (p==8) rgby=buf0&0xffff,rgbi=(rgby==0?0:rgbi);
    else if (p==10) {
      int z=buf0&0xffff;
      if (rgbx && rgby && (z==1 || z==3 || z==4)) f->detectBlock(IMAGE8, rgbi-1, rgbx, rgbx*rgby*z, 512);
      rgbi=0;
    }
  }

  // Detect .tiff file header (2/8/24 bit color, not compressed).
  if (buf1==0x49492a00 && n>i+(int)bswap(buf0)) {
    long savedpos=ftell(in);
    fseek(in, start+i+bswap(buf0)-7, SEEK_SET);

    // read directory
    int dirsize=getc(in);
    int tifx=0,tify=0,tifz=0,tifzb=0,tifc=0,tifofs=0,tifofval=0,b[12];
    if (getc(in)==0) {
      for (int i=0; i<dirsize; i++) {
        for (int j=0; j<12; j++) b[j]=getc(in);
        if (b[11]==EOF) break;
        int tag=b[0]+(b[1]<<8);
        int tagfmt=b[2]+(b[3]<<8);
        int taglen=b[4]+(b[5]<<8)+(b[6]<<16)+(b[7]<<24);
        int tagval=b[8]+(b[9]<<8)+(b[10]<<16)+(b[11]<<24);
        if (tagfmt==3||tagfmt==4) {
          if (tag==256) tifx=tagval;
          else if (tag==257) tify=tagval;
          else if (tag==258) tifzb=taglen==1?tagval:8; // bits per component
          else if (tag==259) tifc=tagval; // 1 = no compression
          else if (tag==273 && tagfmt==4) tifofs=tagval,tifofval=(taglen<=1);
          else if (tag==277) tifz=tagval; // components per pixel
        }
      }
    }
    if (tifx && tify && tifzb && (tifz==1 || tifz==3) && (tifc==1) && (tifofs && tifofs+i<n)) {
      if (!tifofval) {
        fseek(in, start+i+tifofs-7, SEEK_SET);
        for (int j=0; j<4; j++) b[j]=getc(in);
        tifofs=b[0]+(b[1]<<8)+(b[2]<<16)+(b[3]<<24);
      }
      if (tifofs && tifofs<(1<<18) && tifofs+i<n) {
        if (tifz==1 && tifzb==1) f->detectBlock(IMAGE1, i-7, ((tifx-1)>>3)+1, (((tifx-1)>>3)+1)*tify, tifofs);
        else if (tifz==1 && tifzb==8) f->detectBlock(IMAGE8, i-7, tifx, tifx*tify, tifofs);
        else if (tifz==3 && tifzb==8) f->detectBlock(IMAGE24, i-7, tifx*3, tifx*3*tify, tifofs);
      }
    }
    fseek(in, savedpos, SEEK_SET);
  }

  // Detect .tga image (8-bit 256 colors or 24-bit uncompressed)
  if (buf1==0x00010100 && buf0==0x00000118) tga=i,tgax=tgay,tgaz=8,tgat=1;
  else if (buf1==0x00000200 && buf0==0x00000000) tga=i,tgax=tgay,tgaz=24,tgat=2;
  else if (buf1==0x00000300 && buf0==0x00000000) tga=i,tgax=tgay,tgaz=8,tgat=3;
  if (tga) {
    if (i-tga==8) tga=(buf1==0?tga:0),tgax=(bswap(buf0)&0xffff),tgay=(bswap(buf0)>>16);
    else if (i-tga==10) {
      if (tgaz==(int)((buf0&0xffff)>>8) && tgax && tgay) {
        if (tgat==1){
          fseek(in, start+tga+11, SEEK_SET);
          f->detectBlock(IsGrayscalePalette(in)?IMAGE8GRAY:IMAGE8, tga-7, tgax, tgax*tgay, 18+256*3);
        }
        else if (tgat==2) f->detectBlock(IMAGE24, tga-7, tgax*3, tgax*3*tgay, 18);
        else if (tgat==3) f->detectBlock(IMAGE8, tga-7, tgax, tgax*tgay, 18);
      }
      tga=0;
    }
  }

  return false;
}
