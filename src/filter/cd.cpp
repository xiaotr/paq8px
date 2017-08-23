#include "cd.h"

// Function ecc_compute(), edc_compute() and eccedc_init() taken from
// ** UNECM - Decoder for ECM (Error Code Modeler) format.
// ** Version 1.0
// ** Copyright (C) 2002 Neill Corlett

/* LUTs used for computing ECC/EDC */
static U8 ecc_f_lut[256];
static U8 ecc_b_lut[256];
static U32 edc_lut[256];
static int luts_init=0;

void eccedc_init(void) {
  if (luts_init) return;
  U32 i, j, edc;
  for(i = 0; i < 256; i++) {
    j = (i << 1) ^ (i & 0x80 ? 0x11D : 0);
    ecc_f_lut[i] = j;
    ecc_b_lut[i ^ j] = i;
    edc = i;
    for(j = 0; j < 8; j++) edc = (edc >> 1) ^ (edc & 1 ? 0xD8018001 : 0);
    edc_lut[i] = edc;
  }
  luts_init=1;
}

void ecc_compute(U8 *src, U32 major_count, U32 minor_count, U32 major_mult, U32 minor_inc, U8 *dest) {
  U32 size = major_count * minor_count;
  U32 major, minor;
  for(major = 0; major < major_count; major++) {
    U32 index = (major >> 1) * major_mult + (major & 1);
    U8 ecc_a = 0;
    U8 ecc_b = 0;
    for(minor = 0; minor < minor_count; minor++) {
      U8 temp = src[index];
      index += minor_inc;
      if(index >= size) index -= size;
      ecc_a ^= temp;
      ecc_b ^= temp;
      ecc_a = ecc_f_lut[ecc_a];
    }
    ecc_a = ecc_b_lut[ecc_f_lut[ecc_a] ^ ecc_b];
    dest[major              ] = ecc_a;
    dest[major              ] = ecc_a;
    dest[major + major_count] = ecc_a ^ ecc_b;
  }
}

U32 edc_compute(const U8  *src, int size) {
  U32 edc = 0;
  while(size--) edc = (edc >> 8) ^ edc_lut[(edc ^ (*src++)) & 0xFF];
  return edc;
}

int expand_cd_sector(U8 *data, int a, int test) {
  U8 d2[2352];
  eccedc_init();
  d2[0]=d2[11]=0;
  for (int i=1; i<11; i++) d2[i]=255;
  int mode=(data[15]!=1?2:1);
  int form=(data[15]==3?2:1);
  if (a==-1) for (int i=12; i<15; i++) d2[i]=data[i]; else {
    int c1=(a&15)+((a>>4)&15)*10;
    int c2=((a>>8)&15)+((a>>12)&15)*10;
    int c3=((a>>16)&15)+((a>>20)&15)*10;
    c1=(c1+1)%75;
    if (c1==0) {
      c2=(c2+1)%60;
      if (c2==0) c3++;
    }
    d2[12]=(c3%10)+16*(c3/10);
    d2[13]=(c2%10)+16*(c2/10);
    d2[14]=(c1%10)+16*(c1/10);
  }
  d2[15]=mode;
  if (mode==2) for (int i=16; i<24; i++) d2[i]=data[i-4*(i>=20)];
  if (form==1) {
    if (mode==2) {
      d2[1]=d2[12],d2[2]=d2[13],d2[3]=d2[14];
      d2[12]=d2[13]=d2[14]=d2[15]=0;
    } else {
      for(int i=2068; i<2076; i++) d2[i]=0;
    }
    for (int i=16+8*(mode==2); i<2064+8*(mode==2); i++) d2[i]=data[i];
    U32 edc=edc_compute(d2+16*(mode==2), 2064-8*(mode==2));
    for (int i=0; i<4; i++) d2[2064+8*(mode==2)+i]=(edc>>(8*i))&0xff;
    ecc_compute(d2+12, 86, 24,  2, 86, d2+2076);
    ecc_compute(d2+12, 52, 43, 86, 88, d2+2248);
    if (mode==2) {
      d2[12]=d2[1],d2[13]=d2[2],d2[14]=d2[3],d2[15]=2;
      d2[1]=d2[2]=d2[3]=255;
    }
  }
  for (int i=0; i<2352; i++) if (d2[i]!=data[i] && test) form=2;
  if (form==2) {
    for (int i=24; i<2348; i++) d2[i]=data[i];
    U32 edc=edc_compute(d2+16, 2332);
    for (int i=0; i<4; i++) d2[2348+i]=(edc>>(8*i))&0xff;
  }
  for (int i=0; i<2352; i++) if (d2[i]!=data[i] && test) return 0; else data[i]=d2[i];
  return mode+form-1;
}



int encode_cd(FILE* in, FILE* out, int len, int info, int &hdrsize) {
  (void)(hdrsize);
  const int BLOCK=2352;
  U8 blk[BLOCK];
  fputc((len%BLOCK)>>8,out);
  fputc(len%BLOCK,out);
  for (int offset=0; offset<len; offset+=BLOCK) {
    if (offset+BLOCK > len) {
      fread(&blk[0], 1, len-offset, in);
      fwrite(&blk[0], 1, len-offset, out);
    } else {
      fread(&blk[0], 1, BLOCK, in);
      if (info==3) blk[15]=3;
      if (offset==0) fwrite(&blk[12], 1, 4+4*(blk[15]!=1), out);
      fwrite(&blk[16+8*(blk[15]!=1)], 1, 2048+276*(info==3), out);
      if (offset+BLOCK*2 > len && blk[15]!=1) fwrite(&blk[16], 1, 4, out);
    }
  }
  return 1;
}

int decode_cd(Encoder& en, FILE *in, int size, int info, FILE *out, FMode mode, int &diffFound) {
  (void)(en);
  (void)(info);
  const int BLOCK=2352;
  U8 blk[BLOCK];
  long i=0, i2=0;
  int a=-1, bsize=0, q=fgetc(in);
  q=(q<<8)+fgetc(in);
  size-=2;
  while (i<size) {
    if (size-i==q) {
      fread(blk, q, 1, in);
      fwrite(blk, q, 1, out);
      i+=q;
      i2+=q;
    } else if (i==0) {
      fread(blk+12, 4, 1, in);
      if (blk[15]!=1) fread(blk+16, 4, 1, in);
      bsize=2048+(blk[15]==3)*276;
      i+=4*(blk[15]!=1)+4;
    } else {
      a=(blk[12]<<16)+(blk[13]<<8)+blk[14];
    }
    fread(blk+16+(blk[15]!=1)*8, bsize, 1, in);
    i+=bsize;
    if (bsize>2048) blk[15]=3;
    if (blk[15]!=1 && size-q-i==4) {
      fread(blk+16, 4, 1, in);
      i+=4;
    }
    expand_cd_sector(blk, a, 0);
    if (mode==FDECOMPRESS) fwrite(blk, BLOCK, 1, out);
    else if (mode==FCOMPARE) for (int j=0; j<BLOCK; ++j) if (blk[j]!=getc(out) && !diffFound) diffFound=i2+j+1;
    i2+=BLOCK;
  }
  return i2;
}



bool DetectCD::detect(int c, int i, long start, int n, FILE *in, Filetype type) {
  buf1=buf1<<8|buf0>>24;
  buf0=buf0<<8|c;
  if (buf1==0x00ffffff && buf0==0xffffffff && !cdi) cdi=i,cda=-1,cdm=0;
  if (cdi && i>cdi) {
    const int p=(i-cdi)%2352;
    if (p==8 && (buf1!=0xffffff00 || ((buf0&0xff)!=1 && (buf0&0xff)!=2))) cdi=0;
    else if (p==16 && i+2336<n) {
      U8 data[2352];
      long savedpos=ftell(in);
      fseek(in, start+i-23, SEEK_SET);
      fread(data, 1, 2352, in);
      fseek(in, savedpos, SEEK_SET);
      int t=expand_cd_sector(data, cda, 1);
      if (t!=cdm) cdm=t*(i-cdi<2352);
      if (cdm && cda!=10 && (cdm==1 || buf0==buf1)) {
        if (type!=CD) return f->detectBlock(CD, cdi-7, 0, cdm);
        cda=(data[12]<<16)+(data[13]<<8)+data[14];
        if (cdm!=1 && i-cdi>2352 && buf0!=cdf) cda=10;
        if (cdm!=1) cdf=buf0;
      } else cdi=0;
    }
    if (!cdi && type==CD) return f->detectBlock(DEFAULT, i-p-7);
  }
  return false;
}