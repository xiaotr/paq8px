#include "gif.h"

int encode_gif(FILE* in, FILE* out, int len, int info, int &hdrsize) {
  (void)(info);
  int codesize=fgetc(in),diffpos=0,clearpos=0,bsize=0;
  int beginin=ftell(in),beginout=ftell(out);
  U8 output[4096];
  hdrsize=6;
  fputc(hdrsize>>8, out);
  fputc(hdrsize&255, out);
  fputc(bsize, out);
  fputc(clearpos>>8, out);
  fputc(clearpos&255, out);
  fputc(codesize, out);
  for (int phase=0; phase<2; phase++) {
    fseek(in, beginin, SEEK_SET);
    int bits=codesize+1,shift=0,buf=0;
    int blocksize=0,maxcode=(1<<codesize)+1,last=-1,dict[4096];
    bool end=false;
    while ((blocksize=fgetc(in))>0 && ftell(in)-beginin<len && !end) {
      for (int i=0; i<blocksize; i++) {
        buf|=fgetc(in)<<shift;
        shift+=8;
        while (shift>=bits && !end) {
          int code=buf&((1<<bits)-1);
          buf>>=bits;
          shift-=bits;
          if (!bsize && code!=(1<<codesize)) {
            hdrsize+=4; fputc(0, out); fputc(0, out); fputc(0, out); fputc(0, out);
          }
          if (!bsize) bsize=blocksize;
          if (code==(1<<codesize)) {
            if (maxcode>(1<<codesize)+1) {
              if (clearpos && clearpos!=69631-maxcode) return 0;
              clearpos=69631-maxcode;
            }
            bits=codesize+1, maxcode=(1<<codesize)+1, last=-1;
          }
          else if (code==(1<<codesize)+1) end=true;
          else if (code>maxcode+1) return 0;
          else {
            int j=(code<=maxcode?code:last),size=1;
            while (j>=(1<<codesize)) {
              output[4096-(size++)]=dict[j]&255;
              j=dict[j]>>8;
            }
            output[4096-size]=j;
            if (phase==1) fwrite(&output[4096-size], 1, size, out); else diffpos+=size;
            if (code==maxcode+1) { if (phase==1) fputc(j, out); else diffpos++; }
            if (last!=-1) {
              if (++maxcode>=8191) return 0;
              if (maxcode<=4095)
              {
                dict[maxcode]=(last<<8)+j;
                if (phase==0) {
                  bool diff=false;
                  for (int m=(1<<codesize)+2;m<min(maxcode,4095);m++) if (dict[maxcode]==dict[m]) { diff=true; break; }
                  if (diff) {
                    hdrsize+=4;
                    j=diffpos-size-(code==maxcode);
                    fputc((j>>24)&255, out); fputc((j>>16)&255, out); fputc((j>>8)&255, out); fputc(j&255, out);
                    diffpos=size+(code==maxcode);
                  }
                }
              }
              if (maxcode>=((1<<bits)-1) && bits<12) bits++;
            }
            last=code;
          }
        }
      }
    }
  }
  diffpos=ftell(out);
  fseek(out, beginout, SEEK_SET);
  fputc(hdrsize>>8, out);
  fputc(hdrsize&255, out);
  fputc(255-bsize, out);
  fputc((clearpos>>8)&255, out);
  fputc(clearpos&255, out);
  fseek(out, diffpos, SEEK_SET);
  return ftell(in)-beginin==len-1;
}

#define gif_write_block(count) { output[0]=(count);\
if (mode==FDECOMPRESS) fwrite(&output[0], 1, (count)+1, out);\
else if (mode==FCOMPARE) for (int j=0; j<(count)+1; j++) if (output[j]!=getc(out) && !diffFound) diffFound=outsize+j+1;\
outsize+=(count)+1; blocksize=0; }

#define gif_write_code(c) { buf+=(c)<<shift; shift+=bits;\
while (shift>=8) { output[++blocksize]=buf&255; buf>>=8;shift-=8;\
if (blocksize==bsize) gif_write_block(bsize); }}

int decode_gif(Encoder &en, FILE* in, int size, int info, FILE *out, FMode mode, int &diffFound) {
  (void)(en);
  (void)(info);
  int diffcount=fgetc(in), curdiff=0, diffpos[4096];
  diffcount=((diffcount<<8)+fgetc(in)-6)/4;
  int bsize=255-fgetc(in);
  int clearpos=fgetc(in); clearpos=(clearpos<<8)+fgetc(in);
  clearpos=(69631-clearpos)&0xffff;
  int codesize=fgetc(in),bits=codesize+1,shift=0,buf=0,blocksize=0;
  if (diffcount>4096 || clearpos<=(1<<codesize)+2) return 1;
  int maxcode=(1<<codesize)+1,dict[4096],input;
  for (int i=0; i<diffcount; i++) {
    diffpos[i]=fgetc(in);
    diffpos[i]=(diffpos[i]<<8)+fgetc(in);
    diffpos[i]=(diffpos[i]<<8)+fgetc(in);
    diffpos[i]=(diffpos[i]<<8)+fgetc(in);
    if (i>0) diffpos[i]+=diffpos[i-1];
  }
  U8 output[256];
  size-=6+diffcount*4;
  int last=fgetc(in),total=size+1,outsize=1;
  if (mode==FDECOMPRESS) fputc(codesize, out);
  else if (mode==FCOMPARE) if (codesize!=getc(out) && !diffFound) diffFound=1;
  if (diffcount==0 || diffpos[0]!=0) gif_write_code(1<<codesize) else curdiff++;
  while (size-->=0 && (input=fgetc(in))>=0) {
    int code=-1, key=(last<<8)+input;
    for (int i=(1<<codesize)+2; i<=min(maxcode,4095); i++) if (dict[i]==key) code=i;
    if (curdiff<diffcount && total-size>diffpos[curdiff]) curdiff++,code=-1;
    if (code==-1) {
      gif_write_code(last);
      if (maxcode==clearpos) { gif_write_code(1<<codesize); bits=codesize+1, maxcode=(1<<codesize)+1; }
      else
      {
        ++maxcode;
        if (maxcode<=4095) dict[maxcode]=key;
        if (maxcode>=(1<<bits) && bits<12) bits++;
      }
      code=input;
    }
    last=code;
  }
  gif_write_code(last);
  gif_write_code((1<<codesize)+1);
  if (shift>0) {
    output[++blocksize]=buf&255;
    if (blocksize==bsize) gif_write_block(bsize);
  }
  if (blocksize>0) gif_write_block(blocksize);
  if (mode==FDECOMPRESS) fputc(0, out);
  else if (mode==FCOMPARE) if (0!=getc(out) && !diffFound) diffFound=outsize+1;
  return outsize+1;
}

bool DetectGIF::detect(int c, int i, Filetype type, Filetype &nextBlockType) {
  buf1=buf1<<8|buf0>>24;
  buf0=buf0<<8|c;

  // Detect .gif
  if (type==DEFAULT && nextBlockType==GIF && i==0) {
    nextBlockType=DEFAULT;
    if (c==0x2c || c==0x21) gif=2,gifi=2;
  }
  if (!gif && (buf1&0xffff)==0x4749 && (buf0==0x46383961 || buf0==0x46383761)) gif=1,gifi=i+5;
  if (gif) {
    if (gif==1 && i==gifi) gif=2,gifi=i+5+((c&128)?(3*(2<<(c&7))):0);
    if (gif==2 && i==gifi) {
      if ((buf0&0xff0000)==0x210000) gif=5,gifi=i;
      else if ((buf0&0xff0000)==0x2c0000) gif=3,gifi=i;
      else gif=0;
    }
    if (gif==3 && i==gifi+6) gifw=(bswap(buf0)&0xffff);
    if (gif==3 && i==gifi+7) gif=4,gifc=gifb=0,gifa=gifi=i+2+((c&128)?(3*(2<<(c&7))):0);
    if (gif==4 && i==gifi) {
      if (c>0 && gifb && gifc!=gifb) gifw=0;
      if (c>0) gifb=gifc,gifc=c,gifi+=c+1;
      else if (!gifw) gif=2,gifi=i+3;
      else { nextBlockType=GIF; return f->detectBlock(GIF, gifa-1, i-gifa+2, (IMAGE8<<24)+gifw); }
    }
    if (gif==5 && i==gifi) {
      if (c>0) gifi+=c+1; else gif=2,gifi=i+3;
    }
  }
  return false;
}