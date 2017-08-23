#include "audio.h"

bool DetectAUDIO::detect(int c, int i, long start, FILE *in) {
  buf1=buf1<<8|buf0>>24;
  buf0=buf0<<8|c;

  // Detect .wav file header
  if (buf0==0x52494646) wavi=i,wavm=0;
  if (wavi) {
    int p=i-wavi;
    if (p==4) wavsize=bswap(buf0);
    else if (p==8){
      wavtype=(buf0==0x57415645)?1:(buf0==0x7366626B)?2:0;
      if (!wavtype) wavi=0;
    }
    else if (wavtype){
      if (wavtype==1) {
        if (p==16 && (buf1!=0x666d7420 || bswap(buf0)!=16)) wavi=0;
        else if (p==22) wavch=bswap(buf0)&0xffff;
        else if (p==34) wavbps=bswap(buf0)&0xffff;
        else if (p==40+wavm && buf1!=0x64617461) wavm+=bswap(buf0)+8,wavi=(wavm>0xfffff?0:wavi);
        else if (p==40+wavm) {
          int wavd=bswap(buf0);
          if ((wavch==1 || wavch==2) && (wavbps==8 || wavbps==16) && wavd>0 && wavsize>=wavd+36
             && wavd%((wavbps/8)*wavch)==0) f->detectBlock(AUDIO, wavi-3, wavd, wavch+wavbps/4-3, 44+wavm);
          wavi=0;
        }
      }
      else{
        if ((p==16 && buf1!=0x4C495354) || (p==20 && buf0!=0x494E464F))
          wavi=0;
        else if (p>20 && buf1==0x4C495354 && (wavi*=(buf0!=0))){
          wavlen = bswap(buf0);
          wavlist = i;
        }
        else if (wavlist){
          p=i-wavlist;
          if (p==8 && (buf1!=0x73647461 || buf0!=0x736D706C))
            wavi=0;
          else if (p==12){
            int wavd = bswap(buf0);
            if (wavd && (wavd+12)==wavlen)
              f->detectBlock(AUDIO, wavi-3, wavd, 1+16/4-3, (12+wavlist-(wavi-3)+1)&~1);
            wavi=0;
          }
        }
      }
    }
  }

  // Detect .aiff file header
  if (buf0==0x464f524d) aiff=i,aiffs=0; // FORM
  if (aiff) {
    const int p=i-aiff;
    if (p==12 && (buf1!=0x41494646 || buf0!=0x434f4d4d)) aiff=0; // AIFF COMM
    else if (p==24) {
      const int bits=buf0&0xffff, chn=buf1>>16;
      if ((bits==8 || bits==16) && (chn==1 || chn==2)) aiffm=chn+bits/4+1; else aiff=0;
    } else if (p==42+aiffs && buf1!=0x53534e44) aiffs+=(buf0+8)+(buf0&1),aiff=(aiffs>0x400?0:aiff);
    else if (p==42+aiffs) f->detectBlock(AUDIO, aiff-3, buf0-8, aiffm, 54+aiffs);
  }

  // Detect .mod file header
  if ((buf0==0x4d2e4b2e || buf0==0x3643484e || buf0==0x3843484e  // M.K. 6CHN 8CHN
     || buf0==0x464c5434 || buf0==0x464c5438) && (buf1&0xc0c0c0c0)==0 && i>=1083) {
    long savedpos=ftell(in);
    const int chn=((buf0>>24)==0x36?6:(((buf0>>24)==0x38 || (buf0&0xff)==0x38)?8:4));
    int len=0; // total length of samples
    int numpat=1; // number of patterns
    for (int j=0; j<31; j++) {
      fseek(in, start+i-1083+42+j*30, SEEK_SET);
      const int i1=getc(in);
      const int i2=getc(in);
      len+=i1*512+i2*2;
    }
    fseek(in, start+i-131, SEEK_SET);
    for (int j=0; j<128; j++) {
      int x=getc(in);
      if (x+1>numpat) numpat=x+1;
    }
    if (numpat<65) f->detectBlock(AUDIO, i-1083, len, 4, 1084+numpat*256*chn);
    fseek(in, savedpos, SEEK_SET);
  }

  // Detect .s3m file header
  if (buf0==0x1a100000) s3mi=i,s3mno=s3mni=0;
  if (s3mi) {
    const int p=i-s3mi;
    if (p==4) s3mno=bswap(buf0)&0xffff,s3mni=(bswap(buf0)>>16);
    else if (p==16 && (((buf1>>16)&0xff)!=0x13 || buf0!=0x5343524d)) s3mi=0;
    else if (p==16) {
      long savedpos=ftell(in);
      int b[31],sam_start=(1<<16),sam_end=0,ok=1;
      for (int j=0;j<s3mni;j++) {
        fseek(in, start+s3mi-31+0x60+s3mno+j*2, SEEK_SET);
        int i1=getc(in);
        i1+=getc(in)*256;
        fseek(in, start+s3mi-31+i1*16, SEEK_SET);
        i1=getc(in);
        if (i1==1) { // type: sample
          for (int k=0;k<31;k++) b[k]=fgetc(in);
          int len=b[15]+(b[16]<<8);
          int ofs=b[13]+(b[14]<<8);
          if (b[30]>1) ok=0;
          if (ofs*16<sam_start) sam_start=ofs*16;
          if (ofs*16+len>sam_end) sam_end=ofs*16+len;
        }
      }
      if (ok && sam_start<(1<<16)) f->detectBlock(AUDIO, s3mi-31, sam_end-sam_start, 0, sam_start);
      s3mi=0;
      fseek(in, savedpos, SEEK_SET);
    }
  }
  return false;
}