#include "base64.h"

//
// decode/encode base64
//
static const char  table1[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
bool isbase64(unsigned char c) {
    return (isalnum(c) || (c == '+') || (c == '/')|| (c == 10) || (c == 13));
}


inline char valueb(char c){
       const char *p = strchr(table1, c);
       if(p) {
          return p-table1;
       } else {
          return 0;
       }
}

int encode_base64(FILE* in, FILE* out, int len, int info, int &hdrsize) {
  (void)(info);
  (void)(hdrsize);
  int in_len = 0;
  int i = 0;
  int j = 0;
  int b=0;
  int lfp=0;
  int tlf=0;
  char src[4];
  U8 *ptr,*fptr;
  int b64mem=(len>>2)*3+10;
    ptr = (U8*)calloc(b64mem, 1);
    if (!ptr) quit("Out of memory (e_B64)");
    fptr=&ptr[0];
    int olen=5;

  while (b=fgetc(in),in_len++ , ( b != '=') && isbase64(b) && in_len<=len) {
    if (b==13 || b==10) {
       if (lfp==0) lfp=in_len ,tlf=b;
       if (tlf!=b) tlf=0;
       continue;
    }
    src[i++] = b;
    if (i ==4){
          for (j = 0; j <4; j++) src[j] = valueb(src[j]);
          src[0] = (src[0] << 2) + ((src[1] & 0x30) >> 4);
          src[1] = ((src[1] & 0xf) << 4) + ((src[2] & 0x3c) >> 2);
          src[2] = ((src[2] & 0x3) << 6) + src[3];

          fptr[olen++]=src[0];
          fptr[olen++]=src[1];
          fptr[olen++]=src[2];
      i = 0;
    }
  }

  if (i){
    for (j=i;j<4;j++)
      src[j] = 0;

    for (j=0;j<4;j++)
      src[j] = valueb(src[j]);

    src[0] = (src[0] << 2) + ((src[1] & 0x30) >> 4);
    src[1] = ((src[1] & 0xf) << 4) + ((src[2] & 0x3c) >> 2);
    src[2] = ((src[2] & 0x3) << 6) + src[3];

    for (j=0;(j<i-1);j++) {
        fptr[olen++]=src[j];
    }
  }
  fptr[0]=lfp&255; //nl lenght
  fptr[1]=len&255;
  fptr[2]=len>>8&255;
  fptr[3]=len>>16&255;
  if (tlf!=0) {
    if (tlf==10) fptr[4]=128;
    else fptr[4]=64;
  }
  else
      fptr[4]=len>>24&63; //1100 0000
  fwrite(&ptr[0], 1, olen, out);
  free(ptr);
  return 1;
}

int decode_base64(Encoder& en, FILE* in, int size, int info, FILE *out, FMode mode, int &diffFound) {
    (void)(en);
    (void)(size);
    (void)(info);
    U8 inn[3];
    int i,len1=0, len=0, blocksout = 0;
    int fle=0;
    int linesize=0;
    int outlen=0;
    int tlf=0;
    linesize=getc(in);
    outlen=getc(in);
    outlen+=(getc(in)<<8);
    outlen+=(getc(in)<<16);
    tlf=(getc(in));
    outlen+=((tlf&63)<<24);
    U8 *ptr,*fptr;
    ptr = (U8*)calloc((outlen>>2)*4+10, 1);
    if (!ptr) quit("Out of memory (d_B64)");
    fptr=&ptr[0];
    tlf=(tlf&192);
    if (tlf==128)
       tlf=10;        // LF: 10
    else if (tlf==64)
         tlf=13;        // LF: 13
    else
         tlf=0;

    while(fle<outlen){
        len=0;
        for(i=0;i<3;i++){
            inn[i] = getc( in );
            if(!feof(in)){
                len++;
                len1++;
            }
            else {
                inn[i] = 0;
            }
        }
        if(len){
            U8 in0,in1,in2;
            in0=inn[0],in1=inn[1],in2=inn[2];
            fptr[fle++]=(table1[in0>>2]);
            fptr[fle++]=(table1[((in0&0x03)<<4)|((in1&0xf0)>>4)]);
            fptr[fle++]=((len>1?table1[((in1&0x0f)<<2)|((in2&0xc0)>>6)]:'='));
            fptr[fle++]=((len>2?table1[in2&0x3f]:'='));
            blocksout++;
        }
        if(blocksout>=(linesize/4) && linesize!=0){ //no lf if linesize==0
            if( blocksout &&  !feof(in) && fle<=outlen) { //no lf if eof
                if (tlf) fptr[fle++]=(tlf);
                else fptr[fle++]=13,fptr[fle++]=10;
            }
            blocksout = 0;
        }
    }
    //Write out or compare
    if (mode==FDECOMPRESS){
            fwrite(&ptr[0], 1, outlen, out);

        }
    else if (mode==FCOMPARE){
    for(i=0;i<outlen;i++){
        U8 b=fptr[i];


            if (b!=fgetc(out) && !diffFound) diffFound=ftell(out);
        }
    }
    free(ptr);
    return outlen;
}

bool DetectBASE64::detect(int c, int i) {
  buf1=buf1<<8|buf0>>24;
  buf0=buf0<<8|c;

  // Detect base64 encoded data
  if (b64s==0 && buf0==0x73653634 && ((buf1&0xffffff)==0x206261 || (buf1&0xffffff)==0x204261)) b64s=1,b64i=i-6; //' base64' ' Base64'
  if (b64s==0 && ((buf1==0x3b626173 && buf0==0x6536342c) || (buf1==0x215b4344 && buf0==0x4154415b))) b64s=3,b64i=i+1; // ';base64,' '![CDATA['
  if (b64s>0) {
    if (b64s==1 && buf0==0x0d0a0d0a) b64s=((i-b64i>=128)?0:2),b64i=i+1,b64line=0;
    else if (b64s==2 && (buf0&0xffff)==0x0d0a && b64line==0) b64line=i+1-b64i,b64nl=i;
    else if (b64s==2 && (buf0&0xffff)==0x0d0a && b64line>0 && (buf0&0xffffff)!=0x3d0d0a) {
       if (i-b64nl<b64line && buf0!=0x0d0a0d0a) i-=1,b64s=5;
       else if (buf0==0x0d0a0d0a) i-=3,b64s=5;
       else if (i-b64nl==b64line) b64nl=i;
       else b64s=0;
    }
    else if (b64s==2 && (buf0&0xffffff)==0x3d0d0a) i-=1,b64s=5; // '=' or '=='
    else if (b64s==2 && !(isalnum(c) || c=='+' || c=='/' || c==10 || c==13 || c=='=')) b64s=0;
    if (b64line>0 && (b64line<=4 || b64line>255)) b64s=0;
    if (b64s==3 && i>=b64i && !(isalnum(c) || c=='+' || c=='/' || c=='=')) b64s=4;
    if ((b64s==4 && i-b64i>128) || (b64s==5 && i-b64i>512 && i-b64i<(1<<27))) return f->detectBlock(BASE64, b64i, i-b64i);
    if (b64s>3) b64s=0;
  } 
  return false;
}