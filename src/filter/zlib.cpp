#include <zlib.h>
#include "zlib.h"

int parse_zlib_header(int header) {
    switch (header) {
        case 0x2815 : return 0;  case 0x2853 : return 1;  case 0x2891 : return 2;  case 0x28cf : return 3;
        case 0x3811 : return 4;  case 0x384f : return 5;  case 0x388d : return 6;  case 0x38cb : return 7;
        case 0x480d : return 8;  case 0x484b : return 9;  case 0x4889 : return 10; case 0x48c7 : return 11;
        case 0x5809 : return 12; case 0x5847 : return 13; case 0x5885 : return 14; case 0x58c3 : return 15;
        case 0x6805 : return 16; case 0x6843 : return 17; case 0x6881 : return 18; case 0x68de : return 19;
        case 0x7801 : return 20; case 0x785e : return 21; case 0x789c : return 22; case 0x78da : return 23;
    }
    return -1;
}

int zlib_inflateInit(z_streamp strm, int zh) {
    if (zh==-1) return inflateInit2(strm, -MAX_WBITS); else return inflateInit(strm);
}

int encode_zlib(FILE* in, FILE* out, int len, int info, int &hdrsize) {
  (void)(info);
  const int BLOCK=1<<16, LIMIT=128;
  U8 zin[BLOCK*2],zout[BLOCK],zrec[BLOCK*2], diffByte[81*LIMIT];
  int diffPos[81*LIMIT];

  // Step 1 - parse offset type form zlib stream header
  long pos=ftell(in);
  unsigned int h1=fgetc(in), h2=fgetc(in);
  fseek(in, pos, SEEK_SET);
  int zh=parse_zlib_header(h1*256+h2);
  int memlevel,clevel,window=zh==-1?0:MAX_WBITS+10+zh/4,ctype=zh%4;
  int minclevel=window==0?1:ctype==3?7:ctype==2?6:ctype==1?2:1;
  int maxclevel=window==0?9:ctype==3?9:ctype==2?6:ctype==1?5:1;

  // Step 2 - check recompressiblitiy, determine parameters and save differences
  z_stream main_strm, rec_strm[81];
  int diffCount[81], recpos[81], main_ret=Z_STREAM_END;
  main_strm.zalloc=Z_NULL; main_strm.zfree=Z_NULL; main_strm.opaque=Z_NULL;
  main_strm.next_in=Z_NULL; main_strm.avail_in=0;
  if (zlib_inflateInit(&main_strm,zh)!=Z_OK) return false;
  for (int i=0; i<81; i++) {
    memlevel=(i%9)+1;
    clevel=(i/9)+1;
    rec_strm[i].zalloc=Z_NULL; rec_strm[i].zfree=Z_NULL; rec_strm[i].opaque=Z_NULL;
    rec_strm[i].next_in=Z_NULL; rec_strm[i].avail_in=0;
    int ret=deflateInit2(&rec_strm[i], clevel, Z_DEFLATED, window-MAX_WBITS, memlevel, Z_DEFAULT_STRATEGY);
    diffCount[i]=(clevel>=minclevel && clevel<=maxclevel && ret==Z_OK)?0:LIMIT;
    recpos[i]=BLOCK*2;
    diffPos[i*LIMIT]=-1;
    diffByte[i*LIMIT]=0;
  }
  for (int i=0; i<len; i+=BLOCK) {
    unsigned int blsize=min(len-i,BLOCK);
    for (int j=0; j<81; j++) {
      if (diffCount[j]>=LIMIT) continue;
      memmove(&zrec[0], &zrec[BLOCK], BLOCK);
      recpos[j]-=BLOCK;
    }
    memmove(&zin[0], &zin[BLOCK], BLOCK);
    fread(&zin[BLOCK], 1, blsize, in); // Read block from input file

    // Decompress/inflate block
    main_strm.next_in=&zin[BLOCK]; main_strm.avail_in=blsize;
    do {
      main_strm.next_out=&zout[0]; main_strm.avail_out=BLOCK;
      main_ret=inflate(&main_strm, Z_FINISH);

      // Recompress/deflate block with all possible parameters
      for (int j=0; j<81; j++) {
        if (diffCount[j]>=LIMIT) continue;
        rec_strm[j].next_in=&zout[0];  rec_strm[j].avail_in=BLOCK-main_strm.avail_out;
        rec_strm[j].next_out=&zrec[recpos[j]]; rec_strm[j].avail_out=BLOCK*2-recpos[j];
        int ret=deflate(&rec_strm[j], (int)main_strm.total_in == len ? Z_FINISH : Z_NO_FLUSH);
        if (ret!=Z_BUF_ERROR && ret!=Z_STREAM_END && ret!=Z_OK) { diffCount[j]=LIMIT; continue; }

        // Compare
        int end=2*BLOCK-(int)rec_strm[j].avail_out;
        int tail=max(main_ret==Z_STREAM_END ? len-(int)rec_strm[j].total_out : 0,0);
        for (int k=recpos[j]; k<end+tail; k++) {
          if ((k<end && i+k-BLOCK<len && zrec[k]!=zin[k]) || k>=end) {
            if (++diffCount[j]<LIMIT) {
              const int p=j*LIMIT+diffCount[j];
              diffPos[p]=i+k-BLOCK;
              diffByte[p]=zin[k];
            }
          }
        }
        recpos[j]=2*BLOCK-rec_strm[j].avail_out;
      }
    } while (main_strm.avail_out==0 && main_ret==Z_BUF_ERROR);
    if (main_ret!=Z_BUF_ERROR && main_ret!=Z_STREAM_END) break;
  }
  int minCount=LIMIT, index;
  for (int i=80; i>=0; i--) {
    deflateEnd(&rec_strm[i]);
    if (diffCount[i]<minCount) {
      minCount=diffCount[i];
      memlevel=(i%9)+1;
      clevel=(i/9)+1;
      index=i;
    }
  }
  inflateEnd(&main_strm);
  if (minCount==LIMIT) return false;

  // Step 3 - write parameters, differences and precompressed (inflated) data
  fputc(diffCount[index], out);
  fputc(window, out);
  fputc(index, out);
  for (int i=0; i<=diffCount[index]; i++) {
    const int v=i==diffCount[index] ? len-diffPos[index*LIMIT+i]
                                    : diffPos[index*LIMIT+i+1]-diffPos[index*LIMIT+i]-1;
    fputc(v>>24, out); fputc(v>>16, out); fputc(v>>8, out); fputc(v, out);
  }
  for (int i=0; i<diffCount[index]; i++) fputc(diffByte[index*LIMIT+i+1], out);

  fseek(in, pos, SEEK_SET);
  main_strm.zalloc=Z_NULL; main_strm.zfree=Z_NULL; main_strm.opaque=Z_NULL;
  main_strm.next_in=Z_NULL; main_strm.avail_in=0;
  if (zlib_inflateInit(&main_strm,zh)!=Z_OK) return false;
  for (int i=0; i<len; i+=BLOCK) {
    unsigned int blsize=min(len-i,BLOCK);
    fread(&zin[0], 1, blsize, in);
    main_strm.next_in=&zin[0]; main_strm.avail_in=blsize;
    do {
      main_strm.next_out=&zout[0]; main_strm.avail_out=BLOCK;
      main_ret=inflate(&main_strm, Z_FINISH);
      fwrite(&zout[0], 1, BLOCK-main_strm.avail_out, out);
    } while (main_strm.avail_out==0 && main_ret==Z_BUF_ERROR);
    if (main_ret!=Z_BUF_ERROR && main_ret!=Z_STREAM_END) break;
  }
  hdrsize=diffCount[index]*5+7;
  return main_ret==Z_STREAM_END;
}

int decode_zlib(Encoder& en, FILE* in, int size, int info, FILE *out, FMode mode, int &diffFound) {
  (void)(en);
  (void)(info);
  const int BLOCK=1<<16, LIMIT=128;
  U8 zin[BLOCK],zout[BLOCK];
  int diffCount=min(fgetc(in),LIMIT-1);
  int window=fgetc(in)-MAX_WBITS;
  int index=fgetc(in);
  int memlevel=(index%9)+1;
  int clevel=(index/9)+1;
  int len=0;
  int diffPos[LIMIT];
  diffPos[0]=-1;
  for (int i=0; i<=diffCount; i++) {
    int v=fgetc(in)<<24; v|=fgetc(in)<<16; v|=fgetc(in)<<8; v|=fgetc(in);
    if (i==diffCount) len=v+diffPos[i]; else diffPos[i+1]=v+diffPos[i]+1;
  }
  U8 diffByte[LIMIT];
  diffByte[0]=0;
  for (int i=0; i<diffCount; i++) diffByte[i+1]=fgetc(in);
  size-=7+5*diffCount;

  z_stream rec_strm;
  int diffIndex=1,recpos=0;
  rec_strm.zalloc=Z_NULL; rec_strm.zfree=Z_NULL; rec_strm.opaque=Z_NULL;
  rec_strm.next_in=Z_NULL; rec_strm.avail_in=0;
  int ret=deflateInit2(&rec_strm, clevel, Z_DEFLATED, window, memlevel, Z_DEFAULT_STRATEGY);
  if (ret!=Z_OK) return 0;
  for (int i=0; i<size; i+=BLOCK) {
    int blsize=min(size-i,BLOCK);
    fread(&zin[0], 1, blsize, in);
    rec_strm.next_in=&zin[0];  rec_strm.avail_in=blsize;
    do {
      rec_strm.next_out=&zout[0]; rec_strm.avail_out=BLOCK;
      ret=deflate(&rec_strm, i+blsize==size ? Z_FINISH : Z_NO_FLUSH);
      if (ret!=Z_BUF_ERROR && ret!=Z_STREAM_END && ret!=Z_OK) break;
      const int have=min(BLOCK-rec_strm.avail_out,len-recpos);
      while (diffIndex<=diffCount && diffPos[diffIndex]>=recpos && diffPos[diffIndex]<recpos+have) {
        zout[diffPos[diffIndex]-recpos]=diffByte[diffIndex];
        diffIndex++;
      }
      if (mode==FDECOMPRESS) fwrite(&zout[0], 1, have, out);
      else if (mode==FCOMPARE) for (int j=0; j<have; j++) if (zout[j]!=getc(out) && !diffFound) diffFound=recpos+j+1;
      recpos+=have;

    } while (rec_strm.avail_out==0);
  }
  while (diffIndex<=diffCount) {
    if (mode==FDECOMPRESS) fputc(diffByte[diffIndex], out);
    else if (mode==FCOMPARE) if (diffByte[diffIndex]!=getc(out) && !diffFound) diffFound=recpos+1;
    diffIndex++;
    recpos++;
  }
  deflateEnd(&rec_strm);
  return recpos==len ? len : 0;
}


bool DetectZLIB::detect(int c, int i, long start, int n, FILE *in) {
  buf3=buf3<<8|buf2>>24;
  buf2=buf2<<8|buf1>>24;
  buf1=buf1<<8|buf0>>24;
  buf0=buf0<<8|c;

  zbuf[zbufpos]=c;
  zbufpos=(zbufpos+1)%32;
  int zh=parse_zlib_header(((int)zbuf[zbufpos])*256+(int)zbuf[(zbufpos+1)%32]);
  if ((i>=31 && zh!=-1) || zzippos==i) {
    int streamLength=0, ret=0;

    // Quick check possible stream by decompressing first 32 bytes
    z_stream strm;
    strm.zalloc=Z_NULL; strm.zfree=Z_NULL; strm.opaque=Z_NULL;
    strm.next_in=Z_NULL; strm.avail_in=0;
    if (zlib_inflateInit(&strm,zh)==Z_OK) {
      unsigned char tmp[32];
      for (int j=0; j<32; j++) tmp[j]=zbuf[(zbufpos+j)%32];
      strm.next_in=tmp; strm.avail_in=32;
      strm.next_out=zout; strm.avail_out=1<<16;
      ret=inflate(&strm, Z_FINISH);
      ret=(inflateEnd(&strm)==Z_OK && (ret==Z_STREAM_END || ret==Z_BUF_ERROR) && strm.total_in>=16);
    }
    if (ret) {
      // Verify valid stream and determine stream length
      long savedpos=ftell(in);
      strm.zalloc=Z_NULL; strm.zfree=Z_NULL; strm.opaque=Z_NULL;
      strm.next_in=Z_NULL; strm.avail_in=0; strm.total_in=strm.total_out=0;
      if (zlib_inflateInit(&strm,zh)==Z_OK) {
        for (int j=i-31; j<n; j+=1<<16) {
          unsigned int blsize=min(n-j,1<<16);
          fseek(in, start+j, SEEK_SET);
          if (fread(zin, 1, blsize, in)!=blsize) break;
          strm.next_in=zin; strm.avail_in=blsize;
          do {
            strm.next_out=zout; strm.avail_out=1<<16;
            ret=inflate(&strm, Z_FINISH);
          } while (strm.avail_out==0 && ret==Z_BUF_ERROR);
          if (ret==Z_STREAM_END) streamLength=strm.total_in;
          if (ret!=Z_BUF_ERROR) break;
        }
        if (inflateEnd(&strm)!=Z_OK) streamLength=0;
      }
      fseek(in, savedpos, SEEK_SET);
    }
    if (streamLength>0) {
      int info=0;
      if (pdfimw>0 && pdfimh>0) {
        if (pdfimb==8 && (int)strm.total_out==pdfimw*pdfimh) info=(IMAGE8<<24)+pdfimw;
        if (pdfimb==8 && (int)strm.total_out==pdfimw*pdfimh*3) info=(IMAGE24<<24)+pdfimw*3;
        if (pdfimb==4 && (int)strm.total_out==((pdfimw+1)/2)*pdfimh) info=(IMAGE8<<24)+((pdfimw+1)/2);
      }
      return f->detectBlock(ZLIB, i-31, streamLength, info);
    }
  }
  if (zh==-1 && zbuf[zbufpos]=='P' && zbuf[(zbufpos+1)%32]=='K' && zbuf[(zbufpos+2)%32]=='\x3'
    && zbuf[(zbufpos+3)%32]=='\x4' && zbuf[(zbufpos+8)%32]=='\x8' && zbuf[(zbufpos+9)%32]=='\0') {
      int nlen=(int)zbuf[(zbufpos+26)%32]+((int)zbuf[(zbufpos+27)%32])*256
              +(int)zbuf[(zbufpos+28)%32]+((int)zbuf[(zbufpos+29)%32])*256;
      if (nlen<256 && i+30+nlen<n) zzippos=i+30+nlen;
  }
  if (i-pdfimp>1024) pdfim=pdfimw=pdfimh=pdfimb=0;
  if (pdfim>1 && !(isspace(c) || isdigit(c))) pdfim=1;
  if (pdfim==2 && isdigit(c)) pdfimw=pdfimw*10+(c-'0');
  if (pdfim==3 && isdigit(c)) pdfimh=pdfimh*10+(c-'0');
  if (pdfim==4 && isdigit(c)) pdfimb=pdfimb*10+(c-'0');
  if ((buf0&0xffff)==0x3c3c) pdfimp=i,pdfim=1; // <<
  if (pdfim && (buf1&0xffff)==0x2f57 && buf0==0x69647468) pdfim=2,pdfimw=0; // /Width
  if (pdfim && (buf1&0xffffff)==0x2f4865 && buf0==0x69676874) pdfim=3,pdfimh=0; // /Height
  if (pdfim && buf3==0x42697473 && buf2==0x50657243 && buf1==0x6f6d706f
     && buf0==0x6e656e74 && zbuf[(zbufpos+15)%32]=='/') pdfim=4,pdfimb=0; // /BitsPerComponent
  return false;
}
