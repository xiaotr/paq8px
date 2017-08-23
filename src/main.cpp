/*
  paq8px file compressor/archiver. Released on August 23, 2017

  Copyright (C) 2008 Matt Mahoney, Serge Osnach, Alexander Ratushnyak,
  Bill Pettis, Przemyslaw Skibinski, Matthew Fite, wowtiger, Andrew Paterson,
  Jan Ondrus, Andreas Morphis, Pavel L. Holoborodko, Kaido Orav, Simon Berger,
  Neill Corlett, Márcio Pais

  LICENSE

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details at
  Visit <http://www.gnu.org/copyleft/gpl.html>.
*/

#define PROGNAME "paq8px"  // Please change this if you change the program.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <zlib.h>
#define NDEBUG  // remove for debugging (turns on Array bound checks)
#include <assert.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#else
#include <windows.h>
#endif

#include "misc.h"
#include "encoder.h"
#include "filter.h"

#ifndef DEFAULT_OPTION
#define DEFAULT_OPTION 5
#endif

// Compress a file. Split filesize bytes into blocks by type.
// For each block, output
// <type> <size> and call encode_X to convert to type X.
// Test transform and compress.
void compress(const char* filename, long filesize, Encoder& en) {
  assert(en.getMode()==COMPRESS);
  assert(filename && filename[0]);
  FILE *in=fopen(filename, "rb");
  if (!in) perror(filename), quit();
  long start=en.size();
  printf("Block segmentation:\n");
  char blstr[32]="";
  Filter filter;  
  filter.compressRecursive(in, filesize, en, blstr);
  if (in) fclose(in);
  printf("Compressed from %ld to %ld bytes.\n",filesize,en.size()-start);
}

// Try to make a directory, return true if successful
bool makedir(const char* dir) {
#ifdef WINDOWS
  return CreateDirectory(dir, 0)==TRUE;
#else
#ifdef UNIX
  return mkdir(dir, 0777)==0;
#else
  return false;
#endif
#endif
}

// Decompress a file
void decompress(const char* filename, long filesize, Encoder& en) {
  FMode mode=FDECOMPRESS;
  assert(en.getMode()==DECOMPRESS);
  assert(filename && filename[0]);

  // Test if output file exists.  If so, then compare.
  FILE* f=fopen(filename, "rb");
  if (f) mode=FCOMPARE,printf("Comparing");
  else {
    // Create file
    f=fopen(filename, "wb");
    if (!f) {  // Try creating directories in path and try again
      String path(filename);
      for (int i=0; path[i]; ++i) {
        if (path[i]=='/' || path[i]=='\\') {
          char savechar=path[i];
          path[i]=0;
          if (makedir(path.c_str()))
            printf("Created directory %s\n", path.c_str());
          path[i]=savechar;
        }
      }
      f=fopen(filename, "wb");
    }
    if (!f) mode=FDISCARD,printf("Skipping"); else printf("Extracting");
  }
  printf(" %s %ld -> ", filename, filesize);

  // Decompress/Compare
  Filter filter;
  int r=filter.decompressRecursive(f, filesize, en, mode);
  if (mode==FCOMPARE && !r && getc(f)!=EOF) printf("file is longer\n");
  else if (mode==FCOMPARE && r) printf("differ at %d\n",r-1);
  else if (mode==FCOMPARE) printf("identical\n");
  else printf("done   \n");
  if (f) fclose(f);
}

// strings are equal ignoring case?
int equals(const char* a, const char* b) {
  assert(a && b);
  while (*a && *b) {
    int c1=*a;
    if (c1>='A'&&c1<='Z') c1+='a'-'A';
    int c2=*b;
    if (c2>='A'&&c2<='Z') c2+='a'-'A';
    if (c1!=c2) return 0;
    ++a;
    ++b;
  }
  return *a==*b;
}

//////////////////////////// User Interface ////////////////////////////


// int expand(String& archive, String& s, const char* fname, int base) {
// Given file name fname, print its length and base name (beginning
// at fname+base) to archive in format "%ld\t%s\r\n" and append the
// full name (including path) to String s in format "%s\n".  If fname
// is a directory then substitute all of its regular files and recursively
// expand any subdirectories.  Base initially points to the first
// character after the last / in fname, but in subdirectories includes
// the path from the topmost directory.  Return the number of files
// whose names are appended to s and archive.

// Same as expand() except fname is an ordinary file
int putsize(String& archive, String& s, const char* fname, int base) {
  int result=0;
  FILE *f=fopen(fname, "rb");
  if (f) {
    fseek(f, 0, SEEK_END);
    long len=ftell(f);
    if (len>=0) {
      static char blk[24];
      sprintf(blk, "%ld\t", len);
      archive+=blk;
      archive+=(fname+base);
      archive+="\n";
      s+=fname;
      s+="\n";
      ++result;
    }
    fclose(f);
  }
  return result;
}

#ifdef WINDOWS

int expand(String& archive, String& s, const char* fname, int base) {
  int result=0;
  DWORD attr=GetFileAttributes(fname);
  if ((attr != 0xFFFFFFFF) && (attr & FILE_ATTRIBUTE_DIRECTORY)) {
    WIN32_FIND_DATA ffd;
    String fdir(fname);
    fdir+="/*";
    HANDLE h=FindFirstFile(fdir.c_str(), &ffd);
    while (h!=INVALID_HANDLE_VALUE) {
      if (!equals(ffd.cFileName, ".") && !equals(ffd.cFileName, "..")) {
        String d(fname);
        d+="/";
        d+=ffd.cFileName;
        result+=expand(archive, s, d.c_str(), base);
      }
      if (FindNextFile(h, &ffd)!=TRUE) break;
    }
    FindClose(h);
  }
  else // ordinary file
    result=putsize(archive, s, fname, base);
  return result;
}

#else
#ifdef UNIX

int expand(String& archive, String& s, const char* fname, int base) {
  int result=0;
  struct stat sb;
  if (stat(fname, &sb)<0) return 0;

  // If a regular file and readable, get file size
  if (sb.st_mode & S_IFREG && sb.st_mode & 0400)
    result+=putsize(archive, s, fname, base);

  // If a directory with read and execute permission, traverse it
  else if (sb.st_mode & S_IFDIR && sb.st_mode & 0400 && sb.st_mode & 0100) {
    DIR *dirp=opendir(fname);
    if (!dirp) {
      perror("opendir");
      return result;
    }
    dirent *dp;
    while(errno=0, (dp=readdir(dirp))!=0) {
      if (!equals(dp->d_name, ".") && !equals(dp->d_name, "..")) {
        String d(fname);
        d+="/";
        d+=dp->d_name;
        result+=expand(archive, s, d.c_str(), base);
      }
    }
    if (errno) perror("readdir");
    closedir(dirp);
  }
  else printf("%s is not a readable file or directory\n", fname);
  return result;
}

#else  // Not WINDOWS or UNIX, ignore directories

int expand(String& archive, String& s, const char* fname, int base) {
  return putsize(archive, s, fname, base);
}

#endif
#endif


// To compress to file1.paq8px: paq8px [-n] file1 [file2...]
// To decompress: paq8px file1.paq8px [output_dir]
int main(int argc, char** argv) {
  bool pause=argc<=2;  // Pause when done?
  try {

    // Get option
    bool doExtract=false;  // -d option
    bool doList=false;  // -l option
    int level=DEFAULT_OPTION;
    if (argc>1 && argv[1][0]=='-' && argv[1][1] && !argv[1][2]) {
      if (argv[1][1]>='0' && argv[1][1]<='8')
        level=argv[1][1]-'0';
      else if (argv[1][1]=='d')
        doExtract=true;
      else if (argv[1][1]=='l')
        doList=true;
      else
        quit("Valid options are -0 through -8, -d, -l\n");
      --argc;
      ++argv;
      pause=false;
    }

    // Print help message
    if (argc<2) {
      printf(PROGNAME " archiver (C) 2016, Matt Mahoney et al.\n"
        "Free under GPL, http://www.gnu.org/licenses/gpl.txt\n\n"
#ifdef WINDOWS
        "To compress or extract, drop a file or folder on the "
        PROGNAME " icon.\n"
        "The output will be put in the same folder as the input.\n"
        "\n"
        "Or from a command window: "
#endif
        "To compress:\n"
        "  " PROGNAME " -level file               (compresses to file." PROGNAME ")\n"
        "  " PROGNAME " -level archive files...   (creates archive." PROGNAME ")\n"
        "  " PROGNAME " file                      (level -%d, pause when done)\n"
        "level: -0 = store, -1 -2 -3 = faster (uses 35, 48, 59 MB)\n"
        "-4 -5 -6 -7 -8 = smaller (uses 133, 233, 435, 837, 1643 MB)\n"
#if defined(WINDOWS) || defined (UNIX)
        "You may also compress directories.\n"
#endif
        "\n"
        "To extract or compare:\n"
        "  " PROGNAME " -d dir1/archive." PROGNAME "      (extract to dir1)\n"
        "  " PROGNAME " -d dir1/archive." PROGNAME " dir2 (extract to dir2)\n"
        "  " PROGNAME " archive." PROGNAME "              (extract, pause when done)\n"
        "\n"
        "To view contents: " PROGNAME " -l archive." PROGNAME "\n"
        "\n",
        DEFAULT_OPTION);
      quit();
    }

    FILE* archive=0;  // compressed file
    int files=0;  // number of files to compress/decompress
    Array<const char*> fname(1);  // file names (resized to files)
    Array<long> fsize(1);   // file lengths (resized to files)

    // Compress or decompress?  Get archive name
    Mode mode=COMPRESS;
    String archiveName(argv[1]);
    {
      const int prognamesize=strlen(PROGNAME);
      const int arg1size=strlen(argv[1]);
      if (arg1size>prognamesize+1 && argv[1][arg1size-prognamesize-1]=='.'
          && equals(PROGNAME, argv[1]+arg1size-prognamesize)) {
        mode=DECOMPRESS;
      }
      else if (doExtract || doList)
        mode=DECOMPRESS;
      else {
        archiveName+=".";
        archiveName+=PROGNAME;
      }
    }

    // Compress: write archive header, get file names and sizes
    String header_string;
    String filenames;
    if (mode==COMPRESS) {

      // Expand filenames to read later.  Write their base names and sizes
      // to archive.
      int i;
      for (i=1; i<argc; ++i) {
        String name(argv[i]);
        int len=name.size()-1;
        for (int j=0; j<=len; ++j)  // change \ to /
          if (name[j]=='\\') name[j]='/';
        while (len>0 && name[len-1]=='/')  // remove trailing /
          name[--len]=0;
        int base=len-1;
        while (base>=0 && name[base]!='/') --base;  // find last /
        ++base;
        if (base==0 && len>=2 && name[1]==':') base=2;  // chop "C:"
        int expanded=expand(header_string, filenames, name.c_str(), base);
        if (!expanded && (i>1||argc==2))
          printf("%s: not found, skipping...\n", name.c_str());
        files+=expanded;
      }

      // If there is at least one file to compress
      // then create the archive header.
      if (files<1) quit("Nothing to compress\n");
      archive=fopen(archiveName.c_str(), "wb+");
      if (!archive) perror(archiveName.c_str()), quit();
      fprintf(archive, PROGNAME "%c%d", 0, level);
      printf("Creating archive %s with %d file(s)...\n",
        archiveName.c_str(), files);
    }

    // Decompress: open archive for reading and store file names and sizes
    if (mode==DECOMPRESS) {
      archive=fopen(archiveName.c_str(), "rb+");
      if (!archive) perror(archiveName.c_str()), quit();

      // Check for proper format and get option
      String header;
      int len=strlen(PROGNAME)+2, c, i=0;
      header.resize(len+1);
      while (i<len && (c=getc(archive))!=EOF) {
        header[i]=c;
        i++;
      }
      header[i]=0;
      if (strncmp(header.c_str(), PROGNAME "\0", strlen(PROGNAME)+1))
        printf("%s: not a %s file\n", archiveName.c_str(), PROGNAME), quit();
      level=header[strlen(PROGNAME)+1]-'0';
      if (level<0||level>8) level=DEFAULT_OPTION;
    }

    // Set globals according to option
    assert(level>=0 && level<=8);
    Encoder en(mode, archive, level);

    // Compress header
    if (mode==COMPRESS) {
      int len=header_string.size();
      printf("\nFile list (%d bytes)\n", len);
      assert(en.getMode()==COMPRESS);
      long start=en.size();
      en.compress(0); // block type 0
      en.compress(len>>24); en.compress(len>>16); en.compress(len>>8); en.compress(len); // block length
      for (int i=0; i<len; i++) en.compress(header_string[i]);
      printf("Compressed from %d to %ld bytes.\n",len,en.size()-start);
    }

    // Deompress header
    if (mode==DECOMPRESS) {
      if (en.decompress()!=0) printf("%s: header corrupted\n", archiveName.c_str()), quit();
      int len=0;
      len+=en.decompress()<<24;
      len+=en.decompress()<<16;
      len+=en.decompress()<<8;
      len+=en.decompress();
      header_string.resize(len);
      for (int i=0; i<len; i++) {
        header_string[i]=en.decompress();
        if (header_string[i]=='\n') files++;
      }
      if (doList) printf("File list of %s archive:\n%s", archiveName.c_str(), header_string.c_str());
    }

    // Fill fname[files], fsize[files] with input filenames and sizes
    fname.resize(files);
    fsize.resize(files);
    char *p=&header_string[0];
    char* q=&filenames[0];
    for (int i=0; i<files; ++i) {
      assert(p);
      fsize[i]=atol(p);
      assert(fsize[i]>=0);
      while (*p!='\t') ++p; *(p++)='\0';
      fname[i]=mode==COMPRESS?q:p;
      while (*p!='\n') ++p; *(p++)='\0';
      if (mode==COMPRESS) { while (*q!='\n') ++q; *(q++)='\0'; }
    }

    // Compress or decompress files
    assert(fname.size()==files);
    assert(fsize.size()==files);
    long total_size=0;  // sum of file sizes
    for (int i=0; i<files; ++i) total_size+=fsize[i];
    if (mode==COMPRESS) {
      for (int i=0; i<files; ++i) {
        printf("\n%d/%d  Filename: %s (%ld bytes)\n", i+1, files, fname[i], fsize[i]);
        compress(fname[i], fsize[i], en);
      }
      en.flush();
      printf("\nTotal %ld bytes compressed to %ld bytes.\n", total_size, en.size());
    }

    // Decompress files to dir2: paq8px -d dir1/archive.paq8px dir2
    // If there is no dir2, then extract to dir1
    // If there is no dir1, then extract to .
    else if (!doList) {
      assert(argc>=2);
      String dir(argc>2?argv[2]:argv[1]);
      if (argc==2) {  // chop "/archive.paq8px"
        int i;
        for (i=dir.size()-2; i>=0; --i) {
          if (dir[i]=='/' || dir[i]=='\\') {
            dir[i]=0;
            break;
          }
          if (i==1 && dir[i]==':') {  // leave "C:"
            dir[i+1]=0;
            break;
          }
        }
        if (i==-1) dir=".";  // "/" not found
      }
      dir=dir.c_str();
      if (dir[0] && (dir.size()!=3 || dir[1]!=':')) dir+="/";
      for (int i=0; i<files; ++i) {
        String out(dir.c_str());
        out+=fname[i];
        decompress(out.c_str(), fsize[i], en);
      }
    }
    fclose(archive);
    if (!doList) programChecker.print();
  }
  catch(const char* s) {
    if (s) printf("%s\n", s);
  }
  if (pause) {
    printf("\nClose this window or press ENTER to continue...\n");
    getchar();
  }
  return 0;
}
