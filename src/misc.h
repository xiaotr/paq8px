#ifndef MISC_H_INCLUDED
#define MISC_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#else
#include <windows.h>
#endif


#if defined(__x86_64) || defined(_M_X64)
 #define X64
#else
 #undef X64
#endif

// 8, 16, 32 bit unsigned types (adjust as appropriate)
typedef unsigned char  U8;
typedef unsigned short U16;
typedef unsigned int   U32;

// min, max functions
inline int min(int a, int b) {return a<b?a:b;}
inline int max(int a, int b) {return a<b?b:a;}

// Error handler: print message if any, and exit
inline void quit(const char* message=0) {throw message;}

typedef enum {DEFAULT=0, JPEG, HDR, IMAGE1, IMAGE8, IMAGE8GRAY, IMAGE24, IMAGE32, AUDIO, EXE, CD, ZLIB, BASE64, GIF} Filetype;

inline bool hasRecursion(Filetype ft) { return ft==CD || ft==ZLIB || ft==BASE64 || ft==GIF; }
inline bool hasInfo(Filetype ft) { return ft==IMAGE1 || ft==IMAGE8 || ft==IMAGE8GRAY || ft==IMAGE24 || ft==IMAGE32 || ft==AUDIO; }

#define MEM (0x10000<<level)

#define bswap(x) \
+   ((((x) & 0xff000000) >> 24) | \
+    (((x) & 0x00ff0000) >>  8) | \
+    (((x) & 0x0000ff00) <<  8) | \
+    (((x) & 0x000000ff) << 24))

//////////////////////// Program Checker /////////////////////

// Track time and memory used
class ProgramChecker {
  int memused;  // bytes allocated by Array<T> now
  int maxmem;   // most bytes allocated ever
  clock_t start_time;  // in ticks
public:
  void alloc(int n) {  // report memory allocated, may be negative
    memused+=n;
    if (memused>maxmem) maxmem=memused;
  }
  ProgramChecker(): memused(0), maxmem(0) {
    start_time=clock();
    assert(sizeof(U8)==1);
    assert(sizeof(U16)==2);
    assert(sizeof(U32)==4);
    assert(sizeof(short)==2);
    assert(sizeof(int)==4);
  }
  void print() const {  // print time and memory used
    printf("Time %1.2f sec, used %d bytes of memory\n",
      double(clock()-start_time)/CLOCKS_PER_SEC, maxmem);
  }
};

extern ProgramChecker programChecker;


//////////////////////////// Array ////////////////////////////

// Array<T> a(n); creates n elements of T initialized to 0 bits.
// Constructors for T are not called.
// Indexing is bounds checked if assertions are on.
// a.size() returns n.
// a.resize(n) changes size to n, padding with 0 bits or truncating.
// a.push_back(x) appends x and increases size by 1, reserving up to size*2.
// a.pop_back() decreases size by 1, does not free memory.
// Copy and assignment are not supported.

template <class T> class Array {
private:
  int n;     // user size
  int reserved;  // actual size
  char *ptr; // allocated memory, zeroed
  T* data;   // start of n elements of aligned data
  void create(int i);  // create with size i
public:
  explicit Array(int i=0) {create(i);}
  ~Array();
  T& operator[](int i) {
#ifndef NDEBUG
    if (i<0 || i>=n) fprintf(stderr, "%d out of bounds %d\n", i, n), throw 0;
#endif
    return data[i];
  }
  const T& operator[](int i) const {
#ifndef NDEBUG
    if (i<0 || i>=n) fprintf(stderr, "%d out of bounds %d\n", i, n), throw 0;
#endif
    return data[i];
  }
  int size() const {return n;}
  void resize(int i);  // change size to i
  void pop_back() {if (n>0) --n;}  // decrement size
  void push_back(const T& x);  // increment size, append x
private:
  Array(const Array&);  // no copy or assignment
  Array& operator=(const Array&);
};


#ifndef UNIX
typedef unsigned int uint;
typedef unsigned long long qword;

const uint g_PageFlag0 = 0x1000;
const uint g_PageMask0 = 0x1000-1;
extern uint g_PageFlag;
extern uint g_PageMask;

template< class T >
T* VAlloc( qword size ) {
  void* r;
  size *= sizeof(T);
  qword s = (size+g_PageMask) & (~qword(g_PageMask));
  Retry:
//printf( "s=%I64X sizeof(size_t)=%i\n", s, sizeof(size_t) );
#ifndef X64
  if( s>=(qword(1)<<32) ) r=0; else
#endif
  r = VirtualAlloc(0, s, g_PageFlag, PAGE_READWRITE );
//printf( "r=%08X\n", r );
  if( (r==0) && (g_PageMask!=g_PageMask0) ) {
    g_PageFlag = g_PageFlag0;
    g_PageMask = g_PageMask0; 
    s = size;
    goto Retry;
  }
  return (T*)r;
}

template< class T >
void VFree( T* p ) {
  VirtualFree( p, 0, MEM_RELEASE );
}
#endif


template<class T> void Array<T>::resize(int i) {
  if (i<=reserved) {
    n=i;
    return;
  }
  char *saveptr=ptr;
  T *savedata=data;
  int saven=n;
  create(i);
  if (saveptr) {
    if (savedata) {
      memcpy(data, savedata, sizeof(T)*min(i, saven));
      programChecker.alloc(-n*sizeof(T));
    }
    free(saveptr);
  }
}

template<class T> void Array<T>::create(int i) {
  n=reserved=i;
  if (i<=0) {
    data=0;
    ptr=0;
    return;
  }
  const int sz=n*sizeof(T);
  programChecker.alloc(sz);
#ifndef UNIX  
  char* r;

  if( sz>(1<<24) ) {
    r = ptr = VAlloc<char>( sz );
  } else {
    int flag = (sz>=16);
    ptr = (char*)calloc( sz+(flag?15:0), 1 );
    if( !ptr ) throw "Out of memory";
    r = ptr;
    if( flag ) { r = ptr + 15; r -= (r-((char*)0))&15; }
  }

  data = (T*)r; //ptr;
#else
  ptr = (char*)calloc(sz, 1);
  if (!ptr) throw "Out of memory";
  data = (T*)ptr;
#endif
}

template<class T> Array<T>::~Array() {
  programChecker.alloc(-n*sizeof(T));
#ifndef UNIX   
  const int sz=n*sizeof(T);

  if( sz>(1<<24) ) {
    VFree(ptr);
  } else {
    free(ptr);
  }
#else
  free(ptr);
#endif
}

template<class T> void Array<T>::push_back(const T& x) {
  if (n==reserved) {
    int saven=n;
    resize(max(1, n*2));
    n=saven;
  }
  data[n++]=x;
}


/////////////////////////// String /////////////////////////////

// A tiny subset of std::string
// size() includes NUL terminator.

class String: public Array<char> {
public:
  const char* c_str() const {return &(*this)[0];}
  void operator=(const char* s) {
    resize(strlen(s)+1);
    strcpy(&(*this)[0], s);
  }
  void operator+=(const char* s) {
    assert(s);
    pop_back();
    while (*s) push_back(*s++);
    push_back(0);
  }
  String(const char* s=""): Array<char>(1) {
    (*this)+=s;
  }
};



#endif // #ifndef MISC_H_INCLUDED