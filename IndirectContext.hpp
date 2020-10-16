#ifndef PAQ8PX_INDIRECTCONTEXT_HPP
#define PAQ8PX_INDIRECTCONTEXT_HPP

#include "Array.hpp"
#include <cassert>
#include <cstdint>

//间接上下文
template<typename T>
class IndirectContext {
private:
    Array<T> data;
    T *ctx; //ctx是指针
    const uint32_t ctxMask, inputMask, inputBits;

public:
    IndirectContext(const int bitsPerContext, const int inputBits) :
      data(UINT64_C(1) << bitsPerContext), ctx(&data[0]),
      ctxMask((UINT32_C(1) << bitsPerContext) - 1), 
      inputMask((UINT32_C(1) << inputBits) - 1), 
      inputBits(inputBits) {
#ifdef VERBOSE
      printf("Created IndirectContext with bitsPerContext = %d, inputBits = %d\n", bitsPerContext, inputBits);
#endif
      assert(bitsPerContext > 0 && bitsPerContext <= 20);
      assert(inputBits > 0 && inputBits <= 8);
    }

    //将i放到低位
    void operator+=(const uint32_t i) {
      assert(i <= inputMask);
      (*ctx) <<= inputBits;
      (*ctx) |= i;
    };

    //ctxMask是掩码
    void operator=(const uint32_t i) {
      ctx = &data[i & ctxMask];
    }
    
    //()返回cxt指向的值
    auto operator()() -> T & {
      return *ctx;
    };
};

#endif //PAQ8PX_INDIRECTCONTEXT_HPP
