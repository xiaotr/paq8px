#ifndef PAQ8PX_BITCOUNT
#define PAQ8PX_BITCOUNT

#include <cstdint>

//算出v的二进制1的个数
static auto bitCount(uint32_t v) -> uint32_t {
  v -= ((v >> 1U) & 0x55555555U);
  v = ((v >> 2U) & 0x33333333U) + (v & 0x33333333U);
  v = ((v >> 4U) + v) & 0x0f0f0f0fU;
  v = ((v >> 8U) + v) & 0x00ff00ffU;
  v = ((v >> 16U) + v) & 0x0000ffffU;
  return v;
}

#endif