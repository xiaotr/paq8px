#include "Random.hpp"

Random::Random() : table(64) {
  table[0] = 123456789;
  table[1] = 987654321;
  for( uint64_t j = 0; j < 62; j++ ) {
    table[j + 2] = table[j + 1] * 11 + table[j] * 23 / 16;
  }
  i = 0;
}

//返回随机数
auto Random::operator()() -> uint32_t {
  return ++i, table[i & 63U] = table[(i - 24) & 63U] ^ table[(i - 55) & 63U];
}
