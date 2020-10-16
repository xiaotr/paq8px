#ifndef PAQ8PX_HASHTABLE_HPP
#define PAQ8PX_HASHTABLE_HPP

#include "Hash.hpp"
#include "Ilog.hpp"
#include <cassert>
#include <cstdint>
#include <cstring>

//n个条目代表n个上下文
//在每个槽中，第一个字节是使用上下文选择器的上8位的校验和。
//第二个字节是哈希替换的优先级(0 =空)。
//优先级必须由调用者设定(0:最低，255:最高)。
//优先级较低的项目将被替换，以防碰撞。
//对于给定的校验和，只探测另外两个项。
//调用者可以以字节[2..B-1]存储关于给定上下文的任何信息。
//调用者不能修改校验和。
//索引必须是一个乘法散列。

/**
 * A HashTable is an array of n items representing n contexts. Each item is a storage area of @ref B bytes.
 * In every slot the first byte is a checksum using the upper 8 bits of the context selector.
 * The second byte is a priority (0 = empty) for hash replacement.
 * Priorities must be set by the caller (0: lowest, 255: highest).
 * Items with lower priorities will be replaced in case of a collision.
 * Only 2 additional items are probed for the given checksum.
 * The caller can store any information about the given context in bytes [2..B-1].
 * Checksums must not be modified by the caller.
 * The index must be a multiplicative hash.
 * @tparam B the number of bytes per item
 */
template<int B>
class HashTable {
private:
    Array<uint8_t, 64> t; /**< storage area for n items (1 item = b bytes): 0:checksum 1:priority 2:data 3:data  ... b-1:data */
    const int mask; //
    const int hashBits;

public:
    //n是存储区域个数
    /**
     * Creates a hashtable with @ref n slots where n and B must be powers of 2 with n >= B*4, and B >= 2.
     * @param n the number of storage areas
     */
    explicit HashTable(uint64_t n) : t(n), mask(static_cast<int>(n) - 1), hashBits(ilog2(mask + 1)) {
      assert(B >= 2 && isPowerOf2(B));
      assert(n >= B * 4 && isPowerOf2(n));
      assert(n < (1ULL << 31U));
    }

    //B是每个item的byte数
    /**
     * Returns a pointer to the storage area starting from position 1 (e.g. without the checksum)
     * @tparam B the number of bytes per item
     * @param i context selector
     * @return
     */
    inline auto operator[](uint64_t i) -> uint8_t * {
      const uint8_t chk = checksum8(i, hashBits);
      i = finalize64(i, hashBits) * B & mask; /**< index: force bounds */
      //search for the checksum in t
      uint8_t *p = &t[0];
      if( p[i] == chk ) {
        return p + i + 1;
      }
      if( p[i ^ B] == chk ) {
        return p + (i ^ B) + 1;
      }
      if( p[i ^ (B * 2)] == chk ) {
        return p + (i ^ (B * 2)) + 1;
      }
      //未找到
      //not found, let's overwrite the lowest priority element
      if( p[i + 1] > p[(i + 1) ^ B] || p[i + 1] > p[(i + 1) ^ (B * 2)] ) {
        i ^= B;
      }
      if( p[i + 1] > p[(i + 1) ^ B ^ (B * 2)] ) {
        i ^= B ^ (B * 2);
      }
      memset(p + i, 0, B);
      p[i] = chk;
      return p + i + 1;
    };
};

#endif //PAQ8PX_HASHTABLE_HPP
