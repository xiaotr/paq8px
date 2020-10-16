#ifndef PAQ8PX_BUCKET_HPP
#define PAQ8PX_BUCKET_HPP

#include <cstdint>
#include <cstring>
#include "Shared.hpp"
#include "utils.hpp"

/**
 * Hash bucket, 64 bytes
 *
 * Contains bit histories (statistics only) and byte histories (byte values). //包含位历史(仅统计)和字节历史(字节值)。
 * For every byte the bit histories are stored in the "bitState" array  //对于每一个字节，位历史记录都存储在“bitState”数组中
 * (7 bytes) in 3 different buckets:
 * The 1st bucket: 
 *  [bit statistics for bit 0] 
 *  [bit statistics for bit 1 when bit 0 was 0] 
 *  [bit statistics for bit 1 when bit 0 was 1] 
 *  [union of byte run count and the pending flag (0xFF)] //字节运行计数和挂起标志(0xFF)的联合
 *  [byte value: byte in context #1 - belongs to run count]
 *  [byte value: byte in context #2]
 *  [byte value: byte in context #3]
 * The 2nd bucket:
 *  [bit statistics for bit 2]
 *  [bit statistics for bit 3 when bit 2 was 0]
 *  [bit statistics for bit 3 when bit 2 was 1]
 *  [bit statistics for bit 4 when bit 2-3 was 00]
 *  [bit statistics for bit 4 when bit 2-3 was 01]
 *  [bit statistics for bit 4 when bit 2-3 was 10]
 *  [bit statistics for bit 4 when bit 2-3 was 11]
 * The 3rd bucket:
 *  [bit statistics for bit 5]
 *  [bit statistics for bit 6 when bit 5 was 0]
 *  [bit statistics for bit 6 when bit 5 was 1]
 *  [bit statistics for bit 7 when bit 5-6 was 00]
 *  [bit statistics for bit 7 when bit 5-6 was 01]
 *  [bit statistics for bit 7 when bit 5-6 was 10]
 *  [bit statistics for bit 7 when bit 5-6 was 11]
 */

class Bucket {
    uint16_t checksums[7]; /**< byte context checksums */
    uint8_t mostRecentlyUsed; /**< last 2 accesses (0-6) in low, high nibble */
public:
    uint8_t bitState[7][7]; /**< byte context, 3-bit context -> bit history state */
    // bitState[][0] = 1st bit, bitState[][1,2] = 2nd bit, bitState[][3..6] = 3rd bit  第1bit、第2bit在第1bit条件下，第3bit在第1与2bit下
    // bitState[][0] is also a replacement priority, 0 = empty

    /**
     * Find or create hash element matching checksum. 查找或创建匹配校验和的哈希元素。
     * If not found, insert or replace lowest priority (skipping 2 most recent).
     * @param checksum
     * @return
     */

    inline auto find(const uint16_t checksum) -> uint8_t* {
      if( checksums[mostRecentlyUsed & 15] == checksum ) {
        return &bitState[mostRecentlyUsed & 15][0];
      }
      int worst = 0xFFFF, idx = 0;
      for( int i = 0; i < 7; ++i ) {
        if( checksums[i] == checksum ) {
          mostRecentlyUsed = mostRecentlyUsed << 4 | i; //高位是第2近的，低位是最近的
          return &bitState[i][0];
        }
        //不是最近2次访问，且优先级最小
        if( bitState[i][0] < worst && (mostRecentlyUsed & 15) != i && mostRecentlyUsed >> 4 != i ) {
          worst = bitState[i][0];
          idx = i;
        }
      }
      mostRecentlyUsed = 0xF0 | idx;
      checksums[idx] = checksum;
      return (uint8_t *) memset(&bitState[idx][0], 0, 7);
    }
};

#endif //PAQ8PX_BUCKET_HPP
