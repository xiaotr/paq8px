#include "ContextMap.hpp"

//创建hash桶
ContextMap::ContextMap(const Shared* const sh, uint64_t m, const int contexts) : 
  shared(sh), C(contexts), t(m >> 6U), cp(contexts), cp0(contexts), cxt(contexts), chk(contexts), runP(contexts),
  sm(sh, contexts, 256, 1023, StateMap::BitHistory), cn(0),
  mask(uint32_t(t.size() - 1)), hashBits(ilog2(mask + 1)), validFlags(0) {
#ifdef VERBOSE
  printf("Created ContextMap with m = %" PRIu64 ", contexts = %d\n", m, contexts);
#endif
  assert(m >= 64 && isPowerOf2(m));
  static_assert(sizeof(Bucket) == 64, "Size of Bucket should be 64!");
  assert(C <= (int) sizeof(validFlags) * 8); // validFlags is 64 bits - it can't support more than 64 contexts
}

//cn是现存的上下文数
void ContextMap::set(const uint64_t cx) {
  assert(cn >= 0 && cn < C);
  const uint32_t ctx = cxt[cn] = finalize64(cx, hashBits); //cx的hash值
  const uint16_t checksum = chk[cn] = checksum16(cx, hashBits);
  uint8_t* base = cp0[cn] = cp[cn] = t[ctx].find(checksum);
  runP[cn] = base + 3; //runp的位置
  // update pending bit histories for bits 2-7 //更新2-7位暂挂的位历史记录
  //count|d = 10,第一次遇到
  if( base[3] == 2 ) {
    const int c = base[4] + 256;
    //第2bit开始
    uint8_t *p = t[(ctx + (c >> 6U)) & mask].find(checksum);
    p[0] = 1 + ((c >> 5U) & 1U);
    p[1 + ((c >> 5U) & 1U)] = 1 + ((c >> 4U) & 1U);
    p[3 + ((c >> 4U) & 3U)] = 1 + ((c >> 3U) & 1U);
    //第5bit开始
    p = t[(ctx + (c >> 3U)) & mask].find(checksum);
    p[0] = 1 + ((c >> 2U) & 1U);
    p[1 + ((c >> 2U) & 1U)] = 1 + ((c >> 1U) & 1U);
    p[3 + ((c >> 1U) & 3U)] = 1 + (c & 1U);
  }
  cn++;
  validFlags = (validFlags << 1U) + 1; //标识符号
}

void ContextMap::skip() {
  assert(cn >= 0 && cn < C);
  cn++;
  validFlags <<= 1U;
}

void ContextMap::update() {
  INJECT_SHARED_y
  INJECT_SHARED_bpos
  INJECT_SHARED_c1
  INJECT_SHARED_c0
  //对cn个上下文
  for( int i = 0; i < cn; ++i ) {
    //确实存在
    if(((validFlags >> (cn - 1 - i)) & 1) != 0 ) {
      // update bit history state byte
      if( cp[i] != nullptr ) {
        assert(cp[i] >= &t[0].bitState[0][0] && cp[i] <= &t[t.size() - 1].bitState[6][6]);
        assert((uintptr_t(cp[i]) & 63) >= 15);
        StateTable::update(cp[i], y, rnd); //更新bit位历史
      }

      // update context pointers
      if( bpos > 1 && runP[i][0] == 0 ) {
        cp[i] = nullptr;
      } else {
        //对不同的bpos位指向不同位历史
        switch( bpos ) {
          case 1:
          case 3:
          case 6:
            cp[i] = cp0[i] + 1 + (c0 & 1U);
            break;
          case 4:
          case 7:
            cp[i] = cp0[i] + 3 + (c0 & 3U);
            break;
          case 2:
          case 5: {
            const uint16_t checksum = chk[i];
            const uint32_t ctx = cxt[i];
            cp0[i] = cp[i] = t[(ctx + c0) & mask].find(checksum);
            break;
          }
          case 0: {
            // 更新run表
            // update run count of previous context
            if( runP[i][0] == 0 ) { // new context ，之前是空的新加入
              runP[i][0] = 2, runP[i][1] = c1;
            } else if( runP[i][1] != c1 ) { // different byte in context ， 遇到了不同的
              runP[i][0] = 1, runP[i][1] = c1;
            } else if( runP[i][0] < 254 ) { // same byte in context
              runP[i][0] += 2;
            } else if( runP[i][0] == 255 ) { //连续太多了，减半
              runP[i][0] = 128;
            }
            break;
          }
          default:
            assert(false);
        }
      }
    }
  }
  if( bpos == 0 ) {
    cn = 0;
    validFlags = 0;
  }
}

void ContextMap::mix(Mixer &m) {
  shared->GetUpdateBroadcaster()->subscribe(this);//订阅更新
  sm.subscribe();
  INJECT_SHARED_bpos
  INJECT_SHARED_c0
  for( int i = 0; i < cn; ++i ) {
    if(((validFlags >> (cn - 1 - i)) & 1U) != 0 ) {
      // predict from last byte in context
      if((runP[i][1] + 256) >> (8 - bpos) == c0 ) { //如果rup指向的和c0一样，c0为当前bit的
        int rc = runP[i][0]; // count*2, +1 if 2 different bytes seen
        int sign = (runP[i][1] >> (7 - bpos) & 1U) * 2 - 1; // predicted bit + for 1, - for 0
        int c = ilog->log(rc + 1) << (2 + (~rc & 1U));
        m.add(sign * c);
      } else {
        m.add(0); //p=0.5
      }

      //单个bit的预测
      // predict from bit context
      const int s = cp[i] != nullptr ? *cp[i] : 0;
      if( s == 0 ) { //skip context
        sm.skip(i);
        m.add(0);
        m.add(0);
        m.add(0);
        m.add(0);
      } else { //cp非空
        const int p1 = sm.p2(i, s); //statemap的概率
        const int st = stretch(p1) >> 2U;
        const int contextIsYoung = int(s <= 2); //s比2小
        m.add(st >> contextIsYoung);
        m.add((p1 - 2048) >> 3U);
        const int n0 = -!StateTable::next(s, 2); //0的值？
        const int n1 = -!StateTable::next(s, 3); //1的值？
        m.add((n0 | n1) & st); // when both counts are nonzero add(0) otherwise add(st)
        const int p0 = 4095 - p1; //0的概率
        m.add(((p1 & n0) - (p0 & n1)) >> 4U);
      }
    } else { //skipped context
      sm.skip(i);
      m.add(0);
      m.add(0);
      m.add(0);
      m.add(0);
      m.add(0);
    }
  }
}
