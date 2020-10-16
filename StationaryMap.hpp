#ifndef PAQ8PX_STATIONARYMAP_HPP
#define PAQ8PX_STATIONARYMAP_HPP

#include "IPredictor.hpp"
#include "DivisionTable.hpp"
#include "Hash.hpp"
#include "Mixer.hpp"
#include "Stretch.hpp"
#include "UpdateBroadcaster.hpp"

//对静态数据模型有用，直接查找上下文。
//对于建模的每个位元，一个32bit存储一个22位的预测和一个10位的适应速率偏移。
//更低的补偿意味着更快的适应。将在每次出现时增加，直到达到更高的上限。
/**
 * Map for modelling contexts of (nearly-)stationary data.
 * The context is looked up directly. For each bit modelled, a 32bit element stores
 * a 22 bit prediction and a 10 bit adaptation rate offset.
 *
 * - Rate: Initial adaptation rate offset [0..1023]. Lower offsets mean faster adaptation.
 * Will be increased on every occurrence until the higher bound is reached.
 *
 * Uses (2^(BitsOfContext+2))*((2^InputBits)-1) bytes of memory.
 */
class StationaryMap : IPredictor {
public:
    static constexpr int MIXERINPUTS = 2;

private:
    const Shared * const shared;
    Array<uint32_t> data;
    const uint32_t mask, maskBits, stride, bTotal;
    uint32_t context {};
    uint32_t bCount {};
    uint32_t b {};
    uint32_t *cp {};
    int scale;
    const uint16_t limit;
    int *dt;

public:
    //bitsOfContext：每个上下文使用多少位。更高的位被丢弃。
    //inputBits：有多少位[1..8]的输入将为每个上下文建模。必须按这些间隔设置新上下文。
    /**
     * Construct using (2^(BitsOfContext+2))*((2^InputBits)-1) bytes of memory.
     * @param bitsOfContext How many bits to use for each context. Higher bits are discarded.
     * @param inputBits How many bits [1..8] of input are to be modelled for each context. New contexts must be set at those intervals.
     * @param scale
     * @param limit
     */
    StationaryMap(const Shared* const sh, int bitsOfContext, int inputBits, int scale, uint16_t limit);

    /**
     * ctx must be a direct context (no hash)
     * @param ctx
     */
    void setDirect(uint32_t ctx);

    /**
     * ctx must be a hash
     * @param ctx
     */
    void set(uint64_t ctx);
    void setscale(int scale);
    void reset(int rate);
    void update() override;
    void mix(Mixer &m);
};

#endif //PAQ8PX_STATIONARYMAP_HPP
