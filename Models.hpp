#ifndef PAQ8PX_MODELS_HPP
#define PAQ8PX_MODELS_HPP

#include "text/TextModel.hpp"
#include "model/Audio16BitModel.hpp"
#include "model/Audio8BitModel.hpp"
#include "model/CharGroupModel.hpp"
#include "model/DmcForest.hpp"
#include "model/ExeModel.hpp"
#include "model/Image1BitModel.hpp"
#include "model/Image24BitModel.hpp"
#include "model/Image4BitModel.hpp"
#include "model/Image8BitModel.hpp"
#include "model/IndirectModel.hpp"
#include "model/JpegModel.hpp"
#include "model/LinearPredictionModel.hpp"
#include "model/MatchModel.hpp" //匹配模型
#include "model/NestModel.hpp"
#include "model/NormalModel.hpp" //正态模型
#include "model/RecordModel.hpp" //记录模型
#include "model/SparseMatchModel.hpp"
#include "model/SparseModel.hpp" //稀疏模型
#include "model/WordModel.hpp"   //字符模型
#include "model/XMLModel.hpp"
#include "model/DecAlphaModel.hpp"
#include "lstm/LstmModel.hpp"
#include "lstm/LstmFactory.hpp"

/**
 * This is a factory class for lazy object creation for models.
 * Objects created within this class are instantiated on first use and guaranteed to be destroyed.
 */
class Models {
private:
  Shared * const shared;
public:
  explicit Models(Shared* const sh);
  auto normalModel() -> NormalModel &;
  auto dmcForest() -> DmcForest &;
  auto charGroupModel() -> CharGroupModel &;
  auto recordModel() -> RecordModel &;
  auto sparseModel() -> SparseModel &;
  auto matchModel() -> MatchModel &;
  auto sparseMatchModel() -> SparseMatchModel &;
  auto indirectModel() -> IndirectModel &;
  auto textModel() -> TextModel &;
  auto wordModel() -> WordModel &;
  auto nestModel() -> NestModel &;
  auto xmlModel() -> XMLModel &;
  auto exeModel() -> ExeModel &;
  auto linearPredictionModel() -> LinearPredictionModel &;
  auto jpegModel() -> JpegModel &;
  auto image24BitModel() -> Image24BitModel &;
  auto image8BitModel() -> Image8BitModel &;
  auto image4BitModel() -> Image4BitModel &;
  auto image1BitModel() -> Image1BitModel &;
#ifndef DISABLE_AUDIOMODEL
  auto audio8BitModel() -> Audio8BitModel &;
  auto audio16BitModel() -> Audio16BitModel &;
#endif //DISABLE_AUDIOMODEL
  auto lstmModel() -> LstmModel<> &;
  auto decAlphaModel() -> DECAlphaModel &;
};

#endif //PAQ8PX_MODELS_HPP
