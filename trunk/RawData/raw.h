/// \file    raw.h
/// \brief   Collect all the RawData header files together
/// \author  brebel@fnal.gov
/// \version $Id: raw.h,v 1.2 2010/04/15 18:10:18 brebel Exp $
#ifndef RAWDATA_RAW_H
#define RAWDATA_RAW_H
#include "RawData/DAQHeader.h"
#include "RawData/RawDigit.h"
#include "RawData/BeamInfo.h"
#include <vector>

namespace raw{

  void Uncompress(const std::vector<short> adc, std::vector<short> &uncompressed, raw::Compress_t compress);
  void Compress(std::vector<short> &adc,   raw::Compress_t compress);
  void CompressHuffman(std::vector<short> &adc);
  void UncompressHuffman(const std::vector<short> adc, std::vector<short> &uncompressed);

}

#endif
