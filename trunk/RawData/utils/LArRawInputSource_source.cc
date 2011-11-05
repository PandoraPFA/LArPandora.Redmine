// ======================================================================
//
// LArRawInputSource_plugin.cc
//
// ======================================================================

#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/ReaderSource.h"

#include "RawData/utils/LArRawInputDriver.h"

namespace lris {
  typedef art::ReaderSource<LArRawInputDriver> LArRawInputSource;
}

DEFINE_ART_INPUT_SOURCE(lris::LArRawInputSource);
