/// \file  ColorDrawingOptions.h
/// \brief The color scales used by the event display
/// \author messier@indiana.edu
/// \version $Id: ColorScaleTable.h,v 1.1.1.1 2010/11/10 19:44:54 p-novaart Exp $
#ifndef EVD_COLORDRAWINGOPTIONS_H
#define EVD_COLORDRAWINGOPTIONS_H
#ifndef __CINT__
#include "EventDisplayBase/ColorScale.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace evd {
  class ColorDrawingOptions {
  public:

    ColorDrawingOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~ColorDrawingOptions();

    void reconfigure(fhicl::ParameterSet const& pset);

    const evdb::ColorScale& RawQ() const;
    const evdb::ColorScale& CalQ() const;
    const evdb::ColorScale& RawT() const;
    const evdb::ColorScale& CalT() const;

    int    fColorOrGray; ///< 0 = color, 1 = gray
    int    fRawDiv;      ///< number of divisions in raw
    int    fRecoDiv;     ///< number of divisions in raw
    double fRawQLow;     ///< low  edge of ADC values for drawing raw digits
    double fRawQHigh;    ///< high edge of ADC values for drawing raw digits
    double fRecoQLow;    ///< low  edge of ADC values for drawing raw digits
    double fRecoQHigh;   ///< high edge of ADC values for drawing raw digits

  private:
    evdb::ColorScale fColorScaleRaw;
    evdb::ColorScale fColorScaleReco;
    evdb::ColorScale fGrayScaleRaw;
    evdb::ColorScale fGrayScaleReco;
  };
}
#endif // __CINT__
#endif
////////////////////////////////////////////////////////////////////////
