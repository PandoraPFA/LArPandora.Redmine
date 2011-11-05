///
/// \file    HeaderPad.cxx
/// \brief   Drawing pad for time or charge histograms
/// \author  messier@indiana.edu
/// \version $Id: HeaderPad.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
#include "EventDisplay/HeaderPad.h"
#include "EventDisplayBase/evdb.h"
#include "EventDisplay/HeaderDrawer.h"
#include "TPad.h"
#include "TText.h"
using namespace evd;

static const int kTPAD = 0;
static const int kQPAD = 1;
static const int kRAW   = 0;
static const int kCALIB = 1;
static const int kPE =  2;
static const int kTNS = 3;

//......................................................................

HeaderPad::HeaderPad(const char* nm, const char* ti,
		     double x1, double y1,
		     double x2, double y2,
		     const char* /*opt*/) :
  DrawingPad(nm, ti, x1, y1, x2, y2)
{ 
  fView = new evdb::View2D();
}

//......................................................................

HeaderPad::~HeaderPad() 
{
  if (fView!=0) { delete fView; fView = 0; }
}

//......................................................................

void HeaderPad::Draw(const char* /* opt */) 
{
  fView->Clear();

  this->HeaderDraw()->Header(fView);

  this->Pad()->cd();
  fView->Draw();
}

//......................................................................


//////////////////////////////////////////////////////////////////////////
