////////////////////////////////////////////////////////////////////////
//
// HoughLineService class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to find lines (Houghclusters) from clusters 
//  found by DBSCAN after deconvolution and hit finding.
//  The algorithm is based on: 
//  Queisser, A. "Computing the Hough Transform", C/C++ Users Journal 21, 12 (Dec. 2003).
//  Niblack, W. and Petkovic, D. On Improving the Accuracy of the Hough Transform", 
//  Machine Vision and Applications 3, 87 (1990)  
////////////////////////////////////////////////////////////////////////

#include "ClusterFinder/HoughLineService.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>

#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
 
#include "Filters/ChannelFilter.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"


//------------------------------------------------------------------------------
cluster::HoughLineService::HoughLineService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
{
  this->reconfigure(pset);
}

//------------------------------------------------------------------------------
cluster::HoughLineService::~HoughLineService()
{
}

//------------------------------------------------------------------------------
void cluster::HoughLineService::reconfigure(fhicl::ParameterSet const& pset)
{
  fMaxLines            = pset.get< int    >("MaxLines"           );
  fMinHits             = pset.get< int    >("MinHits"            );
  fSaveAccumulator     = pset.get< int    >("SaveAccumulator"    );
  fNumAngleCells       = pset.get< int    >("NumAngleCells"      );
  fMaxDistance         = pset.get< double >("MaxDistance"        );
  fMaxSlope            = pset.get< double >("MaxSlope"           );
  fRhoZeroOutRange     = pset.get< int    >("RhoZeroOutRange"    );
  fThetaZeroOutRange   = pset.get< int    >("ThetaZeroOutRange"  );
  fRhoResolutionFactor = pset.get< int    >("RhoResolutionFactor");
  fPerCluster          = pset.get< int    >("HitsPerCluster"     );
  fMissedHits          = pset.get< int    >("MissedHits"         );

  return;
}

//------------------------------------------------------------------------------
cluster::HoughTransform::HoughTransform()
{  
}

//------------------------------------------------------------------------------
cluster::HoughTransform::~HoughTransform()
{  
}

//------------------------------------------------------------------------------
bool cluster::HoughTransform::AddPoint(int x, int y)
{
  if (x>m_dx || y>m_dy || x<0.0 || y<0.0)
    return false;
  return DoAddPoint(x, y);
}
 
//------------------------------------------------------------------------------
void cluster::HoughTransform::Init(int dx, int dy, int rhores,
				   int numACells)
{
  m_numAngleCells=numACells;
  m_rhoResolutionFactor = rhores;
  m_accum.clear();
  m_accum.resize(m_numAngleCells);
  m_numAccumulated = 0;   
//   m_cosTable.clear();
//   m_sinTable.clear();
  m_cosTable.resize(m_numAngleCells);
  m_sinTable.resize(m_numAngleCells);
  if (dx == m_dx && dy == m_dy)
    return;
  m_dx = dx;
  m_dy = dy;
  m_rowLength = (int)(m_rhoResolutionFactor*2. * sqrt(dx*dx + dy*dy));
  
  int angleIndex;
  double a, angleStep = TMath::Pi()/m_numAngleCells;
  for (a=0.0, angleIndex = 0; angleIndex < m_numAngleCells; ++angleIndex){
    m_cosTable[angleIndex] = cos(a);
    m_sinTable[angleIndex] = sin(a);
    a += angleStep;
  }
}

//------------------------------------------------------------------------------
int cluster::HoughTransform::GetMax(int &xmax, int &ymax)
{
  std::map<int,int>::iterator rhoIter;
  int maxVal = -1;
  for(unsigned int i = 0; i < m_accum.size(); i++){
    for(rhoIter=m_accum[i].begin(); rhoIter!=m_accum[i].end(); ++rhoIter){
      if((*rhoIter).second > maxVal) {
	maxVal = (*rhoIter).second;
	xmax = i;
	ymax = (*rhoIter).first;
      }
    }
  }

  return maxVal;
}

//------------------------------------------------------------------------------
bool cluster::HoughTransform::DoAddPoint(int x, int y)
{
  distCenter = (int)(m_rowLength/2.);
 
  // prime the lastDist variable so our linear fill works below
  lastDist = (int)(distCenter+(m_rhoResolutionFactor*(m_cosTable[0]*x + m_sinTable[0]*y)));
  // loop through all angles a from 0 to 180 degrees
  for (unsigned int a=1; a<m_cosTable.size(); ++a){
    // Calculate the basic line equation dist = cos(a)*x + sin(a)*y.
    // Shift to center of row to cover negative values
    dist = (int)(distCenter+(m_rhoResolutionFactor*(m_cosTable[a]*x + m_sinTable[a]*y)));
    // sanity check to make sure we stay within our row
    if (dist >= 0 && dist<m_rowLength){
      if(lastDist==dist)
	m_accum[a][lastDist]++;
      else{
	// fill in all values in row a, not just a single cell
	stepDir = dist>lastDist ? 1 : -1;
	for (cell=lastDist; cell!=dist; cell+=stepDir){   
	  m_accum[a][cell]++;//maybe add weight of hit here?
	}      
      }
    }      
    lastDist = dist;
  }
  m_numAccumulated++;

  return true;
}

//------------------------------------------------------------------------------
//this method saves a BMP image of the Hough Accumulator, which can be viewed with gimp
void cluster::HoughLineService::HLSSaveBMPFile(const char *fileName, unsigned char *pix, int dx, int dy)
{
  ofstream bmpFile(fileName, std::ios::binary);
  bmpFile.write("B", 1);
  bmpFile.write("M", 1);
  int bitsOffset = 54 +256*4; 
  int size = bitsOffset + dx*dy; //header plus 256 entry LUT plus pixels
  bmpFile.write((const char *)&size, 4);
  int reserved = 0;
  bmpFile.write((const char *)&reserved, 4);
  bmpFile.write((const char *)&bitsOffset, 4);
  int bmiSize = 40;
  bmpFile.write((const char *)&bmiSize, 4);
  bmpFile.write((const char *)&dx, 4);
  bmpFile.write((const char *)&dy, 4);
  short planes = 1;
  bmpFile.write((const char *)&planes, 2);
  short bitCount = 8;
  bmpFile.write((const char *)&bitCount, 2);
  int i, temp = 0;
  for (i=0; i<6; i++)
    bmpFile.write((const char *)&temp, 4);  // zero out optional color info
  // write a linear LUT
  char lutEntry[4]; // blue,green,red
  lutEntry[3] = 0;  // reserved part
  for (i=0; i<256; i++)
    {
      lutEntry[0] = lutEntry[1] = lutEntry[2] = i;
      bmpFile.write(lutEntry, sizeof lutEntry);
    }
  // write the actual pixels
  bmpFile.write((const char *)pix, dx*dy);
}
 
//------------------------------------------------------------------------------
size_t cluster::HoughLineService::Transform(art::PtrVector<recob::Cluster>& clusIn, 
					    std::vector<recob::Cluster>& ccol)
{

  std::vector<int> skip;  

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;
  filter::ChannelFilter chanFilt;
  HoughTransform c;

  extern void SaveBMPFile(const char *f, unsigned char *pix, int dxx, int dyy);
  art::PtrVector<recob::Hit> cHits;
  art::PtrVector<recob::Hit> hit;

  for(size_t t = 0; t < geom->NTPC(); ++t){
    for(unsigned int p = 0; p < geom->Nplanes(t); ++p) {
      art::PtrVector<recob::Cluster>::const_iterator clusterIter = clusIn.begin();
      int clusterID=0;//the unique ID of the cluster
      // This is the loop over clusters. The algorithm searches for lines on a 
      // (DBSCAN) cluster-by-cluster basis. 
      //get the view of the current plane
      geo::View_t view = geom->Plane(p,t).View();

      while(clusterIter!=clusIn.end()) {
	hit.clear();
	cHits.clear();
	if(fPerCluster){
	  if((*clusterIter)->View() == view) hit = (*clusterIter)->Hits();
	}
	else{   
	  while(clusterIter!=clusIn.end()){
	    if( (*clusterIter)->View() == view ){
	      cHits = (*clusterIter)->Hits();
	      if(cHits.size() > 0){
		art::PtrVector<recob::Hit>::const_iterator hitIter = cHits.begin();
		while (hitIter!=cHits.end()){
		  hit.push_back((*hitIter));
		  hitIter++;
		}
	      }
	    }// end if cluster is in correct view
	    clusterIter++;
	  }//end loop over clusters
	}//end if not fPerCluster
	if(hit.size() == 0){ 
	  if(fPerCluster) clusterIter++;
	  continue;
	}
        //factor to make x and y scale the same units
        double xyScale = .001*larprop->DriftVelocity(larprop->Efield(),larprop->Temperature())*detprop->SamplingRate()/geom->WirePitch(0,1,p,t);

	int x, y;
	unsigned int channel, plane, wire, tpc;
	//there must be a better way to find which plane a cluster comes from
	int dx=geom->Nwires(p,t);//number of wires 
	int dy=hit[0]->Wire()->fSignal.size();//number of time samples. 
	skip.clear();
	skip.resize(hit.size());
	std::vector<int> listofxmax;
	std::vector<int> listofymax;  
	std::vector<int> hitTemp;        //indecies ofcandidate hits
	std::vector<int> sequenceHolder; //channels of hits in list
	std::vector<int> currentHits;    //working vector of hits 
	std::vector<int> lastHits;       //best list of hits
	art::PtrVector<recob::Hit> clusterHits;
	double indcolscaling=0.;         //a parameter to account for the different 
                                       //characteristic hit width of induction and collection plane
	double centerofmassx=0;
	double centerofmassy=0;
	double denom = 0; 
	double intercept=0.;
	double slope=0.;
	//this array keeps track of the hits that have already been associated with a line. 
	int xMax=0;
	int yMax=0;
	double rho;
	double theta; 
	int accDx(0), accDy(0);
	
	for (int linenum = 0; linenum < fMaxLines; ++linenum){ 
	  //Init specifies the size of the two-dimensional accumulator 
	  //(based on the arguments, number of wires and number of time samples). 
	  c.Init(dx,dy,fRhoResolutionFactor,fNumAngleCells);
	  //adds all of the hits (that have not yet been associated with a line) to the accumulator
	  
	  for(unsigned int i=0;i < hit.size(); ++i){
	    channel=hit[i]->Wire()->RawDigit()->Channel();
	    geom->ChannelToWire(channel,tpc,plane,wire);
	    if (skip[i] != 1)
	      c.AddPoint(wire,(int)(hit[i]->PeakTime()));
	  }// end loop over hits
 	   
	  //gets the actual two-dimensional size of the accumulator
	  c.GetAccumSize(accDy, accDx);
	  
	  // zeroes out the neighborhood of all previous lines  
	  for(unsigned int i = 0; i < listofxmax.size(); ++i){
	    int yClearStart = listofymax[i]-fRhoZeroOutRange;
	    if (yClearStart < 0) yClearStart = 0;
	    
	    int yClearEnd = listofymax[i]+fRhoZeroOutRange;
	    if (yClearEnd >= accDy) yClearEnd = accDy-1;
	    
	    int xClearStart = listofxmax[i]-fThetaZeroOutRange;
	    if (xClearStart < 0) xClearStart = 0;
	    
	    int xClearEnd = listofxmax[i]+fThetaZeroOutRange;
	    if (xClearEnd >= accDx) xClearEnd = accDx-1;
	    
	    for (y = yClearStart; y <= yClearEnd; ++y){
	      for (x = xClearStart; x <= xClearEnd; ++x){
		c.SetCell(y,x,0);
	      }
	    }
	  }// end loop over size of listxmax
	  
	  //find the weightiest cell in the accumulator.
	  int maxCell = 0;
	  xMax = 0;
	  yMax = 0;
	  maxCell = c.GetMax(yMax,xMax);
	  // break when biggest maximum is smaller than fMinHits
	  if ( maxCell < fMinHits ) 
	    break;
	  
	  //find the center of mass of the 3x3 cell system (with maxCell at the center). 
	  denom=centerofmassx=centerofmassy=0;
	  
	  if(xMax>0 && yMax>0 && xMax+1<accDx && yMax+1<accDy){  
	    for(int i = -1; i < 2; ++i){
	      for(int j = -1; j < 2; ++j){
		denom+=c.GetCell(yMax+i,xMax+j);
		centerofmassx+=j*c.GetCell(yMax+i,xMax+j);
		centerofmassy+=i*c.GetCell(yMax+i,xMax+j);
	      }
	    }
	    centerofmassx/=denom;
	    centerofmassy/=denom;      
	  }
	  else  centerofmassx=centerofmassy=0;
	  
	  //fill the list of cells that have already been found
	  listofxmax.push_back(xMax);
	  listofymax.push_back(yMax);
	  
	  c.GetEquation(yMax+centerofmassy, xMax+centerofmassx, rho, theta);
	  slope=-1./tan(theta);    
	  intercept=(rho/sin(theta));
	  double distance;
	  /// \todo: the collection plane's characteristic hit width's are, 
	  /// \todo: on average, about 5 time samples wider than the induction plane's. 
	  /// \todo: this is hard-coded for now.
	  if(p==0)
	    indcolscaling=5.;
	  else
	    indcolscaling=0.;
	  
	  if(!isinf(slope) && !isnan(slope)){
	    sequenceHolder.clear();
	    hitTemp.clear();
	    for(unsigned int i=0;i<hit.size();i++){
	      channel=hit[i]->Wire()->RawDigit()->Channel();
	      geom->ChannelToWire(channel,tpc,plane,wire);
	      distance=(TMath::Abs(hit[i]->PeakTime()-slope*(double)(wire)-intercept)/(sqrt(pow(xyScale*slope,2)+1)));
	    
	      if(distance < fMaxDistance+((hit[i]->EndTime()-hit[i]->StartTime())/2.)+indcolscaling  && skip[i]!=1){
		hitTemp.push_back(i);
		sequenceHolder.push_back(channel);
	      }
	     
	    }// end loop over hits
	    
	    if(hitTemp.size() < 2) continue;
	    currentHits.clear();  
	    lastHits.clear();
	    int j; 
	    currentHits.push_back(0);
	    for(unsigned int i=0;i+1<sequenceHolder.size();i++) {  
	      j = 1;
	      while((chanFilt.BadChannel(sequenceHolder[i]+j))==true) j++;
	      if(sequenceHolder[i+1]-sequenceHolder[i]<=j+fMissedHits) currentHits.push_back(i+1);
	      else if(currentHits.size() > lastHits.size()) {
		lastHits = currentHits;
		currentHits.clear();
	      }
	      else currentHits.clear();
	    } 
	    
	    if(currentHits.size() > lastHits.size()) lastHits = currentHits;
	    clusterHits.clear();    
	    for(unsigned int i = 0; i < lastHits.size();i++) {
	      clusterHits.push_back(hit[hitTemp[lastHits[i]]]);
	      skip[hitTemp[lastHits[i]]]=1;
	    } 
	    //protection against very steep uncorrelated hits
	    if(TMath::Abs(slope)>fMaxSlope 
	       && TMath::Abs((*clusterHits.begin())->Wire()->RawDigit()->Channel()-
			     clusterHits[clusterHits.size()-1]->Wire()->RawDigit()->Channel())>=0
	       )
	      continue;
	    
	    unsigned int sw = 0;
	    unsigned int ew = 0;
            unsigned int sc = 0;
            unsigned int ec = 0;
	    sc = (*clusterHits.begin())->Wire()->RawDigit()->Channel(); 
	    geom->ChannelToWire(sc,tpc,plane,sw);
	    
	  //there doesn't seem to be a std::vector::back() method to get the last element of the vector here
	    ec = (*(clusterHits.end()-1))->Wire()->RawDigit()->Channel(); 
	    geom->ChannelToWire(ec,tpc,plane,ew);
	  
	    recob::Cluster cluster(clusterHits,
				   sw, 0.,
				   (*clusterHits.begin())->PeakTime(), 0.,
				   ew, 0., 
				   (clusterHits[clusterHits.size()-1])->PeakTime(), 0.,
				   slope, 0., 
				   -999., 0., 
				   clusterID);	      
	    
	    clusterID++;
	    ccol.push_back(cluster);
            //allow double assignment of first and last hits
            for(unsigned int i = 0; i < lastHits.size(); i++) if(skip[hitTemp[lastHits[i]]] ==1) 
            {
	      channel = hit[hitTemp[lastHits[i]]]->Wire()->RawDigit()->Channel();  
              if( channel == sc || channel == ec) skip[i] = 0;
            }
    
              
	  }// end if !isnan

	}// end loop over number of lines found
  
	// saves a bitmap image of the accumulator (useful for debugging), 
	// with scaling based on the maximum cell value
	if(fSaveAccumulator){   
	  unsigned char *outPix = new unsigned char [accDx*accDy];
	  //finds the maximum cell in the accumulator for image scaling
	  int cell, pix=0, maxCell=0;
	  for (y=0; y<accDy; y++) 
            for (x=0; x<accDx; x++)
	      {
		cell = c.GetCell(y,x);
		if (cell > maxCell) maxCell = cell;
	      }
	  for (y=0; y<accDy; y++){
	    for (x=0; x<accDx; x++){ 
	      //scales the pixel weights based on the maximum cell value     
	      if(maxCell>0)
		pix = (int)((1500*c.GetCell(y,x))/maxCell);
	      outPix[y*accDx + x] = pix;
	    }
	  }
	  
	  SaveBMPFile("houghaccum.bmp", outPix, accDx, accDy);
	  delete [] outPix;
	}// end if saving accumulator
	
	hit.clear();
	lastHits.clear();
	if(clusterIter!=clusIn.end()) clusterIter++;
	listofxmax.clear();
	listofymax.clear();
      }//end loop over clusters
      
    }//end loop over planes
  }// end loop over tpcs

  return ccol.size(); 
}

//------------------------------------------------------------------------------
size_t cluster::HoughLineService::Transform(std::vector< art::Ptr<recob::Hit> >& hits,
					    double &slope,
					    double &intercept)
{
  HoughTransform c;

  art::ServiceHandle<geo::Geometry> geom;
  int dx = geom->Nwires(0);               //number of wires 
  int dy = hits[0]->Wire()->fSignal.size();//number of time samples. 

  c.Init(dx,dy,fRhoResolutionFactor,fNumAngleCells);

  unsigned int plane = 0;
  unsigned int wire  = 0;
  unsigned int tpc   = 0;
  for(unsigned int i=0;i < hits.size(); ++i){
    unsigned short channel = hits[i]->Wire()->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);
    c.AddPoint(wire, (int)(hits[i]->PeakTime()));
  }// end loop over hits

  //gets the actual two-dimensional size of the accumulator
  int accDx = 0;
  int accDy = 0;
  c.GetAccumSize(accDy, accDx);

  //find the weightiest cell in the accumulator.
  int xMax = 0;
  int yMax = 0;
  c.GetMax(yMax,xMax);

  //find the center of mass of the 3x3 cell system (with maxCell at the center). 
  double centerofmassx = 0.;
  double centerofmassy = 0.;
  double denom         = 0.;
    
  if(xMax > 0 && yMax > 0 && xMax+1 < accDx && yMax+1 < accDy){  
    for(int i = -1; i < 2; ++i){
      for(int j = -1; j < 2; ++j){
	denom         += c.GetCell(yMax+i,xMax+j);
	centerofmassx += j*c.GetCell(yMax+i,xMax+j);
	centerofmassy += i*c.GetCell(yMax+i,xMax+j);
      }
    }
    centerofmassx /= denom;
    centerofmassy /= denom;      
  }
  else  centerofmassx = centerofmassy = 0;

  double rho   = 0.;
  double theta = 0.;
  c.GetEquation(yMax+centerofmassy, xMax+centerofmassx, rho, theta);
  slope     = -1./tan(theta);    
  intercept = rho/sin(theta);
  
  // \todo could eventually refine this method to throw out hits that are 
  // \todo far from the hough line and refine the fit

  return hits.size();
}


