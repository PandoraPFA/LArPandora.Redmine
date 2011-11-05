////////////////////////////////////////////////////////////////////////
//
// EndPointService class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to find (weak) vertices from hits after deconvolution and hit finding. 
//  A weak end point is a end point that has been found using a dedicated end point finding algorithm only. A 
//  strong end point is a end point that has been found using a dedicated end point finding algorithm and matched 
//  to a crossing of two or more HoughLineFinder lines. The End PointMatch module finds strong vertices.
////////////////////////////////////////////////////////////////////////
/// The algorithm is based on:
///C. Harris and M. Stephens (1988). "A combined corner and edge detector". Proceedings of the 4th Alvey 
///Vision Conference. pp. 147-151.
///B. Morgan (2010). "Interest Point Detection for Reconstruction in High Granularity Tracking Detectors". 
///arXiv:1006.3012v1 [physics.ins-det]
//Thanks to B. Morgan of U. of Warwick for comments and suggestions

#include <iostream>
#include <vector>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ClusterFinder/EndPointService.h"

#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "TMath.h"


#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

//-----------------------------------------------------------------------------
cluster::EndPointService::EndPointService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) :
 
  fTimeBins        (pset.get< int    >("TimeBins")      ),
  fMaxCorners      (pset.get< int    >("MaxCorners")    ),
  fGsigma          (pset.get< double >("Gsigma")        ),
  fWindow          (pset.get< int    >("Window")        ),
  fThreshold       (pset.get< double >("Threshold")     ),
  fSaveVertexMap   (pset.get< int    >("SaveVertexMap") )
{
}

//-----------------------------------------------------------------------------
cluster::EndPointService::~EndPointService()
{
}

//-----------------------------------------------------------------------------
double cluster::EndPointService::Gaussian(int x, int y, double sigma)
{
  double Norm=1./sqrt(2*TMath::Pi()*pow(sigma,2));
  double value=Norm*exp(-(pow(x,2)+pow(y,2))/(2*pow(sigma,2)));
  return value;
}

//-----------------------------------------------------------------------------
double cluster::EndPointService::GaussianDerivativeX(int x,int y)
{
  double Norm=1./(sqrt(2*TMath::Pi())*pow(fGsigma,3));
  double value=Norm*(-x)*exp(-(pow(x,2)+pow(y,2))/(2*pow(fGsigma,2)));
  return value;
}

//-----------------------------------------------------------------------------
double cluster::EndPointService::GaussianDerivativeY(int x,int y)
{
  double Norm=1./(sqrt(2*TMath::Pi())*pow(fGsigma,3));
  double value=Norm*(-y)*exp(-(pow(x,2)+pow(y,2))/(2*pow(fGsigma,2)));
  return value;
}


//-----------------------------------------------------------------------------
//this method saves a BMP image of the vertex map space, which can be viewed with gimp
void cluster::EndPointService::VSSaveBMPFile(const char *fileName, unsigned char *pix, int dx, int dy)
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
      lutEntry[0] =i;
      lutEntry[1] =i+1;
      lutEntry[2] = i+2;
      bmpFile.write(lutEntry, sizeof lutEntry);
    }
  // write the actual pixels
  bmpFile.write((const char *)pix, dx*dy);
}

//......................................................
size_t cluster::EndPointService::EndPoint(art::PtrVector<recob::Cluster>& clusIn, std::vector<recob::EndPoint2D>& vtxcol)
{

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larp;
   
  // // function below does not link, so commenting
  //extern void SaveBMPFile(const char *f, unsigned char *pix, int dxx, int dyy);
  
  
  //Point to a collection of vertices to output.
 
  art::PtrVector<recob::Hit> cHits;
  art::PtrVector<recob::Hit> hit;
  

 

  int flag=0;
  int windex=0;//the wire index to make sure the end point finder does not fall off the edge of the hit map
  int tindex=0;//the time index to make sure the end point finder does not fall off the edge of the hit map
  int n=0; //index of window cell. There are 49 cells in the 7X7 Gaussian and Gaussian derivative windows
  unsigned int numberwires;
  double numbertimesamples;
  double MatrixAAsum,MatrixBBsum,MatrixCCsum;
  std::vector<double> Cornerness2;
  //gaussian window definitions. The cell weights are calculated here to help the algorithm's speed
  double  w[49]={0.};
  double wx[49]={0.};
  double wy[49]={0.};
  int ctr = 0;
  for(int i = -3; i < 4; ++i){
    for(int j = 3; j > -4; --j){
      w[ctr] = Gaussian(i, j, fGsigma);
      wx[ctr] = GaussianDerivativeX(i,j);
      wy[ctr] = GaussianDerivativeY(i,j);
      ++ctr;
    }
  }
  
  unsigned int channel,plane,wire,wire2,tpc;
  for(unsigned int t = 0; t < geom->NTPC(); ++t){
    for(unsigned int p = 0; p < geom->Nplanes(t); p++) {
      art::PtrVector<recob::Hit> vHits;
      art::PtrVector<recob::Cluster>::const_iterator clusterIter = clusIn.begin();
      geo::View_t view = geom->Plane(p,t).View();
      hit.clear();
      cHits.clear();      
      while(clusterIter!= clusIn.end() ) {
	if((*clusterIter)->View() == view){
	  cHits = (*clusterIter)->Hits();
	  if(cHits.size() > 0)
	    //hit.insert(hit.end(),cHits.begin(),cHits.end());
	    for(unsigned int i = 0; i < cHits.size(); i++)
	      hit.push_back(cHits[i]);
	}//end if cluster is in the correct view
	clusterIter++;  
      } 
      if(hit.size() == 0) 
        continue;

       numberwires=geom->Nwires(p,t);
       numbertimesamples=hit[0]->Wire()->fSignal.size();
      std::cout << " --- endpoints check " << numberwires << " " << numbertimesamples << " " << fTimeBins << std::endl;
 
      std::vector < std::vector < double > > MatrixAsum;
      MatrixAsum.resize(numberwires);
      std::vector < std::vector < double > >  MatrixBsum;
      MatrixBsum.resize(numberwires);
      std::vector < std::vector < double > >  hit_map;//the map of hits 
      hit_map.resize(numberwires);
      std::vector < std::vector < int > > hit_loc;//the index of the hit that corresponds to the potential corner
      hit_loc.resize(numberwires);
      std::vector < std::vector < double > >  Cornerness;//the "weight" of a corner
      Cornerness.resize(numberwires);
  
  
      for(unsigned int wi=0;wi < numberwires; wi++)
      {  hit_map[wi].resize(fTimeBins,0);
         hit_loc[wi].resize(fTimeBins,-1);
	 Cornerness[wi].resize(fTimeBins,0);
	 MatrixAsum[wi].resize(fTimeBins,0);
	 MatrixBsum[wi].resize(fTimeBins,0);
/*      
	for(unsigned int timebin=0;timebin < fTimeBins; timebin++){
	  hit_map[wi][timebin]=0.;
	  hit_loc[wi][timebin]=-1;
	  Cornerness[wi][timebin]=0.;
	  MatrixAsum[wi][timebin]=0.;
	  MatrixBsum[wi][timebin]=0.;
	}*/
      }      
      for(unsigned int i=0;i < hit.size(); i++){
	channel=hit[i]->Wire()->RawDigit()->Channel();
	geom->ChannelToWire(channel,tpc,plane,wire);
	//pixelization using a Gaussian
	for(int j=0;j <= (int)(hit[i]->EndTime()-hit[i]->StartTime()+.5); j++)    
	  hit_map[wire][(int)((hit[i]->StartTime()+j)*(fTimeBins/numbertimesamples)+.5)]+=Gaussian((int)(j-((hit[i]->EndTime()-hit[i]->StartTime())/2.)+.5),0,hit[i]->EndTime()-hit[i]->StartTime());      
      }
	
      ////Gaussian derivative convolution  
      for(unsigned int wire=1;wire < numberwires-1; wire++)
      {//std:: cout << "loop for wire " << wire << std::endl;
	for(int timebin=1;timebin < fTimeBins-1; timebin++){
	  MatrixAsum[wire][timebin]=0.;
	  MatrixBsum[wire][timebin]=0.;
	  n=0;
	  for(int i = -3; i <= 3; i++) {
	      windex=wire+i;
	      if(windex < 0 ) windex=0;
	      // this is ok, because the line before makes sure it's not negative
	      else if ((unsigned int)windex >= numberwires) windex=numberwires-1; 
	    for(int j = -3; j <= 3; j++) {
	      tindex=timebin+j;
	      if(tindex<0) tindex=0;
	      else if (tindex>=fTimeBins) tindex=fTimeBins-1;
	      
	      MatrixAsum[wire][timebin]+=wx[n]*hit_map[windex][tindex];  
	      MatrixBsum[wire][timebin]+=wy[n]*hit_map[windex][tindex]; 
	      n++;
	    }
	  }
	 }
      }
      //calculate the cornerness of each pixel while making sure not to fall off the hit map.
      for(unsigned int wire=1;wire < numberwires-1; wire++)
      { //std:: cout << "loop for wire " << wire << std::endl;
	for(int timebin=1;timebin < fTimeBins-1; timebin++){    
	  MatrixAAsum=0;
	  MatrixBBsum=0;
	  MatrixCCsum=0;
	  //Gaussian smoothing convolution
	  n=0;
	  for(int i = -3; i <= 3; i++) {
	     windex=wire+i;
	      if(windex<0) windex=0;
	      else if ((unsigned int)windex >= numberwires) windex=numberwires-1; // this is ok, because the line before makes sure it's not negative
	    for(int j = -3; j <= 3; j++) {
	      tindex=timebin+j;
	      if(tindex<0) tindex=0;
	      else if (tindex>=fTimeBins) tindex=fTimeBins-1;
	     
	      MatrixAAsum+=w[n]*pow(MatrixAsum[windex][tindex],2);  
	      MatrixBBsum+=w[n]*pow(MatrixBsum[windex][tindex],2);                   
	      MatrixCCsum+=w[n]*MatrixAsum[windex][tindex]*MatrixBsum[windex][tindex]; 
	      n++;
	    }
	  }
	  if((MatrixAAsum+MatrixBBsum)>0)		
	    Cornerness[wire][timebin]=(MatrixAAsum*MatrixBBsum-pow(MatrixCCsum,2))/(MatrixAAsum+MatrixBBsum);
	  else
	    Cornerness[wire][timebin]=0;
	  
	  if(Cornerness[wire][timebin]>0){	  
	    for(unsigned int i=0;i < hit.size(); i++){
	      channel=hit[i]->Wire()->RawDigit()->Channel();
	      geom->ChannelToWire(channel,tpc,plane,wire2);	 
	      //make sure the end point candidate coincides with an actual hit.
	      if(wire==wire2 && hit[i]->StartTime()<timebin*(numbertimesamples/fTimeBins) 
		 && hit[i]->EndTime()>timebin*(numbertimesamples/fTimeBins)){ 	       
		//this index keeps track of the hit number
		hit_loc[wire][timebin]=i;
		Cornerness2.push_back(Cornerness[wire][timebin]);
		break;
	      } 	        
	    }	     
	  }	    
	}      
    }  // end wire loop 
      std::sort(Cornerness2.rbegin(),Cornerness2.rend());

      for(int vertexnum=0;vertexnum<fMaxCorners;vertexnum++){
	flag=0;
	for(unsigned int wire=0;wire < numberwires && flag==0; wire++)
	  for(int timebin=0;timebin < fTimeBins && flag==0; timebin++){    
	    if(Cornerness2.size()>(unsigned int)vertexnum)
	      if(Cornerness[wire][timebin]==Cornerness2[vertexnum] 
		 && Cornerness[wire][timebin]>0. 
		 && hit_loc[wire][timebin]>-1){
		flag++;      
		//thresholding
		if(Cornerness2.size())
		  if(Cornerness[wire][timebin]<(fThreshold*Cornerness2[0]))
		    vertexnum=fMaxCorners;
		vHits.push_back(hit[hit_loc[wire][timebin]]);
		recob::EndPoint2D endpoint(vHits);
		endpoint.SetWireNum(wire);
		endpoint.SetDriftTime(hit[hit_loc[wire][timebin]]->PeakTime());
		//weak vertices are given endpoint id=1
		endpoint.SetID(1);
		endpoint.SetStrength(Cornerness[wire][timebin]);	  
		vtxcol.push_back(endpoint);
		vHits.clear();
		// non-maximal suppression on a square window. The wire coordinate units are 
		// converted to time ticks so that the window is truly square. 
		// Note that there are 1/0.0743=13.46 time samples per 4.0 mm (wire pitch in ArgoNeuT), 
		// assuming a 1.5 mm/us drift velocity for a 500 V/cm E-field 
		
   
	      //introduce detector independent version - note: change .198 to detector ind. value once this it's sorted out.
		double drifttick=larp->DriftVelocity(larp->Efield(),larp->Temperature())*(.198);
		double wirepitch=geom->WirePitch(0,1,0);
		double corrfactor=drifttick/wirepitch;
				
		for(int wireout=(int)wire-(int)((fWindow*(numbertimesamples/fTimeBins)*corrfactor)+.5);
		    wireout <= (int)wire+(int)((fWindow*(numbertimesamples/fTimeBins)*corrfactor)+.5) ; wireout++)
		  for(int timebinout=timebin-fWindow;timebinout <= timebin+fWindow; timebinout++)
		    if(sqrt(pow(wire-wireout,2)+pow(timebin-timebinout,2))<fWindow)//circular window 
		      Cornerness[wireout][timebinout]=0;	  
	      }     
	  }
      }
      Cornerness2.clear();
      hit.clear();
      if(clusterIter!=clusIn.end()) clusterIter++;

      if((int)p==fSaveVertexMap){ 
	unsigned char *outPix = new unsigned char [fTimeBins*numberwires];
	//finds the maximum cell in the map for image scaling
	int cell, pix=0, maxCell=0;
	//int xmaxx, ymaxx;
	for (int y=0; y<fTimeBins; y++)
	  for (unsigned int x=0; x<numberwires; x++){
	    cell = (int)(hit_map[x][y]*1000);
	    if (cell > maxCell){
	      maxCell = cell;
	    //  xmaxx=x;
	    //  ymaxx=y;
	    }
	  }
	
	for (int y=0; y<fTimeBins; y++)
	  for (unsigned int x=0; x<numberwires; x++){ 
	    //scales the pixel weights based on the maximum cell value     
	    if(maxCell>0)
	      pix = (int)((1500000*hit_map[x][y])/maxCell);
	      outPix[y*numberwires + x] = pix;
	  }
	  // add 3x3 pixel squares to with the harris vertex finders to the .bmp file

	  for(unsigned int ii=0;ii<vtxcol.size();ii++)
	      {
		if(vtxcol[ii].View()==(unsigned int)view)
		{
		   pix = (int)(255);
		   outPix[(int)(vtxcol[ii].DriftTime()*(fTimeBins/numbertimesamples))*numberwires + vtxcol[ii].WireNum()] = pix;
		   outPix[(int)(vtxcol[ii].DriftTime()*(fTimeBins/numbertimesamples))*numberwires + vtxcol[ii].WireNum()-1] = pix;
		   outPix[(int)(vtxcol[ii].DriftTime()*(fTimeBins/numbertimesamples))*numberwires + vtxcol[ii].WireNum()+1] = pix;
		   outPix[(int)((vtxcol[ii].DriftTime()*(fTimeBins/numbertimesamples))-1)*numberwires + vtxcol[ii].WireNum()] = pix;
		   outPix[(int)((vtxcol[ii].DriftTime()*(fTimeBins/numbertimesamples))+1)*numberwires + vtxcol[ii].WireNum()] = pix;
		   outPix[(int)((vtxcol[ii].DriftTime()*(fTimeBins/numbertimesamples))-1)*numberwires + vtxcol[ii].WireNum()-1] = pix;
		   outPix[(int)((vtxcol[ii].DriftTime()*(fTimeBins/numbertimesamples))-1)*numberwires + vtxcol[ii].WireNum()-2] = pix;
		   outPix[(int)((vtxcol[ii].DriftTime()*(fTimeBins/numbertimesamples))+1)*numberwires + vtxcol[ii].WireNum()+1] = pix;
		   outPix[(int)((vtxcol[ii].DriftTime()*(fTimeBins/numbertimesamples))+1)*numberwires + vtxcol[ii].WireNum()+2] = pix;
		}
	      }
	     
	VSSaveBMPFile(Form("harrisvertexmap_%d_%d.bmp",(*clusIn.begin())->ID(),p), outPix, numberwires, fTimeBins);
	delete [] outPix;
      }   
    } // end loop over planes
  }//end loop over tpcs
    
  return vtxcol.size();   
}

