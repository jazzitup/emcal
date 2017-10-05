//Create a 2-D histogram from an image.
//Author: Olivier Couet

#include <vector>
#include "commonUtility.h"
#include <TASImage.h>
#include <TEllipse.h>
double zeroInt = -1;

bool isDone ( TH2D* hh=0) ;
bool isNeighbor( vector<int>& px, vector<int>& py, int xx=0, int yy=0) ;
double getWgtCenter( vector<int>& x, vector<double>& e) ;
double getSum( vector<double>& e) ;

void defectFinder()
{
  double seedThr = 0.6;
  double bkgThr = 0.3;
  int searchRange = 50;
  float rad = 6;

  gStyle->SetPalette(52);


  TASImage image("inputPics/defectExample.png");

  
  UInt_t yPixels = image.GetHeight();
  UInt_t xPixels = image.GetWidth();
  UInt_t *argb   = image.GetArgbArray();

  TH2D* h = new TH2D("h","",xPixels,.5,xPixels+1,yPixels+1,.5,yPixels);
  TH1D* h1d = new TH1D("h1d","1D histogram",256,0,1);
  TH1D* henergy = new TH1D("henergy",";energy of clusters",300,0,300);
  TH2D* hzero = (TH2D*)h->Clone("hzero");
  
   h1d->Sumw2();
   for (int row=0; row<xPixels; ++row) {
     for (int col=0; col<yPixels; ++col) {
       int index = col*xPixels+row;
	//	 cout << "argb[index] = " << argb[index] << endl;
         float grey = float(argb[index]&0xff)/256;
         h->SetBinContent(row+1,yPixels-col,grey);
	 h1d->Fill(grey);
      }
   }
   
   TCanvas* c1 = new TCanvas("c1","",400,400);
   h->Draw("colz");
   TCanvas* c2 = new TCanvas("c2","",400,400);
   h1d->Draw();


   

}

bool isDone ( TH2D* hh ) { 
  bool flag = true;
  for ( int ix=1 ; ix<=hh->GetNbinsX() ; ix++) {
    for ( int iy=1 ; iy<=hh->GetNbinsY() ; iy++) {
      if ( hh->GetBinContent(ix,iy) != zeroInt ) {
	flag = false; 
	break;
      }
    }
  }
  return flag;
}

bool isNeighbor( vector<int>& px, vector<int>& py, int xx, int yy) {
  bool flag = false;
  for ( int ix = 0 ; ix< (int) px.size() ; ix++) { 
    for ( int iy = 0 ; iy< (int) py.size() ; iy++) { 
      if (  ( px[ix] == xx)  && (  abs( py[iy] - yy ) <= 1 ) )	{ 
	//	cout << "found! :  " <<  px[ix] << ", " << xx << " and " << py[iy] << ", " << yy << endl;
	return true;
      }
      if (  (py[iy]==yy)  && (  abs( px[ix] - xx ) <= 1 ) )  {
	//	cout << "found! :  " <<  px[ix] << ", " << xx << " and " << py[iy] << ", " << yy << endl;
	return true;
      }
      
    }}
  //  cout << "No they are not neighbor" << endl;
  return false;
}
double getWgtCenter( vector<int>& x, vector<double>& e) { 
  double sum=0;
  double totEnergy = 0 ;
  for ( int ii = 0 ; ii < x.size() ; ii++) {
    sum = sum + x[ii]*e[ii];
    totEnergy = totEnergy + e[ii] ;
  }
  return sum/totEnergy;
}

double getSum( vector<double>& e) { 
  double sum = 0;
  for ( int ii = 0 ; ii < e.size() ; ii++) {
    sum = sum + e[ii];
  }
  
  return sum;
}
