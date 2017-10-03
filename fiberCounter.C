//Create a 2-D histogram from an image.
//Author: Olivier Couet

#include <vector>
#include "commonUtility.h"
#include <TASImage.h>
#include <TEllipse.h>
double zeroInt = -1;

bool isDone ( TH2D* hh=0) ;
double getWgtCenter( vector<int>& x, vector<double>& e) ;
//void findNeighbors(TH2F* hInd, vector<int>& px, vector<int>& py);
double getSum( vector<double>& e) ;

void fiberCounter()
{
  short seedThr = 0.9 * 256;
  short bkgThr = 0.5 * 256;
  const int searchRange = 50;
  float rad = 10;
  
  const int maxX = 3000;
  const int maxY = 3000;
  gStyle->SetPalette(52);
  
  //  TASImage image("inputPics/kodak2_sept11/block2.jpg");
  //  TASImage image("inputPics/anablesPics/100_0183_trimmed.JPG");
  TASImage image("/Users/yongsunkim/uiucAnalysis/emcal/inputPics/anablesPics/100_0183_trimmed_small-1.JPG");
  
  UInt_t yPixels = image.GetHeight();
  UInt_t xPixels = image.GetWidth();
  cout << " xPixels = " << xPixels<<endl;
  cout << " yPixels = " << yPixels<<endl;
  
  if ( (maxX<xPixels) || ( maxY <yPixels) )   {
    cout <<"  (maxX<xPixels) || ( maxY <yPixels) ! " << endl;
    return;
  }
  UInt_t *argb   = image.GetArgbArray();
  
  TH2D* hOriginal = new TH2D("hOriginal","",xPixels,.5,xPixels+1,yPixels+1,.5,yPixels);
  TH2D* h = hOriginal;


  TH1D* h1d = new TH1D("h1d","1D histogram",256,0,256);
  TH1D* henergy = new TH1D("henergy",";energy of clusters",300,0,300);
  
  h1d->Sumw2();
  for (int row=0; row<xPixels; ++row) {
    for (int col=0; col<yPixels; ++col) {
      int index = col*xPixels+row;
      short grey = argb[index]&0xff ;
      h->SetBinContent(row+1,yPixels-col,grey);
      h1d->Fill(grey);
    }
  }
  
  TCanvas* c0 = new TCanvas("c0","",800,400);
  c0->Divide(2,1);
  c0->cd(1);
  h1d->Draw();
  gPad->SetLogy();
  c0->cd(2);
  h->SetAxisRange(0,256,"z");
  h->Draw("colz");
  gPad->SetRightMargin(0.2);
  
  c0->SaveAs("histo-1d.gif");

  
  TH2D* h3 = (TH2D*)h->Clone("h3");
  
  std::vector<int> px;
  std::vector<int> py;
  //  std::vector<short> inten;
     
  TCanvas * ctemp = new TCanvas("ctemp","",500,500);
  
  // Clean up backgrounds ;
  short arrInten[maxX+1][maxY+1];
  cout << "here 1 " << endl;
  short arrFlag[maxX+1][maxY+1];
  int kIn=1;  int kOut=-1;   int kUnde=0;
  int nXbins = h3->GetNbinsX() ; 
  int nYbins = h3->GetNbinsY() ; 

  for ( int ix0=1 ; ix0<=nXbins ; ix0++) {
    for ( int iy0=1 ; iy0<=nYbins ; iy0++) {
      short val0 = h3->GetBinContent(ix0,iy0);
      
      arrInten[ix0][iy0] = val0;
  cout << "here 2 " << endl;
      arrFlag[ix0][iy0] = kUnde;	
      if (val0 < bkgThr )  {
	h3->SetBinContent(ix0, iy0, zeroInt);
	arrInten[ix0][iy0] = 0;
	cout << "here 3 " << endl;
	arrFlag[ix0][iy0] =  kOut;
      }
    }
  }
  
  int nClst = 0;
  
  TCanvas* c2 = new TCanvas("c2","",400,400);
  for ( int ix0=1 ; ix0<=nXbins ; ix0++) {
    for ( int iy0=1 ; iy0<=nYbins ; iy0++) {
      short val0 = arrInten[ix0][iy0];
      
      if ( val0 < seedThr ) 
	continue;
      
      // Found a seed!     (ix0, iy0, val0) are the seed! 
      nClst++;
      
      // Begin Clustering 
      //      inten.clear();
      px.clear();
      py.clear();
      
      //      inten.push_back(val0);
      px.push_back(ix0);   
      py.push_back(iy0);
  cout << "here 4 " << endl;
      arrFlag[ix0][iy0] = kIn;
      
      bool completeFlag = false;
      int nIter = 0;
      while ( completeFlag==false)  
	 {
	   nIter++;
	   completeFlag = true;
	   for ( int vi = 0 ; vi<= (int)(px.size()) ; vi++) 
	     {
	       int xCand = px[vi];     int yCand = py[vi]; 
	       //short intenCand = arrInten[xCand][yCand];   // intensity below background threshold is already out 
	       
	         cout << "here 5 " << endl;

	       if ( (xCand < nXbins) && (arrFlag[xCand+1][yCand]==kUnde) ) {	       // right end 
		 completeFlag=false; 
		 px.push_back( xCand+1 ) ;
		 py.push_back( yCand   ) ;
		 arrFlag[xCand+1][yCand]=kIn; 
	       }
	       cout << "here 6 " << endl;
  
	       if ( (xCand > 1    ) && (arrFlag[xCand-1][yCand]==kUnde) ) {	       // Left end 
		 completeFlag=false; 
		 px.push_back( xCand-1 ) ;
		 py.push_back( yCand   ) ;
		 arrFlag[xCand-1][yCand]=kIn; 
	       }
	       cout << "here 7 " << endl;
	       if ( (yCand < nYbins) && (arrFlag[xCand][yCand+1]==kUnde) ) { 	       // top end 
		 completeFlag=false; 
		 px.push_back( xCand ) ;
		 py.push_back( yCand+1   ) ;
		 arrFlag[xCand][yCand+1]=kIn; 
	       }
	       cout << "here 8 " << endl;
	       if ( (yCand > 1    ) && (arrFlag[xCand][yCand-1]==kUnde) ) {  	       // bottom end 
		 completeFlag=false; 
		 px.push_back( xCand ) ;
		 py.push_back( yCand-1   ) ;
		 arrFlag[xCand][yCand-1]=kIn; 
	       }
	       cout << "here 9 " << endl;
	     }
	   cout << "number of hits in "<<nClst<<"th cluster: " << px.size() << ",  nIteration = "<<nIter<<endl;  
	 }
      
      return;
    }}
  // end of clustering 
  
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
