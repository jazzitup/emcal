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

void fiberCounter()
{
  double seedThr = 0.5;
  double bkgThr = 0.3;
  int searchRange = 30;
  float rad = 6;

  gStyle->SetPalette(52);


  //  TASImage image("inputPics/picTakenAtNPL_piece1.png");
  //  TASImage image("inputPics/picTakenAtNPL.png");
  //  TASImage image("inputPics/exmample_defacts.png");
  //  TASImage image("inputPics/exmample_goodClusters.png");
  //  TASImage image("inputPics/pic_aug22_samll.png");
  //  TASImage image("inputPics/pic_aug22.png");
  TASImage image("inputPics/kodak/timmed_100_0015.png");

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
   



   TCanvas* c0 = new TCanvas("c0","",400,400);
   //   h1d->SetAxisRange(0.5,1,"X");
   h1d->Draw();
   gPad->SetLogy();
   c0->SaveAs("histo-1d.gif");

   TCanvas* c5 = new TCanvas("c5","",400,400);
   hzero->Reset();
   for ( int ix0=1 ; ix0<=h->GetNbinsX() ; ix0++) {
     for ( int iy0=1 ; iy0<=h->GetNbinsY() ; iy0++) {
       double val0 = h->GetBinContent(ix0,iy0);
       //       if (val0 == 0 )
       //       if (val0 < 0.05 ) 
       if ( (val0 >0.1 ) && ( val0 < 0.2) )
         hzero->SetBinContent(ix0, iy0, 1);
     }}

   hzero->Draw("colz");




   TCanvas* c1 = new TCanvas("c1","",800,400);
   c1->Divide(2,1);
   c1->cd(1);
   h->SetAxisRange(0,1,"z");
   //   h->SetAxisRange(300,400,"X");    h->SetAxisRange(300,400,"Y");
   h->Draw("colz");
   gPad->SetRightMargin(0.2);
   
   c1->cd(2);
   TH2D* h2 = (TH2D*)h->Clone("h2");
   h2->Reset();
   for (int row=0; row<xPixels; ++row) {
     for (int col=0; col<yPixels; ++col) {
       double colVal = h->GetBinContent(row+1,yPixels-col);
       if ((colVal >=.2) && ( colVal<1.1))
	 h2->SetBinContent(row+1,yPixels-col,1);
      }
   }

  
   cout << "total number of pixels = " << (xPixels+1)*(yPixels+1) << endl;
   //   gStyle->SetPalette(53);
   h2->SetAxisRange(0,1,"z");
   //   h2->DrawCopy("colz");

   c1->SaveAs("histo-2d.gif");
   //   return;

   TH2D* h3 = (TH2D*)h->Clone("h3");
   
   int px0; 
   int py0;
   double inten0;
   std::vector<int> px;
   std::vector<int> py;
   std::vector<double> inten;

   std::vector<int> wgtx;
   std::vector<int> wgty;
   std::vector<double> energy;


   
   // Remove backgrounds ;
   for ( int ix0=1 ; ix0<=h3->GetNbinsX() ; ix0++) {
     for ( int iy0=1 ; iy0<=h3->GetNbinsY() ; iy0++) {

       double val0 = h3->GetBinContent(ix0,iy0);
       if (val0 < bkgThr ) 
	 h3->SetBinContent(ix0, iy0, zeroInt);
     }}
   

   //   return;
   
   int nClst = 0;
   
   TCanvas* c2 = new TCanvas("c2","",400,400);
   
   for ( int ix0=1 ; ix0<=h3->GetNbinsX() ; ix0++) {
     for ( int iy0=1 ; iy0<=h3->GetNbinsY() ; iy0++) {
       
       double val0 = h3->GetBinContent(ix0,iy0);
       /*       if (val0 == zeroInt ) 
       continue;
       */
       if ( val0 < seedThr ) 
	 continue;
       

       // Found the quasi-seed
       px.clear();   py.clear();   inten.clear();
       
       px.push_back(ix0); 
       py.push_back(iy0);
       inten.push_back(val0);
       h3->SetBinContent(ix0,iy0,zeroInt);
       
       // Iterate NxN around it.  
       int nIteration =  3;
       for ( int iter = 1 ; iter<= nIteration ; iter++) {
	 for ( int ix= ix0 - searchRange ; ix<= ix0 + searchRange ; ix++) {
	   for ( int iy= iy0 - searchRange ; iy<= iy0 + searchRange ; iy++) {
	     
	     if (  ( ix < 1 ) ||  ( iy < 1 ) || ( ix > h3->GetNbinsX() ) ||  ( iy > h3->GetNbinsY() ) ) 
	       continue;

	     double val = h3->GetBinContent(ix,iy);
	     if (val  == zeroInt )
	       continue;
	     if ( isNeighbor(px, py, ix,iy) )  {
	       px.push_back(ix);
	       py.push_back(iy);
	       inten.push_back(val);
	       h3->SetBinContent(ix,iy,zeroInt);
	     }
	     
	   }} // end of second iter
       }  // Repeat by nIteration times
       
       cout << "Iteration " << nClst << "size of the cluster = " << inten.size() << endl;

       wgtx.push_back(getWgtCenter( px, inten) );
       wgty.push_back(getWgtCenter( py, inten) );
       energy.push_back( getSum(inten) );
       henergy->Fill( getSum(inten) ) ;

       //       cout << " x,y,totE = " << wgtx[nClst] << ", " <<  wgty[nClst] <<", " <<  energy[nClst] <<", " << endl;
       nClst++;
       

       //       if (nClst == 5)        	 return;

       h3->SetAxisRange(0,2,"z");
       h3->Draw("colz");
       h3->SetXTitle(Form("Fiber counts : %d",nClst));

       //       c2->SaveAs(Form("figs/iter_%d.png",nClst));
       px.clear();   py.clear();   inten.clear();
       //       return;
     }} // end of first iter
   
   
   cout << "number of clusters = " << nClst << endl;
   
   
   TCanvas* c3 = new TCanvas("c3","",400,400);
   h->Draw("colz");
   cout << "n = " << wgtx.size() << endl;
   for ( int ii = 0 ; ii<wgtx.size() ; ii++) { 
     TEllipse *el3 = new TEllipse( wgtx[ii], wgty[ii], rad);
     el3->SetLineColor(kYellow);
     el3->SetLineWidth(2);
     el3->SetFillStyle(0);
     el3->Draw();
   }
   //   c3->Update();
   c3->SaveAs("result.gif");

   TCanvas* c4 = new TCanvas("c4","",400,400);
   handsomeTH1(henergy);
   henergy->Draw();
   c4->SaveAs("energyDist.gif");

   TFile* hout = new TFile("result.root","recreate");
   h->Write();
   h1d->Write();
   henergy->Write();
   hout->Close();
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
