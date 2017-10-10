//Create a 2-D histogram from an image.
//Author: Olivier Couet

#include <vector>
#include "commonUtility.h"
#include <TASImage.h>
#include <TEllipse.h>
double zeroInt = 1;

bool isDone ( TH2D* hh=0) ;
double getWgtCenter( vector<int>& x, vector<double>& e) ;
//void findNeighbors(TH2F* hInd, vector<int>& px, vector<int>& py);
double getSum( vector<double>& e) ;

void fiberCounter(int num=1)
{
  short seedThr =100;
  short bkgThr = 110;
  const int searchRange = 50;
  //  float rad = 10;
  
  const int maxX = 5000;
  const int maxY = 5000;
  gStyle->SetPalette(52);
  
//  TString infName = Form("inputPics/anablesPics/Oct3_%d.JPG",num);
  TString infName = "inputPics/anablesPics/both_ends/DBN_61-BL_22-WG.JPG";
//  TString infName = "inputPics/anablesPics/both_ends/DBN_61-BL_22-NG.JPG";
  TASImage image(infName);
  //  TASImage image("/Users/yongsunkim/uiucAnalysis/emcal/inputPics/anablesPics/100_0183_trimmed_small-1.JPG");
  
  UInt_t yPixels = image.GetHeight();
  UInt_t xPixels = image.GetWidth();
  cout << " xPixels = " << xPixels<<endl;
  cout << " yPixels = " << yPixels<<endl;
  
  if ( (maxX<xPixels) || ( maxY <yPixels) )   {
    cout <<"  (maxX<xPixels) || ( maxY <yPixels) ! " << endl;
    return;
  }
  
  TH1D* hNfib = new TH1D("hNfib",";number of clusters;Entries",4000,0,4000);


  UInt_t *argb   = image.GetArgbArray();
  
  TH2D* hOriginal = new TH2D("hOriginal","",xPixels,.5,xPixels+1,yPixels+1,.5,yPixels);
  TH2D* h = hOriginal;


  TH1D* h1d = new TH1D("h1d","1D histogram",256,0,256);
  TH1D* henergy = new TH1D(Form("henergy_%d",num),";energy of clusters",300,0,300000);
  TH1D* henergyNorm = new TH1D(Form("henergyNorm_%d",num),";Self-normalized intensity",300,0,3);
  TH1D* hMeanE = (TH1D*)henergy->Clone("hMeanE");
  hMeanE->SetXTitle("<cluster intensity>");

  TH1D* hRMSE = new TH1D("hRMSE",";RMS;Entries",1000,0,100000);
  TH1D* hRMSNorm = new TH1D("hRMSNorm",";RMS normalized by mean;Entries",1000,0,1);
  
  h1d->Sumw2();
  for (int row=0; row<xPixels; ++row) {
    for (int col=0; col<yPixels; ++col) {
      int index = col*xPixels+row;
      short grey = argb[index]&0xff ;
      h->SetBinContent(row+1,yPixels-col,grey);
      h1d->Fill(grey);
    }
  }
  
  TCanvas* c0 = new TCanvas("c0","",1200,600);
  c0->Divide(2,1);
  c0->cd(1);
  h1d->SetAxisRange(0.5,2e6,"Y");
  h1d->Draw();
  jumSun(seedThr,0.5,seedThr,2e6);
  jumSun(bkgThr,0.5,bkgThr,2e6);
  gPad->SetLogy();
  c0->cd(2);
  h->SetAxisRange(0,256,"z");
  h->Draw("colz");
  gPad->SetRightMargin(0.2);
  
  c0->SaveAs("histo-1d.gif");

  
  TH2D* h3 = (TH2D*)h->Clone("h3");
  
  TH2D* h4 = (TH2D*)h->Clone("h4");  
  for ( int ix0=1 ; ix0<=h4->GetNbinsX() ; ix0++) {
    for ( int iy0=1 ; iy0<=h4->GetNbinsY() ; iy0++) {
      h4->SetBinContent(ix0,iy0,0.0001);
    }}
  
  std::vector<int> px;
  std::vector<int> py;
  std::vector<int> pinten;
  
  std::vector<int> vClstE;  
  std::vector<float> vClstX;  
  std::vector<float> vClstY;  


  //  std::vector<short> inten;
     

  
  // Clean up backgrounds ;
  short arrInten[maxX+1][maxY+1];
  short arrFlag[maxX+1][maxY+1];
  int kIn=1;  int kOut=-1;   int kUnde=0;
  int nXbins = h3->GetNbinsX() ; 
  int nYbins = h3->GetNbinsY() ; 

  for ( int ix0=1 ; ix0<=nXbins ; ix0++) {
    for ( int iy0=1 ; iy0<=nYbins ; iy0++) {
      short val0 = h3->GetBinContent(ix0,iy0);
      
      arrInten[ix0][iy0] = val0;
      arrFlag[ix0][iy0] = kUnde;	
      if (val0 < bkgThr )  {
	h3->SetBinContent(ix0, iy0, zeroInt);
	arrInten[ix0][iy0] = 0;
	arrFlag[ix0][iy0] =  kOut;
      }
    }
  }
  
  int nClst = 0;
  
  for ( int ix0=1 ; ix0<=nXbins ; ix0++) {
    for ( int iy0=1 ; iy0<=nYbins ; iy0++) {
      short val0 = arrInten[ix0][iy0];
      
      if ( val0 < seedThr ) 
	continue;
      if ( arrFlag[ix0][iy0] != kUnde ) // if it is already included in other clusters.
	continue;

      cout << " found new seed!   " ;
      cout <<"(x,y,intensity) = " <<  ix0<<", "<<iy0<<", "<<val0<<endl;


      // Found a seed!     (ix0, iy0, val0) are the seed! 
      nClst++;
      
      // Begin Clustering 
      //      inten.clear();
      px.clear();
      py.clear();
      pinten.clear();

      //      inten.push_back(val0);
      px.push_back(ix0);   
      py.push_back(iy0);
      pinten.push_back( val0 ) ;

      arrFlag[ix0][iy0] = kIn;
      
      bool completeFlag = false;
      int nIter = 0;
      int countedHits=0;
      int itrBegin=0;  // <== this number will evolve over the iteration :  itrBegin = itrEnd;
      int itrEnd;

      while ( completeFlag==false)  
	{
	  //	  cout << " px.size() = " << px.size() << endl;
	  nIter++;
	  completeFlag = true;
	  
	  itrEnd = int(px.size());
	  //	  cout << " itrBegin = " << itrBegin ;
	  //	  cout << ", itrEnd = " << itrEnd<<endl;
	  for ( int vi = itrBegin ; vi < itrEnd ; vi++) {
	    
	    int xCand = px[vi];     int yCand = py[vi];      short intenCand = pinten[vi];
	    
	    if ( (xCand < nXbins) && (arrFlag[xCand+1][yCand]==kUnde) ) {	       // right end 
	      completeFlag=false; 
	      countedHits++;
	      px.push_back( xCand+1 ) ;
	      py.push_back( yCand   ) ;
	      pinten.push_back( arrInten[xCand+1][yCand]   ) ;
	      arrFlag[xCand+1][yCand]=kIn; 
	    }
	    
	    if ( (xCand > 1    ) && (arrFlag[xCand-1][yCand]==kUnde) ) {	       // Left end 
	      completeFlag=false; 
	      countedHits++;
	      px.push_back( xCand-1 ) ;
	      py.push_back( yCand   ) ;
	      pinten.push_back( arrInten[xCand-1][yCand]   ) ;
	      arrFlag[xCand-1][yCand]=kIn; 
	    }
	    if ( (yCand < nYbins) && (arrFlag[xCand][yCand+1]==kUnde) ) { 	       // top end 
	      completeFlag=false; 
	      countedHits++;
	      px.push_back( xCand ) ;
	      py.push_back( yCand+1   ) ;
	      pinten.push_back( arrInten[xCand][yCand+1]   ) ;
	      arrFlag[xCand][yCand+1]=kIn; 
	    }
	    if ( (yCand > 1    ) && (arrFlag[xCand][yCand-1]==kUnde) ) {  	       // bottom end 
	      completeFlag=false; 
	      countedHits++;
	      px.push_back( xCand ) ;
	      py.push_back( yCand-1   ) ;
	      pinten.push_back( arrInten[xCand][yCand-1]  ) ;
	      arrFlag[xCand][yCand-1]=kIn; 
	    }
	  }
	  itrBegin = itrEnd;


	  /*	  if ( nClst==2 ) {
            TCanvas* ctemp = new TCanvas("ctemp","", 400,400);
            TH2D* htemp = new TH2D("htemp","",50, px[0]-15, px[0]+35, 50, py[0]-25, py[0]+25);
            for ( int ii=0 ; ii<px.size(); ii++)   {
              int theBin = htemp->FindBin( px[ii], py[ii]) ;
              htemp->SetBinContent( theBin, pinten[ii] );
            }
	    htemp->SetAxisRange(200,256,"Z");
	    htemp->Draw("colz");
	    drawText(Form("Iter. %d",nIter),0.2,0.8);
	    ctemp->SaveAs(Form("clustering_itr%d.png",nIter));
	      
	    }*/





	}
      cout << "number of hits in "<<nClst<<"th cluster: " << px.size() << ",  nIteration = "<<nIter<<endl;  
      int sumEnergy = 0;
      for ( int vi=0 ; vi < px.size() ; vi++) {
	sumEnergy = sumEnergy + pinten[vi];
      }
      cout << " total energy = " << sumEnergy << endl;

      float xmean=0;
      float ymean=0;
      for ( int vi=0 ; vi < px.size() ; vi++) { 
	xmean = xmean + px[vi]*pinten[vi] ;
	ymean = ymean + py[vi]*pinten[vi] ;
      }
      xmean = xmean / sumEnergy ;
      ymean = ymean / sumEnergy ;

      henergy->Fill(sumEnergy);
      vClstE.push_back(sumEnergy);
      vClstX.push_back(xmean);
      vClstY.push_back(ymean);

      //      for ( int vi=0 ; vi < px.size() ; vi++) {
      //	h4->SetBinContent( px[vi],  py[vi], pinten[vi]);
      //      }
      //      TCanvas* c11 = new TCanvas("c11","",400,400);
      //      h4->Draw("colz");
      //      c11->SaveAs(Form("Cluster_%d.png",nClst));
      
    }}
  
  for ( int ci = 0 ; ci< vClstE.size() ; ci++) {
    henergyNorm->Fill ( vClstE[ci] / henergy->GetMean() ) ;
  }
	  
  
  hNfib->Fill(nClst);
  hMeanE->Fill( henergy->GetMean() );
  hRMSE->Fill( henergy->GetRMS() );
  hRMSNorm->Fill(  henergy->GetRMS() / henergy->GetMean() ) ;
  
  TCanvas* c15 = new TCanvas("c15","",600,600);
  h->SetAxisRange(0,256,"z");
  h->Draw("colz");
  gPad->SetRightMargin(0.2);
  for ( int ci = 0 ; ci< vClstX.size() ; ci++) {
    TEllipse *el3 = new TEllipse( vClstX[ci], vClstY[ci], 5);
    el3->SetLineColor(kRed);
    el3->SetLineWidth(2);
    el3->SetFillStyle(0);
    el3->Draw();
  }

  // end of clustering 
  TCanvas* c2 = new TCanvas("c2","",400,400);
  handsomeTH1(henergy,1);
  henergy->SetXTitle("Intensity of each cluster");
  henergy->SetYTitle("Entries");
  henergy->Draw();
  
  TFile* fout = new TFile(Form("%s_outputHistograms.root",infName.Data()),"RECREATE");
  henergy->Write();
  henergyNorm->Write();
  hNfib->Write();
  hMeanE->Write();
  hRMSE->Write();
  hRMSNorm->Write();
  fout->Close();
  
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
