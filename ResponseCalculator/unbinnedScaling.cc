#include<iomanip> // provides setprecision
#include<iostream>
#include<sstream>
#include<string>

// ROOT
#include<TCanvas.h>
#include<TEfficiency.h>
#include<TFile.h>
#include<TGraphAsymmErrors.h>
#include<TH3F.h>
#include<TLine.h>
#include<TMath.h>
#include<TProfile2D.h>
#include<TROOT.h>
#include<TChain.h>

// user incuded files
#include "Style.h"

using namespace std;


TH3F closure3d( TChain& tree, const TH3F& scales3d ) {
  auto closure = *((TH3F*)scales3d.Clone());
  closure.Reset();

  float r, eta, e;
  tree.SetBranchAddress("r",&r);
  tree.SetBranchAddress("e",&e);
  tree.SetBranchAddress("eta",&eta);

  for( int i=0; i<tree.GetEntries(); i++ ) {
    tree.GetEntry(i);
    if( r <= 0.85 ) continue;
    auto bin = scales3d.FindFixBin( e, eta, r );
    auto scale = scales3d.GetBinContent( bin );
    closure.Fill( e, eta, r*scale );
  }

  return closure;
}

void drawClosure( const TH3F& fullh3, const TH3F& fasth3, const TH3F& modih3 ) {
  gStyle->SetOptStat(0);
  for( int xbin=1; xbin< fasth3.GetNbinsX()+1; ++xbin ) { // E_gen
    for( int ybin=1; ybin< fasth3.GetNbinsY()+1; ++ybin ) { // eta_gen
      // Create the 1d histograms
      auto h1_full = fullh3.ProjectionZ( "full", xbin, xbin, ybin, ybin );
      auto h1_fast = fasth3.ProjectionZ( "fast", xbin, xbin, ybin, ybin );
      auto h1_modi = modih3.ProjectionZ( "mod", xbin, xbin, ybin, ybin );
      h1_full->SetLineColor(1);
      h1_fast->SetLineColor(2);
      h1_modi->SetLineColor( kBlue );
      h1_modi->SetLineWidth(2);

      // scale to unity for drawing
      h1_full->Scale( 1./h1_full->GetEntries() );
      h1_fast->Scale( 1./h1_fast->GetEntries() );
      h1_modi->Scale( 1./h1_modi->GetEntries() );

      // set minimum xaxis
      //h1_modi->GetXaxis()->SetRangeUser( 0.8, 1.05 );

      // Rebin
      //h1_modi->Rebin(40);
      //h1_full->Rebin(40);
      //h1_fast->Rebin(40);

      h1_modi->Draw();
      h1_full->Draw("same");
      h1_fast->Draw("same");
      gPad->SaveAs( (std::string("plots/checker_")+to_string(xbin)+"vs"+to_string(ybin)+".pdf").c_str() );

    }
  }
}

TH3F calculateResponseFromTree( TChain& fasttree, TChain& fulltree, TH3F h ) {
  float rfast, efast,etafast;
  fasttree.SetBranchAddress("r",&rfast);
  fasttree.SetBranchAddress("e",&efast);
  fasttree.SetBranchAddress("eta",&etafast);

  float rfull;
  fulltree.SetBranchAddress("r",&rfull);

  for( int i=0; i<100; i++ ) {
    fasttree.GetEntry(i);
    for( int j=0; j<100; j++ ) {
      fulltree.GetEntry(j);
      h.Fill( efast, etafast, rfull/rfast );
    }
  }

  return h;

}

std::vector<float> getVectorFromTree( TChain& tree, bool sort=true ) {
  std::vector<float> vec;
  vec.reserve( tree.GetEntries() );
  float r;
  tree.SetBranchAddress("r",&r);
  //for( int i=0; i<100000; i++ ) {
  for( int i=0; i<tree.GetEntries(); i++ ) {
    tree.GetEntry(i);
    if( r<0.3) continue;
    vec.push_back( r );
  }
  if( sort ) {
    std::sort( vec.begin(), vec.end() );
  }

  return vec;
}



TH3F calculateResponseUnbinned( TChain& fasttree, TChain& fulltree, TH3F h ) {

  auto fastvector = getVectorFromTree( fasttree );
  auto fullvector = getVectorFromTree( fulltree );

  printf( "fast %f to %f\n", fastvector[0], fastvector[fastvector.size()-1] );
  printf( "full %f to %f\n", fullvector[0], fullvector[fullvector.size()-1] );

  TProfile profile("profile", "title", h.GetZaxis()->GetNbins(), h.GetZaxis()->GetXmin(), h.GetZaxis()->GetXmax() );
  for( unsigned i=0; i<fastvector.size(); i++ ) {
    auto fa = fastvector[i];
    auto jRel = 1.*i/fastvector.size()*fullvector.size();
    int j = (int) jRel;
    auto fu = fullvector[j]*(jRel-j)+fullvector[j+1]*(j-jRel+1);
    profile.Fill( fa, fu/fa );
  }

  for( int i=0; i<profile.GetNbinsX()+2;i++ ) {
    h.SetBinContent( 1, 1, i, profile.GetBinContent(i) );
  }


  profile.SetMinimum( 0.95 );
  profile.SetMaximum( 1.05 );
  profile.Draw();
  gPad->SaveAs("scaleProfile.pdf");

  return h;

}

int main( int argc, char** argv ) {
//  string fastname = "../../CMSSW/CMSSW_7_3_0/src/Analyzer/ECALScaleFactorCalculator/3d_fast.root";
//  string fullname = "../../CMSSW/CMSSW_7_3_0/src/Analyzer/ECALScaleFactorCalculator/3d_full.root";
  string fastname = "../3d_fast.root";
  string fullname = "../3d_full.root";


  TChain fasttree("ecalScaleFactorCalculator/responseTree");
  fasttree.AddFile( fastname.c_str() );
  TChain fulltree("ecalScaleFactorCalculator/responseTree");
  fulltree.AddFile( fullname.c_str() );

  //TH3F h3default("responseVsEVsEta", ";E_{gen};#eta_{gen};E/E_{gen}", 100, 5, 1005, 400, 0, 3.2, 2000, 0, 1.05 );
  TH3F h3default("responseVsEVsEta", ";E_{gen};#eta_{gen};E/E_{gen}", 1, -1, 1e6, 1, -1, 5, 1000, 0.8, 1.01 );

  auto fasth3 = (TH3F*) h3default.Clone("fasth3");
  auto fullh3 = (TH3F*) h3default.Clone("fullh3");

  fasttree.Draw("r:eta:e>>fasth3", "1");
  fulltree.Draw("r:eta:e>>fullh3", "1");

  cout << "Calculate scale" << endl;
  auto scales3d = calculateResponseUnbinned( fasttree, fulltree, h3default );

  cout << "Apply scale" << endl;
  auto closureh3d = closure3d( fasttree, scales3d );

  drawClosure( *fullh3, *fasth3, closureh3d );


}





































