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
#include<TVector3.h>

// user incuded files
#include "Style.h"

using namespace std;

TChain* getChain( std::string const & filename, std::string const & treename ) {
    TChain* ch = new TChain( treename.c_str() );
    ch->AddFile( filename.c_str() );
    return ch;
}


template <class HIST>
HIST getHist( std::string const & filename, std::string const & histname ) {
    // Reads a histogram from a file

    TFile file( filename.c_str() );
    if( file.IsZombie() ) {
        std::cerr << "ERROR: Could not open file " << filename << std::endl;
        exit(1);
    }

    TObject* obj = file.Get( histname.c_str() );
    if( obj == 0 ) {
        std::cerr << "ERROR: Could not extract " << histname
            << " histogram from " << filename << std::endl;
        exit(1);
    }
    // Before closing the file, the histogram has to be written on memory
    gROOT->cd();
    HIST hist( *((HIST*)obj) );
    file.Close();

    return hist;
}

TGraphAsymmErrors modifyScale( const TGraphAsymmErrors& origScale, double mean ) {
    // Resplaces the points with too large uncertainty with the mean and removes all uncertainties.
    // This is done to prevent statistical fluctuations.

    auto modScale = (TGraphAsymmErrors*) origScale.Clone();

    for( auto i=0; i<modScale->GetN(); i++) {
        double x, y;
        modScale->GetPoint(i, x, y );
        double errorUp = modScale->GetErrorYhigh(i);
        double errorDn = modScale->GetErrorYlow(i);
        if(
            (errorUp+errorDn)/2 > std::max( 0.05, std::abs( mean - 1 ) )// uncertainty larger than correction
        ) {
            modScale->SetPoint( i, x, mean );
        }

        // The uncertainties will not be used and are therefore set to 0
        modScale->SetPointEYhigh( i, 0 );
        modScale->SetPointEYlow ( i, 0 );
    }
    return *modScale;
}


TH1D getSimplifiedScale( const TH1D& h1_fast, const TH1D& h1_full ) {
    auto h1_scale = *((TH1D*)h1_fast.Clone());

    float intFast = h1_fast.Integral();
    float intFull = h1_full.Integral();

    for( auto binFast=h1_fast.GetNbinsX()+1; binFast>0; binFast-- ) {
        auto fastInt = h1_fast.Integral( binFast, -1 )/intFast;
        int binFull;
        for( binFull=h1_fast.GetNbinsX()+1; binFull>0; binFull-- ) {
            auto fullInt = h1_full.Integral( binFull, -1 )/intFull;
            if( fullInt > fastInt ) break;
        }
        h1_scale.SetBinContent( binFast, h1_fast.GetXaxis()->GetBinCenter(binFull)/h1_fast.GetXaxis()->GetBinCenter(binFast) );

    }



    return h1_scale;
}

TGraphAsymmErrors getScaleWithUncertainties( const TH1D& h1_fast, const TH1D& h1_full ) {
    /* For each fastsim energy, calculate the are from -inf to the energy.
     * Search then the fullsim energy, which corresponds to the same area.
     * The statistical uncertanity of fullsim, fastsim and the binning uncertainty is taken into account.
     */

    // This object will be returned
    TGraphAsymmErrors out = TGraphAsymmErrors();
    // The x-title will be inherited from the input histo, the y-title is newly set
    out.SetTitle( (std::string(";")+h1_fast.GetXaxis()->GetTitle()+";Scale    ").c_str() );

    // Unweighted histograms are assumed
    // A cast to int is done, to bypass numerical (un)precission
    if( round(h1_fast.GetEntries()) != round(h1_fast.GetEffectiveEntries()) ||
            round(h1_full.GetEntries()) != round(h1_full.GetEffectiveEntries()) ) {
        std::cerr << "Please provide unweighted histograms" << std::endl;
        return out;
    }

    // Calculate confidence level ( approx 0.683 )
    double alpha = TMath::Erf( 1./TMath::Sqrt2() );

    // Used for calculating uncertainties
    TEfficiency eff = TEfficiency();

    int entriesFast = h1_fast.GetEntries();
    int entriesFull = h1_full.GetEntries();

    int nBinsFast = h1_fast.GetNbinsX()+2;
    int nBinsFull = h1_fast.GetNbinsX()+2;

    // To improve acess time, the entries of h1_full will be saved as vector
    std::vector<int> cumulativeFull;
    cumulativeFull.reserve( nBinsFull );
    cumulativeFull.push_back( int(h1_full.GetBinContent(0)) );
    for( int i=1; i<nBinsFull; ++i ) {
        cumulativeFull.push_back( cumulativeFull.back() + int(h1_full.GetBinContent(i)) );
    }

    int summedEntriesFast = 0;
    for( int binFast=1; binFast<nBinsFast; ++binFast ) {
        summedEntriesFast += h1_fast.GetBinContent( binFast );
        //double areaFast = summedEntriesFast/entriesFast; // not needed, since the mean is calculated as (up+down)/2

        // Take into account the statistical precission of h1
        double areaFastUp = eff.ClopperPearson( entriesFast, summedEntriesFast, alpha, true );
        double areaFastDn = eff.ClopperPearson( entriesFast, summedEntriesFast, alpha, false );

        // Number of entries in the fullsim corresponding to these bondaries
        int entriesFull_StatFastDn = floor( areaFastDn * entriesFull );
        int entriesFull_StatFastUp = ceil ( areaFastUp * entriesFull );


        // Find the bins in hFull, which corresponds to the same area from hFast
        // This is the propagation of the statistical uncertainty of hFast to hFull
        int binFull_StatFastDn = -1;
        int binFull_StatFastUp = -1;
        for( int binFull=nBinsFull-1; 0<=binFull; --binFull ) {
            if( cumulativeFull.at(binFull) >= entriesFull_StatFastDn ) {
                binFull_StatFastDn = binFull;
            }
        }
        for( int binFull=0; binFull<nBinsFull-1; ++binFull ) {
            if( cumulativeFull.at(binFull) <= entriesFull_StatFastUp ) {
                binFull_StatFastUp = binFull;
            }
        }

        // Take the middle. This is somehow arbitrary, but feel free to envolve a better method
        int binFullMiddle = (binFull_StatFastUp + binFull_StatFastDn)/2;

        // Compute the statistical precission of hFull
        double areaFullUp = eff.ClopperPearson( entriesFull, cumulativeFull.at(binFullMiddle), alpha, true );
        double areaFullDn = eff.ClopperPearson( entriesFull, cumulativeFull.at(binFullMiddle), alpha, false );

        // Number of entries in the fullsim corresponding to these bondaries
        int entriesFull_StatFullDn = floor( areaFullDn * entriesFull );
        int entriesFull_StatFullUp = ceil ( areaFullUp * entriesFull );
        int binFull_StatFullDn = -1;
        int binFull_StatFullUp = -1;
        for( int binFull=nBinsFull-1; 0<=binFull; --binFull ) {
            if( cumulativeFull.at(binFull) >= entriesFull_StatFullDn ) {
                binFull_StatFullDn = binFull;
            }
        }
        for( int binFull=0; binFull<nBinsFull-1; ++binFull ) {
            if( cumulativeFull.at(binFull) <= entriesFull_StatFullUp ) {
                binFull_StatFullUp = binFull;
            }
        }

        // Calculate the scale
        double efull = h1_fast.GetBinCenter( binFullMiddle );
        double efast = h1_fast.GetBinCenter( binFast );
        double scale = efull / efast;

        // add binning uncertanity, assuming equidistant binning
        double binningUncertSquared = pow( h1_fast.GetBinWidth(1), 2 ) / 12; // assume uniform distribution
        binningUncertSquared = 0; // binning uncertainty are taken into account by including/excluding the border-bins

        // now combine binning uncertainty, and stat h1 and stat h2
        double errorUp = sqrt( pow(h1_fast.GetBinCenter( binFull_StatFastUp )-efull, 2 ) + pow(h1_fast.GetBinCenter( binFull_StatFullUp )-efull, 2 ) + binningUncertSquared ) / efast;
        double errorDn = sqrt( pow(h1_fast.GetBinCenter( binFull_StatFastDn )-efull, 2 ) + pow(h1_fast.GetBinCenter( binFull_StatFullDn )-efull, 2 ) + binningUncertSquared ) / efast;

        out.SetPoint( binFast, efast, scale );
        out.SetPointError( binFast, 0, 0, errorDn, errorUp );
    }

    return out;

}

TH1F graphToHisto( const TGraph& gr ) {
    auto h = gr.GetHistogram();
    h->Reset( "ICESM" );
    for( int i=0;i<gr.GetN();i++) {
        double x,y;
        gr.GetPoint(i,x,y);
        h->SetBinContent( h->FindBin(x), y );
    }
    return *h;
}

void applyScale( TH1D h1_fast, TH1D h1_full, const TGraphAsymmErrors& scale, const std::string& savename ) {

    auto fast_clone = (TH1D*)h1_fast.Clone();

    auto closure = (TH1D*)h1_fast.Clone();
    closure->Reset( "ICESM" );

    auto h_scale = graphToHisto( scale );

    for (int i=0;i<fast_clone->GetNbinsX()+2;i++) {
        double e = fast_clone->GetBinCenter(i);
        double thisScale = h_scale.GetBinContent( h_scale.FindBin( e )-1 ); // why -1 ?
        double newE = e*thisScale;
        int newBin = fast_clone->FindBin( newE );
        closure->SetBinContent( newBin, fast_clone->GetBinContent(i) + closure->GetBinContent( newBin ) );
    }
    h1_fast = *closure;

    // clone inputs, since we modify them for drawing
    h1_fast.Scale( 1./h1_fast.Integral() );
    h1_full.Scale( 1./h1_full.Integral() );
    h1_full.SetLineColor(1);
    h1_fast.SetLineColor(2);
    auto minBin = std::min( h1_fast.FindFirstBinAbove(), h1_full.FindFirstBinAbove() );
    auto maxBin = std::max( h1_fast.FindLastBinAbove(),  h1_full.FindLastBinAbove()  );

    h1_fast.GetXaxis()->SetRange( minBin, maxBin );
    h1_full.GetXaxis()->SetRange( minBin, maxBin );

    h1_fast.SetMaximum( 1.05*std::max( h1_fast.GetMaximum(), h1_full.GetMaximum() ) );

    TCanvas* can = new TCanvas();
    can->cd();
    h1_fast.Draw("hist");
    h1_full.Draw("hist same");

    can->SaveAs( (std::string("plots/")+savename+"_closure.pdf").c_str() );
}

void drawAll( TH1D h1_fast, TH1D h1_full, TGraphAsymmErrors scale, TGraphAsymmErrors corrScale, const std::string& savename ) {
    // The inputs are cloned, so we can modify them


    // "Closure test"
    //applyScale( h1_fast, h1_full, corrScale, savename );

    h1_fast.Scale( 1./h1_fast.GetEntries() );
    h1_full.Scale( 1./h1_full.GetEntries() );

    h1_full.SetLineColor(1);
    h1_fast.SetLineColor(2);

    auto minBin = std::min( h1_fast.FindFirstBinAbove(), h1_full.FindFirstBinAbove() );
    auto maxBin = std::max( h1_fast.FindLastBinAbove(),  h1_full.FindLastBinAbove()  );
    minBin=0;
    maxBin=h1_fast.FindBin(1.05);

    h1_fast.GetXaxis()->SetRange( minBin, maxBin );
    h1_full.GetXaxis()->SetRange( minBin, maxBin );

    scale.GetXaxis()->SetLimits( h1_fast.GetBinLowEdge( minBin ),
        h1_fast.GetBinLowEdge( maxBin+1 ));

    scale.SetMaximum(1.1);
    scale.SetMinimum(0.9);
    //scale.GetYaxis()->SetNdivisions( 2, 0, 2 );
    scale.SetMarkerStyle( kFullDotMedium );
    corrScale.SetLineColor(6);

    h1_fast.SetMaximum( 1.05*std::max( h1_fast.GetMaximum(), h1_full.GetMaximum() ) );

    TCanvas* can = new TCanvas();
    can->cd();
    h1_fast.Draw("hist");
    h1_full.Draw("hist same");

    newRatioPad(0.6);
    scale.Draw("apz0");
    corrScale.Draw("same l");
    TLine* oneLine = new TLine();
    oneLine->SetLineStyle(2);
    oneLine->DrawLine( h1_fast.GetXaxis()->GetBinLowEdge( minBin ), 1, h1_fast.GetXaxis()->GetBinLowEdge( maxBin+1 ), 1 );

    can->SaveAs( (std::string("plots/")+savename+".pdf").c_str() );

}

std::string getBinLabel( const TH3F& h3, int xbin, int ybin ) {

    // Get name (bin ranges) for this specific histogram
    // e.g. E=20GeV, 0.2 < eta < 0.4
    std::ostringstream sname;
    sname << "E^{gen} = " << int(h3.GetXaxis()->GetBinCenter( xbin )) << " GeV, "
          << std::setprecision(3) << h3.GetYaxis()->GetBinLowEdge( ybin )
          << " #leq |#eta^{gen}| < "
          << std::setprecision(3) << h3.GetYaxis()->GetBinUpEdge( ybin );
    std::string name = sname.str();
    return name;

}

TH3F calculateResponse( const TH3F& h3_fast, const TH3F& h3_full ) {

    // This is the output histogram
    auto h3_scale = *((TH3F*)h3_fast.Clone("responseVsEVsEta"));
    h3_scale.Reset();

    // For each E_gen, eta_gen bin, extract the one dimensional E_sim/E_gen histogram
    for( int xbin=1; xbin< h3_fast.GetNbinsX()+1; ++xbin ) { // E_gen
        for( int ybin=1; ybin< h3_fast.GetNbinsY()+1; ++ybin ) { // eta_gen

            auto name = getBinLabel( h3_fast, xbin, ybin );

            // Create the 1d histograms
            auto h1_fast = h3_fast.ProjectionZ( (name+" fast").c_str(), xbin, xbin, ybin, ybin );
            auto h1_full = h3_full.ProjectionZ( (name+" full").c_str(), xbin, xbin, ybin, ybin );

            if( !h1_fast->GetEntries() || !h1_full->GetEntries() ) continue;

            name += ";E/E_{gen}; Normalized Entries  ";
            h1_fast->SetTitle( name.c_str() );
            h1_full->SetTitle( name.c_str() );

/*
            auto h1_scale = getSimplifiedScale( *h1_fast, *h1_full );
            for( int i=1;i<h1_scale.GetNbinsX()+1;i++){
                h3_scale.SetBinContent( xbin, ybin, i, h1_scale.GetBinContent(i) );
            }
*/

            auto scale = getScaleWithUncertainties( *h1_fast, *h1_full );
            auto corrScale = modifyScale( scale, h1_full->GetMean()/h1_fast->GetMean() );

            // Push back the scale into the output histogram
            for( auto i=0; i<corrScale.GetN(); i++) {
                double x,y;
                corrScale.GetPoint(i,x,y);
                // i+1 corresponds to the z-bin
                h3_scale.SetBinContent( xbin, ybin, i+1, y );
            }

            std::string savename = std::to_string(xbin) + "and" + std::to_string(ybin);
            drawAll( (TH1D)(*h1_fast), (TH1D)(*h1_full), scale, corrScale, savename );
        }
    }

    return h3_scale;
}

TH2D drawMeanResponse( const TH3F& h3_fast, const TH3F& h3_full ) {
    // Analysis function. The result of this function is NOT used for scaling,
    // but for analysis purpose only. The mean scale (E_sim/E_gen) is plotted
    // as a function of E_gen and eta_gen

    auto h2_fast = h3_fast.Project3DProfile( "xy" )->ProjectionXY("pxy");
    auto h2_full = h3_full.Project3DProfile( "xy2" )->ProjectionXY("pxy2");

    h2_full->Divide( h2_fast );

    h2_full->SetTitle("");
    h2_full->SetXTitle( h3_fast.GetXaxis()->GetTitle() );
    h2_full->SetYTitle( h3_fast.GetYaxis()->GetTitle() );
    h2_full->SetZTitle( "Mean scale" );

    h2_full->SetMaximum(1.05);
    h2_full->SetMinimum(0.95);

    h2_full->GetYaxis()->SetRangeUser(5, 105);

    setStyle2d();
    TCanvas* can = new TCanvas();
    can->cd();

    h2_full->Draw("colz");

    can->SaveAs( "meanResponse.pdf" );
    return *h2_full;

}

TH3F meanResponseAsH3( const TH3F& h3_fast, const TH3F& h3_full ) {
    /* I'm to lazy to rewrite the cmssw code, so just put the th2 into a th3
     * the z-axis of courese does not make any sence
     */
    auto h2 = drawMeanResponse( h3_fast, h3_full );

    // This is the output histogram
    auto h3_scale = *((TH3F*)h3_fast.Clone());
    h3_scale.Reset();

    // For each E_gen, eta_gen bin, extract the one dimensional E_sim/E_gen histogram
    for( int xbin=1; xbin< h3_scale.GetNbinsX()+1; ++xbin ) { // E_gen
        float genE = h3_fast.GetXaxis()->GetBinCenter( xbin );
        for( int ybin=1; ybin< h3_scale.GetNbinsY()+1; ++ybin ) { // eta_gen
            float genEta = h3_fast.GetYaxis()->GetBinCenter( ybin );
            float content = h2.GetBinContent( h2.FindFixBin( genEta, genE ) );
            for( int zbin=1; zbin< h3_scale.GetNbinsZ()+1; ++zbin ) { // energy
                h3_scale.SetBinContent( xbin, ybin, zbin, content );
            }
        }
    }

    return h3_scale;

}

TH3F fill3dHist_simple( TChain& chain ) {
    //auto h = TH3F("responseVsEVsEta", ";E_{gen};#eta_{gen};E/E_{gen}", 100, 5, 1005, 400, 0, 3.2, 100, 0.9, 1.05 );
    auto h = TH3F("responseVsEVsEta", ";E_{gen};#eta_{gen};E/E_{gen}", 100, 5, 1005, 400, 0, 3.2, 2000, 0, 1.05 );
    h.Rebin3D(1, 10, 1 );

    float e, eta, r;

    chain.SetBranchAddress( "e", &e );
    chain.SetBranchAddress( "r", &r );
    chain.SetBranchAddress( "eta", &eta );

    for( int i=0; i<chain.GetEntries(); i++ ) {
        chain.GetEntry(i);
        h.Fill( e, eta, r, 1 );
    }
    return h;
}

TH3F fill3dHist( TChain& chain ) {
    auto h = TH3F("responseVsEVsEta", ";E_{gen};#eta_{gen};E/E_{gen}", 100, 5, 1005, 400, 0, 3.2, 100, 0.9, 1.05 );
    //auto h = TH3F("responseVsEVsEta", ";E_{gen};#eta_{gen};E/E_{gen}", 100, 5, 1005, 400, 0, 3.2, 2000, 0, 1.05 );
    h.Rebin3D(1, 10, 1 );
    TVector3* genVec = 0;
    TVector3* simVec = 0;
    chain.SetBranchAddress( "genVec", &genVec );
    chain.SetBranchAddress( "hitVec", &simVec );
    for( int i=0; i<chain.GetEntries(); i++ ) {
        chain.GetEntry(i);
        float genE = genVec->Mag();
        float genEta = genVec->Eta();
        float oldRes = simVec->Mag() / genVec->Mag();
        h.Fill( genE, genEta, oldRes );
    }
    return h;
}



int main_old( int argc, char** argv ) {
    setStyle();

    if ( argc < 3 ) {
        std::cerr << "Usage: " << argv[0] << " fastsim.root fullsim.root" << std::endl;
    }
    std::string filenameFast = argv[1];
    std::string filenameFull = argv[2];


    auto fastTree = getChain( argv[1], "SimTreeProducer/SimTree" );
    auto fullTree = getChain( argv[2], "SimTreeProducer/SimTree" );
    const char* cutString = "hitVec.Mag()/genVec.Mag() > 0.9 &&  genVec.Eta()>1.5";
    //const char* cutString = "genVec.Eta()>0.0035 && genVec.Eta()<1.475";
    //const char* cutString = "genVec.Eta()>0.0035 && genVec.Eta()<1.475";
    fastTree = (TChain*)fastTree->CopyTree( cutString );
    fullTree = (TChain*)fullTree->CopyTree( cutString );


    auto h3d_scale = calculateResponse( fill3dHist( *fastTree ), fill3dHist( *fullTree ) );
    //auto h3d_scale = getHist<TH3F>( "scaleECALFastsim.root", "responseVsEVsEta" );

    TVector3* genVec = 0;
    TVector3* simVec = 0;

    fastTree->SetBranchAddress( "genVec", &genVec );
    fastTree->SetBranchAddress( "hitVec", &simVec );
    TH1F h1_fastRes("h1_fastRes", "", 100, 0.9, 1.01 );
    TH1F h1_fullRes("h1_fullRes", "", 100, 0.9, 1.01 );


    for( int i=0; i<fastTree->GetEntries(); i++ ) {
        fastTree->GetEntry(i);
        float genE = genVec->Mag();
        float genEta = genVec->Eta();
        float oldRes = simVec->Mag() / genVec->Mag();
        float scale = h3d_scale.GetBinContent( h3d_scale.FindFixBin( genE, genEta, oldRes ) );
        float newRes = oldRes * scale;
        h1_fastRes.Fill( newRes );
    }
    auto fullTree2 = getChain( argv[2], "SimTreeProducer/SimTree" );
    fullTree2 = (TChain*)fullTree2->CopyTree( cutString );
    fullTree2->Draw("hitVec.Mag()/genVec.Mag()>>h1_fullRes", "1", "goff" );

    TCanvas c1;
    h1_fullRes.Draw();
    h1_fastRes.SetLineColor(2);
    h1_fastRes.Scale( h1_fullRes.Integral() / h1_fastRes.Integral() );
    h1_fastRes.Draw("same");
    gPad->SaveAs("test.pdf");


    return 0;

}

int main( int argc, char** argv ) {
    setStyle();

    auto fastTree = getChain( "../../CMSSW/CMSSW_7_3_0/src/Analyzer/ECALScaleFactorCalculator/3d_fast.root_E30.root", "ecalScaleFactorCalculator/responseTree" );
    auto fullTree = getChain( "../../CMSSW/CMSSW_7_3_0/src/Analyzer/ECALScaleFactorCalculator/3d_full.root_E30.root", "ecalScaleFactorCalculator/responseTree" );
//    auto h3d_scale = calculateResponse( fill3dHist_simple( *fastTree ), fill3dHist_simple( *fullTree ) );
    auto h3d_scale = getHist<TH3F>( "scaleECALFastsim.root", "responseVsEVsEta" );




    float e, eta, r;

    fastTree->SetBranchAddress( "e", &e );
    fastTree->SetBranchAddress( "r", &r );
    fastTree->SetBranchAddress( "eta", &eta );
    TH1F h1_fastRes("h1_fastRes", "", 100, 0.9, 1.01 );
    TH1F h1_fullRes("h1_fullRes", "", 100, 0.9, 1.01 );


    for( int i=0; i<fastTree->GetEntries(); i++ ) {
        fastTree->GetEntry(i);
        float scale = h3d_scale.GetBinContent( h3d_scale.FindFixBin( e,eta,r ) );
        float newRes = r * scale;
        h1_fastRes.Fill( newRes );
    }
    fullTree->Draw("r>>h1_fullRes", "1", "goff" );

    TCanvas c1;
    h1_fullRes.Draw();
    h1_fastRes.SetLineColor(2);
    h1_fastRes.Scale( h1_fullRes.Integral() / h1_fastRes.Integral() );
    h1_fastRes.Draw("same");
    gPad->SaveAs("test.pdf");

    return 0;

}


