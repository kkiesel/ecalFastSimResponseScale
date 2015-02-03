#include<iostream>
#include<string>
#include<sstream>
#include<iomanip> // provides setprecision

// ROOT
#include<TCanvas.h>
#include<TEfficiency.h>
#include<TFile.h>
#include<TGraphAsymmErrors.h>
#include<TH3F.h>
#include<TLine.h>
#include<TMath.h>
#include<TProfile2D.h>

// user incuded files
#include "Style.h"

using namespace std;

template <class HIST>
HIST getHisto( std::string const & filename, std::string const & histname ) {
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
    //HIST hist( *((HIST*)obj));
    HIST hist = *((HIST*)obj);
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
            (errorUp+errorDn)/2 > std::abs( mean - 1 ) // uncertainty larger than correction
        ) {
            modScale->SetPoint( i, x, mean );
        }
        /*
        if( y + modScale->GetErrorYhigh(i) > 1 and y - modScale->GetErrorYlow(i) < 1 ) { // not significant
            // don't do this if this is a unity-crossing, eg. if the scale is closer to one than the mean
            if( std::abs( y - 1) > std::abs(mean - y) ) {
                cout << x << " : " << std::abs( y -1 ) << " \t" << (modScale->GetErrorYhigh(i)+modScale->GetErrorYlow(i))/2 << endl;
                if( std::abs(  mean -1 ) < (modScale->GetErrorYhigh(i)+modScale->GetErrorYlow(i))/2 ) { // low error
                    modScale->SetPoint( i, x, mean );
                }
            }
        }
        */

        // The uncertainties will not be used and are therefore set to 0
        modScale->SetPointEYhigh( i, 0 );
        modScale->SetPointEYlow ( i, 0 );
    }
    return *modScale;
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
    cout << "Entries Fast = " << entriesFast << endl;
    cout << "Entries Full = " << entriesFull << endl;

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

        bool debug = binFast == 375 || binFast == 1290; // pick an example bin
        if(debug) cout << "Calculate scale for fastsim bin " << binFast << " corresponding to E/E_gen=" << h1_fast.GetBinCenter(binFast) << endl;

        summedEntriesFast += h1_fast.GetBinContent( binFast );
        if(debug) cout << "Fastsim entries up to this bin (integral) = " << summedEntriesFast << endl;
        //double areaFast = summedEntriesFast/entriesFast; // not needed, since the mean is calculated as (up+down)/2

        // Take into account the statistical precission of h1
        double areaFastUp = eff.ClopperPearson( entriesFast, summedEntriesFast, alpha, true );
        double areaFastDn = eff.ClopperPearson( entriesFast, summedEntriesFast, alpha, false );
        if(debug) cout << "This corresponds to an normalised area ( acceptance) of something in between " << areaFastDn << " to " << areaFastUp << endl;

        // Number of entries in the fullsim corresponding to these bondaries
        int entriesFull_StatFastDn = floor( areaFastDn * entriesFull );
        int entriesFull_StatFastUp = ceil ( areaFastUp * entriesFull );
        if(debug) cout << "To aquire the same acceptance (area) in fullsim, " << entriesFull_StatFastDn << " to " << entriesFull_StatFastUp << " fullsim entries are needed" <<  endl;


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
        if(debug) cout << "The lowest possible bin with this requirement is " << binFull_StatFastDn << " (" << h1_full.GetBinCenter( binFull_StatFastDn ) << ") and the highest possible bin is " << binFull_StatFastUp << " (" << h1_full.GetBinCenter( binFull_StatFastUp ) << ")" << endl;

        // Take the middle. This is somehow arbitrary, but feel free to envolve a better method
        int binFullMiddle = (binFull_StatFastUp + binFull_StatFastDn)/2;
        if(debug) cout << "As a simplification, take the mean bin " << binFullMiddle << " (" << h1_full.GetBinCenter( binFullMiddle ) << ")" << endl;

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
        if(debug) cout << "The statistical uncertainty of fullsim corresponds to an region of " << binFull_StatFullDn << " (" << h1_full.GetBinCenter( binFull_StatFullDn ) << ")" << endl;

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

        if(debug) cout << "This yields a scale of " << scale << " - " << errorDn << " + " << errorUp << endl;

        out.SetPoint( binFast, efast, scale );
        out.SetPointError( binFast, 0, 0, errorDn, errorUp );
    }

    return out;

}

TH1F graphToHisto( const TGraphAsymmErrors& gr ) {
    auto h = gr.GetHistogram();
    h->Reset( "ICESM" );
    for( int i=0;i<gr.GetN();i++) {
        double x,y;
        gr.GetPoint(i,x,y);
        h->SetBinContent( h->FindBin(x), y );
    }
    return *h;
}

void applyScale( TH1D h1_fast, TH1D h1_full, TGraphAsymmErrors scale, const std::string& savename ) {

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

    applyScale( h1_fast, h1_full, corrScale, savename );

    // clone inputs, since we modify them for drawing
    h1_fast.Scale( 1./h1_fast.GetEntries() );
    h1_full.Scale( 1./h1_full.GetEntries() );

    h1_full.SetLineColor(1);
    h1_fast.SetLineColor(2);

    auto minBin = std::min( h1_fast.FindFirstBinAbove(), h1_full.FindFirstBinAbove() );
    auto maxBin = std::max( h1_fast.FindLastBinAbove(),  h1_full.FindLastBinAbove()  );
    h1_fast.GetXaxis()->SetRange( minBin, maxBin );
    h1_full.GetXaxis()->SetRange( minBin, maxBin );

    scale.GetXaxis()->SetLimits( h1_fast.GetBinLowEdge( minBin ),
        h1_fast.GetBinLowEdge( maxBin+1 ));

    scale.SetMaximum(1.05);
    scale.SetMinimum(0.95);
    scale.GetYaxis()->SetNdivisions( 2, 0, 2 );
    corrScale.SetLineColor(6);

    h1_fast.SetMaximum( 1.05*std::max( h1_fast.GetMaximum(), h1_full.GetMaximum() ) );

    TCanvas* can = new TCanvas();
    can->cd();
    h1_fast.Draw("hist");
    h1_full.Draw("hist same");

    newRatioPad();
    scale.Draw("a*0");
    corrScale.Draw("same l");
    TLine* oneLine = new TLine();
    oneLine->SetLineStyle(2);
    oneLine->DrawLine( h1_fast.GetXaxis()->GetBinLowEdge( minBin ), 1, h1_fast.GetXaxis()->GetBinLowEdge( maxBin+1 ), 1 );

    can->SaveAs( (std::string("plots/")+savename+".pdf").c_str() );

}


TH3F calculateResponse( const TH3F& h3_fast, const TH3F& h3_full ) {

    // This is the output histogram
    auto h3_scale = *((TH3F*)h3_fast.Clone());
    h3_scale.Reset();


    // For each E_gen, eta_gen bin, extract the one dimensional E_sim/E_gen histogram
    for( int xbin=1; xbin<h3_fast.GetNbinsX()+1; ++xbin ) { // E_gen
        for( int ybin=1; ybin<h3_fast.GetNbinsY()+1; ++ybin ) { // eta_gen
            xbin = 2;
            ybin=2;

            // Get name (bin ranges) for this specific histogram
            // e.g. E=20GeV, 0.2 < eta < 0.4
            std::ostringstream sname;
            sname << "E^{gen} = " << int(h3_fast.GetXaxis()->GetBinCenter( xbin )) << " GeV"
                  << ", "
                  << std::setprecision(3) << h3_fast.GetYaxis()->GetBinLowEdge( ybin )
                  << " #leq |#eta^{gen}| < "
                  << std::setprecision(3) << h3_fast.GetYaxis()->GetBinUpEdge( ybin );
            std::string name = sname.str();

            // Create the 1d histograms
            auto h1_fast = h3_fast.ProjectionZ( (name+" fast").c_str(), xbin, xbin, ybin, ybin );
            auto h1_full = h3_full.ProjectionZ( (name+" full").c_str(), xbin, xbin, ybin, ybin );

            // Here the binning of E_sim/E_gen can be adjusted
//            h1_fast->Rebin(500);
//            h1_full->Rebin(500);

            if( !h1_fast->GetEntries() || !h1_full->GetEntries() ) continue;

            name += ";E/E_{gen}; Normalized Entries  ";
            h1_fast->SetTitle( name.c_str() );
            h1_full->SetTitle( name.c_str() );
            auto scale = getScaleWithUncertainties( *h1_fast, *h1_full );
            auto corrScale = modifyScale( scale, h1_full->GetMean()/h1_fast->GetMean() );

            // Push back the scale into the output histogram
            for( auto i=0; i<corrScale.GetN(); i++) {
                double x,y;
                corrScale.GetPoint(i,x,y);

                int bin = h3_scale.FindBin( // in three dimensions
                    h3_fast.GetXaxis()->GetBinCenter( xbin ),
                    h3_fast.GetYaxis()->GetBinCenter( ybin ),
                    x );

                h3_scale.SetBinContent( bin, y );
            }

            std::string savename = std::to_string(xbin) + "and" + std::to_string(ybin);
           drawAll( (TH1D)(*h1_fast), (TH1D)(*h1_full), scale, corrScale, savename );
    exit(0);
        }
    }

    return h3_scale;
}

void drawMeanResponse( const TH3F& h3_fast, const TH3F& h3_full ) {
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

    auto st = setStyle2d();
    st->SetPalette(1);
    TCanvas* can = new TCanvas();
    can->cd();
    h2_full->Draw("colz");

    can->SaveAs( "meanResponse.pdf" );
    exit(0);

}

int main( int argc, char** argv ) {

    if ( argc < 3 ) {
        std::cerr << "Usage: " << argv[0] << " fastsim.root fullsim.root" << std::endl;
    }
    std::string filenameFast = argv[1];
    std::string filenameFull = argv[2];

    std::string histname = "ecalScaleFactorCalculator/energyVsEVsEta";

    setStyle();
    //gErrorIgnoreLevel = kWarning;

    //drawMeanResponse( getHisto<TH3F>( filenameFast, histname ), getHisto<TH3F>( filenameFull, histname ) ); // analysis function
    auto h = calculateResponse( getHisto<TH3F>( filenameFast, histname ), getHisto<TH3F>( filenameFull, histname ) );


    TFile file( "scaleECALFastsim.root", "recreate" );
    file.cd();
    h.Write();
    file.Close();
}


