#include <TColor.h>
#include <TStyle.h>
#include <TPad.h>

TStyle* setStyle() {
    TStyle *myStyle = new TStyle("knutStyle","Style defined myself");

    myStyle->SetCanvasDefH(600);
    myStyle->SetCanvasDefW(600);

    myStyle->SetCanvasColor( kWhite );
    myStyle->SetPadColor( kWhite );
    myStyle->SetCanvasBorderMode(0);
    myStyle->SetPadBorderMode(0);

    myStyle->SetPadTickX(1);
    myStyle->SetPadTickY(1);

    myStyle->SetOptStat(0);

    myStyle->SetTitleBorderSize(1);
    myStyle->SetTitleFillColor( kWhite );

    myStyle->SetPadTopMargin(0.05);
    myStyle->SetPadBottomMargin(0.09);
    myStyle->SetPadLeftMargin(0.14);
    myStyle->SetPadRightMargin(0.02);

    myStyle->SetTitleYOffset(1.6);

    myStyle->SetPalette(56);
    myStyle->SetNumberContours(999);

    myStyle->SetLegendBorderSize(0);
    myStyle->SetLegendFillColor( kWhite );

    myStyle->cd();
    return myStyle;
}

TStyle* setStyle2d() {
    auto myStyle = setStyle();
    myStyle->SetOptLogx(0);
    myStyle->SetOptLogy(0);
    myStyle->SetPadRightMargin(0.18);

    int NRGBs = 3, NCont = 999;
    gStyle->SetNumberContours(NCont);
    Double_t stops[NRGBs] = { 0.00, 0.50, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 1.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 1.00, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 1.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

    myStyle->cd();
    return myStyle;
}


TPad* newRatioPad( double csf=0.2 ) {
    // csf: the ratio in which the pad is splitted

    // Delete label and title of all histograms in the current pad
    TIter next( gPad->GetListOfPrimitives() );
    TObject* obj;
    while( ( obj = next() ) ) {
        try {

            auto xaxis = ((TH1*)obj)->GetXaxis();
            xaxis->SetLabelSize(0);
            xaxis->SetLabelColor(0);
            xaxis->SetLabelOffset(1000);
            xaxis->SetTitle("");
            xaxis->SetTitleColor(0);
            xaxis->SetTitleSize(0);

        } catch ( ... ) { }
    }

    gPad->SetBottomMargin( csf + (1-csf)*gPad->GetBottomMargin() - csf*gPad->GetTopMargin() );
    TPad* rPad = new TPad( "rPad", "ratio", 0, 0, 1, 1 );
    rPad->SetTopMargin( (1-csf) - (1-csf)*rPad->GetBottomMargin() + csf*rPad->GetTopMargin() );
    rPad->SetFillStyle(3955);
    rPad->Draw();
    rPad->cd();
    rPad->SetLogy(0);

    return rPad;
}

