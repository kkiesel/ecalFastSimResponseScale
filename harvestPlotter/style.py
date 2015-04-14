import ROOT


def defaultStyle():
    st = ROOT.TStyle("defaultStyle", "Knut's owns style" )

    st.SetCanvasColor( ROOT.kWhite )
    st.SetCanvasDefH(600)
    st.SetCanvasDefW(600)

    st.SetPadTickX( 1 )
    st.SetPadTickY( 1 )

    st.SetPadColor( ROOT.kWhite )

    # Margins:
    st.SetPadTopMargin(0.06)
    st.SetPadBottomMargin(0.12)
    st.SetPadLeftMargin(0.16)
    st.SetPadRightMargin(0.04)

    st.SetTitleFillColor( ROOT.kWhite )
    st.SetTitleBorderSize( 0 )

    st.SetTitleOffset( 1.1, "x" )
    st.SetTitleOffset( 1.6, "y" )

    st.SetStatBorderSize(1)
    st.SetStatColor(0)

    st.SetLegendBorderSize(0)
    st.SetLegendFillColor( ROOT.kWhite )

    st.SetOptStat(0)

    textSize = 0.05
    st.SetLabelSize( textSize, "xyz" )
    st.SetTitleSize( textSize, "xyz" )

    st.SetTextFont( st.GetLabelFont() )
    st.SetTextSize( st.GetLabelSize() )

    st.SetNdivisions( 505, "xyz" )

    st.cd()
    return st

defaultStyle()

# not style, but similar
ROOT.gROOT.SetBatch()
ROOT.TH1.SetDefaultSumw2()
ROOT.gROOT.ForceStyle()
