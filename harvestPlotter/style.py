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
    st.SetPadTopMargin(0.05)
    #st.SetPadBottomMargin(0.13)
    st.SetPadLeftMargin(0.13)
    st.SetPadRightMargin(0.03)

    st.SetTitleFillColor( ROOT.kWhite )
    st.SetTitleBorderSize( 0 )

    st.SetStatBorderSize(1)
    st.SetStatColor(0)

    st.SetLegendBorderSize(0)
    st.SetLegendFillColor( ROOT.kWhite )


    font = 42
    textSize = 0.035

    st.SetTextFont( font )
    st.SetTextSize( textSize )




    st.cd()
    return st

defaultStyle()

# not style, but similar
ROOT.gROOT.SetBatch()
ROOT.TH1.SetDefaultSumw2()

