import ROOT
import ratio
import style
import math
ROOT.gROOT.SetBatch()


def readHisto( filename, histoname ):
    f = ROOT.TFile( filename )
    h = f.Get( histoname )
    if not h: print "No such histo found"
    h = ROOT.gROOT.CloneObject( h )
    return h

def getObjectNamesRec( filename, path ):
    file1 = ROOT.TFile( filename )
    outList = []
    tempDir = file1.GetDirectory( path )
    for element in tempDir.GetListOfKeys():
        if isinstance( element.ReadObj(), ROOT.TH1 ):
            outList.append( path+"/"+element.GetName() )
        elif isinstance( element.ReadObj(), ROOT.TDirectory ):
            outList.extend( getObjectNamesRec( file1, path+"/"+element.GetName() ) )
        elif isinstance( element.ReadObj(), ROOT.TObjString ):
            pass
        else:
            print "Do not know what to do with", element.GetName(), type(element.ReadObj())

    file1.Close()
    return outList

def getFirstDifferentPosition( names ):
    # Finds the position in the string where the strings in 'names' differ
    for i, letters in enumerate(zip(*names)):
        if list(letters) != [letters[0]]*len(letters):
            return i

def getAbsFromHist( origHist ):
    origNbins = origHist.GetNbinsX()
    origXmin = origHist.GetBinLowEdge(1)
    origXmax = origHist.GetBinLowEdge(origHist.GetNbinsX()+1)
    if origNbins%2:
        print "cant handle histos with odd number of bins"
        return origHist
    if origXmin + origXmax:
        print "cant handle assymetric histograms"
        return origHist

    h = ROOT.TH1F( "", origHist.GetTitle(), origNbins/2, 0, origXmax )

    for origBin in range( origNbins+2 ):
        newBin = int(abs(origBin - (origNbins+1.)/2)) + 1
        h.SetBinContent( newBin, h.GetBinContent(newBin) + origHist.GetBinContent(origBin) )
        h.SetBinError( newBin, math.sqrt( h.GetBinError(newBin)**2 + origHist.GetBinError(origBin)**2 ) )

    h.Scale(0.5)
    return h


def drawTwo( hfull, hfast, name ):

    c = ROOT.TCanvas()
    if not hfull.GetEntries() and not hfast.GetEntries(): return

    for h in hfull, hfast:
        h.Sumw2()
        h.Scale( 1./h.GetEntries() )

    if isinstance( hfull, ROOT.TH2 ):
        hfull = hfull.ProfileX("xfull")
        hfast = hfast.ProfileX("xfast")

    #hfull = getAbsFromHist( hfull )
    #hfast = getAbsFromHist( hfast )
    hfull.SetLineColor(1)
    hfast.SetLineColor(2)
    hfast.SetFillColor(2)
    hfast.SetFillStyle(3004)

    histoMaximum = max([ h.GetMaximum() for h in hfull,hfast ])
    histoMinimum = min([ h.GetMinimum() for h in hfull,hfast ])
    for h in hfull, hfast:
        if histoMinimum is 0:
            h.SetMinimum( 0 )
        h.SetMaximum( 1.01*histoMaximum )
        if "h_ele_PoPtrueVsPt" in name:
            h.SetMinimum(0.9)
        h.SetTitleOffset( 1.6, "y" )
        h.SetTitleOffset( 1, "x" )
        h.SetMarkerSize(0)
        h.SetTitleSize( ROOT.gStyle.GetTextSize(), "xyz" )
        h.SetTitleFont( ROOT.gStyle.GetTextFont(), "xyz" )
        h.SetLabelSize( ROOT.gStyle.GetTextSize(), "xyz" )
        h.SetLabelFont( ROOT.gStyle.GetTextFont(), "xyz" )

    hfull.SetName("FullSim")
    hfast.SetName("FastSim")

    #first draw to get stat boxes
    ROOT.gStyle.SetOptStat(0)

    if "h_scl_et" in name:
        for h in hfull, hfast:
            h.SetMinimum(0)
            h.SetMaximum(0.068)
            h.SetTitle(";Supercluster E_{T} (GeV);Normalised Entries   ")
    if "h_ele_mee_os" in name:
        for h in hfull, hfast:
            h.SetMinimum(0)
            h.SetMaximum(0.158)
            h.SetTitle(";m_{e^{+}e^{-}} (GeV);Normalised Entries   ")
    if "h_ele_chargedHadronIso" in name:
        for h in hfull, hfast:
            h.SetMinimum(0)
            h.SetMaximum(0.072)
            h.GetXaxis().SetRangeUser(0,9)
            h.SetTitle(";Charged Hadron Isolation (GeV);Normalised Entries   ")
    if "r9All" in name:
        for h in hfull, hfast:
            h.GetXaxis().SetRangeUser( 0.6, 1 )
            h.SetTitle(";r9;Normalised Entries   ")
    if "sigmaIetaIetaAll" in name:
        for h in hfull, hfast:
            h.GetXaxis().SetRangeUser( 0, 0.04 )
            h.SetTitle(";#sigma_{i#etai#eta};Normalised Entries   ")
    if "eResBarrel" in name:
        for h in hfull, hfast:
            h.SetTitle(";Barrel E/E_{gen};Normalised Entries   ")
    if "eResEndcap" in name:
        for h in hfull, hfast:
            h.SetTitle(";Endcap E/E_{gen};Normalised Entries   ")

    hfull.Draw("e")
    hfast.Draw("e same")

    label = ROOT.TLatex()
    label.SetTextFont( ROOT.gStyle.GetTextFont() )
    label.SetTextSize( ROOT.gStyle.GetTextSize() )
    sample = "H#rightarrow #gamma#gamma"
    saveName = "Hgg"
    if "/h_" in name:
        sample = "Z#rightarrow ee"
        saveName ="Zee"
    label.DrawLatexNDC( .05, .96, "CMS Private Work    %s   13TeV           #color[1]{FullSim} #color[2]{FastSim}"%sample )


    #r = ratio.Ratio( "h1/h2", h1, h2 )
    #r.draw(0.5,1.5)

    ROOT.gPad.GetCanvas().SaveAs( "plots%sTwo/%s.pdf"%(saveName,name.split("/")[-1]) )

def compareTwoHistograms( fullname, fastname, path ):

    names = getObjectNamesRec( fastname, path )
    names = [
        path+"/h_scl_et",
        path+"/h_ele_mee_os",
        path+"/h_ele_chargedHadronIso"
    ]
    names = [
        path+"/r9All",
        path+"/sigmaIetaIetaAll",
        path+"/eResBarrel",
        path+"/eResEndcap"
    ]

    for name in names:
        hfull = readHisto( fullname, name )
        hfast = readHisto( fastname, name )
        drawTwo( hfull, hfast, name )




if __name__ == "__main__":
    fullname = "../../CMSSW/CMSSW_7_3_0/src/harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValZEE_13__official_FullSim.root"
    fastname = "../../CMSSW/CMSSW_7_3_0/src/harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValZEE_13__official_FastSim.root"
    path = "DQMData/Run 1/EgammaV/Run summary/ElectronMcSignalValidator"

    fullname = "../../CMSSW/CMSSW_7_3_0/src/harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValH130GGgluonfusion_13__official_FullSim.root"
    fastname = "../../CMSSW/CMSSW_7_3_0/src/harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValH130GGgluonfusion_13__official_FastSim.root"
    path = "DQMData/Run 1/EgammaV/Run summary/pfPhotonValidator/Photons"



    compareTwoHistograms( fullname, fastname, path )
