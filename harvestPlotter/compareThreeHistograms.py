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


def drawThree( hfull, hfast, hmod, name ):
    c = ROOT.TCanvas()
    if not hfull.GetEntries() or not hfast.GetEntries() or not hmod.GetEntries(): return

    for h in hfull, hfast, hmod:
        h.Sumw2()
        h.Scale( 1./h.GetEntries() )

    if isinstance( hfull, ROOT.TH2 ):
        hfull = hfull.ProfileX("xfull")
        hfast = hfast.ProfileX("xfast")
        hmod = hmod.ProfileX("xmod")

    #hfull = getAbsFromHist( hfull )
    #hfast = getAbsFromHist( hfast )
    #hfmod = getAbsFromHist( hmod )


    histoMaximum = max([ h.GetMaximum() for h in hfull,hfast,hmod ])
    for h in hfull, hfast, hmod:
        h.SetMaximum( 1.02*histoMaximum )
        if "h_ele_PoPtrueVsPt" in name:
            h.SetMinimum(0.9)
        if "h_scl_EoEtrue" in name:
            h.GetXaxis().SetRangeUser( 0.7, 1.2 )



    hfull.SetName("FullSim")
    hfast.SetName("FastSim")
    hmod.SetName("FastSim+ResponseScaling")

    #first draw to get stat boxes
    hfull.Draw()
    ROOT.gPad.Update()
    statsfull = hfull.GetListOfFunctions().FindObject("stats").Clone("stats1")
    hfast.Draw()
    ROOT.gPad.Update()
    statsfast = hfast.GetListOfFunctions().FindObject("stats").Clone("stats2")
    hmod.Draw()
    ROOT.gPad.Update()
    statsmod = hmod.GetListOfFunctions().FindObject("stats").Clone("stats3")

    hfull.SetLineColor(1)

    hfast.SetLineColor(2)
    statsfast.SetTextColor(2)
    statsfast.SetY1NDC(.62)
    statsfast.SetY2NDC(.82)

    hmod.SetLineColor(4)
    statsmod.SetTextColor(4)
    statsmod.SetY1NDC(.4)
    statsmod.SetY2NDC(.6)

    hfull.Draw("e")
    hfast.Draw("e same")
    hmod.Draw("e same")
    hfull.Draw("e same h")
    statsfast.Draw()
    statsmod.Draw()


    label = ROOT.TLatex()
    label.SetTextFont(42)
    label.DrawLatexNDC( .1, .9, "Z#rightarrow ee    13TeV" )

    #r = ratio.Ratio( "h1/h2", h1, h2 )
    #r.draw(0.5,1.5)

    ROOT.gPad.GetCanvas().SaveAs( "plots/%s.pdf"%name.split("/")[-1] )

def compareThreeHistograms( fullname, fastname, modname, path ):

    names = getObjectNamesRec( fastname, path )
    names = filter(lambda x: "scl_EoEtrue_" in x,  names )

    for name in names:
        hfull = readHisto( fullname, name )
        hfast = readHisto( fastname, name )
        hmod = readHisto( modname, name )
        drawThree( hfull, hfast, hmod, name )




if __name__ == "__main__":
    fastname = "../../CMSSW/CMSSW_7_3_0/src/testRunTheMatrix/fast/DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO_fast.root"
    fullname = "../../CMSSW/CMSSW_7_3_0/src/testRunTheMatrix/fast/DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO_full.root"
    modname  = "../../CMSSW/CMSSW_7_3_0/src/testRunTheMatrix/fast/DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO_mod1.root"
    path = "DQMData/Run 1/EgammaV/Run summary/ElectronMcSignalValidator"

    #fastname = "../../CMSSW/CMSSW_7_3_0/src/Analyzer/SimTreeWriter/fastsim_muchStat.root"
    #fullname = "../../CMSSW/CMSSW_7_3_0/src/Analyzer/SimTreeWriter/fullsim_muchStat.root"
    #modname = "../../CMSSW/CMSSW_7_3_0/src/Analyzer/SimTreeWriter/fastsim_val.root"
    #path = "SimTreeProducer"

    #fastname = "../../CMSSW/CMSSW_7_3_0/src/Analyzer/SimTreeWriter/closure_fast.root"
    #fullname = "../../CMSSW/CMSSW_7_3_0/src/Analyzer/SimTreeWriter/closure_full.root"
    #modname = "../../CMSSW/CMSSW_7_3_0/src/Analyzer/SimTreeWriter/closure_fast_validation.root"
    #path = "SimTreeProducer"


    compareThreeHistograms( fullname, fastname, modname, path )
