#!/usr/bin/env python2
# -*- coding: utf-8 -*-


import ConfigParser
import ROOT
import math
from random import randint
from sys import maxint

# private libs
import ratio
import style
import multiplot

ROOT.gROOT.ForceStyle()

def randomName():
    """
    Generate a random string. This function is useful to give ROOT objects
    different names to avoid overwriting.
    """
    return "%x"%(randint(0, maxint))

def readHist( filename, histoname ):
    f = ROOT.TFile( filename )
    h = f.Get( histoname )
    h = ROOT.gROOT.CloneObject( h )
    h.drawOption_ = ""
    return h

def getObjectNames( filename, path ):
    f = ROOT.TFile( filename )
    tmpDir = f.GetDirectory( path )

    outList = []
    for element in tmpDir.GetListOfKeys():
        if isinstance( element.ReadObj(), ROOT.TH1 ):
            outList.append( element.GetName() )
        elif isinstance( element.ReadObj(), ROOT.TObjString ):
            pass
        else:
            print "do not know what to do with", element.GetName()

    return outList

def getNgenDQM( filename ):
    h = readHist( filename, "DQMData/Run 1/Generator/Run summary/GenParticles/nEvt" )
    nGen = h.GetBinContent( 1 )
    return nGen

def getValAndError( val, err, sig=2 ):
    from math import floor, log10
    digit = sig - int(floor(log10(err))) - 1
    return "{} #pm {}".format( round(val,digit), round(err,digit) )

def getOwnStatBox( h, x1, y1 ):

    text = ROOT.TLatex()
    text.SetTextColor( h.GetLineColor() )

    mean = h.GetMean()
    mean_err = h.GetMeanError()

    text.DrawLatexNDC( x1, y1, "#mu = "+getValAndError( mean, mean_err ) )

def absHistWeighted( origHist ):
    origNbins = origHist.GetNbinsX()
    origXmin = origHist.GetBinLowEdge(1)
    origXmax = origHist.GetBinLowEdge(origHist.GetNbinsX()+1)
    if origXmin + origXmax:
        print "cant handle assymetric histograms"
        return origHist

    h = ROOT.TH1F( "", origHist.GetTitle(), int(math.ceil(origNbins/2.)), 0, origXmax )

    for origBin in range( origNbins+2 ):
        newBin = int(abs(origBin - (origNbins+1.)/2)) + 1

        c1 = origHist.GetBinContent( origBin )
        e1 = origHist.GetBinError( origBin )
        c2 = h.GetBinContent( newBin )
        e2 = h.GetBinError( newBin )

        if e1 and e2:
            h.SetBinContent( newBin, ( c1*e1**-2 + c2*e2**-2 )/(e1**-2 + e2**-2) )
            h.SetBinError( newBin, 1./math.sqrt( e1**-2 + e2**-2 ) )

        else:
            h.SetBinContent( newBin, origHist.GetBinContent(origBin) )
            h.SetBinError( newBin, origHist.GetBinError(origBin) )

    return h


def draw( file1, file2, path, name, config ):

    c = ROOT.TCanvas()

    if "H130GG" in file1:
        processTex = "H#rightarrow#gamma#gamma   13TeV"
        processName = "Hgg"
    if "ZEE" in file1:
        processTex = "Z#rightarrowee  13TeV"
        processName = "Zee"
    processTex += "          FullSim #color[2]{FastSim}"

    h1 = readHist( file1, "{}/{}".format( path, name ) )
    h2 = readHist( file2, "{}/{}".format( path, name ) )
    if not round( h1.GetEntries() ) or not round( h2.GetEntries() ): return
    if not round( h1.Integral() ) or not round( h2.Integral() ): return

    if isinstance( h1, ROOT.TH2 ):
        h1.Sumw2()
        h2.Sumw2()
        h1 = h1.ProfileX( randomName() )
        h2 = h2.ProfileX( randomName() )

    if "VsEta" in name:
        h1 = absHistWeighted( h1 )
        h2 = absHistWeighted( h2 )

    h2.SetLineColor(2)

    h1.SetName("FullSim")
    h2.SetName("FastSim")

    for h in h1, h2:
        # constumize histograms
        xmin = h.GetXaxis().GetXmin()
        xmax = h.GetXaxis().GetXmax()
        if config.has_option( name, "xmin" ): xmin = config.getfloat( name, "xmin")
        if config.has_option( name, "xmax" ): xmax = config.getfloat( name, "xmax" )
        h.GetXaxis().SetRangeUser( xmin, xmax )

        if config.has_option( name, "title" ):
            h.SetTitle( config.get( name, "title" )+ "  " )
        else:
            h.SetXTitle( h.GetTitle() +"    "+ h.GetXaxis().GetTitle() )
            h.SetTitle( "" )

        h.drawOption_ = "hist e"
        h.Sumw2()
        h.SetMarkerSize(0)

    for h in h1, h2:
        if not isinstance( h, ROOT.TProfile ) and not "VsEta" in name:
            h.Scale( 1./h.GetEntries() )

    m = multiplot.Multiplot()
    m.add( h1 )
    m.add( h2 )
    m.Draw()

    #getOwnStatBox( h1, .6,.9 )
    #getOwnStatBox( h2, .6,.85 )


    label = ROOT.TLatex()
    label.DrawLatexNDC( .01, .96, "#font[61]{CMS} Private Work   "+processTex )

    r = ratio.Ratio( "Full/Fast  ", h1, h2 )
    r.draw(0.5,1.5)
    ROOT.gPad.GetCanvas().SaveAs("plots/%s_%s.pdf"%(processName, name ))

def compareHistograms( file1, file2, path ):
    names = getObjectNames( file1, path )
    names = [ "pEResVsEtaAll" ]

    configuration = ConfigParser.SafeConfigParser()
    configuration.read( "histoDefinitions.cfg" )


    for name in names:
        draw( file1, file2, path, name, configuration )


if __name__ == "__main__":

    compareHistograms(
        "../../DQM_V0001_R000000001__CMSSW_7_3_0__RelValH130GGgluonfusion_13__official_FullSim.root",
        "../../DQM_V0001_R000000001__CMSSW_7_3_0__RelValH130GGgluonfusion_13__official_FastSim.root",
        "DQMData/Run 1/EgammaV/Run summary/pfPhotonValidator/Photons",
    )
    #compareHistograms(
    #    "../../DQM_V0001_R000000001__CMSSW_7_3_0__RelValH130GGgluonfusion_13__official_FullSim.root",
    #    "../../DQM_V0001_R000000001__CMSSW_7_3_0__RelValH130GGgluonfusion_13__official_FastSim.root",
    #    "DQMData/Run 1/EgammaV/Run summary/pfPhotonValidator/Efficiencies",
    #)
    #compareHistograms(
    #    "../../DQM_V0001_R000000001__CMSSW_7_3_0__RelValZEE_13__official_FullSim.root",
    #    "../../DQM_V0001_R000000001__CMSSW_7_3_0__RelValZEE_13__official_FastSim.root",
    #    "DQMData/Run 1/EgammaV/Run summary/ElectronMcSignalValidator",
    #)
