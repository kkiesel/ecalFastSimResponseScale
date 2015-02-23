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
    if not h.GetSumw2N():
        h.Sumw2()
    h.drawOption_ = ""
    return h

def getObjectNames( filename, path ):
    f = ROOT.TFile( filename )
    tmpDir = f.GetDirectory( path )

    outList = []
    for element in tmpDir.GetListOfKeys():
        obj = element.ReadObj()
        if isinstance( obj, ROOT.TH1 ):
            outList.append( element.GetName() )
        elif isinstance( obj, ROOT.TObjString ) or isinstance( obj, ROOT.TTree ):
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
        if origXmin:
            print "cant handle assymetric histograms"
        # else: print "already symetric?"
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


def draw( files, path, name, config ):

    c = ROOT.TCanvas()

    processTex = ""
    processName = files[0].split("/")[6]
    if "H130GG" in files[0]:
        processTex = "H#rightarrow#gamma#gamma   13TeV"
        processName = "Hgg"
    if "ZEE" in files[0]:
        processTex = "Z#rightarrowee  13TeV"
        processName = "Zee"
    processTex += "          FullSim #color[2]{FastSim}"
    if len(files) > 2:
        processTex += " #color[4]{FastSim+Mod}"

    hists = []
    for file in files:
        h = readHist( file, "{}/{}".format( path, name ) )
        if not round( h.GetEntries()) or not round( h.Integral() ): return

        if isinstance( h, ROOT.TH2 ):
            h = h.ProfileX( randomName() )
        if "VsEta" in name:
            h = absHistWeighted( h )

        hists.append( h )

    hists[0].SetName("FullSim")
    if len(hists) > 1:
        hists[1].SetName("FastSim")
        hists[1].SetLineColor(2)
    if len(hists) > 2:
        hists[2].SetName("FastSim+Mod")
        hists[2].SetLineColor( ROOT.kBlue )

    for h in hists:
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
        h.SetMarkerSize(0)

    for h in hists:
        if not isinstance( h, ROOT.TProfile ) and not "VsEta" in name:
            h.Scale( 1./h.GetEntries() )

    m = multiplot.Multiplot()
    for h in hists:
        m.add( h )
    m.Draw()

    #for i, h in enumerate(hists):
    #    getOwnStatBox( h, .6,.9-0.05*i )


    label = ROOT.TLatex()
    label.DrawLatexNDC( .01, .96, "#font[61]{CMS} Private Work   "+processTex )

    r = ratio.Ratio( "Full/Fast  ", hists[0], hists[1] )
    if len(hists)>2:
        r = ratio.Ratio( "Full/Mod  ", hists[0], hists[2] )
    r.draw( )
    ROOT.gPad.GetCanvas().SaveAs("plots/%s_%s.pdf"%(processName, name ))

def compareHistograms( files, path="SimTreeProducer" ):
    names = getObjectNames( files[0], path )

    configuration = ConfigParser.SafeConfigParser()
    configuration.read( "histoDefinitions.cfg" )


    for name in names:
        draw( files, path, name, configuration )


if __name__ == "__main__":

    ##compareHistograms( [ "../../CMSSW/CMSSW_7_3_0/src/Analyzer/SimTreeWriter/fullsim_muchStat.root", "../../CMSSW/CMSSW_7_3_0/src/Analyzer/SimTreeWriter/fastsim_muchStat.root", "../../CMSSW/CMSSW_7_3_0/src/Analyzer/SimTreeWriter/fastsim_val.root" ]  )
    ##compareHistograms( [ "../../CMSSW/CMSSW_7_3_0/src/validateScaling/closure/fullsim.root", "../../CMSSW/CMSSW_7_3_0/src/validateScaling/closure/fastsim.root", "../../CMSSW/CMSSW_7_3_0/src/validateScaling/closure/fastsim_mod.root" ] )

    compareHistograms( [ "../../CMSSW/CMSSW_7_3_0/src/validateScaling/closure/fullsim.root", "../../CMSSW/CMSSW_7_3_0/src/validateScaling/closure/fastsim.root", "../../CMSSW/CMSSW_7_3_0/src/validateScaling/closure/fastsim_mod.root" ] )

    compareHistograms( [ "../../CMSSW/CMSSW_7_3_0/src/validateScaling/closureAllE/fullsim.root", "../../CMSSW/CMSSW_7_3_0/src/validateScaling/closureAllE/fastsim.root", "../../CMSSW/CMSSW_7_3_0/src/validateScaling/closureAllE/fastsim_mod.root" ] )

    compareHistograms( [ "../../CMSSW/CMSSW_7_3_0/src/validateScaling/enableVtxSmearing/fullsim.root", "../../CMSSW/CMSSW_7_3_0/src/validateScaling/enableVtxSmearing/fastsim.root", "../../CMSSW/CMSSW_7_3_0/src/validateScaling/enableVtxSmearing/fastsim_mod.root" ] )

    compareHistograms( [ "../../CMSSW/CMSSW_7_3_0/src/validateScaling/tracker/fullsim.root", "../../CMSSW/CMSSW_7_3_0/src/validateScaling/tracker/fastsim.root", "../../CMSSW/CMSSW_7_3_0/src/validateScaling/tracker/fastsim_mod.root" ] )

