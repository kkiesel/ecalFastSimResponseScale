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
    if not h:
        print "Histogram {} not found in file {}".format(histoname, filename)
        return
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
    if "closure" in files[0]:
        processTex = "Single e^{#minus}"
    processTex += "          FullSim #color[2]{FastSim}"
    if len(files) > 2:
        processTex += " #color[4]{FastSim+Mod}"

    hists = []
    for file in files:
        h = readHist( file, "{}/{}".format( path, name ) )
        if not h:
            return
        if not round( h.GetEntries()) or not round( h.Integral() ):
            print "no entries in {}".format( name )
            return

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
    label.DrawLatexNDC( .01, .96, "#font[61]{CMS} #lower[-0.2]{#it{#scale[0.6]{#splitline{Private Work}{Simulation}}}}  "+processTex )

    if len(hists) == 2:
        r = ratio.Ratio( "Full/Fast  ", hists[0], hists[1] )
    else:
        r = ratio.Ratio( "Full/Mod  ", hists[0], hists[2] )
    r.draw()
    if len(hists)>2:
        processName = "modIncl_DQMscale_"+processName

    if name == "h_ele_PoPtrueVsEta":
        for bin in range( 1, r.ratio.GetNbinsX()+1 ):
            x = r.ratio.GetBinLowEdge( bin+1 )
            y = r.ratio.GetBinContent( bin )
            ey1 = r.ratio.GetBinError( bin )
            ey2 = r.totalUncert.GetBinError(bin)
            ey = ROOT.TMath.Sqrt( ey1**2 + ey2**2 )
            #if abs(y-1) < ey: y = 1.
            print "else if( genEta < %s ) { scale = %s; }"%(x,y)

    ROOT.gPad.GetCanvas().SaveAs("plots/%s_%s.pdf"%(processName, name ))

def compareHistograms( files, path="SimTreeProducer" ):
    names = getObjectNames( files[0], path )

    configuration = ConfigParser.SafeConfigParser()
    configuration.read( "histoDefinitions.cfg" )


    for name in names:
        draw( files, path, name, configuration )


if __name__ == "__main__":

    cmsswPath = "../../CMSSW/CMSSW_7_3_0/src/"

    ##compareHistograms( [ cmsswPath+"Analyzer/SimTreeWriter/fullsim_muchStat.root", cmsswPath+"Analyzer/SimTreeWriter/fastsim_muchStat.root", cmsswPath+"Analyzer/SimTreeWriter/fastsim_val.root" ]  )
    ##compareHistograms( [ cmsswPath+"validateScaling/closure/fullsim.root", cmsswPath+"validateScaling/closure/fastsim.root", cmsswPath+"validateScaling/closure/fastsim_mod.root" ] )
    #compareHistograms( [ cmsswPath+"validateScaling/enableVtxSmearing/fullsim.root", cmsswPath+"validateScaling/enableVtxSmearing/fastsim.root", cmsswPath+"validateScaling/enableVtxSmearing/fastsim_mod.root" ] )
    #compareHistograms( [ cmsswPath+"validateScaling/tracker/fullsim.root", cmsswPath+"validateScaling/tracker/fastsim.root", cmsswPath+"validateScaling/tracker/fastsim_mod.root" ] )


    # closure single E
    #compareHistograms( [ cmsswPath+"validateScaling/closure/fullsim.root", cmsswPath+"validateScaling/closure/fastsim.root", cmsswPath+"validateScaling/closure/fastsim_mod.root" ] )

    # closure full E
    #compareHistograms( [ cmsswPath+"validateScaling/closureAllE/fullsim.root", cmsswPath+"validateScaling/closureAllE/fastsim.root", cmsswPath+"validateScaling/closureAllE/fastsim_mod.root" ] )

    # Hgg
    #compareHistograms( [ cmsswPath+"harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValH130GGgluonfusion_13__official_FullSim.root", cmsswPath+"harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValH130GGgluonfusion_13__official_FastSim.root" ], "DQMData/Run 1/EgammaV/Run summary/PhotonValidator/Photons" )

    # Hgg, incl mod
    #compareHistograms( [ cmsswPath+"harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValH130GGgluonfusion_13__official_FullSim.root", cmsswPath+"harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValH130GGgluonfusion_13__official_FastSim.root", cmsswPath+"testRunTheMatrix/Hgg/DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root" ], "DQMData/Run 1/EgammaV/Run summary/PhotonValidator/Photons" )

    # Zee
    #compareHistograms( [ cmsswPath+"harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValZEE_13__official_FullSim.root", cmsswPath+"harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValZEE_13__official_FastSim.root" ], "DQMData/Run 1/EgammaV/Run summary/ElectronMcSignalValidator" )

    # Zee, incl mod
    #compareHistograms( [ cmsswPath+"harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValZEE_13__official_FullSim.root", cmsswPath+"harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValZEE_13__official_FastSim.root", cmsswPath+"testRunTheMatrix/fast/DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO_mod1.root" ], "DQMData/Run 1/EgammaV/Run summary/ElectronMcSignalValidator" )

    # Zee, incl dqm scaling mod
    #compareHistograms( [ cmsswPath+"harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValZEE_13__official_FullSim.root", cmsswPath+"harvest/DQM_V0001_R000000001__CMSSW_7_3_0__RelValZEE_13__official_FastSim.root", cmsswPath+"testRunTheMatrix/fast2/DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root" ], "DQMData/Run 1/EgammaV/Run summary/ElectronMcSignalValidator" )

    # different steps
    compareHistograms( [ cmsswPath+"Analyzer/ResponseOnDifferentLevels/full.root", cmsswPath+"Analyzer/ResponseOnDifferentLevels/fast.root" ], "ana" )

