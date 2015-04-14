#!/usr/bin/env python2

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(54)

#ROOT.TH1.AddDirectory( ROOT.kFALSE )

# x - y - z == E - eta - resp

def getH( fname, oname ):
    f = ROOT.TFile( fname )
    obj = f.Get( oname )
    #obj.SetDirectory(0)
    obj = ROOT.gROOT.CloneObject( obj )
    return obj

def main():
    fastFileName = "../3d_fast.root"
    fullFileName = "../3d_full.root"
    path = "ecalScaleFactorCalculator/"

    for hname in [
#        "sEtaEtaVsEVsEta_woGaps",
#        "sEtaPhiVsEVsEta_woGaps",
        "sPhiPhiVsEVsEta_woGaps"
        ]:



        hfast3d = getH( fastFileName, path+hname )
        hfull3d = getH( fullFileName, path+hname )
        var = hfast3d.GetZaxis().GetTitle()
        hfast3d.Scale( hfull3d.GetEntries() / hfast3d.GetEntries() )

        hfast1d = hfast3d.ProjectionZ("fast")
        hfull1d = hfull3d.ProjectionZ("full")

        hfast1d.GetXaxis().SetRangeUser(0, 0.05 )
        hfast1d.SetLineColor(2)

        hfast1d.Draw()
        hfull1d.Draw("same")
        ROOT.gPad.SaveAs("plots/%s/total.pdf"%hname)


        # E - eta plot

        hfast2d = hfast3d.Project3D( "xyfast" )
        hfull2d = hfull3d.Project3D( "xyfull" )

        for h in hfast2d, hfull2d:
            h.Rebin2D( 1, 1 )

        hfull2d.Divide( hfast2d )
        hfull2d.GetYaxis().SetRangeUser(10,100)

        hfull2d.SetMaximum(2)
        hfull2d.SetTitle( "#LT%s#GT_{full} / #LT%s#GT_{fast}"%(var,var))

        hfull2d.Draw("colz")
        ROOT.gPad.SaveAs("plots/%s/E-Eta.pdf"%hname)

        for energy in [ 30, 60, 90 ]:
            for etaRange in [ (0.1,0.3), (.9, 1.3), (2, 2.4 ) ]:
                for h in hfast3d, hfull3d:
                    h.GetXaxis().SetRangeUser( energy, energy )
                    h.GetYaxis().SetRangeUser( *etaRange )
                hfast1d = hfast3d.Project3D( "zfast")
                hfull1d = hfull3d.Project3D( "zfull")

                hfast1d.GetXaxis().SetRangeUser(0.005, 0.035 )
                hfast1d.Scale( hfull1d.Integral() / hfast1d.Integral() )
                hfast1d.SetTitle( "E=%sGeV, %s<#eta<%s"%(energy, etaRange[0], etaRange[1]) )
                hfast1d.SetLineColor(2)
                hfast1d.SetMaximum( max( hfast1d.GetMaximum(), hfull1d.GetMaximum())+ROOT.gStyle.GetHistTopMargin() )

                hfast1d.Draw()
                hfull1d.Draw("same")
                ROOT.gPad.SaveAs("plots/%s/e%s_eta%s.pdf"%(hname,energy, etaRange[0]))






        del hfast3d, hfull3d


if __name__ == "__main__":
    main()
