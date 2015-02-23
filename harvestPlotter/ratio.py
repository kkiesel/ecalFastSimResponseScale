import ROOT
from math import sqrt

def randomName():
    """
    Generate a random string. This function is useful to give ROOT objects
    different names to avoid overwriting.
    """
    from random import randint
    from sys import maxint
    return "%x"%(randint(0, maxint))


class Ratio:
    def __init__( self, title, numerator, denominator, sysHisto=None ):
        if isinstance( numerator, ROOT.TProfile ):
            numerator = numerator.ProjectionX( randomName() )
        if isinstance( denominator, ROOT.TProfile ):
            denominator = denominator.ProjectionX( randomName() )

        self.title = title
        self.numerator = numerator
        self.denominator = denominator
        self.sysHisto = sysHisto
        self.ratio = numerator.Clone( randomName() )
        self.ratioSys = denominator.Clone( randomName() )
        self.totalUncert = denominator.Clone( randomName() )
        self.allowUnsymmetricYaxis_ = False

    def calculateRatio( self ):
        for bin in range(self.numerator.GetNbinsX()+2):
            if self.denominator.GetBinContent(bin):
                self.ratio.SetBinContent( bin, self.numerator.GetBinContent(bin) / self.denominator.GetBinContent(bin) )
                self.ratio.SetBinError(   bin, self.numerator.GetBinError(bin)   / self.denominator.GetBinContent(bin) )
                if self.sysHisto:
                    self.ratioSys.SetBinContent( bin, self.sysHisto.GetBinContent(bin) / self.denominator.GetBinContent(bin) )
                    self.ratioSys.SetBinError(   bin, self.sysHisto.GetBinError(bin)   / self.denominator.GetBinContent(bin) )
                    combinedError = sqrt( self.sysHisto.GetBinError(bin)**2 + self.denominator.GetBinError(bin)**2 )
                else:
                    self.ratioSys.SetBinContent( bin, 1 )
                    self.ratioSys.SetBinError( bin, 0 )
                    combinedError = self.denominator.GetBinError(bin)
                self.totalUncert.SetBinError( bin, combinedError/self.denominator.GetBinContent(bin) )
                self.totalUncert.SetBinContent( bin, 1 )
            elif self.numerator.GetBinContent(bin):
                self.ratio.SetBinError( bin, 1.*self.numerator.GetBinError(bin)/self.numerator.GetBinContent(bin) )
                self.ratio.SetBinContent(bin,0)

    def allowUnsymmetricYaxis( self, allow=True ):
        self.allowUnsymmetricYaxis_ = allow

    def getYrange( self ):
        # If no minimum or maximum is specified, choose a minimum from 0 to .5
        # and a maximum from 1.5 to 50
        yMin = 0
        minimum = self.ratio.GetBinContent(self.ratio.GetMaximumBin())
        for bin in range( self.ratio.GetNbinsX()+2 ):
            minInBin = self.ratio.GetBinContent(bin)
            if minInBin < minimum and minInBin > 0:
                minimum = minInBin
        yMin = minimum*.95

        from math import ceil
        yMax = min( max(1.5, ceil(self.ratio.GetBinContent(self.ratio.GetMaximumBin()))), 50 )
        yMax = self.ratio.GetBinContent(self.ratio.GetMaximumBin())*1.05


        yValues = [ self.ratio.GetBinContent(bin) for bin in range( self.ratio.GetNbinsX()+2 ) ]
        yValues = filter( lambda a: a != 0, yValues )
        if self.allowUnsymmetricYaxis_:
            return yMin, yMax
        else:
            yValuesAbsDiff = [ abs(x-1) for x in yValues ]
            yValuesAbsDiff.sort()
            maxYDiff = yValuesAbsDiff[-1]
            y = 0.5
            if maxYDiff < 0.05:
                y = 0.05
            if maxYDiff < 0.01:
                y = 0.01

            return 1-y, 1+y


    def draw( self ):
        self.calculateRatio()

        yMin, yMax = self.getYrange()

        # Set ratio properties
        for hist in [ self.ratio, self.ratioSys, self.totalUncert ]:
            hist.GetYaxis().SetNdivisions(2, 0, 2)
            hist.SetTitleOffset(.9, "Y")
            hist.SetYTitle( self.title )
            hist.SetMinimum( yMin )
            hist.SetMaximum( yMax )

        self.ratioSys.SetFillStyle(3254)
        self.ratioSys.SetMarkerSize(0)
        self.ratioSys.SetFillColor(46)

        self.totalUncert.SetFillStyle(3002)
        self.totalUncert.SetMarkerSize(0)
        self.totalUncert.SetFillColor(2)
        self.totalUncert.SetLineColor(0)


        # Delete label and title of all histograms in the current pad
        for ding in ROOT.gPad.GetListOfPrimitives():
            if isinstance( ding, ROOT.TH1 ) or isinstance( ding, ROOT.THStack ):
                xaxis = ding.GetXaxis()
                xaxis.SetNdivisions(0,50,0, True)
                xaxis.SetLabelSize(0)
                xaxis.SetLabelColor(0)
                xaxis.SetLabelOffset(1000)
                xaxis.SetTitle("")
                xaxis.SetTitleColor(0)
                xaxis.SetTitleSize(0)


        csf = 0.2 # the ratio in which the pad is splitted
        ROOT.gPad.SetBottomMargin( csf + (1-csf)*ROOT.gPad.GetBottomMargin() - csf*ROOT.gPad.GetTopMargin() )
        rPad = ROOT.TPad( "rPad", "ratio", 0, 0, 1, 1 )
        rPad.SetTopMargin( (1-csf) - (1-csf)*rPad.GetBottomMargin() + csf*rPad.GetTopMargin() )
        rPad.SetFillStyle(3955)
        rPad.Draw()
        rPad.cd()
        rPad.SetLogy(0)

        self.totalUncert.Draw("e2")
        if self.sysHisto:
            self.ratioSys.Draw("e2 same")
        self.ratio.Draw("e0 same")
        if yMin < 1 and yMax > 1:
            oneLine = ROOT.TLine()
            oneLine.SetLineStyle(2)
            axis = self.ratio.GetXaxis()
            oneLine.DrawLine( self.ratio.GetBinLowEdge( axis.GetFirst() ), 1.0,
                self.ratio.GetBinLowEdge( axis.GetLast() ), 1.0 )


