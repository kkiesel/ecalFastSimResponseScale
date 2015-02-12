from math import sqrt

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import Styles
import argparse
from sys import stdout
from prettifyFunctions import *

# To use user defined help message, sys.arv has to be sent to python and not
# to TApplication.

Styles.tdrStyle()
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libTreeObjects.so")

def rebin2D( oldHist, xList, yList ):
	"""Rebins a TH2F.
	oldHist: input TH2F
	xList: list which mark the bin edges of the new xaxis
	yList: list which mark the bin edges of the new yaxis
	"""
	import ROOT
	import array
	newHist = ROOT.TH2F( oldHist.GetName()+randomName(), "%s;%s;%s"%(oldHist.GetTitle(),oldHist.GetXaxis().GetTitle(),oldHist.GetYaxis().GetTitle()), \
			len(xList)-1, array.array('d', xList), \
			len(yList)-1, array.array('d', yList) )
	newHist.Sumw2()
	for i in range(oldHist.GetNbinsX()+2):
		for j in range(oldHist.GetNbinsY()+2):
			# newContent = newContent + oldContent
			# newError = sqrt( newError**2 + oldError**2 )
			newBin = newHist.FindBin( oldHist.GetXaxis().GetBinCenter(i), oldHist.GetYaxis().GetBinCenter(j) )
			newHist.SetBinContent( newBin, oldHist.GetBinContent(i,j)+ newHist.GetBinContent(newBin) )
			newHist.SetBinError( newBin, sqrt(newHist.GetBinError(newBin)**2 + oldHist.GetBinError( i, j )**2) )
	return newHist

def createHistoFromTree2D(tree, variable, weight, nBinsX=[], nBinsY=[] ):
	from ROOT import TH2F
	import array
	name = randomName()
	if isinstance( nBinsX, list ):
		xBins = array.array( 'd', nBinsX )
		if isinstance( nBinsY, list ):
			yBins = array.array( 'd', nBinsY )
			result = TH2F( name, variable, len(nBinsX)-1, xBins, len(nBinsY)-1, yBins )
			result.Sumw2()
			tree.Draw("%s>>%s"%(variable, name), weight, "goff")
			#result.Scale(1, "width")
	yLabel, xLabel = variable.split(":")
	result.SetTitle(";%s;%s"%( getAxisTitle( xLabel ), getAxisTitle( yLabel ) ) )
	return result

def createHistoFromTree(tree, variable, weight="", nBins=20, firstBin=None, lastBin=None ):
	"""
	tree: tree to create histo from
	variable: variable to plot (must be a branch of the tree)
	weight: weights to apply (e.g. "var1*(var2 > 15)" will use weights from var1 and cut on var2 > 15
	nBins, firstBin, lastBin: number of bins, first bin and last bin (same as in TH1F constructor)
	nBins: if nBins is a list, and to a int, a user binned plot will be generated
	returns: histogram
	"""
	from ROOT import TH1F
	name = randomName()
	if isinstance( nBins, list ):
		import array
		xBins = array.array('d', nBins )
		result = TH1F(name, variable, len(nBins)-1, xBins)
		result.Sumw2()
		tree.Draw("%s>>%s"%(variable, name), weight, "goff")
		result.Scale(1,"width")
	elif firstBin==None and lastBin==None:
		import ROOT
		tree.Draw("%s>>%s(%s,,)"%(variable,name,nBins), weight, "goff")
		result = ROOT.gDirectory.Get( name )
		if isinstance( result, ROOT.TTree ):
			print "Warning, no entries"
			return ROOT.TH1F()
		result.Sumw2() # applying the errors here is perhaps not entirely correct
	else:
		result = TH1F(name, variable, nBins, firstBin, lastBin)
		result.Sumw2()
		tree.Draw("%s>>%s"%(variable, name), weight, "goff")
	return result

def readTree( filename, treename = "susyTree" ):
	"""
	filename: name of file containing the tree
	treename: name of the tree
	returns: TChain Object
	"""
	import ROOT
	tree = ROOT.TChain( treename )
	tree.AddFile( filename )
	return tree

def readHisto( filename, histoname="eventNumbers" ):
	"""
	filename: name of file containing the histogram
	histoname: name of the histogram
	returns: Object with the given name, should be a histogram
	"""
	import ROOT
	import os
	if not os.path.isfile(filename):
		print( "File %s does not exist"%filename)
	f = ROOT.TFile( filename )
	# with this command, the histo is cloned to root and not deleted after end
	# of function
	ROOT.gROOT.cd()
	histo = f.Get( histoname )
	# the name +Clone is only temporaly, since for TH1::Clone a different name is expected
	if not isinstance( histo, ROOT.TH1 ):
		print 'Histogram "%s" not found in file "%s".'%(histoname, filename)
		return
	if not histo.GetSumw2():
		histo.Sumw2()
	histo.SetName(histoname+"Clone")
	histo = histo.Clone( histoname )
	return histo

def getAxisBinningFromHisto( axis ):
	result = []
	xbins = axis.GetXbins()
	if xbins.GetSize(): # variable binning
		for i in range( xbins.GetSize() ):
			result.append( xbins[i] )
	else:
		result = [ axis.GetBinLowEdge(1) ]
		width = axis.GetBinWidth(1)
		for i in range( 1, axis.GetNbins()+1 ):
			result.append( result[-1]+width )
	return result

def appendFlowBin2D( oldHist, firstBinX=0, lastBinX=0, firstBinY=0, lastBinY=0 ):
	# this is still buggy. If the over/underflow bin is too small, it is not drawn.
	oldXList = getAxisBinningFromHisto( oldHist.GetXaxis() )
	oldYList = getAxisBinningFromHisto( oldHist.GetYaxis() )
	if firstBinX:
		oldXList.insert( 0, oldXList[0] - firstBinX )
	if lastBinX:
		oldXList.append( oldXList[-1] + lastBinX )
	if firstBinY:
		oldYList.insert( 0, oldYList[0] - firstBinY )
	if lastBinY:
		oldYList.append( oldYList[-1] + lastBinY )
	import array
	newX = array.array("d", oldXList )
	newY = array.array("d", oldYList )
	import ROOT
	title = oldHist.GetTitle()+";"+oldHist.GetXaxis().GetTitle()+";"+oldHist.GetYaxis().GetTitle()
	newHist = ROOT.TH2D( randomName(), title, len(newX)-1, newX, len(newY)-1, newY )
	for i in range(oldHist.GetNbinsX()+2):
		for j in range(oldHist.GetNbinsY()+2):
			newBin = newHist.FindBin( oldHist.GetXaxis().GetBinCenter(i), oldHist.GetYaxis().GetBinCenter(j) )
			newHist.SetBinContent( newBin, oldHist.GetBinContent(i,j) )
			newHist.SetBinError( newBin, oldHist.GetBinError(i,j) )
	return newHist

def appendOverflowBin( oldHist, overflow ):
	"""Append the overflow bin to a histogram.
	oldHist: input histogram
	overflow: size of the overflow bin, which should be drawn
	"""
	import ROOT
	import array
	oldArray = oldHist.GetXaxis().GetXbins()
	newArray = array.array('d', [0]*(oldArray.GetSize()+1))
	for i in range( oldArray.GetSize() ):
		newArray[i] = oldArray[i]
	newArray[-1] = newArray[-2] + overflow
	newHist = ROOT.TH1F( randomName(), "", len(newArray)-1, newArray )

	for bin in range(1, len(newArray)):
		newHist.SetBinContent(bin, oldHist.GetBinContent(bin ))
		newHist.SetBinError(bin, oldHist.GetBinError(bin ))
	return newHist

def getXMinXMax( histo_list ):
	"""Find minimum and maximum on the x-axis for a list of histograms.
	histo_lists: list of histograms
	"""
	from sys import maxint
	mini = maxint
	maxi = -maxint
	for histo in histo_list:
		mini = min( mini, histo.GetBinLowEdge(1) )
		lastBin = histo.GetNbinsX()
		maxi = max( maxi, histo.GetBinLowEdge(lastBin)+histo.GetBinWidth(lastBin) )
	return mini, maxi

def getHisto( tree, plot, cut="1", overflow=0, weight="weight", color=1, nBins=None, firstBin=None, lastBin=None, fillEmptyBins=False, appendOverflowBin=True ):
	"""Creates a histogram and apply the axis settings
	cut: cutstring applied to the tree
	overflow: size of the overflow bin, if overflow>0
	"""
	if nBins == None:
		nBins = 20
		label, unit, binning = readAxisConf( plot )
	else:
		label, unit, binning = readAxisConf( plot )
		binning = nBins


	if firstBin != None and lastBin != None:
		if "Length$(" in plot or "nVertex" == plot:
			firstBin -= .5
			lastBin += .5
			nBins = int(lastBin-firstBin)

	if "AddVariables"+tree.GetName() in tree.GetFile().GetListOfKeys():
		tree.AddFriend( "AddVariables"+tree.GetName(), tree.GetFile().GetName() )

	if binning:
		histo = createHistoFromTree( tree, plot, "%s*(%s)"%(weight, cut), nBins=binning)
	else:
		histo = createHistoFromTree( tree, plot, "%s*(%s)"%(weight, cut), nBins=nBins, firstBin=firstBin, lastBin=lastBin )
	if overflow > 0:
		histo = appendOverflowBin(histo, overflow)

	# for syncronnization with christian
	#fillEmptyBins = False

	if fillEmptyBins:
		# If a neigbour of a empty bin is filled, the error of the bin will be
		# set to the poisson error for 0 times the weight.
		poissonZeroError = 1.14787446444
		weightH = createHistoFromTree( tree, "weight", "weight", 100 )
		weight = weightH.GetMean()

		for bin in range(1, histo.GetNbinsX()+2):
			# if the bin left or right is not empty but the bin itself, set the error
			if not histo.GetBinContent( bin ) and ( histo.GetBinContent( bin-1 ) or histo.GetBinContent( bin+1 ) ) and histo.GetBinWidth(bin):
				histo.SetBinError( bin, poissonZeroError*weight / histo.GetBinWidth(bin) )
			pass

	if appendOverflowBin:
		lastBin = histo.GetNbinsX()
		histo.SetBinError( lastBin, sqrt(histo.GetBinError(lastBin)**2+histo.GetBinError(lastBin+1)**2) )
		histo.SetBinContent( lastBin, histo.GetBinContent(lastBin) + histo.GetBinContent(lastBin+1) )
		histo.SetBinContent( lastBin+1,0 )
		histo.SetBinError( lastBin+1, 0 )
	histo.SetLineColor( color )
	histo.SetMarkerColor( color )
	histo.SetLineWidth(2)

	histo.SetTitle( getHistoTitle( histo, plot, label, unit, binning ) )
	return histo

def tryAddFriend( tree, nameAppendix="AddVariables" ):
	newTreename = tree.GetName()+nameAppendix
	file = tree.GetFile()
	if newTreename in file.GetListOfKeys():
		tree.AddFriend( newTreename, file.GetName() )

	return tree


def getHists( filenames, plot="met", cut="1", treeName="photonTree" ):
	endHist = None
	for filename in filenames:
		tree = readTree( filename, treeName )
		tree = tryAddFriend( tree )

		hist = getHisto( tree, plot, color=1, fillEmptyBins=not ("PhotonHad" in filename), cut=cut )
		if endHist:
			endHist.Add( hist )
		else:
			endHist = hist

	return endHist

def setRelativeUncertainty( hist, uncert ):
	hist = hist.Clone( randomName() )
	for bin in range( hist.GetNbinsX()+2 ):
		hist.SetBinError( bin, uncert * hist.GetBinContent(bin) )
	return hist


def getQCDErrorHisto( tree, plot, cut="1", overflow=0, nBins=20, firstBin=None, lastBin=None ):
	"""Applies w_qcd +- w_qcd_error for qcd error propagation.
	The returned histo's content is the mean, the error are the shifts up and down.
	"""
	if firstBin != None and lastBin != None:
		if "Length$(" in plot or "nVertex" == plot:
			firstBin -= .5
			lastBin += .5
			nBins = int(lastBin-firstBin)

	label, unit, binning = readAxisConf( plot )
	if binning:
		histoUp = createHistoFromTree( tree, plot, "(weight*(w_qcd+w_qcd_error))*(%s)"%(cut), nBins=binning)
		histoDown = createHistoFromTree( tree, plot, "(weight*(w_qcd-w_qcd_error))*(%s)"%(cut), nBins=binning)
	else:
		histoUp = createHistoFromTree( tree, plot, "(weight*(w_qcd+w_qcd_error))*(%s)"%(cut), nBins=nBins, firstBin=firstBin, lastBin=lastBin )
		histoDown = createHistoFromTree( tree, plot, "(weight*(w_qcd-w_qcd_error))*(%s)"%(cut), nBins=nBins, firstBin=firstBin, lastBin=lastBin )
	if overflow > 0:
		histoUp = appendOverflowBin(histoUp, overflow)
		histoDown = appendOverflowBin(histoDown, overflow)

	outHisto = histoUp.Clone( randomName() )
	for bin in range( histoUp.GetNbinsX()+1 ):
		up = histoUp.GetBinContent(bin)
		down = histoDown.GetBinContent(bin)
		outHisto.SetBinContent( bin, (up+down)/2 )
		outHisto.SetBinError( bin, (up-down)/2 )
	return outHisto

def getAxisTitle( plot ):
	objectReplacement = {
			"photons": "#gamma",
			"electrons": "e",
			"muons": "#mu"
		}
	variableReplacement = {
			"phi":"#phi",
			"eta":"#eta",
			"pt":"p_{T ",
			"ptJet":"p_{T* ",
			"sigmaIetaIeta":"#sigma_{i#etai#eta",
			"hadTowOverEm":"H/E",
			"chargedIso": "Iso^{#pm}",
			"neutralIso": "Iso^{0}",
			"photonIso": "Iso^{#gamma}"
		}
	import re
	objVarExpr = "([a-zA-Z]+)\[{0,1}(\d*)\]{0,1}\.([_a-zA-Z]+)" # matches eg photon.pt
	if "Length$(" in plot:
		obj, nObj, var = re.match( "Length\$\(%s\)"%objVarExpr, plot ).groups()
		obj = reduce(lambda x, y: x.replace(y, objectReplacement[y]), objectReplacement, obj )
		label = "N_{%s}"%obj
	elif "." in plot and re.match( objVarExpr, plot ):
		obj, nObj, var = re.match( objVarExpr, plot ).groups()
		var = reduce(lambda x, y: x.replace(y, variableReplacement[y]), variableReplacement, var )
		obj = reduce(lambda x, y: x.replace(y, objectReplacement[y]), objectReplacement, obj )
		if "_{" not in var:
			var += "_{"
		if nObj:
			label = "%s%s.%s}"%(var,int(nObj)+1,obj)
		else:
			label = "%s%s}"%(var,obj)
	elif plot == "ht":
		label = "H_{T}"
	else:
		label = plot
	return label

def getHistoTitle( histo, plot, label, unit, binning ):
	ytitle = "Entries"
	if not label:
		label = getAxisTitle( plot )
	if binning:
		if unit:
			ytitle+= " / GeV"
		if unit:
			label+= " (%s)"%unit
	else:
		if histo.GetBinWidth(1) != 1:
			ytitle+= " / {}".format(roundToSignificantDigits(histo.GetBinWidth(1),2))
		if unit:
			label+= " (%s)"%unit
			ytitle+= " %s"%unit
	return ";%s;%s"%( label, ytitle )

def addHistos( histos, scales=None ):
	"""
	add several histos with different scales to one single histo
	histos: list of histograms
	scales: list of scales (of same size as histos)

	return: single histogram
	"""
	if not scales:
		scales = [1]*len(histos)
	else:
		if len(histos) != len(scales):
			print "Histos and scales have to have same dimension"
	sumHist = histos[0].Clone( randomName() )
	sumHist.Scale( scales[0] )
	for i, h in enumerate(histos[1:]):
		sumHist.Add( h, scales[i+1] )
	return sumHist

def divideHistos( numerator, denominator, bayes=False ):
	"""Divide numerator by denominator."""
	option = "B" if bayes else ""
	resultHisto = numerator.Clone( randomName() )
	resultHisto.Divide( numerator, denominator, 1,1, option )
	return resultHisto

def applyFakeRateEWK( histo, fakeRate=None, fakeRateError=None ):
	"""Apply the hard-coded fake rate and error to a histogram. The prediction
	is returned."""

	if not fakeRate:
		fakeRate = 0.0084
	if not fakeRateError:
		fakeRateError = 0.0008 # syst and stat

	# correct fake rate, if it is estimated with yutaros method
	fakeRateError = fakeRateError / (1-fakeRate)**2
	fakeRate = fakeRate / ( 1 - fakeRate )

	for i in range( histo.GetNbinsX() +1 ):
		content = histo.GetBinContent(i)
		histo.SetBinContent( i,fakeRate * content )
		# sigma_{ef} = sqrt( ( e*e_f )**2 + ( e_e*f )**2 )
		histo.SetBinError( i, sqrt( (fakeRateError*content)**2 + (fakeRate*histo.GetBinError(i))**2 ) )
	return histo

def getIntAndError( h, bin ):
	err = ROOT.Double()
	con = h.IntegralAndError( bin, bin, err, "width" )
	return con, err


