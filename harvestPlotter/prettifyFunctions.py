"""Here functions are collected which prettify the plots and the way of coding."""

def isValidFile( fileName ):
	"""This function can be used for argparse to check if the file exists."""
	import argparse
	import os
	if not os.path.isfile( fileName ):
		raise argparse.ArgumentTypeError( "File %s does not exist"%fileName )
	else:
		return fileName

def getDatasetAbbr( fileName, slim=True ):
	"""Parses the abbrevation for a sample from a root fileName."""
	import re

	datamatch = re.match("PhotonHad([ABCD]{0,1})_V.*", fileName.split("/")[-1])
	if datamatch:
		return "Data %s"%datamatch.groups()[0]

	prefix = "slim" if slim else ""
	match = re.match("%s(.*)_V.*_tree.root"%prefix, fileName.split("/")[-1] )
	if match:
		return match.groups()[0]
	else:
		import os
		return os.path.basename(fileName).replace(".root", "").replace("slim", "" )

def replaceListElementsBy( list_, search, replacement ):
	if len([x for x in list_ if x in search ]) == len(search):
		return [x for x in list_ if x not in search ]+[replacement]
	return list_

def mergeDatasetAbbr( datasetAbbrs ):
	datasetAbbrs = replaceListElementsBy( datasetAbbrs, ["GJets_200_400", "GJets_400_inf"], "GJets" )
	datasetAbbrs = replaceListElementsBy( datasetAbbrs, ["QCD_250_500", "QCD_500_1000", "QCD_1000_inf"], "QCD" )
	datasetAbbrs = replaceListElementsBy( datasetAbbrs, ["GJets", "QCD"], "AllQCD" )
	datasetAbbrs = replaceListElementsBy( datasetAbbrs, ["WJets_250_300", "WJets_300_400", "WJets_400_inf"], "W" )
	datasetAbbrs = replaceListElementsBy( datasetAbbrs, ["WGamma_50_130", "WGamma_130_inf"], "WGamma" )
	datasetAbbrs = replaceListElementsBy( datasetAbbrs, ["ZGammaNuNu", "ZGammaLL"], "ZGamma2" )
	datasetAbbrs = replaceListElementsBy( datasetAbbrs, ["Data A", "Data B", "Data C", "Data D"], "Data" )
	return datasetAbbrs


def datasetToLatex( datasetAbbr ):
	"""Translates the dataset name to a TLatex name"""
	sets = { "AllQCD": "#gammaJet+QCD",
			"GJets": "#gamma#text{Jets}",
			"TTJets": "t#bar{t}",
			"WJets": "W",
			"QCD": "QCD",
			"WGamma": "#gamma W",
			"TTGamma": "#gamma t#bar{t}",
			"ZGamma": "#gamma#text{Z}",
			"ZGammaLL": "#gammaZ#rightarrow#gammall",
			"ZGammaNuNu": "#gamma#text{Z}",
			"ZGammaNuNu3_400_inf": "Z"
			}
	try:
		return sets[datasetAbbr]
	except:
		# find out if it is a signal sample
		from re import match
		matchResult = match("([WB])_(\d{3,4})_(\d{3,4})_375", datasetAbbr)
		if matchResult:
			flavor, mg, mq = matchResult.groups()
			return "%sino %s,%s"%(flavor,mg,mq)
		else:
			return datasetAbbr

def createDatasetLabel( datasetAbbr ):
	import ROOT
	"""Creates sample info which then can be printed."""
	datasetLabel = ROOT.TPaveText(.4,.94,.6,1, "ndc")
	datasetLabel.SetFillColor(0)
	datasetLabel.SetBorderSize(0)
	datasetLabel.AddText( datasetToLatex(datasetAbbr) )
	return datasetLabel

def randomName():
	"""
	Generate a random string. This function is useful to give ROOT objects
	different names to avoid overwriting.
	"""
	from random import randint
	from sys import maxint
	return "%x"%(randint(0, maxint))

def roundToSignificantDigits(x, sig=2):
	"""Round number to 'sig' significant digits. If the number is large enough,
	just print the integer.
	"""
	from math import log10, floor
	if x >= 10**(sig-1):
		return int(round(x))
	if x>0:
		return round(x, sig-int(floor(log10(x)))-1)
	elif x<0:
		return round(-x, sig-int(floor(log10(-x)))-1)
	return x

def myLegend( x1, y1, x2=0,y2=0 ):
	import ROOT
	if x2==0 or y2==0:
		style = ROOT.gROOT.GetStyle("tdrStyle")
		x2 = 1 - style.GetPadRightMargin()
		y2 = 1 - style.GetPadTopMargin()
	leg = ROOT.TLegend(x1,y1,x2,y2)
	leg.SetFillColor(0)
	leg.SetBorderSize(0)
	return leg

configurationFileName="axis.cfg"
import ConfigParser
configuration = ConfigParser.SafeConfigParser()
configuration.read( configurationFileName )

def readAxisConf( plot, configurationFileName="axis.cfg" ):
	"""Read the configuration file for the axis.
	returns the label, the unit and the binning as list if avaible
	"""
	#import ConfigParser
	#configuration = ConfigParser.SafeConfigParser()
	#configuration.read( configurationFileName )
	#brackets are identified as sections, so they have to be deleted
	plot = plot.replace("[","").replace("]","")
	if not configuration.has_section( plot ):
		return "","",""
	label = configuration.get( plot, "label" )
	unit = configuration.get( plot, "unit" )
	binning = configuration.get( plot, "binning" )
	if binning:
		binning = map(float, binning.split(" "))
	else:
		binning = []
	return label, unit, binning

def getLumiWeight( datasetAbbr, nGenEvents, integratedLumi=19800, configName="dataset.cfg" ):
	import ConfigParser
	datasetConf = ConfigParser.SafeConfigParser()
	datasetConf.read( configName )
	if datasetConf.has_option( datasetAbbr, "crosssection" ):
		crosssection = datasetConf.getfloat( datasetAbbr, "crosssection" )
	else:
		import re
		# Signal sample?
		fileNameMatch = re.match("(.)_(\d+)_(\d+)_(\d+)", datasetAbbr )
		if fileNameMatch:
			nlsp, mg, mq, mNlsp = fileNameMatch.groups()
			xSectionFile = "/home/home4/institut_1b/kiesel/infos/Spectra_gsq_%s_8TeV.xsec"%nlsp
			f = open( xSectionFile )
			text = f.readlines()
			f.close()

			floatMatch = "\-{0,1}\d\.\d+e[+-]\d{2}"
			crosssection = None
			for t in text:
				matches = re.match(" (\d+)\s+%s\s+%s\s+(\d+)\s+(\d+) LO: (%s) \+ (%s) \- (%s) NLO: (%s) \+ (%s) \- (%s)\n"%(mg, mq, floatMatch,floatMatch,floatMatch,floatMatch,floatMatch,floatMatch), t )
				if matches:
					crosssection = float(matches.groups()[6])
			if not crosssection:
				raise NameError( "Signal not in xsection file %s."%xSectionFile )
		else:
			raise NameError( "Configuration for %s not found"%datasetAbbr )

	return 1. * integratedLumi * crosssection / nGenEvents

def manipulateSaveName( saveName ):
	"""Replace some charakters, so root nor unix have problems to read them."""
	#saveName = saveName.replace("/","VS")
	saveName = saveName.replace(" ","_")
	unallowedCharacters = ["{","}","(",")","#","|",".","[","]","/","$"]
	for char in unallowedCharacters:
		saveName = saveName.replace( char, "" )
	return saveName

def SaveAs( can, name, folder="plots", endings=["pdf"] ):
	"""Save ROOT.TCanvas in specified folder with a cleaned plot name."""
	for ending in endings:
		can.SaveAs( folder+"/"+manipulateSaveName( name )+"."+ending )

def SavePad( name, folder="plots", endings=["pdf"] ):
	import ROOT
	SaveAs( ROOT.gPad.GetCanvas(), name, folder, endings )

class PlotCaption:
	"""Creates the superscription for each plot, eg
	'19fb^{-1} sqrt{s)=8TeV #geq1#gamma #geq2jets'
	"""
	def __init__( self, x0=.96, y0=.965, treeName="photonTree", analysisInfo=True, option="ndc", signal=False, control=False ):
		import ROOT
		self.x0 = x0
		self.text = ROOT.TLatex( x0, y0, "" )
		self.text.SetTextSize(0.03)
		self.text.SetNDC()
		self.treeName = treeName
		if analysisInfo:
			self.addAnalysisInfo()
		if signal:
			self.signalCut()
		if control:
			self.controlCut()

	def addAnalysisInfo( self, lumi=19800, e=8 ):
		photonNameAppendix = ""
		if self.treeName == "photonTree":
			photonNameAppendix = "_{tight}"
		if self.treeName == "photonJetTree":
			photonNameAppendix = "_{loose}"
		elif self.treeName == "photonElectronTree":
			photonNameAppendix = "_{pixel}"
		defaultcuts = "#geq1#gamma%s,#geq2jets"%photonNameAppendix
		self.text.SetText( self.text.GetX(), self.text.GetY(), "%.1ffb^{-1} #sqrt{s}=%sTeV %s"%( lumi/1000., e, defaultcuts ) )

	def appendEnd( self, string ):
		newText = self.text.GetTitle() + string
		self.text.Clear()
		self.text.SetText( self.text.GetX(), self.text.GetY(), newText )

	def appendFront( self, string ):
		newText = string + self.text.GetTitle()
		self.text.Clear()
		self.text.SetText( self.text.GetX(), self.text.GetY(), newText )

	def controlCut( self ):
		self.appendEnd(",#slash{E}_{T}<100GeV")

	def signalCut( self ):
		self.appendEnd(",#slash{E}_{T}#geq100GeV")

	def Draw( self ):
		import ROOT
		shiftNDC = self.text.GetXsize() / ( ROOT.gPad.GetX2() - ROOT.gPad.GetX1() )
		self.text.SetX( self.x0-shiftNDC )
		self.text.Draw()

def readSignalPdfUncertainty( filename ):
	"""Read xsection and other informations for various signal MC from a file
	found at https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA3PrivateSignalMC2012
	The syntax is printf("%i %i %i %i %i %f %f\n",nevents,mgluino,msquark,mbino,mwino,xsecpdferrs,acceppdferrs
	returns list[ (point1, point2) ] = ( xsecpdferrs,acceppdferrs )
	"""
	f = open( filename )
	text = f.readlines()
	f.close()

	info = {}

	import re
	for t in text:
		matches = re.match("(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)", t ).groups()

		info[ (int(matches[2]), int(matches[1]) ) ] = ( float(matches[5]), float(matches[6]) )

	return info

def readSignalXSection( filename ):
	"""Read xsection and other informations for various signal MC from a file
	found at https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA3PrivateSignalMC2012
	The syntax is printf("%6.0f %6.0f %6.0f %6.0f %6.0f LO: %9.3e + %9.3e - %9.3e NLO: %9.3e + %9.3e - %9.3e\n", $nevents, $msquark, $mgluino, $mbino, $mwwino, $loxsec,$lohierr,$loloerr,$nloxsec,$nlohierr,$nloloerr);

	returns list[ (point1, point2) ] = ( sigmaNLO, errorUp, errorDown )
	"""
	f = open( filename )
	text = f.readlines()
	f.close()

	info = {}

	floatMatch = "\-{0,1}\d\.\d+e[+-]\d{2}"
	import re
	for t in text:
		matches = re.match(" (\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+) LO: (%s) \+ (%s) \- (%s) NLO: (%s) \+ (%s) \- (%s)\n"%(floatMatch,floatMatch,floatMatch,floatMatch,floatMatch,floatMatch), t ).groups()

		info[ (int(matches[1]), int(matches[2]) ) ] = ( float(matches[8]), float(matches[9]), float(matches[10]) )

	return info

def readSModelXsection( filename ):
	import re
	f = open( filename )

	xSections = {}

	for line in f.readlines():
		if line.startswith("#"): continue
		m, xsec, uncert = line.split(" ")
		xSections[int(m)] = (float(xsec), float(uncert) )
	f.close()

	return xSections

def getSaveNameFromDatasets( filenames ):
	return "".join( mergeDatasetAbbr( [ getDatasetAbbr(x) for x in filenames ] ) )

def shortName( filenames ):
	# just an alias
	return getSaveNameFromDatasets( filenames )

def correctTiksPlot( filename, scale1=True ):
	# TTexDump is not perfect now, so replace necessary things, like 10^{-2}
	#  in Math mode
	import fileinput, sys,re
	for line in fileinput.FileInput( filename, inplace=1 ):

		# each object with scale=0 is deleted
		line = re.sub( ".*scale=0[^\.].*", "", line )

		# expressions with 10^x will be in math mode now
		line = re.sub( "[^\$](10\^\{-?[0-9]+\})[^\$]", "{$\\1$}", line )

		# '%' without a leading \ will be replaced by \%
		line = re.sub( "(.*[^\\\\])%(.*)", "\\1\%\\2", line )

		# remove double '$$'
		line = line.replace( '$$', '$' )

		if scale1:
			# scale all objects to 1
			line = re.sub( "(.*scale=)(\d*\.\d*)([, ]+.*)", "\\1 1 \\3", line )

		# sys.stdout is redirected to the file
		sys.stdout.write(line)


class Infix:
	# definition of an Infix operator class
	# this recipe also works in jython
	# calling sequence for the infix is either:
	#  x |op| y
	# or:
	# x <<op>> y
	def __init__(self, function):
		self.function = function
	def __ror__(self, other):
		return Infix(lambda x, self=self, other=other: self.function(other, x))
	def __or__(self, other):
		return self.function(other)
	def __rlshift__(self, other):
		return Infix(lambda x, self=self, other=other: self.function(other, x))
	def __rshift__(self, other):
		return self.function(other)
	def __call__(self, value1, value2):
		return self.function(value1, value2)


# quadratic addition, usage: 3 |qPlus| 4
from math import sqrt
qPlus = Infix( lambda x,y: sqrt(x**2+y**2) )
qMinus = Infix( lambda x,y: sqrt(x**2-y**2) )

def integralAndError( h, binx1=0, binx2=-1, option="" ):
	import ROOT
	err = ROOT.Double(0)
	integral = h.IntegralAndError( binx1, binx2, err, option )
	return integral, err

def asimovSignificance( s, b, b_uncert=0 ):
	if not b:
		print "Background zero, can't calculate significance"
		return 0
	if b_uncert:
		return sqrt( 2 * ( (s+b)*log( (s+b)*(b+b_uncert**2)/(b**2+(s+b)*b_uncert**2) ) - b**2/b_uncert**2 * log( 1+ b_uncert**2*s/b/(b+b_uncert**2) ) ) )
	return sqrt( 2* ( (s+b)*log(1+s/b) - s ) )

def getPdfUncertaintiesAcceptance( inFilename ):
    out = {}
    with open( inFilename ) as f:
        for line in f:
            try:
                m1, m2, uncert = line.split(" ")
            except:
                m1, m2, uncert = 0, 0,0
            m1 = int(m1)
            m2 = int(m2)
            uncert = float(uncert)
            out[(m1,m2)] = uncert
    return out





