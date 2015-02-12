import ROOT
from treeFunctions import myLegend

class Multihisto:
	def __init__( self ):
		self.orderByIntegral = True
		self.histos = []
		self.histosToStack = []
		self.minimum = None
		self.maximum = None
		self.leg = myLegend(.7,.7,.93,.92)

	def setMinimum( self, m ):
		self.minimum = m

	def setMaximum( self, m ):
		self.maximum = m

	def addHisto( self, singleHisto, label=None, toStack=False, draw="hist" ):
		if toStack:
			self.histosToStack.append( ( singleHisto, label, draw ) )
		else:
			self.histos.append( ( singleHisto, label, draw ) )

	def stackHistos( self ):
		if not self.histosToStack:
			return
		if self.orderByIntegral:
			# sort histograms first
			self.histosToStack.sort( key=lambda x: x[0].Integral() )

		stack = ROOT.THStack()
		stack.SetTitle( ";%s;%s"%(self.histosToStack[0][0].GetXaxis().GetTitle(),self.histosToStack[0][0].GetYaxis().GetTitle()) )
		for h in self.histosToStack:
			h[0].SetFillColor( h[0].GetLineColor() )
			h[0].SetLineColor(1)
			stack.Add( h[0] )

		self.stack = stack
		return stack

	def GetMinimum( self, histos ):
		values = []

		for hist, label, draw in histos + self.histosToStack:
			if not isinstance( hist, ROOT.THStack ):
				values.append( hist.GetMinimum(0) )

		return min( values )

	def GetMaximum( self, histos ):
		return max( [ x[0].GetMaximum() for x in histos ] )

	def fillLegend( self ):
		for hist, label, draw in self.histos:
			if label:
				legendOption = "l"
				if "p" in draw:
					legendOption = "pe"
				elif "hist" in draw:
					legendOption = "l"
				elif draw == "e2":
					legendOption = "f"
				if label == "Simulation":
					legendOption = "lep"
				if label == "Prediction":
					legendOption = "lfe"
				self.leg.AddEntry( hist, label, legendOption )
		for hist, label, draw in reversed(self.histosToStack):
			if label:
				self.leg.AddEntry( hist, label, "f" )

	def Draw( self  ):
		if not self.histos and not self.histosToStack:
			print "No histogram added"
			return

		stack = self.stackHistos()
		self.stack = stack

		histosToDraw = self.histos
		if stack:
			histosToDraw = [(stack, "", "hist" )] + self.histos

		maximum = self.GetMaximum( histosToDraw )
		minimum = self.GetMinimum( histosToDraw )
		if not ROOT.gPad or ROOT.gPad.GetLogy():
			maximum = 2.5*maximum
			minimum = 0.3*minimum
		else:
			maximum = maximum + (maximum-minimum)*.1
			minimum = minimum - (maximum-minimum)*.1

		if self.maximum != None:
			maximum = self.maximum
		if self.minimum != None:
			minimum = self.minimum

		histosToDraw[0][0].SetMaximum( maximum )
		histosToDraw[0][0].SetMinimum( minimum )

		histosToDraw[0][0].Draw( histosToDraw[0][2] )
		for hist, label, draw in histosToDraw[1:]:
			hist.Draw("same %s"%draw)


		self.fillLegend()
		if self.leg.GetNRows():
			self.leg.Draw()

