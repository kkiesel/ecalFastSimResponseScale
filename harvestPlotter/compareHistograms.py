import ROOT
import ratio
import style
ROOT.gROOT.SetBatch()


def readHisto( filename, histoname ):
    f = ROOT.TFile( filename )
    h = f.Get( histoname )
    h = ROOT.gROOT.CloneObject( h )
    return h

def getObjectNamesRec( file1, path ):
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

    return outList

def getFirstDifferentPosition( names ):
    # Finds the position in the string where the strings in 'names' differ
    for i, letters in enumerate(zip(*names)):
        if list(letters) != [letters[0]]*len(letters):
            return i


def draw( h1, h2, savename ):
    c = ROOT.TCanvas()
    if not h1.GetEntries() and not h2.GetEntries(): return

    h2.Scale(2.5)

    for h in h1,h2:
        h.SetMarkerSize(0)
    h2.SetLineColor(2)
    h1.SetMaximum( max( h1.GetMaximum(), h2.GetMaximum() ) )

    h1.Draw("hist e")
    h2.Draw("same e")
    r = ratio.Ratio( "h1/h2", h1, h2 )
    r.draw(0.5,1.5)
    ROOT.gPad.GetCanvas().SaveAs( "plots/%s.pdf"%savename)

def compareDQMHistograms( filename1, filename2, path, name ):

    file1 = ROOT.TFile( filename1 )
    names = getObjectNamesRec( file1, path )

    diffPos = getFirstDifferentPosition( names )

    for name in names:
        #if "summary/PhotonValidator/Background/phoBkgDEta" not in name: continue
        h1 = readHisto( filename1, name )
        h2 = readHisto( filename2, name )
        draw( h1, h2, name[diffPos:].replace("/","_") )




if __name__ == "__main__":
    filename1 = "fastsim_dqm.root"
    filename2 = "fullsim_dqm.root"
    path = "DQMData/Run 1/EcalBarrel/Run summary/EBClusterTask"
    path = "DQMData/Run 1/EgammaV/Run summary/PhotonValidator"

    compareDQMHistograms( filename1, filename2, path, "test" )
