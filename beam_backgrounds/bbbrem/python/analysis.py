
import ROOT
import functions

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


def main():

    hist_ele_E = ROOT.TH1D("hist_ele_E", "Electron energy", 100, 0, 50)

    evts = functions.load_evts()
    for evt in evts:
        hist_ele_E.Fill(evt.ele.E())

    functions.plot_hist(hist_ele_E, "./", title="Electron energy", xMin=0, xMax=50, xLabel="Electron energy (GeV)", yLabel="Events")



if __name__ == "__main__":
    main()