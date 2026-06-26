
import sys,array,ROOT,math,os,copy

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

sys.path.insert(0, f'{os.path.dirname(os.path.realpath(__file__))}/../../python')
import plotter

def recoil_plot(outName, files, hists, labels, colors, xMin, xMax, yMin, yMax, xTitle, yTitle, rebin=1, logy=False, norm=False, legLabel="", topRight="", legposx=.2):

    n = len(files) + (0 if legLabel=="" else 1)
    leg = ROOT.TLegend(legposx, 0.9-0.05*n, 0.95, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(1)
    leg.SetTextSize(0.035)
    leg.SetMargin(0.15)
    if legLabel != "":
        leg.SetHeader(legLabel)

    cfg = {

        'logy'              : logy,
        'logx'              : False,

        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,

        'xtitle'            : xTitle,
        'ytitle'            : yTitle,

        'topRight'          : topRight, 
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}",
    }

    plotter.cfg = cfg
    canvas = plotter.canvas()
    dummy = plotter.dummy()
    dummy.Draw("HIST")

    for i,h in enumerate(hists):
        lumi = lumi_240 if "ecm240" in hists[i] else lumi_365
        fIn = ROOT.TFile(files[i])
        hist = copy.deepcopy(fIn.Get(hists[i]))
        hist.Scale(lumi)
        if norm:
            x1, x2 = hist.FindBin(xMin), hist.FindBin(xMax)
            norm_ = hist.Integral(x1, x2)
            #norm_ = hist.Integral()
            hist.Scale(1./norm_)

        hist.Rebin(rebin)
        hist.SetLineColor(colors[i])
        hist.SetLineWidth(2)
        leg.AddEntry(hist, labels[i], "L")
        hist.Draw("SAME HIST")
        #print(hist.Integral())


    leg.Draw("SAME") 
    plotter.aux()
    canvas.SetGrid()  
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs(f"{outDir}/{outName}.png")
    canvas.SaveAs(f"{outDir}/{outName}.pdf")
    canvas.Close()



if __name__ == "__main__":

    tag = "test"
    lumi_240 = 10.8
    lumi_365 = 3


    topRight = "#sqrt{{s}} = 240 GeV, {} ab^{{#minus1}}".format(lumi_240)
    outDir = "/home/submit/jaeyserm/public_html/fccee/h_mass/plots_compare/"
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+1]
    colors = [
        ROOT.kRed+1,
        ROOT.kBlue-4,
        ROOT.kGreen+3,
        #ROOT.kOrange-3,
        ROOT.kBlack
    ]

    files = ["output/h_mass/histmaker/ecm240/wzp6_ee_mumuH_ecm240.root", "output/h_mass/histmaker/ecm240/wzp6_ee_eeH_ecm240.root", "output/h_mass/histmaker/ecm365/wzp6_ee_mumuH_ecm365.root"]
    hists = ["mumu_zll_recoil_m", "ee_zll_recoil_m", "mumu_zll_recoil_m"]
    labels = ["Z(#mu^{#plus}#mu^{#minus})H, 240 GeV", "Z(e^{#plus}e^{#minus})H, 240 GeV", "Z(#mu^{#plus}#mu^{#minus})H, 365 GeV"]
    recoil_plot("240_365_mumu_ee", files, hists, labels, colors, 122, 132, 0, 0.07, "Recoil mass (GeV)", "Events (normalized)", rebin=10, legLabel="", norm=True, topRight="")

    files = ["output/h_mass/histmaker/ecm240/wzp6_ee_mumuH_noBES_ecm240.root", "output/h_mass/histmaker/ecm240/wzp6_ee_mumuH_ecm240.root", "output/h_mass/histmaker/ecm365/wzp6_ee_mumuH_ecm365.root"]
    hists = ["mumu_zll_recoil_m", "mumu_zll_recoil_m", "mumu_zll_recoil_m"]
    labels = ["240 GeV, no beam energy spread", "240 GeV", "365 GeV"]
    recoil_plot("240_365_mumu_bes_nobes", files, hists, labels, colors, 122, 132, 0, 0.15, "Recoil mass (GeV)", "Events (normalized)", rebin=10, legLabel="Z(#mu^{#plus}#mu^{#minus})H", norm=True, topRight="")

    files = ["output/h_mass/histmaker/ecm240/wzp6_ee_mumuH_noBES_ecm240.root", "output/h_mass/histmaker/ecm240/wzp6_ee_mumuH_ecm240.root", "output/h_mass/histmaker/ecm240/wzp6_ee_eeH_ecm240.root", "output/h_mass/histmaker/ecm365/wzp6_ee_mumuH_ecm365.root"]
    hists = ["mumu_zll_recoil_m", "mumu_zll_recoil_m", "ee_zll_recoil_m", "mumu_zll_recoil_m"]
    labels = ["Z(#mu^{#plus}#mu^{#minus})H, 240 GeV, no beam energy spread", "Z(#mu^{#plus}#mu^{#minus})H, 240 GeV", "Z(e^{#plus}e^{#minus})H, 240 GeV", "Z(#mu^{#plus}#mu^{#minus})H, 365 GeV"]
    recoil_plot("240_365_mumu_ee_bes_nobes", files, hists, labels, colors, 122, 132, 0, 0.15, "Recoil mass (GeV)", "Events (normalized)", rebin=10, legLabel="", norm=True, topRight="")


    files = ["output/h_mass/histmaker/ecm240/wzp6_ee_mumuH_ecm240.root", "output/h_mass/histmaker/ecm240/gen/wzp6_ee_mumuH_ecm240.root", "output/h_mass/histmaker/ecm240/wzp6_ee_mumuH_ecm240_3T.root", "output/h_mass/histmaker/ecm240/wzp6_ee_mumuH_ecm240_SiTracking.root"]
    hists = ["mumu_zll_recoil_m", "mumu_zll_recoil_m", "mumu_zll_recoil_m", "mumu_zll_recoil_m"]
    labels = ["Default", "Perfect resolution", "Magnetic field 3T", "Silicon tracker"]
    recoil_plot("240_mumu_variations", files, hists, labels, colors, 122, 132, 0, 0.08, "Recoil mass (GeV)", "Events (normalized)", rebin=10, legLabel="Z(#mu^{#plus}#mu^{#minus})H, 240 GeV", norm=True, topRight="", legposx=0.55)

    #files = ["output/h_mh/histmaker_test/wz3p6_ee_mumuH_ecm240.root", "output/h_mh/histmaker_test/wz3p6_ee_mumuH_ecm240_WithWires.root"]
    #hists = ["mumu_zll_recoil", "mumu_zll_recoil"]
    #labels = ["IDEA FSR", "IDEA wire material"]
    #recoil_plot("IDEA_comparison", files, hists, labels, colors, 122, 132, 0, 0.008, "Recoil (GeV)", "Events", rebin=10, legLabel="Muon final state Z(#mu^{#plus}#mu^{#minus})H", norm=True, topRight=topRight)

