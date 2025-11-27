

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

if False:



    cut = 0.99975
    hName, fOut = "cosThetaMiss_cut5", "costhetamiss_zgamma" # cosTheta_miss_new cosThetaMiss_
    xMin, xMax, xTitle = 0.9, 1, "cos(theta miss)"
    rebin = 1

    #cut = 0
    #hName, fOut = "missingEnergy_", "missingenergy"
    #xMin, xMax, xTitle = 0, 100, "Missing energy (GeV)"
    #rebin = 1


    file = ROOT.TFile.Open("output/h_zh/histmaker/wzp6_ee_mumu_ecm240.root")
    h1 = file.Get(hName)
    h1.Scale(1./h1.Integral())

    cut_bin = h1.GetXaxis().FindBin(cut)
    left = h1.Integral(0, cut_bin)
    right = h1.Integral(cut_bin, 1)

    print(left, right)

    h1.SetLineColor(ROOT.kBlack)
    h1.SetLineWidth(2)
    h1.Rebin(rebin)

    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
    #canvas.SetLogy()

    h1.Draw("HIST")
    h1.GetXaxis().SetRangeUser(xMin, xMax)
    h1.GetXaxis().SetTitle(xTitle)
    h1.GetYaxis().SetTitle("Entries")


    canvas.SaveAs(f"/home/submit/jaeyserm/public_html/fccee/h_zh/costhetamiss/{fOut}.png")

    file.Close()

if True:

    decay_pdgids = [4, 5,  15, 21, 22, 23, 24]
    #decay_pdgids = [21, 22, 23, 24]
    decay_names = [r"cc", r"bb", r"#tau#tau", r"gg", r"#gamma#gamma", r"ZZ", r"WW"]
    decay_names_py = ["cc", "bb", "tautau", "gg", "gaga", "ZZ", "WW"]
    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen, ROOT.kOrange, ROOT.kCyan, ROOT.kGray, ROOT.kGreen+1, ROOT.kGreen+2]

    hists = []



    cut = 0.98
    hName, fOut = "cosTheta_miss_new", "costhetamiss" # cosTheta_miss_new cosThetaMiss_
    xMin, xMax, xTitle = 0.9, 1, "cos(theta miss)"
    rebin = 10

    cut = 0.99975
    hName, fOut = "cosTheta_miss_new", "costhetamiss_ILCCUT" # cosTheta_miss_new cosThetaMiss_
    xMin, xMax, xTitle = 0.85, 1, "cos(theta miss)"
    rebin = 1


    cut = 10
    hName, fOut = "visibleEnergy_new", "visibleEnergy" # cosTheta_miss_new cosThetaMiss_
    xMin, xMax, xTitle = 0, 150, "Visible energy"
    rebin = 10


    #cut = 0
    #hName, fOut = "missingEnergy_", "missingenergy"
    #xMin, xMax, xTitle = 0, 100, "Missing energy (GeV)"
    #rebin = 1


    file = ROOT.TFile.Open("output/h_zh/histmaker_newcostheta_ilcsel/wzp6_ee_mumuH_ecm240.root")
    h2 = file.Get(hName)  # Replace with the actual name of your 2D histogram

    for i,pdgid in enumerate(decay_pdgids):
        h1 = h2.ProjectionX(f"y{pdgid}", pdgid+1, pdgid+1)
        
        h1.SetName(decay_names_py[i])
        h1.Scale(1./h1.Integral())

        cut_bin = h1.GetXaxis().FindBin(cut)
        left = h1.Integral(0, cut_bin)
        right = h1.Integral(cut_bin, 1)

        h1.SetTitle(f"{decay_names[i]} ({left:.3f}/{right:.3f})")
        h1.SetLineColor(colors[i])
        h1.SetLineWidth(2)
        h1.Rebin(rebin)
        hists.append(h1)



    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
    #canvas.SetLogy()

    legend = ROOT.TLegend(0.2, 0.5, 0.5, 0.85)
    legend.SetHeader(f"Cut {cut}: yields below/yields above")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

        

    for i, hist in enumerate(hists):

        if i == 0:
            hist.Draw("HIST")
            hist.GetXaxis().SetRangeUser(xMin, xMax)
            hist.GetXaxis().SetTitle(xTitle)
            hist.GetYaxis().SetTitle("Entries")

        else:
            hist.Draw("SAME HIST")
        legend.AddEntry(hist, hist.GetTitle(), "L")

    legend.Draw()
    canvas.SaveAs(f"/home/submit/jaeyserm/public_html/fccee/h_zh/costhetamiss/{fOut}.png")

    file.Close()