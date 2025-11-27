
import sys,array,ROOT,math,os,copy
import argparse
import ctypes

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

parser = argparse.ArgumentParser()
parser.add_argument("--flavor", type=str, help="Flavor (mumu or ee)", default="mumu")
parser.add_argument("--type", type=str, help="Run type (mass or xsec)", choices=["mass", "xsec"], default="xsec")
args = parser.parse_args()



def makePlot(g_pulls, h_pulls, avg=0, avg_err=0):

    canvas = ROOT.TCanvas("c", "c", 800, 800)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.1)
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.05)
    canvas.SetFillStyle(4000) # transparency?
    canvas.SetGrid(1, 0)
    canvas.SetTickx(1)

    xTitle = "Selection efficiency Z(#mu^{#plus}#mu^{#minus})H  (%)"
    if flavor == "ee": xTitle = "Selection efficiency Z(e^{#plus}e^{#minus})H  (%)"


    h_pulls.GetXaxis().SetTitleSize(0.04)
    h_pulls.GetXaxis().SetLabelSize(0.035)
    h_pulls.GetXaxis().SetTitle(xTitle)
    h_pulls.GetXaxis().SetTitleOffset(1)
    h_pulls.GetYaxis().SetLabelSize(0.055)
    h_pulls.GetYaxis().SetTickLength(0)
    h_pulls.GetYaxis().LabelsOption('v')
    h_pulls.SetNdivisions(506, 'XYZ')
    h_pulls.Draw("HIST 0")
   

    maxx = 9
    line = ROOT.TLine(avg, 0, avg, maxx)
    line.SetLineColor(ROOT.kGray)
    line.SetLineWidth(2)
    #line.Draw("SAME")
    

    
    shade = ROOT.TGraph()
    shade.SetPoint(0, avg*0.999, 0)
    shade.SetPoint(1, avg*1.001, 0)
    shade.SetPoint(2, avg*1.001, maxx)
    shade.SetPoint(3, avg*0.999, maxx)
    shade.SetPoint(4, avg*0.999, 0)
    #shade.SetFillStyle(3013)
    shade.SetFillColor(16)
    shade.SetFillColorAlpha(16, 0.35);
    shade.Draw("SAME F")

    g_pulls.SetMarkerSize(1.2)
    g_pulls.SetMarkerStyle(20)
    g_pulls.SetLineWidth(2)
    g_pulls.Draw('P0 SAME')
    
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(30) # 0 special vertical aligment with subscripts
    latex.DrawLatex(0.95, 0.925, "#sqrt{s} = 240 GeV, 10.8 ab^{#minus1}")

    latex.SetTextAlign(13)
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.15, 0.96, "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}")
    
    latex.SetTextAlign(13)
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.SetTextColor(ROOT.kGray+1)
    latex.DrawLatex(0.2, 0.17, f"{avg:.2f} #pm {avg_err:.2f} %")
    
    canvas.SaveAs(f"{outDir}/decay_mode_independence.png")
    canvas.SaveAs(f"{outDir}/decay_mode_independence.pdf")
    


if __name__ == "__main__":

    flavor = args.flavor
    inputDir = "output/h_zh/histmaker_newcostheta/"
    outDir = "/home/submit/jaeyserm/public_html/fccee/h_zh/costhetamiss/"
    proc = f"wzp6_ee_{flavor}H_ecm240"


    # cuts with cos(theta)
    cuts = ["_cut0", "_cut1", "_cut2", "_cut3", "_cut4", "_cut5", "_cut6"]
    cut_labels = [r"No selection", r"$\geq 1 \mu$ + ISO", r"$\geq 2 \mu$ + OS", r"$86 < m_{\mu^{+}\mu^{-}} < 96$", r"$20 < p_{\mu^{+}\mu^{-}} < 70$", r"$120 < m_{rec} < 140$", r"$\cos(\theta_{miss}) < 0.98$"]

    # cuts without cos(theta)
    #cuts = ["_cut0", "_cut1", "_cut2", "_cut3", "_cut4", "_cut6"]
    #cut_labels = [r"No selection", r"$\geq 1 \mu$ + ISO", r"$\geq 2 \mu$ + OS", r"$86 < m_{\mu^{+}\mu^{-}} < 96$", r"$20 < p_{\mu^{+}\mu^{-}} < 70$", r"$120 < m_{rec} < 140$"]

    xMin, xMax = 70, 75
    if flavor == "ee":
        xMin, xMax = 60, 70
        cut_labels = [x.replace("\mu", "e") for x in cut_labels]


    decay_pdgids = [4, 5, 13, 15, 21, 22, 23, 24]
    decay_names = [r"cc", r"bb", r"$\mu\mu$", r"$\tau\tau$", r"gg", r"$\gamma\gamma$", r"ZZ", r"WW"]
    decay_names_py = ["cc", "bb", "mumu", "tautau", "gg", "gaga", "ZZ", "WW"]
    decay_names_tex = ["cc", "bb", "#mu#mu", "#tau#tau", "gg", "#gamma#gamma", "ZZ", "WW"]
    fIn = ROOT.TFile(f"{inputDir}/{proc}.root")


    h_pulls = ROOT.TH2F("pulls", "pulls", (xMax-xMin)*10, xMin, xMax, len(decay_pdgids)+1, 0, len(decay_pdgids)+1)
    g_pulls = ROOT.TGraphErrors(len(decay_pdgids)+1)

    f = open(f"{outDir}/decay_mode_independence.txt", 'w')
    sys.stdout = f

    print("BRANCHING RATIOS (%)")
    formatted_row = '{:<10} {:<25} {:<25} {:<25} {:<25} {:<25} {:<25} {:<25} {:<25} {:<25}'
    print(formatted_row.format(*(["Cut"]+decay_names_py+["Total"])))
    print(formatted_row.format(*(["----------"]+["-----------------------"]*9)))
    for i, cut in enumerate(cuts):
        row = [f"Cut {i}"]
        h = fIn.Get(f"higgs_decay{cut}")
        y = h.Integral()
        br_tot, br_tot_err = 0, 0
        for pdg in decay_pdgids:
            br, br_err = h.GetBinContent(pdg+1)*100./y, h.GetBinError(pdg+1)*100./y
            row.append("%.3e +/- %.3e" % (br, br_err))
            br_tot += br
            br_tot_err += br_err**2
        row.append("%.3e +/- %.3e" % (br_tot, br_tot_err**0.5))
        print(formatted_row.format(*row))


    print("\n\n\n")
    print("SELECTION EFFICIENCIES (%)") 
    formatted_row = '{:<10} {:<25} {:<25} {:<25} {:<25} {:<25} {:<25} {:<25} {:<25} {:<25}'
    print(formatted_row.format(*(["Cut"]+decay_names_py+["Average"])))
    print(formatted_row.format(*(["----------"]+["-----------------------"]*9)))
    ip = 0 # tgraph point counter
    for i,cut in enumerate(cuts):
        row = ["Cut %d"%i]

        idx = 0 if i == 0 else i-1
        idx = 0
        h_ref = fIn.Get(f"higgs_decay{cuts[idx]}")
        h_cut = fIn.Get(f"higgs_decay{cut}")

        y_ref_tot_err = ctypes.c_double(1.)
        y_ref_tot = h_ref.IntegralAndError(0, h_ref.GetNbinsX() + 1, y_ref_tot_err)
        y_ref_tot_err = y_ref_tot_err.value

        y_cut_tot_err = ctypes.c_double(1.)
        y_cut_tot = h_cut.IntegralAndError(0, h_cut.GetNbinsX() + 1, y_cut_tot_err)
        y_cut_tot_err = y_cut_tot_err.value

        sel_eff_tot = y_cut_tot / y_ref_tot
        sel_eff_tot_err = sel_eff_tot * ( (y_cut_tot_err/y_cut_tot)**2 + (y_ref_tot_err/y_ref_tot)**2)**0.5


        if i == len(cuts)-1:

            g_pulls.SetPoint(ip, sel_eff_tot*100., float(ip) + 0.5)
            #g_pulls.SetPointError(ip, sel_eff_tot_err*100., sel_eff_tot_err*100., 0., 0.)
            g_pulls.SetPointError(ip, sel_eff_tot_err*100., 0.)
            h_pulls.GetYaxis().SetBinLabel(ip + 1, "Average")
            ip += 1

        for k, pdg in enumerate(decay_pdgids):

            y_ref = h_ref.GetBinContent(pdg+1)
            y_cut = h_cut.GetBinContent(pdg+1)
            y_ref_err = h_ref.GetBinError(pdg+1)
            y_cut_err = h_cut.GetBinError(pdg+1)
            #y_ref_err = y_ref**0.5
            #y_cut_err = y_cut**0.5
            sel_eff = y_cut / y_ref
            sel_eff_err = sel_eff * ( (y_cut_err/y_cut)**2 + (y_ref_err/y_ref)**2)**0.5
            if sel_eff_err == 0.1: sel_eff_err = 0.0999
            row.append("%.3e +/- %.3e" % (sel_eff*100., sel_eff_err*100.))

            if i == len(cuts)-1:
                g_pulls.SetPoint(ip, sel_eff*100., float(ip) + 0.5)
                g_pulls.SetPointError(ip, sel_eff_err*100., 0.)
                h_pulls.GetYaxis().SetBinLabel(ip + 1, decay_names_tex[k])
                ip += 1

        row.append("%.3e +/- %.3e" % (sel_eff_tot*100., sel_eff_tot_err*100.)) # average
        print(formatted_row.format(*row))

    makePlot(g_pulls, h_pulls, avg=sel_eff_tot*100., avg_err=sel_eff_tot_err*100.)

    fIn.Close()