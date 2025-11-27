

import ROOT
import copy

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


def plot(hName, xMin, xMax, xTitle, rebin=1, cuts=[0, 500]):

    modes = ['wzp6_ee_qqH_ecm240', 'wzp6_ee_qqH_Hss_ecm240', 'wzp6_ee_qqH_Hcc_ecm240', 'wzp6_ee_qqH_Hbb_ecm240', 'wzp6_ee_qqH_Hgg_ecm240', 'wzp6_ee_qqH_HZZ_ecm240', 'wzp6_ee_qqH_HWW_ecm240','wzp6_ee_qqH_Haa_ecm240', 'wzp6_ee_qqH_HZa_ecm240', 'wzp6_ee_qqH_Hmumu_ecm240']
    labels = ['inclusive', 'ss', 'cc', 'bb', 'gg', 'ZZ', 'WW', '#gamma#gamma', 'Z#gamma', '#mu#mu']

    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen, ROOT.kOrange, ROOT.kCyan, ROOT.kGray, ROOT.kGreen+1, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kRed+2]

    hists = []
    for i,mode in enumerate(modes):
        file = ROOT.TFile.Open(f"output/h_zh_hadronic/histmaker/{mode}.root")
        h1 = copy.deepcopy(file.Get(hName))
        h1.Scale(1./h1.Integral())

        cut_bin_left = h1.GetXaxis().FindBin(cuts[0])
        cut_bin_right = h1.GetXaxis().FindBin(cuts[1])
        cut_yield = h1.Integral(cut_bin_left, cut_bin_right)

        #h1.SetTitle(f"{labels[i]}")
        h1.SetTitle(f"{labels[i]} [{cuts[0]:.0f}, {cuts[1]:.0f}]={cut_yield:.3f}")
        h1.SetLineColor(colors[i])
        h1.SetLineWidth(2)
        h1.Rebin(rebin)
        hists.append(h1)
        




    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
    #canvas.SetLogy()

    legend = ROOT.TLegend(0.15, 0.5, 0.5, 0.85)
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
    canvas.SaveAs(f"/home/submit/jaeyserm/public_html/fccee/h_zh_hadronic/plots/{hName}.png")

    file.Close()

if __name__ == '__main__':


    #plot("njets", 0, 10, "njets")
    #plot("zqq_m", 0, 150, "zqq_m")
    #plot("zqq_recoil_m", 0, 200, "zqq_recoil_m")
    
    plot("zqq_m_best", 0, 150, "zqq_m", cuts=[65, 115])
    plot("zqq_recoil_m_best", 0, 200, "zqq_recoil_m", cuts=[110, 140])