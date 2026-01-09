
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


class RadiativeBhabha:

    def __init__(self, q2_1, q2_2, q2_3, q2_4, p2_1, p2_2, p2_3, p2_4, qk_1, qk_2, qk_3, qk_4, weight, sqrtmt, t, d1, d2, s1):

        self.ele = ROOT.Math.PxPyPzEVector(float(p2_1), float(p2_2), float(p2_3), float(p2_4))
        self.pos = ROOT.Math.PxPyPzEVector(float(q2_1), float(q2_2), float(q2_3), float(q2_4))
        self.gam = ROOT.Math.PxPyPzEVector(float(qk_1), float(qk_2), float(qk_3), float(qk_4))
        self.weight = float(weight)
        self.sqrtmt = float(sqrtmt)
        self.t = float(t)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.s1 = float(s1)


def plot_hist(h, outdir, title="", xMin=-1, xMax=-1, yMin=-1, yMax=-1, xLabel="", yLabel="Events", logY=False):

    c = ROOT.TCanvas("c", "c", 800, 600)
    if logY:
        c.SetLogy(True)
    
    h.Draw("HIST")
    if title != "":
        h.SetTitle(title)
    if xLabel != "":
        h.GetXaxis().SetTitle(xLabel)
    if yLabel != "":
        h.GetYaxis().SetTitle(yLabel)

    if xMin != -1 and xMax != -1:
        h.GetXaxis().SetRangeUser(xMin, xMax)
    if yMin != -1 and yMax != -1:
        h.SetMinimum(yMin)
        h.SetMaximum(yMax)

    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.Update()

    c.SaveAs(f"{outdir}/{h.GetName()}.png")
    c.SaveAs(f"{outdir}/{h.GetName()}.pdf")



def load_evts(fInName='bbrems_ntuple.tfs'):
    evts = []

    with open('bbrems_ntuple.tfs', 'r') as file:
        for line in file:
            line = line.strip()
            if line[0] == '@' or line[0] == '*':
                continue
            q2_1, q2_2, q2_3, q2_4, p2_1, p2_2, p2_3, p2_4, qk_1, qk_2, qk_3, qk_4, weight, sqrtmt, t, d1, d2, s1 = line.split()
            rb = RadiativeBhabha(q2_1, q2_2, q2_3, q2_4, p2_1, p2_2, p2_3, p2_4, qk_1, qk_2, qk_3, qk_4, weight, sqrtmt, t, d1, d2, s1)
            evts.append(rb)
    return evts