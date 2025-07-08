
import sys,argparse
import os
import ROOT
import math
from collections import OrderedDict

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

# function to read the geometry defined in the TrackCovariance module
def read_delphes_card(delphes_card):
    geometry = []
    detectors = []
    with open(delphes_card) as f:
        isGeom = False
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if not line.isspace() and line.replace(' ', '')[0] == "#":
                continue
            
            if "set" in line and "DetectorGeometry" in line:
                isGeom = True
                continue
            if not isGeom:
                continue
            if "}" in line:
                isGeom = False
                continue
            tmp = line.split()
            cfg = [int(tmp[0]), str(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5]), float(tmp[6])]
            if not str(tmp[1]) in detectors:
                detectors.append(str(tmp[1]))
            geometry.append(cfg)

    return geometry, detectors

def main(args, detector_groups, colors):

    geometry, detectors = read_delphes_card(args.input)
    #print(detectors)

    # make an histogram for each material
    hists = {}
    for m in detectors:
        h = ROOT.TH1D(m, "", angleBinsN, angleMin, angleMax)
        hists[m] = h

    # loop over the geometry and calculate material budget for each material
    for xBin in range(1, angleBinsN+1):
        if angleDef == "theta":
            theta = hists[detectors[0]].GetBinCenter(xBin)
        else:
            theta = math.degrees(math.acos(hists[detectors[0]].GetBinCenter(xBin)))
        sint = math.sin(math.radians(theta))
        cost = math.cos(math.radians(theta))
        tant = math.tan(math.radians(theta))
        for ig, g in enumerate(geometry):
            if g[1] == 'MAG':
                continue
            x0 = 0

            # barrel
            if g[0] == 1:
                r = g[4]
                z_max = g[3]
                z = r / tant
                if z < z_max:
                    x0 = 100.*g[5]/g[6]/sint

            # endcap
            if g[0] == 2:
                r_min = g[2]
                r_max = g[3]
                z = g[4]
                r = z * tant
                if z < 0: continue # only one quadrant, assume symmetric
                if r < r_max and r > r_min:
                    x0 = 100.*g[5]/g[6]/cost

            hists[g[1]].SetBinContent(xBin, hists[g[1]].GetBinContent(xBin) + x0)


    # plotting
    c = ROOT.TCanvas("", "", 800, 800)
    c.SetTopMargin(0.055)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.11)

    dummy = ROOT.TH1D("dummy", "", angleBinsN, angleMin, angleMax)
    dummy.GetXaxis().SetTitle("#theta (deg)" if angleDef == "theta" else "cos(#theta)")
    dummy.GetYaxis().SetTitle("Material budget x/X_{0} (%)")

    dummy.GetXaxis().SetTitleFont(43)
    dummy.GetXaxis().SetTitleSize(32)
    dummy.GetXaxis().SetLabelFont(43)
    dummy.GetXaxis().SetLabelSize(28)

    dummy.GetYaxis().SetTitleFont(43)
    dummy.GetYaxis().SetTitleSize(32)
    dummy.GetYaxis().SetLabelFont(43)
    dummy.GetYaxis().SetLabelSize(28)

    if angleDef == "theta":
        leg = ROOT.TLegend(0.5, 0.9-0.04*(len(detector_groups)+1), 0.9, 0.9)
    else:
        leg = ROOT.TLegend(0.2, 0.9-0.04*(len(detector_groups)+1), 0.55, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)
    leg.SetHeader(args.title)

    st = ROOT.THStack()
    st.SetName("stack")
    hists_to_write = []
    for i,m in enumerate(detector_groups.keys()):
        hists_merge = [hists[x] for x in detector_groups[m]]
        hist = hists_merge[0]
        for h in hists_merge[1:]: hist.Add(h)

        hist.SetName(m)
        hist.SetFillColor(colors[i])
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineWidth(1)
        hist.SetLineStyle(1)
        st.Add(hist)
        leg.AddEntry(hist, m, "F")
        hists_to_write.append(hist)

    dummy.GetYaxis().SetRangeUser(0, args.ymax if args.ymax > 0 else st.GetStack().Last().GetMaximum()*1.05)
    dummy.Draw("HIST")
    st.Draw("SAME HIST")
    leg.Draw()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()
    c.SaveAs(f"{args.outDir}/delphes_x0.png")
    c.SaveAs(f"{args.outDir}/delphes_x0.pdf")

    fOut = ROOT.TFile(f"{args.outDir}/delphes_x0.root", "RECREATE")
    for h in hists_to_write:
        h.Write()
    fOut.Close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="delphes_card_IDEA.tcl", type=str, help="Input Delphes card (.tcl)")
    parser.add_argument("-o", "--outDir", default="/home/submit/jaeyserm/public_html/fccee/material_budget/", type=str, help="Output filename (will produce .png, .pdf, .root)")
    parser.add_argument("-t", "--title", type=str, default="IDEA CLD silicon tracker (Delphes)", help="Name of the detector")
    parser.add_argument("-y", "--ymax", type=int, default=-1, help="Maximum y-axis")
    args = parser.parse_args()

    angleMin = 0
    angleMax = 0.99
    angleDef = "cosTheta"
    angleBinning = 0.01
    angleBinsN = (int)((angleMax - angleMin) / angleBinning) # number of bins


    # define detector groups
    if "CLD" in args.input or "SiTracking" in args.input:
        detector_groups = OrderedDict()
        detector_groups['Beam pipe'] = ['PIPE']
        #detector_groups['Vertex detector'] = ['VTXLOW', 'VTXDSK'] # new CLD
        detector_groups['Vertex detector'] = ['VTX'] # old CLD
        detector_groups['Inner tracker'] = ['ITK', 'ITKDSK']
        detector_groups['Outer tracker'] = ['OTK', 'OTKDSK']


    if "IDEA" in args.input:
        detector_groups = OrderedDict()
        detector_groups['Beampipe'] = ['PIPE']
        detector_groups['Vertex'] = ['VTXLOW', 'VTXHIGH', 'VTXDSK']
        #detector_groups['Drift chamber'] = ['DCHCANI', 'DCH', 'DCHWALL', 'DCHCANO']
        #detector_groups['Silicon wrapper'] = ['FSILWRP', 'BSILWRP']
        #detector_groups['Preshower'] = ['BPRESH', 'FPRESH']
        #detector_groups['FRAD'] = ['FRAD']

    colors = [ROOT.kRed, ROOT.kOrange, ROOT.kBlue, ROOT.kGreen, ROOT.kGray+1, ROOT.kMagenta]
    main(args, detector_groups, colors)