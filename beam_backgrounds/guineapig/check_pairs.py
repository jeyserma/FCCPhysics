import os
import glob

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptTitle(0)

# Directory with .pair files
directory = "/ceph/submit/data/group/fcc/ee/detector/guineapig//FCCee_Z_4IP_04may23_FCCee_Z128_noTracking/"
directory = "/ceph/submit/data/group/fcc/ee/detector/guineapig//FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_halfWindow/"
directory = "/ceph/submit/data/group/fcc/ee/detector/guineapig//FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_reducedWindow/"


directory = "/ceph/submit/data/group/fcc/ee/detector/guineapig//FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_halfZWindow/"
directory = "/ceph/submit/data/group/fcc/ee/detector/guineapig//FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_halfYWindow/"
directory = "/ceph/submit/data/group/fcc/ee/detector/guineapig//FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_doubleZWindow/"

directory = "/ceph/submit/data/group/fcc/ee/detector/guineapig//FCCee_Z_4IP_04may23_noXing_FCCee_Z128_noTracking_reducedWindow/"
directory = "/ceph/submit/data/group/fcc/ee/detector/guineapig//FCCee_Z_4IP_04may23_noXing_FCCee_Z128_noTracking/"


directory = "/ceph/submit/data/group/fcc/ee/detector/guineapig//FCCee_Z_4IP_04may23_FCCee_Z256_noTracking_reducedWindow/"
directory = "/ceph/submit/data/group/fcc/ee/detector/guineapig//FCCee_Z_4IP_04may23_FCCee_Z256_noTracking/"


directory = "/ceph/submit/data/group/fcc/ee/detector/guineapig_devTrackingWindow//FCCee_Z_4IP_04may23_FCCee_Z128/"

name = os.path.basename(os.path.normpath(directory))
outDir = f"/home/submit/jaeyserm/public_html/fccee/pairs_check/{name}/"

os.system(f"mkdir -p {outDir}")
os.system(f"cp {outDir}../index.php {outDir}")


# base window
doProdKin = True # true/false using production kinematics (pairs0.dat)
wx = 1.32555
wy = 1000.*0.00186
wz = 25.4
n = 500
extraWindow = 1.1

doProdKinSuffix = '0' if doProdKin else ''

# Create ROOT histograms
hist_E = ROOT.TH1F("hist_e", "Energy distribution; Energy (MeV); Entries", 1000, 0, 10)
hist_x = ROOT.TH1F("hist_x", "x distribution; x (mm); Entries", n, -wx*extraWindow, wx*extraWindow)
hist_y = ROOT.TH1F("hist_y", "y distribution; y (um); Entries", n, -wy*extraWindow, wy*extraWindow)
hist_z = ROOT.TH1F("hist_z", "z distribution; z (mm); Entries", n, -wz*extraWindow, wz*extraWindow)


for filepath in glob.glob(f"{directory}/output{doProdKinSuffix}_*.pairs"):
    print(filepath)
    with open(filepath, "r") as file:
        for line in file:
            if line.startswith("E") or line.strip() == "":
                continue  # Skip header or empty lines
            try:
                E, px, py, pz, x, y, z, v, v,v = map(float, line.split())
                #print(x/1000/1000)
                hist_x.Fill(x/1000/1000) # nm to mm
                hist_y.Fill(y/1000) # nm to um
                hist_z.Fill(z/1000/1000) # nm to mm
                hist_E.Fill(E*1000) # to MeV
            except ValueError:
                print(f"Skipping line in {filepath}: {line.strip()}")

hist_x.Scale(1./hist_x.Integral())
hist_y.Scale(1./hist_y.Integral())
hist_z.Scale(1./hist_z.Integral())
hist_E.Scale(1./hist_E.Integral())

# Draw histograms
c1 = ROOT.TCanvas("c1", "Histograms", 1000, 1000)
c1.SetTopMargin(0.055)
c1.SetRightMargin(0.05)
c1.SetLeftMargin(0.12)
c1.SetBottomMargin(0.11)

hist_x.Draw("HIST")
c1.SaveAs(f"{outDir}/hist_x{doProdKinSuffix}.png")

c1.Clear()
hist_y.Draw("HIST")
c1.SaveAs(f"{outDir}/hist_y{doProdKinSuffix}.png")

c1.Clear()
hist_z.Draw("HIST")
c1.SaveAs(f"{outDir}/hist_z{doProdKinSuffix}.png")

c1.Clear()
c1.SetLogy()
hist_E.Draw("HIST")
c1.SaveAs(f"{outDir}/hist_E{doProdKinSuffix}.png")


print(hist_x.Integral())
print(hist_y.Integral())
print(hist_z.Integral())


fOut = ROOT.TFile(f"{outDir}/output{doProdKinSuffix}.root", "RECREATE")
hist_x.Write()
hist_y.Write()
hist_z.Write()
hist_E.Write()
fOut.Close()

