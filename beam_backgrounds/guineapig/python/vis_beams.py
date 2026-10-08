import os, sys
import numpy as np
import math
from scipy.constants import c, micro, nano, pi, milli, micro

from PIL import Image
import glob

import functions

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.TH1.SetDefaultSumw2(False)




# take from acc.dat
sigmaz = 15200.00*micro
sigmax = 8837.42*nano
sigmay = 38.34*nano
npart = 21.80e10
n0 = npart / (sigmax * sigmay * sigmaz * (2.*pi)**(3./2.))





ROOT.EnableImplicitMT(4) # use all cores
#ROOT.DisableImplicitMT() # single core

# load libraries
ROOT.gInterpreter.Declare('#include "python/vis_beams.h"')




# time
n_iterations = 256
iterations = range(n_iterations)
# box
nx = 128*4
ny = 128*4
nz = 128*4
Lx = 150*sigmax*1
Ly = 100*sigmay*1
Lz = 4*sigmaz*1
gridx = np.linspace(-0.5*Lx, 0.5*Lx, nx+1)
gridy = np.linspace(-0.5*Ly, 0.5*Ly, ny+1)
gridz = np.linspace(-0.5*Lz, 0.5*Lz, nz+1)
dx = gridx[1]-gridx[0]
dy = gridy[1]-gridy[0]
dz = gridz[1]-gridz[0]
nmacropart = 1e5
w0 = npart / nmacropart # weight

def rm_h(h):
    nx = h.GetNbinsX()
    ny = h.GetNbinsY()

    minZ, maxZ = 1e12, -1e12
    for ix in range(1, nx+1):
        for iy in range(1, ny+1):
            v = h.GetBinContent(ix, iy)
            if v == 0:
                h.SetBinContent(ix, iy, -1e12)
            if v > maxZ:
                maxZ = v
            if v < minZ:
                minZ = v
    return minZ, maxZ

def one_step_gp(n):
    """
    inputs
        n: current timestep
        gp_dir: simulation folder
    outputs
        H_zx, H_zy: density of the beams integrated along y and x, resp
    """
    global w0, dx, dy, dz, gridx, gridy, gridz

    # get beams' data, columns are:
    # Particle Energy [GeV] | x [um] | y [um] | z [um] | x' [urad] | y' [urad]
    data1 = np.loadtxt(os.path.join(INPUT_DIR, 'b1.%d' % n))
    data2 = np.loadtxt(os.path.join(INPUT_DIR, 'b2.%d' % n))

    # get the number of macroparticles in the beams
    N1 = np.shape(data1)[0]
    N2 = np.shape(data2)[0]

    # stack the data together and convert to SI
    E_data = np.hstack((data1[:,0], data2[:,0]))
    x_data = np.hstack((data1[:,1], data2[:,1]))*micro
    y_data = np.hstack((data1[:,2], data2[:,2]))*micro
    z_data = np.hstack((data1[:,3], data2[:,3]))*micro
    ax_data = np.hstack((data1[:,4], data2[:,4]))
    ay_data = np.hstack((data1[:,5], data2[:,5]))

    print(E_data, ax_data)

    weights = np.ones(N1+N2) * w0
    weights[-N2:] *= -1.0 # assign negative weights to the second beam

    H_zx = ROOT.TH2D("H_zx", "", nz+1, -0.5*Lz, 0.5*Lz, nx+1, -0.5*Lx, 0.5*Lx)
    for i in range(x_data.size):
        H_zx.Fill(z_data[i], x_data[i], weights[i])
    #H_zx.Scale(1./(dz*dx))
    minZ, maxZ = rm_h(H_zx)
    H_zx.SetMinimum(minZ)
    H_zx.SetMaximum(maxZ)
    return H_zx




def plot_one_step(step_number=1):
    H_zx = one_step_gp(step_number)


    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kBird)

    c = ROOT.TCanvas("c", "", 800, 800)
    c.SetBottomMargin(0.10)
    c.SetLeftMargin(0.15)
    c.SetGrid()
    H_zx.SetTitle(f"STEP {step_number:03d}")
    H_zx.GetXaxis().SetTitle("z (mm)")
    H_zx.GetYaxis().SetTitle("x (mm)")
    H_zx.Draw("COL")
    c.SaveAs(f"{OUTPUT_DIR}/raw/step_{step_number:03d}.png")


def make_gif():
    files = sorted(glob.glob(f"{OUTPUT_DIR}/raw/*.png"))
    images = [Image.open(f) for f in files]

    images[0].save(
        f"{OUTPUT_DIR}/crossing.gif",
        save_all=True,
        append_images=images[1:],
        duration=0.01, # milliseconds per frame
        loop=0
    )


#for step in iterations:
#    plot_one_step(step_number=step)

#make_gif()






def main():

    ACC = "FCCee_Z_GHC_V25p1"
    PAR = "CFG_DEF_64_4Z_XYADJ"

    label = "GHC V25.1 (FSR) 91.2 GeV"

    tag = f"{ACC}_{PAR}"
    input_dir = f"/ceph/submit/data/group/fcc/ee/beam_backgrounds/guineapig/visualization/{ACC}_{PAR}/"
    output_dir = f"/home/submit/jaeyserm/public_html/fccee/guineapig/visualization_edm/{ACC}_{PAR}"
    index_php = f"/home/submit/jaeyserm/public_html/fccee/guineapig/visualization_edm/index.php"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        os.system(f"cp {index_php} {output_dir}")

    fInName = f"{input_dir}/output.root"
    reader = functions.GuineaPigReader(fInName)
    
    cut_x = reader.get_metdata('cut_x') # nm
    cut_y = reader.get_metdata('cut_y') # nm
    cut_z = reader.get_metdata('cut_z') # um
    min_x_ = reader.get_metdata('min_x_') # nm
    max_x_ = reader.get_metdata('max_x_') # nm
    min_y_ = reader.get_metdata('min_y_') # nm
    max_y_ = reader.get_metdata('max_y_') # nm
    min_z_ = reader.get_metdata('min_z_') # nm
    max_z_ = reader.get_metdata('max_z_') # nm
    n_x = reader.get_metdata('n_x')
    n_y = reader.get_metdata('n_y')
    n_z = reader.get_metdata('n_z')
    npart = reader.get_metdata('beam1_pars_particles')
    sigmax = reader.get_metdata('beam1_pars_sigma_x') # nm
    sigmay = reader.get_metdata('beam1_pars_sigma_y') # nm
    sigmaz = reader.get_metdata('beam1_pars_sigma_z') # um
    betay = reader.get_metdata('beam1_pars_beta_y') # um
    n0 = npart / (sigmax * sigmay * sigmaz * (2.*pi)**(3./2.))
    nmacropart = 1e5
    w0 = npart / nmacropart # weight




    min_x_calc_, max_x_calc_ = -(n_x-2)/(n_x)*cut_x, (n_x-2)/(n_x)*cut_x
    min_y_calc_, max_y_calc_ = -(n_y-2)/(n_y)*cut_y, (n_y-2)/(n_y)*cut_y
    min_z_calc_, max_z_calc_ = -cut_z*1000, cut_z*1000

    print("****Inner grid dimensions")
    print(f"n_x={n_x}")
    print(f"n_y={n_y}")
    print(f"n_z={n_z}")

    print(f"sigmax={sigmax} nm")
    print(f"sigmay={sigmay} nm")
    print(f"sigmaz={sigmaz} um")
    print(f"betay={betay} um")

    print(f"cut_x={cut_x} nm")
    print(f"cut_y={cut_y} nm")
    print(f"cut_z={cut_z} um")

    print(f"min/max x={min_x_}/{max_x_} nm")
    print(f"min/max y={min_y_}/{max_y_} nm")
    print(f"min/max z={min_z_}/{max_z_} nm")

    print(f"min/max calc x={min_x_calc_}/{max_x_calc_} nm")
    print(f"min/max calc y={min_y_calc_}/{max_y_calc_} nm")
    print(f"min/max calc z={min_z_calc_}/{max_z_calc_} nm")

    infl = 1.2
    bins_x = (int(n_x), infl*min_x_/1e3, infl*max_x_/1e3) # um
    bins_y = (int(n_y), infl*min_y_/1e0, infl*max_y_/1e0) # nm
    bins_z = (int(n_z), infl*min_z_/1e6, infl*max_z_/1e6) # mm
    bins_t = (int(n_z), 0, n_z*2)


    print(bins_x)
    print(bins_y)
    print(bins_z)
    print(bins_t)
    #quit()

    ## FILLING: https://gitlab.cern.ch/jaeyserm/guinea-pig/-/blob/master/src/EDMwriterCPP.cc#L147-148

    if True:

        print("**** Create RDF")
        df = ROOT.RDataFrame("events", fInName)
        df = df.Define("beam1_slice", "getSlice(Beam1Slice)")
        df = df.Define("beam1_x", "getPos(Beam1Slice, 0)") # um
        df = df.Define("beam1_y", "getPos(Beam1Slice, 1)*1e3") # nm
        df = df.Define("beam1_z", "getPos(Beam1Slice, 2)*1e-3") # mm
        df = df.Define("beam1_px", "getMom(Beam1Slice, 0)/1e3") # MeV
        df = df.Define("beam1_py", "getMom(Beam1Slice, 1)/1e3") # MeV
        df = df.Define("beam1_pz", "getMom(Beam1Slice, 2)/1e6") # GeV
      
        df = df.Define("beam2_slice", "getSlice(Beam2Slice)")
        df = df.Define("beam2_x", "getPos(Beam2Slice, 0)") # um
        df = df.Define("beam2_y", "getPos(Beam2Slice, 1)*1e3") # nm
        df = df.Define("beam2_z", "getPos(Beam2Slice, 2)*1e-3") # mm
        df = df.Define("beam2_px", "getMom(Beam2Slice, 0)/1e3") # MeV
        df = df.Define("beam2_py", "getMom(Beam2Slice, 1)/1e3") # MeV
        df = df.Define("beam2_pz", "getMom(Beam2Slice, 2)/1e6") # GeV



        #df = df.Define("xing", "0.03")
        #df = df.Define("cosx", "cos(xing)")
        #df = df.Define("sinx", "sin(xing)")
        #df = df.Define("beam1_x_rot", "beam1_x*cosx + beam1_z*sinx")
        #df = df.Define("beam1_z_rot", "-beam1_x*sinx + beam1_z*cosx")

        df = df.Define("t", "ROOT::VecOps::Concatenate(beam1_slice, beam2_slice)")
        df = df.Define("x", "ROOT::VecOps::Concatenate(beam1_x, beam2_x)")
        df = df.Define("y", "ROOT::VecOps::Concatenate(beam1_y, beam2_y)")
        df = df.Define("z", "ROOT::VecOps::Concatenate(beam1_z, beam2_z)")

        #df = df.Define("vx", "ROOT::VecOps::Concatenate(beam1_vx, beam2_vx)")
        #df = df.Define("vy", "ROOT::VecOps::Concatenate(beam1_vy, beam2_vy)")
        #df = df.Define("vz", "ROOT::VecOps::Concatenate(beam1_vz, beam2_vz)")     
        #df = df.Define("test", "print(y)")
        #df = df.Filter("test")

        hists = []
        hists.append(df.Histo3D(("zx", "", *(bins_z+bins_x+bins_t)), "z", "x", "t"))
        hists.append(df.Histo3D(("zy", "", *(bins_z+bins_y+bins_t)), "z", "y", "t"))
        hists.append(df.Histo3D(("xy", "", *(bins_x+bins_y+bins_t)), "x", "y", "t"))
        hists.append(df.Histo3D(("zx1", "", *(bins_z+bins_x+bins_t)), "beam1_z", "beam1_x", "beam1_slice"))
        hists.append(df.Histo3D(("zy1", "", *(bins_z+bins_y+bins_t)), "beam1_z", "beam1_y", "beam1_slice"))
        hists.append(df.Histo3D(("xy1", "", *(bins_x+bins_y+bins_t)), "beam1_x", "beam1_y", "beam1_slice"))
        hists.append(df.Histo3D(("zx2", "", *(bins_z+bins_x+bins_t)), "beam2_z", "beam2_x", "beam2_slice"))
        hists.append(df.Histo3D(("zy2", "", *(bins_z+bins_y+bins_t)), "beam2_z", "beam2_y", "beam2_slice"))
        hists.append(df.Histo3D(("xy2", "", *(bins_x+bins_y+bins_t)), "beam2_x", "beam2_y", "beam2_slice"))
        hists.append(df.Histo2D(("x1", "", *(bins_x+bins_t)),  "beam1_x", "beam1_slice"))
        hists.append(df.Histo2D(("x2", "", *(bins_x+bins_t)),  "beam2_x", "beam2_slice"))
        hists.append(df.Histo2D(("y1", "", *(bins_y+bins_t)),  "beam1_y", "beam1_slice"))
        hists.append(df.Histo2D(("y2", "", *(bins_y+bins_t)),  "beam2_y", "beam2_slice"))
        hists.append(df.Histo2D(("z1", "", *(bins_z+bins_t)),  "beam1_z", "beam1_slice"))
        hists.append(df.Histo2D(("z2", "", *(bins_z+bins_t)),  "beam2_z", "beam2_slice"))

        hists.append(df.Histo2D(("px1", "", *((20000, -100, 100)+bins_t)),  "beam1_px", "beam1_slice"))
        hists.append(df.Histo2D(("py1", "", *((20000, -100, 100)+bins_t)),  "beam1_py", "beam1_slice"))
        hists.append(df.Histo2D(("pz1", "", *((10000, 40, 50)+bins_t)),  "beam1_pz", "beam1_slice"))

        hists.append(df.Histo2D(("px2", "", *((20000, -100, 100)+bins_t)),  "beam2_px", "beam2_slice"))
        hists.append(df.Histo2D(("py2", "", *((20000, -100, 100)+bins_t)),  "beam2_py", "beam2_slice"))
        hists.append(df.Histo2D(("pz2", "", *((10000, 40, 50)+bins_t)),  "beam2_pz", "beam2_slice"))

        print("**** RunGraphs")
        ROOT.RDF.RunGraphs(hists)
        #h_zx, h_zy, h_xy, h_x, h_y, h_z = hists[0], hists[1], hists[2], hists[3], hists[5], hists[7]

        print("**** Save histograms")
        fOut = ROOT.TFile(f"{tag}.root", "RECREATE")
        for h in hists:
            print(" write", h.GetName())
            h.Write()
        fOut.Close()


    fIn = ROOT.TFile(f"{tag}.root")
    h_zx1 = fIn.Get("zx1")
    h_zy1 = fIn.Get("zy1")
    h_xy1 = fIn.Get("xy1")
    h_zx2 = fIn.Get("zx2")
    h_zy2 = fIn.Get("zy2")
    h_xy2 = fIn.Get("xy2")

    h_x1 = fIn.Get("x1")
    h_y1 = fIn.Get("y1")
    h_z1 = fIn.Get("z1")
    h_x2 = fIn.Get("x2")
    h_y2 = fIn.Get("y2")
    h_z2 = fIn.Get("z2")


    h_px1 = fIn.Get("px1")
    h_py1 = fIn.Get("py1")
    h_pz1 = fIn.Get("pz1")
    h_px2 = fIn.Get("px2")
    h_py2 = fIn.Get("py2")
    h_pz2 = fIn.Get("pz2")

    print("**** Plotting")
    c = ROOT.TCanvas("c", "", 800, 800)
    c.SetBottomMargin(0.10)
    c.SetLeftMargin(0.15)
    c.SetGrid()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kBird)

    # hourglass envelopes (betay must be in mm)
    g_up, g_dn = hourglass_graphs(sigmay, betay/1e3, min_z_/1e6, max_z_/1e6, 4)


    dir_zx = f"{output_dir}/zx/"
    dir_zy = f"{output_dir}/zy/"
    dir_xy = f"{output_dir}/xy/"
    dir_x = f"{output_dir}/x/"
    dir_y = f"{output_dir}/y/"
    dir_z = f"{output_dir}/z/"
    dir_px = f"{output_dir}/px/"
    dir_py = f"{output_dir}/py/"
    dir_pz = f"{output_dir}/pz/"
    if not os.path.exists(dir_zx):
        os.makedirs(dir_zx)
        os.system(f"cp {index_php} {dir_zx}")
    if not os.path.exists(dir_zy):
        os.makedirs(dir_zy)
        os.system(f"cp {index_php} {dir_zy}")
    if not os.path.exists(dir_xy):
        os.makedirs(dir_xy)
        os.system(f"cp {index_php} {dir_xy}")
    if not os.path.exists(dir_x):
        os.makedirs(dir_x)
        os.system(f"cp {index_php} {dir_x}")
    if not os.path.exists(dir_y):
        os.makedirs(dir_y)
        os.system(f"cp {index_php} {dir_y}")
    if not os.path.exists(dir_z):
        os.makedirs(dir_z)
        os.system(f"cp {index_php} {dir_z}")
    if not os.path.exists(dir_px):
        os.makedirs(dir_px)
        os.system(f"cp {index_php} {dir_px}")
    if not os.path.exists(dir_py):
        os.makedirs(dir_py)
        os.system(f"cp {index_php} {dir_py}")
    if not os.path.exists(dir_pz):
        os.makedirs(dir_pz)
        os.system(f"cp {index_php} {dir_pz}")

    beam1_color = ROOT.TColor.GetColor("#0072B2")  # blue
    beam2_color = ROOT.TColor.GetColor("#D55E00")  # orange/vermilion

    nsteps = int(n_z)+1
    for step in range(1, nsteps):
        print(step)

        ################################################################################
        h_zx1.GetZaxis().SetRange(step, step)
        h1 = h_zx1.Project3D("yx")
        h1.SetName(f"zx_{step:03d}")
        h1.SetTitle(f"Step {step:03d}/{nsteps-1:03d}")
        h1.GetXaxis().SetTitle("z (mm)")
        h1.GetYaxis().SetTitle("x (#mum)")

        h_zx2.GetZaxis().SetRange(step, step)
        h2 = h_zx2.Project3D("yx")

        h1.SetMarkerColor(beam1_color)
        h2.SetMarkerColor(beam2_color)

        h1.SetMarkerStyle(ROOT.kDot)
        h2.SetMarkerStyle(ROOT.kDot)

        h1.Draw("SCAT=0.5")
        h2.Draw("SCAT=0.5 SAME")

        box = ROOT.TBox(min_z_/1e6, min_x_/1e3, max_z_/1e6, max_x_/1e3)
        box.SetFillStyle(0)
        box.SetLineWidth(2)
        box.SetLineColor(ROOT.kBlack)
        box.Draw("SAME")

        header_2d(label)

        c.SaveAs(f"{dir_zx}/step_{step:03d}.png")
        c.SaveAs(f"{dir_zx}/step_{step:03d}.pdf")
        c.Clear()

        ################################################################################
        h_zy1.GetZaxis().SetRange(step, step)
        h1 = h_zy1.Project3D("yx")
        h1.SetName(f"zy_{step:03d}")
        h1.SetTitle(f"Step {step:03d}/{nsteps-1:03d}")
        h1.GetXaxis().SetTitle("z (mm)")
        h1.GetYaxis().SetTitle("y (nm)")

        h_zy2.GetZaxis().SetRange(step, step)
        h2 = h_zy2.Project3D("yx")

        h1.SetMarkerColor(beam1_color)
        h2.SetMarkerColor(beam2_color)

        h1.SetMarkerStyle(ROOT.kDot)
        h2.SetMarkerStyle(ROOT.kDot)

        h1.Draw("SCAT=0.5")
        h2.Draw("SCAT=0.5 SAME")


        # draw hourglass
        g_up.Draw("L SAME")
        g_dn.Draw("L SAME")

        box = ROOT.TBox(min_z_/1e6, min_y_, max_z_/1e6, max_y_)
        box.SetFillStyle(0)
        box.SetLineWidth(2)
        box.SetLineColor(ROOT.kBlack)
        box.Draw("SAME")

        header_2d(label)

        c.SaveAs(f"{dir_zy}/step_{step:03d}.png")
        c.SaveAs(f"{dir_zy}/step_{step:03d}.pdf")
        c.Clear()

        ################################################################################
        h_xy1.GetZaxis().SetRange(step, step)
        h1 = h_xy1.Project3D("yx")
        h1.SetName(f"xy_{step:03d}")
        h1.SetTitle(f"Step {step:03d}/{nsteps-1:03d}")
        h1.GetXaxis().SetTitle("x (#mum)")
        h1.GetYaxis().SetTitle("y (nm)")

        h_xy2.GetZaxis().SetRange(step, step)
        h2 = h_xy2.Project3D("yx")

        h1.SetMarkerColor(beam1_color)
        h2.SetMarkerColor(beam2_color)

        h1.SetMarkerStyle(ROOT.kDot)
        h2.SetMarkerStyle(ROOT.kDot)

        

        h1.Draw("SCAT=0.5")
        h2.Draw("SCAT=0.5 SAME")


        box = ROOT.TBox(min_x_/1e3, min_y_, max_x_/1e3, max_y_)
        box.SetFillStyle(0)
        box.SetLineWidth(2)
        box.SetLineColor(ROOT.kBlack)
        box.Draw("SAME")

        header_2d(label)

        c.SaveAs(f"{dir_xy}/step_{step:03d}.png")
        c.SaveAs(f"{dir_xy}/step_{step:03d}.pdf")
        c.Clear()

        ################################################################################
        h = h_x1.ProjectionX(f"x_{step:03d}", step, step)
        h.SetTitle(f"STEP {step:03d}")
        h.GetXaxis().SetTitle("x (um)")
        h.GetYaxis().SetTitle("Counts")
        h.Draw("HIST")
        c.SaveAs(f"{dir_x}/step_{step:03d}.png")
        c.SaveAs(f"{dir_x}/step_{step:03d}.pdf")
        c.Clear()

        h = h_y1.ProjectionX(f"y_{step:03d}", step, step)
        h.SetTitle(f"STEP {step:03d}")
        h.GetXaxis().SetTitle("y (nm)")
        h.GetYaxis().SetTitle("Counts")
        h.Draw("HIST")
        c.SaveAs(f"{dir_y}/step_{step:03d}.png")
        c.SaveAs(f"{dir_y}/step_{step:03d}.pdf")
        c.Clear()

        h = h_z1.ProjectionX(f"z_{step:03d}", step, step)
        h.SetTitle(f"STEP {step:03d}")
        h.GetXaxis().SetTitle("z (mm)")
        h.GetYaxis().SetTitle("Counts")
        h.Draw("HIST")
        c.SaveAs(f"{dir_z}/step_{step:03d}.png")
        c.SaveAs(f"{dir_z}/step_{step:03d}.pdf")
        c.Clear()

        ## x momenta
        xmin, xmax, rebin = -20, 20, 10
        h1 = h_px1.ProjectionX(f"px1_{step:03d}", step, step)
        h1.SetTitle(f"STEP {step:03d}")
        h1.GetXaxis().SetTitle("p_{x} (MeV)")
        h1.GetYaxis().SetTitle("Counts")
        h1.Rebin(rebin)
        h1.GetXaxis().SetRangeUser(xmin, xmax)
        h1.Draw("HIST")

        h2 = h_px2.ProjectionX(f"px2_{step:03d}", step, step)
        h2.SetLineColor(ROOT.kRed)
        h2.Rebin(rebin)
        h1.GetXaxis().SetRangeUser(xmin, xmax)
        h2.Draw("SAME HIST")

        fgaus1 = ROOT.TF1("fgaus1", "gaus", h1.GetXaxis().GetXmin(), h1.GetXaxis().GetXmax())
        h1.Fit(fgaus1, "R")
        mean1 = fgaus1.GetParameter(1)
        mean_err1 = fgaus1.GetParError(1)

        fgaus2 = ROOT.TF1("fgaus2", "gaus", h1.GetXaxis().GetXmin(), h1.GetXaxis().GetXmax())
        h1.Fit(fgaus2, "R")
        mean2 = fgaus2.GetParameter(1)
        mean_err2 = fgaus2.GetParError(1)

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.DrawLatex(0.60, 0.85, f"#mu = {mean1:.3f} #pm {mean_err1:.3f}")
        latex.DrawLatex(0.60, 0.80, f"#mu = {mean2:.3f} #pm {mean_err2:.3f}")       

        c.SaveAs(f"{dir_px}/step_{step:03d}.png")
        c.SaveAs(f"{dir_px}/step_{step:03d}.pdf")
        c.Clear()


        ## y momenta
        xmin, xmax, rebin = -20, 20, 10
        h1 = h_py1.ProjectionX(f"py1_{step:03d}", step, step)
        h1.SetTitle(f"STEP {step:03d}")
        h1.GetXaxis().SetTitle("p_{y} (MeV)")
        h1.GetYaxis().SetTitle("Counts")
        h1.Rebin(rebin)
        h1.GetXaxis().SetRangeUser(xmin, xmax)
        h1.Draw("HIST")

        h2 = h_py2.ProjectionX(f"py2_{step:03d}", step, step)
        h2.SetLineColor(ROOT.kRed)
        h2.Rebin(rebin)
        h1.GetXaxis().SetRangeUser(xmin, xmax)
        h2.Draw("SAME HIST")

        fgaus1 = ROOT.TF1("fgaus1", "gaus", h1.GetXaxis().GetXmin(), h1.GetXaxis().GetXmax())
        h1.Fit(fgaus1, "R")
        mean1 = fgaus1.GetParameter(1)
        mean_err1 = fgaus1.GetParError(1)

        fgaus2 = ROOT.TF1("fgaus2", "gaus", h1.GetXaxis().GetXmin(), h1.GetXaxis().GetXmax())
        h1.Fit(fgaus2, "R")
        mean2 = fgaus2.GetParameter(1)
        mean_err2 = fgaus2.GetParError(1)

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.DrawLatex(0.60, 0.85, f"#mu = {mean1:.3f} #pm {mean_err1:.3f}")
        latex.DrawLatex(0.60, 0.80, f"#mu = {mean2:.3f} #pm {mean_err2:.3f}")       

        c.SaveAs(f"{dir_py}/step_{step:03d}.png")
        c.SaveAs(f"{dir_py}/step_{step:03d}.pdf")
        c.Clear()


        ## z momenta
        h1 = h_pz1.ProjectionX(f"pz1_{step:03d}", step, step)
        h1.SetTitle(f"STEP {step:03d}")
        h1.GetXaxis().SetTitle("pz (GeV)")
        h1.GetYaxis().SetTitle("Counts")
        h1.Draw("HIST")
        h2 = h_pz2.ProjectionX(f"pz2_{step:03d}", step, step)
        h2.SetLineColor(ROOT.kRed)
        h2.Draw("SAME HIST")

        c.SaveAs(f"{dir_pz}/step_{step:03d}.png")
        c.SaveAs(f"{dir_pz}/step_{step:03d}.pdf")
        c.Clear()

        #quit()

def header_2d(label):
    latex = ROOT.TLatex()
    latex.SetTextSize(0.035)
    latex.SetTextFont(42)
    latex.SetTextAlign(13)
    latex.DrawLatexNDC(0.15, 0.93, "#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}")
    
    latex = ROOT.TLatex()
    latex.SetTextSize(0.035)
    latex.SetTextFont(42)
    latex.SetTextAlign(33)
    latex.DrawLatexNDC(0.90, 0.93, label)


def hourglass_graphs(sigma_y_star, beta_y_star, zmin, zmax, nSigma):

    nPoints=2000

    g_up = ROOT.TGraph(nPoints)
    g_dn = ROOT.TGraph(nPoints)

    for i in range(nPoints):
        z = zmin + (zmax - zmin) * i / (nPoints - 1.0)
        sigma_y_z = sigma_y_star * math.sqrt(1.0 + (z*z)/(beta_y_star*beta_y_star))
        g_up.SetPoint(i, z, +nSigma * sigma_y_z)
        g_dn.SetPoint(i, z, -nSigma * sigma_y_z)

    # style
    for g in (g_up, g_dn):
        g.SetLineColor(ROOT.kBlack)
        g.SetLineWidth(2)
        g.SetLineStyle(3)

    return g_up, g_dn


if __name__ == "__main__":
    main()
