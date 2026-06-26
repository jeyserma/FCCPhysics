
import sys,copy,array,os,subprocess
import ROOT
import numpy as np
import argparse
import math

sys.path.insert(0, f'{os.path.dirname(os.path.realpath(__file__))}/../../../python')
import plotter
import functions

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


parser = argparse.ArgumentParser()
parser.add_argument("--mode", type=str, help="Detector mode", choices=["IDEA", "IDEA_MC", "IDEA_3T", "CLD", "CLD_FullSim", "IDEA_noBES", "IDEA_2E", "IDEA_BES6pct"], default="IDEA")
parser.add_argument("--lumi", type=float, help="Luminosity scale", default=-1.0)
parser.add_argument("--statOnly", action='store_true', help="Stat-only fit")
parser.add_argument("--ecm", type=int, help="Center-of-mass", choices=[-1, 240, 365], default=-1)
parser.add_argument("--tag", type=str, help="Analysis tag", default="")
parser.add_argument('--combination', action='store_true', help="Do full combination")
args = parser.parse_args()

def findCrossing(xv, yv, left=True, flip=125, cross=1.):

    closestPoint, idx = 1e9, -1
    for i in range(0, len(xv)):
        if left and xv[i] > flip: continue
        if not left and xv[i] < flip: continue
        dy = abs(yv[i]-cross)
        if dy < closestPoint: 
            closestPoint = dy
            idx = i
    # find correct indices around crossing
    if left: 
        if yv[idx] > cross: idx_ = idx+1
        else: idx_ = idx-1
    else:
        if yv[idx] > cross: idx_ = idx-1
        else: idx_ = idx+1

    # do interpolation
    omega = (yv[idx]-yv[idx_])/(xv[idx]-xv[idx_])
    return (cross-yv[idx])/omega + xv[idx]

def analyzeMass(runDir, outDir, xMin=-1, xMax=-1, yMin=0, yMax=2, label="label"):

    if not os.path.exists(outDir): os.makedirs(outDir)

    fIn = ROOT.TFile(f"{runDir}/higgsCombinemass.MultiDimFit.mH125.root", "READ")
    t = fIn.Get("limit")

    str_out = ""

    xv, yv = [], []
    for i in range(0, t.GetEntries()):

        t.GetEntry(i)
        
        if t.quantileExpected < -1.5: continue
        if t.deltaNLL > 1000: continue
        if t.deltaNLL > 20: continue
        xv.append(getattr(t, "MH"))
        yv.append(t.deltaNLL*2.)

    xv, yv = zip(*sorted(zip(xv, yv)))
    g = ROOT.TGraph(len(xv), array.array('d', xv), array.array('d', yv))

    # bestfit = minimum
    mass = 1e9
    for i in range(g.GetN()):
        if g.GetY()[i] == 0.: mass = g.GetX()[i]

    # extract uncertainties at crossing = 1
    unc_m = findCrossing(xv, yv, left=True, flip=mass)
    unc_p = findCrossing(xv, yv, left=False, flip=mass)
    unc = 0.5*(abs(mass-unc_m) + abs(unc_p-mass))

    ########### PLOTTING ###########
    cfg = {

        'logy'              : False,
        'logx'              : False,

        'xmin'              : min(xv) if xMin < 0 else xMin,
        'xmax'              : max(xv) if xMax < 0 else xMax,
        'ymin'              : yMin,
        'ymax'              : yMax , # max(yv)

        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "-2#DeltaNLL",

        'topRight'          : topRight, 
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Internal}}",
        }

    plotter.cfg = cfg

    canvas = plotter.canvas()
    canvas.SetGrid()
    dummy = plotter.dummy()

    dummy.GetXaxis().SetNdivisions(507)
    dummy.Draw("HIST")

    g.SetMarkerStyle(20)
    g.SetMarkerColor(ROOT.kRed)
    g.SetMarkerSize(1)
    g.SetLineColor(ROOT.kRed)
    g.SetLineWidth(2)
    g.Draw("SAME LP")

    line = ROOT.TLine(float(cfg['xmin']), 1, float(cfg['xmax']), 1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(2)
    line.Draw("SAME")

    leg = ROOT.TLegend(.20, 0.825, 0.90, .9)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.035)
    leg.SetMargin(0.15)
    leg.SetBorderSize(1)
    leg.AddEntry(g, "%s, #delta(m_{h}) = %.2f MeV" % (label, unc*1000.), "LP")
    leg.Draw()

    plotter.aux()
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/mass{suffix}.png")
    canvas.SaveAs(f"{outDir}/mass{suffix}.pdf")


    # write values to text file
    str_out = f"{unc_m} {unc_p} {unc} {mass}\n"
    for i in range(0, len(xv)): str_out += f"{xv[i]} {yv[i]}\n"
    tFile = open(f"{outDir}/mass{suffix}.txt", "w")
    tFile.write(str_out)
    tFile.close()
    tFile = open(f"{runDir}/mass{suffix}.txt", "w")
    tFile.write(str_out)
    tFile.close()


def doFit_mass(runDir, mhMin=124.99, mhMax=125.01, npoints=50, combineOptions = ""):

    cmd = f"combine -M MultiDimFit -t -1 --setParameterRanges MH={mhMin},{mhMax} --points={npoints} --algo=grid ws.root --expectSignal=1 -m 125 --redefineSignalPOIs MH --X-rtd TMCSO_AdaptivePseudoAsimov -v 10 --X-rtd ADDNLL_CBNLL=0 -n mass {combineOptions}"
    cmd_tot = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c '{cmd}'"
    subprocess.call(cmd_tot, shell=True, cwd=runDir)


def doFitDiagnostics_mass(runDir, mhMin=124.99, mhMax=125.01, combineOptions = ""):

    #cmd = f"combine -M FitDiagnostics -t -1 --setParameterRanges MH={mhMin},{mhMax} ws.root --expectSignal=1 -v 10  -m 125 --redefineSignalPOIs MH --floatParameters MH --X-rtd TMCSO_AdaptivePseudoAsimov --X-rtd ADDNLL_CBNLL=0 -n mass {combineOptions}"
    cmd = f"combine -M MultiDimFit --algo singles -t -1 --setParameterRanges MH={mhMin},{mhMax} ws.root --expectSignal=1 -v 10  -m 125 --redefineSignalPOIs MH --floatParameters MH --X-rtd TMCSO_AdaptivePseudoAsimov --X-rtd ADDNLL_CBNLL=0 -n mass {combineOptions}"

    cmd_tot = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c '{cmd}'"
    subprocess.call(cmd_tot, shell=True, cwd=runDir)

    # get the uncertainty
    fIn = ROOT.TFile(f"{runDir}/higgsCombinemass.MultiDimFit.mH125.root")
    fIn.ls()
    tt = fIn.Get("limit")
    vals = []
    for i in range(tt.GetEntries()):
        tt.GetEntry(i)
        vals.append(float(tt.MH))

    vals = sorted(vals)
    lo   = vals[0]
    best = vals[1]
    hi   = vals[2]

    err_down = best - lo
    err_up   = hi - best
    err_avg  = 0.5 * (err_up + err_down)
    fIn.Close()
    return err_avg



def plotMultiple(tags, labels, fOut, xMin=-1, xMax=-1, yMin=0, yMax=2, legLabel="", forceStat=[], legMargin=0.15):

    global suffix
    best_mass, best_xsec = [], []
    unc_mass, unc_xsec = [], []
    g_mass, g_xsec = [], []
    statOnly = False
    if len(forceStat) == 0:
        statOnly = True
        forceStat = [False]*len(tags)

    for i, tag in enumerate(tags):
        xv, yv = [], []
        fIn = open("%s/mass%s.txt" % (tag, suffix+"_stat" if forceStat[i] else suffix), "r")
        for i,line in enumerate(fIn.readlines()):

            line = line.rstrip()
            if i == 0: 
                best_mass.append(float(line.split(" ")[3]))
                unc_mass.append(float(line.split(" ")[2]))
            else:
                xv.append(float(line.split(" ")[0]))
                yv.append(float(line.split(" ")[1]))

        g = ROOT.TGraph(len(xv), array.array('d', xv), array.array('d', yv))    
        g_mass.append(g)


    ########### PLOTTING ###########
    cfg = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,

        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "-2#DeltaNLL",

        'topRight'          : topRight,
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}",
        }

    plotter.cfg = cfg

    canvas = plotter.canvas()
    canvas.SetGrid()
    dummy = plotter.dummy()

    dummy.GetXaxis().SetNdivisions(507)  
    dummy.Draw("HIST")

    n = len(g_mass) + (0 if legLabel=="" else 1)
    leg = ROOT.TLegend(.20, 0.9-n*0.05, 0.90, .9)
    leg.SetBorderSize(0)
    #leg.SetFillStyle(0) 
    leg.SetTextSize(0.03)
    leg.SetMargin(legMargin)
    leg.SetBorderSize(1)
    if legLabel != "":
        leg.SetHeader(legLabel)

    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+1]
    for i,g in enumerate(g_mass):
        g.SetMarkerStyle(20)
        g.SetMarkerColor(colors[i])
        g.SetMarkerSize(1)
        g.SetLineColor(colors[i])
        g.SetLineWidth(4)
        g.Draw("SAME L")
        leg.AddEntry(g, "%s #delta(m_{h}) = %.2f MeV" % (labels[i], unc_mass[i]*1000.), "L")

    leg.Draw()
    line = ROOT.TLine(float(cfg['xmin']), 1, float(cfg['xmax']), 1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(2)
    line.Draw("SAME")

    plotter.aux()
    canvas.Modify()
    canvas.Update()
    canvas.Draw()

    canvas.SaveAs("%s%s.png" % (fOut, suffix))
    canvas.SaveAs("%s%s.pdf" % (fOut, suffix))

def plotMultiple_xsec(tags, labels, fOut, xMin=-1, xMax=-1, yMin=0, yMax=2):

    best_xsec, unc_xsec, g_xsec = [], [], []
    for tag in tags:
        xv, yv = [], []
        fIn = open("%s/xsec.txt" % (tag), "r")
        for i,line in enumerate(fIn.readlines()):

            line = line.rstrip()
            if i == 0: 
                best_xsec.append(float(line.split(" ")[3]))
                unc_xsec.append(float(line.split(" ")[2]))
            else:
                xv.append(float(line.split(" ")[0]))
                yv.append(float(line.split(" ")[1]))

        g = ROOT.TGraph(len(xv), array.array('d', xv), array.array('d', yv))    
        g_xsec.append(g)

    cfg = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,

        'xtitle'            : "#sigma(ZH#rightarrowl^{#plus}l^{#minus})/#sigma_{ref}",
        'ytitle'            : "-2#DeltaNLL",

        'topRight'          : topRight, 
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}",
        }

    plotter.cfg = cfg

    canvas = plotter.canvas()
    canvas.SetGrid()
    dummy = plotter.dummy()

    dummy.GetXaxis().SetNdivisions(507)
    dummy.Draw("HIST")

    totEntries = len(g_xsec)
    leg = ROOT.TLegend(.20, 0.9-totEntries*0.05, 0.90, .9)
    leg.SetBorderSize(0)
    #leg.SetFillStyle(0) 
    leg.SetTextSize(0.03)
    leg.SetMargin(0.15)
    leg.SetBorderSize(1)

    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+1]
    for i,g in enumerate(g_xsec):
        g.SetMarkerStyle(20)
        g.SetMarkerColor(colors[i])
        g.SetMarkerSize(1)
        g.SetLineColor(colors[i])
        g.SetLineWidth(4)
        g.Draw("SAME L")
        leg.AddEntry(g, "%s #delta(#sigma) = %.2f %%" % (labels[i], unc_xsec[i]*100.), "L")

    leg.Draw()
    line = ROOT.TLine(float(cfg['xmin']), 1, float(cfg['xmax']), 1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(2)
    line.Draw("SAME")

    plotter.aux()
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s.png" % (fOut))
    canvas.SaveAs("%s.pdf" % (fOut))




def breakDown(outDir_):
    def getUnc(tag, type_):

        xv, yv = [], []
        fIn = open("%s/%s%s.txt" % (outDir_, type_, tag), "r")
        for i,line in enumerate(fIn.readlines()):

            line = line.rstrip()
            if i == 0: 
                best = float(line.split(" ")[3])
                unc = float(line.split(" ")[2])
                break
                
        if type_ == "mass": unc*= 1000. # convert to MeV
        if type_ == "xsec": unc*= 100. # convert to %
        return best, unc

    ############# mass
    canvas = ROOT.TCanvas("c", "c", 1000, 1000)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.1)
    canvas.SetLeftMargin(0.25)
    canvas.SetRightMargin(0.05)
    canvas.SetFillStyle(4000) # transparency?
    canvas.SetGrid(1, 0)
    canvas.SetTickx(1)


    xMin, xMax = -3, 3
    xTitle = "#sigma_{syst.}(m_{h}) (MeV)"

    ref = ""
    best_ref, unc_ref = getUnc("_stat", "mass")
    params = ["_BES", "_SQRTS", "_LEPSCALE_MU", "_LEPSCALE_EL", ""]
    labels = ["BES", "#sqrt{s} #pm 2 MeV", "Muon scale (~10^{-5})", "El. scale (~10^{-5})", "Syst. combined"]

    n_params = len(params)
    h_pulls = ROOT.TH2F("pulls", "pulls", 6, xMin, xMax, n_params, 0, n_params)
    g_pulls = ROOT.TGraphAsymmErrors(n_params)

    i = n_params
    for p in range(n_params):
        i -= 1
        best, unc = getUnc(params[p], "mass")
        #unc = math.sqrt(unc_ref**2 - unc**2)
        print(unc_ref, unc)
        unc = math.sqrt(unc**2 - unc_ref**2)
        g_pulls.SetPoint(i, 0, float(i) + 0.5)
        g_pulls.SetPointError(i, unc, unc, 0., 0.)
        h_pulls.GetYaxis().SetBinLabel(i + 1, "#splitline{%s}{(%.2f MeV)}" % (labels[p], unc))


    h_pulls.GetXaxis().SetTitleSize(0.04)
    h_pulls.GetXaxis().SetLabelSize(0.03)
    h_pulls.GetXaxis().SetTitle(xTitle)
    h_pulls.GetXaxis().SetTitleOffset(1)
    h_pulls.GetYaxis().SetLabelSize(0.045)
    h_pulls.GetYaxis().SetTickLength(0)
    h_pulls.GetYaxis().LabelsOption('v')
    h_pulls.SetNdivisions(506, 'XYZ')
    h_pulls.Draw("HIST")


    g_pulls.SetMarkerSize(0.8)
    g_pulls.SetMarkerStyle(20)
    g_pulls.SetLineWidth(2)
    g_pulls.Draw('P SAME')
    
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(30) # 0 special vertical aligment with subscripts
    latex.DrawLatex(0.95, 0.925, topRight)

    latex.SetTextAlign(13)
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.25, 0.96, "#bf{FCCee} #scale[0.7]{#it{Simulation}}")

    canvas.SaveAs("%s/mass_breakdown_impacts.png" % (outDir_))
    canvas.SaveAs("%s/mass_breakdown_impacts.pdf" % (outDir_))
    del canvas, g_pulls, h_pulls

    #################################
    ## lln plots

    ref = ""
    best_ref, unc_ref = getUnc("_stat", "mass")
    params = ["_stat", "_BES", "_LEPSCALE", "_SQRTS", ""]
    labels = ["Stat. only", "Beam energy spread", "Lepton scale", "Center-of-mass energy", "Stat. + syst. combined"]

    best_mass, unc_mass, g_mass = [], [], []

    tags = [f"{outDir_}/mass{p}.txt" for p in params]
    for i, tag in enumerate(tags):
        xv, yv = [], []
        fIn = open(tag, "r")
        for i,line in enumerate(fIn.readlines()):

            line = line.rstrip()
            if i == 0: 
                best_mass.append(float(line.split(" ")[3]))
                unc_mass.append(float(line.split(" ")[2]))
            else:
                xv.append(float(line.split(" ")[0]))
                yv.append(float(line.split(" ")[1]))

        g = ROOT.TGraph(len(xv), array.array('d', xv), array.array('d', yv))    
        g_mass.append(g)


    ########### PLOTTING ###########
    cfg = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : 124.995,
        'xmax'              : 125.005,
        'ymin'              : 0,
        'ymax'              : 2,

        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "-2#DeltaNLL",

        'topRight'          : topRight,
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}",
        }

    plotter.cfg = cfg

    canvas = plotter.canvas()
    canvas.SetGrid()
    dummy = plotter.dummy()

    dummy.GetXaxis().SetNdivisions(507)  
    dummy.Draw("HIST")

    legLabel = ""
    legMargin = 0.1
    n = len(g_mass) + (0 if legLabel=="" else 1)
    leg = ROOT.TLegend(.20, 0.9-n*0.05, 0.90, .9)
    leg.SetBorderSize(0)
    #leg.SetFillStyle(0) 
    leg.SetTextSize(0.03)
    leg.SetMargin(legMargin)
    leg.SetBorderSize(1)
    if legLabel != "":
        leg.SetHeader(legLabel)

    
    line = ROOT.TLine(float(cfg['xmin']), 1, float(cfg['xmax']), 1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(2)
    line.Draw("SAME")

    colors = [
        ROOT.kRed+1,     # stat only
        ROOT.kBlue-4,    # beam energy spread
        ROOT.kOrange-3,  # lepton scale
        ROOT.kGreen+3,   # center-of-mass energy
        ROOT.kBlack      # stat + syst combined
    ]

    widths = [4, 2, 2, 2, 4]
    styles = [1, 2, 2, 2, 1]
    for i,g in enumerate(g_mass):
        g.SetMarkerStyle(styles[i])
        g.SetMarkerColor(colors[i])
        g.SetMarkerSize(1)
        g.SetLineColor(colors[i])
        g.SetLineWidth(widths[i])
        g.Draw("SAME L")
        leg.AddEntry(g, "%s #delta(m_{h}) = %.2f MeV" % (labels[i], unc_mass[i]*1000.), "L")

    leg.Draw()

    plotter.aux()
    canvas.Modify()
    canvas.Update()
    canvas.Draw()

    canvas.SaveAs("%s/mass_breakdown.png" % (outDir_))
    canvas.SaveAs("%s/mass_breakdown.pdf" % (outDir_))


def text2workspace(runDir):

    cmd = "text2workspace.py datacard.txt -o ws.root  -v 10 --X-allow-no-background"
    cmd_tot = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c '{cmd}'"
    subprocess.call(cmd_tot, shell=True, cwd=runDir)

def combineCards(runDir, input_=[]):

    if not os.path.exists(runDir): os.makedirs(runDir)
    
    input_ = ["%s/%s" % (os.getcwd(), i) for i in input_]
    cards = ' '.join(input_)

    cmd = "combineCards.py   --force-shape %s > datacard.txt" % cards
    cmd_tot = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c '{cmd}'"
    subprocess.call(cmd_tot, shell=True, cwd=runDir)
    text2workspace(runDir)


if __name__ == "__main__":
    mode = args.mode # detector mode
    ecm = str(args.ecm)
    tag = args.tag
    lumiScale = args.lumi
    lumiStr = str(int(args.lumi)) if args.lumi.is_integer() else str(args.lumi)
    lumiLabel = lumiStr.replace(".", "p")

    topRight = f"#sqrt{{s}} = {ecm} GeV, {lumiStr} ab^{{#minus1}}"
    topLeft = f"#bf{{FCC-ee}} #scale[0.7]{{#it{{Internal}}}}"

    combineDir = f"output/h_mass/combine/{tag}/{mode}/lumi{lumiLabel}/" # baseline combine dir
    outDir = f"/work/submit/jaeyserm/public_html/fccee/h_mass/{tag}/combine/{mode}/lumi{lumiLabel}/"



    def mHrange(mh_err):
        if 1.5*mh_err > 0.05: # bound to 50 MeV
            return 124.95, 125.05
        return 125-1.5*mh_err, 125+1.5*mh_err

    combineOptions = ""
    freezeParameters = []
    setParameters = []
    suffix=""

    doSyst = not args.statOnly
    if not doSyst:
        suffix = "_stat"
        #freezeParameters.extend(["BES", "ISR", "SQRTS", "LEPSCALE_MU", "LEPSCALE_EL"])
        freezeParameters.extend(["BES_ecm240", "SQRTS_ecm240", "LEPSCALE_MU_ecm240", "LEPSCALE_EL_ecm240"])
        freezeParameters.extend(["bkg_norm_mumu_ecm240", "bkg_norm_ee_ecm240", "bkg_norm_mumu_ecm365", "bkg_norm_ee_ecm365"])
        freezeParameters.extend(["BES_ecm365", "SQRTS_ecm365", "LEPSCALE_MU_ecm365", "LEPSCALE_EL_ecm365"])




    freezeBkg = False
    if freezeBkg:
        suffix = "_freezeBkg%s"%suffix
        freezeParameters.extend(["bkg_norm_mumu_ecm240", "bkg_norm_ee_ecm240", "bkg_norm_mumu_ecm365", "bkg_norm_ee_ecm365"])


    noBkg = False
    if noBkg:
        suffix = "_noBkg%s"%suffix
        freezeParameters.extend(["bkg_norm_mumu_ecm240", "bkg_norm_ee_ecm240", "bkg_norm_mumu_ecm365", "bkg_norm_ee_ecm365", "r"]) # r to be frozen to avoid issues with no-bkg fit
        setParameters.extend(["bkg_norm_mumu_ecm240=0", "bkg_norm_ee_ecm240=0", "bkg_norm_mumu_ecm365=0", "bkg_norm_ee_ecm365=0"])


    ## systematic variations, unfreeze them
    systs = ["BES_ecm240", "SQRTS_ecm240", "LEPSCALE_MU_ecm240", "LEPSCALE_EL_ecm240", "BES_ecm365", "SQRTS_ecm365", "LEPSCALE_MU_ecm365", "LEPSCALE_EL_ecm365"]

    #suffix+="_BES"
    #systs.remove("BES_ecm240")
    #systs.remove("BES_ecm365")

    #suffix+="_SQRTS"
    #systs.remove("SQRTS_ecm240")
    #systs.remove("SQRTS_ecm365")

    #suffix+="_LEPSCALE_MU"
    #systs.remove("LEPSCALE_MU_ecm240")
    #systs.remove("LEPSCALE_MU_ecm365")

    #suffix+="_LEPSCALE_EL"
    #systs.remove("LEPSCALE_EL_ecm240")
    #systs.remove("LEPSCALE_EL_ecm365")

    suffix+="_LEPSCALE"
    systs.remove("LEPSCALE_EL_ecm240")
    systs.remove("LEPSCALE_EL_ecm365")
    systs.remove("LEPSCALE_MU_ecm240")
    systs.remove("LEPSCALE_MU_ecm365")

    freezeParameters.extend(systs)


    doBreakDown = True
    if doBreakDown:
        breakDown(outDir + "/mumu_ee_combined_categorized_ecm240")
        quit()
        breakDown("mumu_ee_combined_inclusive")
        breakDown("mumu_cat0")
        breakDown("mumu_cat1")
        breakDown("mumu_cat2")
        breakDown("mumu_cat3")
        breakDown("mumu_combined")
        breakDown("ee_cat0")
        breakDown("ee_cat1")
        breakDown("ee_cat2")
        breakDown("ee_cat3")
        breakDown("ee_combined")
        quit()

    doSummary = False
    if doSummary:
        outDir__ = "/work/submit/jaeyserm/public_html/fccee/higgs_mass_xsec/%s/combine/" % args.tag
        outDir   = "/work/submit/jaeyserm/public_html/fccee/higgs_mass_xsec/%s/combine/summaryPlots/" % args.tag
        plotMultiple(["%s/IDEA/lumi10p8/mumu_combined_ecm240/"%outDir__, "%s/IDEA_MC/lumi10p8/mumu_combined_ecm240/"%outDir__, "%s/IDEA_3T/lumi10p8/mumu_combined_ecm240/"%outDir__, "%s/CLD/lumi10p8/mumu_combined_ecm240/"%outDir__], ["IDEA", "IDEA perfect resolution", "IDEA 3T", "IDEA CLD silicon tracker"], "%s/IDEA_IDEAL_2T_3T_CLD_mumu"%outDir, xMin=124.99, xMax=125.01, legLabel="Muon final state Z(#mu^{#plus}#mu^{#minus})H (stat. + syst.)")
        plotMultiple(["%s/IDEA/lumi10p8/mumu_combined_ecm240/"%outDir__, "%s/IDEA_MC/lumi10p8/mumu_combined_ecm240/"%outDir__, "%s/IDEA_3T/lumi10p8/mumu_combined_ecm240/"%outDir__, "%s/CLD/lumi10p8/mumu_combined_ecm240/"%outDir__], ["IDEA", "IDEA perfect resolution", "IDEA 3T", "IDEA CLD silicon tracker"], "%s/IDEA_IDEAL_2T_3T_CLD_mumu_stat"%outDir, xMin=124.99, xMax=125.01, legLabel="Muon final state Z(#mu^{#plus}#mu^{#minus})H (stat. only)", forceStat=[True, True, True, True])
        plotMultiple(["%s/IDEA/lumi10p8/mumu_ee_combined_categorized_ecm240/"%outDir__, "%s/IDEA/lumi10p8/mumu_ee_combined_categorized_ecm240/"%outDir__], ["Statistical", "Statistical+systematic"], "%s/IDEA_stat_syst"%outDir, xMin=124.995, xMax=125.005, legLabel="Combined muon and electron final states", forceStat=[True, False])
        quit()


    ##################################
    if len(freezeParameters) > 0:
        combineOptions += " --freezeParameters " + ",".join(freezeParameters)
    if len(setParameters) > 0:
        combineOptions += " --setParameters " + ",".join(setParameters)


    ############### MUON
    if True:


        tag_cat0, label,  = "mumu_cat0_ecm%s"%ecm, "#mu^{#plus}#mu^{#minus}, inclusive"
        rMin, rMax = 0.98, 1.02
        mh_err = doFitDiagnostics_mass("%s/%s" % (combineDir, tag_cat0), mhMin=124.97, mhMax=125.03, combineOptions=combineOptions)
        mhMin, mhMax = mHrange(mh_err)
        doFit_mass("%s/%s" % (combineDir, tag_cat0), mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
        analyzeMass("%s/%s" % (combineDir, tag_cat0), "%s/%s/" % (outDir, tag_cat0), label=label, xMin=mhMin, xMax=mhMax)

        tag_cat1, label = "mumu_cat1_ecm%s"%ecm, "#mu^{#plus}#mu^{#minus}, central-central"
        if ecm == "240": # at 365 it is beyond 50 MeV
            mh_err = doFitDiagnostics_mass("%s/%s" % (combineDir, tag_cat1), mhMin=124.95, mhMax=125.05, combineOptions=combineOptions)
            mhMin, mhMax = mHrange(mh_err)
            doFit_mass("%s/%s" % (combineDir, tag_cat1), mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
            analyzeMass("%s/%s" % (combineDir, tag_cat1), "%s/%s/" % (outDir, tag_cat1), label=label, xMin=mhMin, xMax=mhMax)

        tag_cat2, label = "mumu_cat2_ecm%s"%ecm, "#mu^{#plus}#mu^{#minus}, central-forward"
        
        if ecm == "240": # at 365 it is beyond 50 MeV
            mh_err = doFitDiagnostics_mass("%s/%s" % (combineDir, tag_cat2), mhMin=124.95, mhMax=125.05, combineOptions=combineOptions)
            mhMin, mhMax = mHrange(mh_err)
            doFit_mass("%s/%s" % (combineDir, tag_cat2), mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
            analyzeMass("%s/%s" % (combineDir, tag_cat2), "%s/%s/" % (outDir, tag_cat2), label=label, xMin=mhMin, xMax=mhMax)

        
        tag_cat3, label = "mumu_cat3_ecm%s"%ecm, "#mu^{#plus}#mu^{#minus}, forward-forward"
        
        if ecm == "240": # at 365 it is beyond 50 MeV
            mh_err = doFitDiagnostics_mass("%s/%s" % (combineDir, tag_cat3), mhMin=124.95, mhMax=125.05, combineOptions=combineOptions)
            mhMin, mhMax = mHrange(mh_err)
            doFit_mass("%s/%s" % (combineDir, tag_cat3), mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
            analyzeMass("%s/%s" % (combineDir, tag_cat3), "%s/%s/" % (outDir, tag_cat3), label=label, xMin=mhMin, xMax=mhMax)
        
        tag_combined, label = "mumu_combined_ecm%s"%ecm, "#mu^{#plus}#mu^{#minus}, combined"
        combineCards("%s/%s" % (combineDir, tag_combined), ["%s/%s/datacard.txt"%(combineDir,tag_cat1), "%s/%s/datacard.txt"%(combineDir,tag_cat2), "%s/%s/datacard.txt"%(combineDir,tag_cat3)])
        mh_err = doFitDiagnostics_mass("%s/%s" % (combineDir, tag_combined), mhMin=124.95, mhMax=125.05, combineOptions=combineOptions)
        mhMin, mhMax = mHrange(mh_err)
        doFit_mass("%s/%s" % (combineDir, tag_combined), mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
        analyzeMass("%s/%s" % (combineDir, tag_combined), "%s/%s/" % (outDir, tag_combined), label=label, xMin=mhMin, xMax=mhMax)


    ############### ELECTRON
    if True:

        tag_cat0, label,  = "ee_cat0_ecm%s"%ecm, "e^{#plus}e^{#minus}, inclusive"
        rMin, rMax = 0.98, 1.02
        mh_err = doFitDiagnostics_mass("%s/%s" % (combineDir, tag_cat0), mhMin=124.95, mhMax=125.05, combineOptions=combineOptions)
        #if ecm == 240: # at 365 it is beyond 50 MeV
        mhMin, mhMax = mHrange(mh_err)
        doFit_mass("%s/%s" % (combineDir, tag_cat0), mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
        analyzeMass("%s/%s" % (combineDir, tag_cat0), "%s/%s/" % (outDir, tag_cat0), label=label, xMin=mhMin, xMax=mhMax)

        tag_cat1, label = "ee_cat1_ecm%s"%ecm, "e^{#plus}e^{#minus}, central-central"
        if ecm == "240": # at 365 it is beyond 50 MeV
            mh_err = doFitDiagnostics_mass("%s/%s" % (combineDir, tag_cat1), mhMin=124.95, mhMax=125.05, combineOptions=combineOptions)
            mhMin, mhMax = mHrange(mh_err)
            doFit_mass("%s/%s" % (combineDir, tag_cat1), mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
            analyzeMass("%s/%s" % (combineDir, tag_cat1), "%s/%s/" % (outDir, tag_cat1), label=label, xMin=mhMin, xMax=mhMax)

        tag_cat2, label = "ee_cat2_ecm%s"%ecm, "e^{#plus}e^{#minus}, central-forward"

        if ecm == "240": # at 365 it is beyond 50 MeV
            mh_err = doFitDiagnostics_mass("%s/%s" % (combineDir, tag_cat2), mhMin=124.95, mhMax=125.05, combineOptions=combineOptions)
            mhMin, mhMax = mHrange(mh_err)
            doFit_mass("%s/%s" % (combineDir, tag_cat2), mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
            analyzeMass("%s/%s" % (combineDir, tag_cat2), "%s/%s/" % (outDir, tag_cat2), label=label, xMin=mhMin, xMax=mhMax)

        tag_cat3, label = "ee_cat3_ecm%s"%ecm, "e^{#plus}e^{#minus}, forward-forward"
        if ecm == "240": # at 365 it is beyond 50 MeV
            mh_err = doFitDiagnostics_mass("%s/%s" % (combineDir, tag_cat3), mhMin=124.95, mhMax=125.05, combineOptions=combineOptions)
            mhMin, mhMax = mHrange(mh_err)
            doFit_mass("%s/%s" % (combineDir, tag_cat3), mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
            analyzeMass("%s/%s" % (combineDir, tag_cat3), "%s/%s/" % (outDir, tag_cat3), label=label, xMin=mhMin, xMax=mhMax)

        tag_combined, label = "ee_combined_ecm%s"%ecm, "e^{#plus}e^{#minus}, combined"
        combineCards("%s/%s" % (combineDir, tag_combined), ["%s/%s/datacard.txt"%(combineDir,tag_cat1), "%s/%s/datacard.txt"%(combineDir,tag_cat2), "%s/%s/datacard.txt"%(combineDir,tag_cat3)])
        mh_err = doFitDiagnostics_mass("%s/%s" % (combineDir, tag_combined), mhMin=124.95, mhMax=125.05, combineOptions=combineOptions)
        mhMin, mhMax = mHrange(mh_err)
        doFit_mass("%s/%s" % (combineDir, tag_combined), mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
        analyzeMass("%s/%s" % (combineDir, tag_combined), "%s/%s/" % (outDir, tag_combined), label=label, xMin=mhMin, xMax=mhMax)


    ############### MUON+ELECTRON
    if True:
        if args.ecm < 0 and args.lumi < 0:
            print("ERROR: ecm and lumi to be defined")
            quit()

        # check if lumi and ecm are defined
        tag, label = "mumu_ee_combined_inclusive_ecm%s"%ecm, "#mu^{#plus}#mu^{#minus}+e^{#plus}e^{#minus}, inclusive"
        combineCards("%s/%s" % (combineDir, tag), ["%s/mumu_cat0_ecm%s/datacard.txt"%(combineDir, ecm), "%s/ee_cat0_ecm%s/datacard.txt"%(combineDir, ecm)])
        
        mh_err = doFitDiagnostics_mass("%s/%s" % (combineDir, tag), mhMin=124.95, mhMax=125.05, combineOptions=combineOptions)
        mhMin, mhMax = mHrange(mh_err)
        doFit_mass("%s/%s" % (combineDir, tag), mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
        analyzeMass("%s/%s" % (combineDir, tag), "%s/%s/" % (outDir, tag), label=label, xMin=mhMin, xMax=mhMax)
        #quit()

        tag, label = "mumu_ee_combined_categorized_ecm%s"%ecm, "#mu^{#plus}#mu^{#minus}+e^{#plus}e^{#minus}, categorized"
        combineCards("%s/%s" % (combineDir, tag), ["/%s/mumu_combined_ecm%s/datacard.txt"%(combineDir, ecm), "/%s/ee_combined_ecm%s/datacard.txt"%(combineDir, ecm)])
        mh_err = doFitDiagnostics_mass("%s/%s" % (combineDir, tag), mhMin=124.95, mhMax=125.05, combineOptions=combineOptions)
        mhMin, mhMax = mHrange(mh_err)
        doFit_mass("%s/%s" % (combineDir, tag), mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
        analyzeMass("%s/%s" % (combineDir, tag), "%s/%s/" % (outDir, tag), label=label, xMin=mhMin, xMax=mhMax)

        '''
        if ecm == "240":
            plotMultiple(["%s/mumu_cat0/"%outDir, "%s/ee_cat0/"%outDir, "%s/mumu_ee_combined_inclusive/"%outDir], ["#mu^{#plus}#mu^{#minus}, inclusive", "e^{#plus}e^{#minus}, inclusive", "#mu^{#plus}#mu^{#minus} + e^{#plus}e^{#minus}, inclusive"], "%s/mumu_ee_inclusive"%outDir, xMin=124.99, xMax=125.01)
            plotMultiple(["%s/mumu_combined/"%outDir, "%s/ee_combined/"%outDir, "%s/mumu_ee_combined_categorized/"%outDir], ["#mu^{#plus}#mu^{#minus}, categorized", "e^{#plus}e^{#minus}, categorized", "#mu^{#plus}#mu^{#minus} + e^{#plus}e^{#minus}, categorized"], "%s/mumu_ee_categorized"%outDir, xMin=124.99, xMax=125.01)
            plotMultiple(["%s/mumu_cat0/"%outDir, "%s/mumu_combined/"%outDir, "%s/ee_cat0/"%outDir, "%s/ee_combined/"%outDir], ["#mu^{#plus}#mu^{#minus}, inclusive", "#mu^{#plus}#mu^{#minus}, categorized", "e^{#plus}e^{#minus}, inclusive", "e^{#plus}e^{#minus}, categorized"], "%s/mumu_ee_inclusive_categorized"%outDir, xMin=124.99, xMax=125.01)
        '''

        plotMultiple(["%s/mumu_combined_ecm%s/"%(outDir,ecm), "%s/ee_combined_ecm%s/"%(outDir,ecm), "%s/mumu_ee_combined_categorized_ecm%s/"%(outDir,ecm)], ["#mu^{#plus}#mu^{#minus}", "e^{#plus}e^{#minus}", "#mu^{#plus}#mu^{#minus} + e^{#plus}e^{#minus}"], "%s/mumu_ee_combined_categorized_ecm%s"%(outDir,ecm), xMin=124.99, xMax=125.01)
        quit()
        #plotMultiple(["%s/mumu_combined_ecm%s/"%(outDir, ecm), "%s/ee_combined_ecm%s/"%(outDir, ecm), "%s/mumu_ee_combined_categorized_ecm%s/"%(outDir, ecm)], ["#mu^{#plus}#mu^{#minus}, categorized", "e^{#plus}e^{#minus}, categorized", "#mu^{#plus}#mu^{#minus} + e^{#plus}e^{#minus}, categorized"], "%s/mumu_ee_categorized"%outDir, xMin=124.95, xMax=125.05)


    # 240+365
    if args.combination:
    
        lumi_240, lumi_365 = "10p8", "3"
        fit_tag = "mumu_ee_combined_categorized"
        label = "#mu^{#plus}#mu^{#minus}+e^{#plus}e^{#minus}, categorized"
        #fit_tag = "mumu_ee_combined_inclusive"
        #label = "#mu^{#plus}#mu^{#minus}+e^{#plus}e^{#minus}, inclusive"

        #fit_tag = "mumu_combined"
        #label = "#mu^{#plus}#mu^{#minus}, categorized"
       
        #fit_tag = "ee_combined"
        #label = "e^{#plus}e^{#minus}, categorized"

        #fit_tag = "mumu_cat0"
        #label = "#mu^{#plus}#mu^{#minus}, inclusive"
        

        tag = "combined_ecm_%s_%s" % (lumi_240, lumi_365)
        topRight = "#sqrt{s} = 240/365 GeV, %s/%s ab^{#minus1}" % (lumi_240.replace("p", "."), lumi_365.replace("p", "."))

        base_dir = "combine/higgs_mass_%s/%s/" % (args.tag, mode)
        combineDir = "%s/%s/%s/" % (base_dir, tag, fit_tag)
        combineDir_240 = "%s/lumi%s/%s_ecm240" % (base_dir, lumi_240, fit_tag)
        combineDir_365 = "%s/lumi%s/%s_ecm365" % (base_dir, lumi_365, fit_tag)

        outDir_ = "/work/submit/jaeyserm/public_html/fccee/higgs_mass_xsec/%s/combine/%s/%s/%s/" % (args.tag, mode, tag, fit_tag)
        outDir_240 = "/work/submit/jaeyserm/public_html/fccee/higgs_mass_xsec/%s/combine/%s/lumi%s/%s_ecm240/" % (args.tag, mode, lumi_240, fit_tag)
        outDir_365 = "/work/submit/jaeyserm/public_html/fccee/higgs_mass_xsec/%s/combine/%s/lumi%s/%s_ecm365/" % (args.tag, mode, lumi_365, fit_tag)

        doBreakDown = True
        if doBreakDown:
            breakDown(outDir_)
            quit()

        combineCards(combineDir, [combineDir_240+"/datacard.txt", combineDir_365+"/datacard.txt"])
        #quit()
        mh_err = doFitDiagnostics_mass(combineDir, mhMin=124.95, mhMax=125.05, combineOptions=combineOptions)
        #quit()
        mhMin, mhMax = mHrange(mh_err)
        doFit_mass(combineDir, mhMin=mhMin, mhMax=mhMax, npoints=50, combineOptions=combineOptions)
        analyzeMass(combineDir, outDir_, label=label, xMin=mhMin, xMax=mhMax)
        plotMultiple([outDir_240, outDir_365, outDir_], ["#sqrt{s} = 240 GeV", "#sqrt{s} = 365 GeV", "Combination"], outDir_, outName="comparison", xMin=124.98, xMax=125.02)

    # Fast-Fullsim comparison
    if False:
        outName = "/work/submit/jaeyserm/public_html/fccee/higgs_mass_xsec/%s/combine/comparison" % (args.tag)
        outDir_fastsim = "/work/submit/jaeyserm/public_html/fccee/higgs_mass_xsec/%s/combine/CLD/lumi10p8/mumu_combined_ecm240/" % (args.tag)
        outDir_fullsim = "/work/submit/jaeyserm/public_html/fccee/higgs_mass_xsec/%s/combine/CLD_FullSim/lumi10p8/mumu_combined_ecm240/" % (args.tag)
        plotMultiple([outDir_fastsim, outDir_fullsim], ["IDEA CLD silicon tracker (Delphes)", "CLD FullSim"], outName, xMin=124.99, xMax=125.01, legMargin=0.075)