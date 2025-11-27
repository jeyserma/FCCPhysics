import sys,array,ROOT,math,os,copy
import argparse
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

parser = argparse.ArgumentParser()
parser.add_argument("--cats", type=str, help="Category (qq, mumu, ee)")
parser.add_argument("--ecm", type=int, help="Center-of-mass energy", default=240)
parser.add_argument("--run", help="Run combine", action='store_true')
parser.add_argument("--bbb", help="Enable bin-by-bin statistics (BB)", action='store_true')
parser.add_argument("--target", type=str, help="Target pseudodata", default="bb")
parser.add_argument("--pert", type=float, help="Target pseudodata size", default=1.0)
parser.add_argument("--freezeBackgrounds", help="Freeze backgrounds", action='store_true')
parser.add_argument("--floatBackgrounds", help="Float backgrounds", action='store_true')
args = parser.parse_args()

def getMetaInfo(proc):
    if "Hinv" in proc:
        proc = proc.replace("wzp6", "wz3p6")
    fIn = ROOT.TFile(f"{inputDir}/{proc}.root")
    xsec = fIn.Get("crossSection").GetVal()
    
    if "HZZ" in proc: # HZZ contains invisible, remove xsec
        xsec_inv = getMetaInfo(proc.replace("wzp6", "wz3p6").replace("HZZ", "Hinv"))
        print("REMOVE INV FROM ZZ XSEC", proc, xsec, xsec-xsec_inv)
        xsec = xsec - xsec_inv
    return xsec

def removeNegativeBins(hist):
    totNeg, tot = 0., 0.
    if "TH1" in hist.ClassName():
        pass
    elif "TH2" in hist.ClassName():
        pass

    elif "TH3" in hist.ClassName():
        nbinsX = hist.GetNbinsX()
        nbinsY = hist.GetNbinsY()
        nbinsZ = hist.GetNbinsZ()
        for x in range(1, nbinsX + 1):
            for y in range(1, nbinsY + 1):
                for z in range(1, nbinsZ + 1):
                    content = hist.GetBinContent(x, y, z)
                    error = hist.GetBinError(x, y, z)  # Retrieve bin error
                    tot += content
                    if content < 0:
                        #print("WARNING: NEGATIVE BIN CONTENT", content, hist.GetName())
                        totNeg += content
                        hist.SetBinContent(x, y, z, 0)
                        hist.SetBinError(x, y, z, 0)
    if totNeg != 0:
        print(f"WARNING: TOTAL {tot}, NEGATIVE {totNeg}, FRACTION {totNeg/tot}")
    return hist

def getSingleHist(hName, proc):
    if "Hinv" in proc:
        proc = proc.replace("wzp6", "wz3p6")
    fIn = ROOT.TFile(f"{inputDir}/{proc}.root")
    h = fIn.Get(hName)
    h.SetDirectory(0)
    fIn.Close()
    '''
    if "HWW" in proc: h.Scale(72.733/72.69)
    if "HZZ" in proc: h.Scale(72.733/72.84)
    if "HZa" in proc: h.Scale(72.733/72.90)
    if "Haa" in proc: h.Scale(72.733/72.66)
    if "Hinv" in proc: h.Scale(72.733/72.84)
    if "Hbb" in proc: h.Scale(72.733/72.61)
    if "Hcc" in proc: h.Scale(72.733/72.67)
    if "Hgg" in proc: h.Scale(72.733/72.74)
    if "Hss" in proc: h.Scale(72.733/72.73)
    if "Hmumu" in proc: h.Scale(72.733/72.66)
    if "Htautau" in proc: h.Scale(72.733/72.73)
    '''
    #if "HZZ" in proc:
    #    print("REMOVE INV FROM ZZ", proc)
    #    h_inv = getSingleHist(hName, proc.replace('wzp6', 'wz3p6').replace('HZZ', 'Hinv'))
    #    h.Add(h_inv, -1)
    h = removeNegativeBins(h) ## need to remove potential negative bins here at single process level (otherwise issues with pseudodata etc.)
    return h

def getHists(hName, procName):
    hist = None
    procs = procs_cfg[procName]
    for k,proc in enumerate(procs):
        h = getSingleHist(hName, proc)
        if hist == None:
            hist = h
        else:
            hist.Add(h)
        #fIn.Close()
    if procName in proc_scales:
        hist.Scale(proc_scales[procName])
        print(f"SCALE {procName} with factor {proc_scales[procName]}")
    return hist


def smooth_histogram_gaussian_kernel(hist, sigma, fOut):

    n_bins = hist.GetNbinsX()
    smoothed_values = []

    # Loop over all bins
    for i in range(1, n_bins + 1):
        bin_center = hist.GetBinCenter(i)
        kernel_sum = 0.0
        value_sum = 0.0

        # Apply Gaussian kernel to surrounding bins
        for j in range(1, n_bins + 1):
            neighbor_center = hist.GetBinCenter(j)
            weight = np.exp(-0.5 * ((bin_center - neighbor_center) / sigma) ** 2)
            value_sum += hist.GetBinContent(j) * weight
            kernel_sum += weight

        # Normalize the value by the kernel sum
        smoothed_values.append(value_sum / kernel_sum)

    smoothed_hist = hist.Clone(f"{hist.GetName()}_smoothed")
    smoothed_hist.Reset()  # Clear bin contents
    for i, value in enumerate(smoothed_values, start=1):
        smoothed_hist.SetBinContent(i, value)


    canvas = ROOT.TCanvas("", "", 800, 800)

    # Draw the original histogram
    hist.SetLineColor(ROOT.kBlue)
    hist.SetLineWidth(2)
    hist.Draw("HIST")

    smoothed_hist.SetLineColor(ROOT.kRed)
    smoothed_hist.SetLineWidth(2)
    smoothed_hist.Draw("HIST SAME")

    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    legend.AddEntry(hist, "Original", "l")
    legend.AddEntry(smoothed_hist, "Smoothed", "l")
    legend.Draw()

    canvas.SaveAs(f"{plot_dir}/{fOut}.png")
    canvas.SaveAs(f"{plot_dir}/{fOut}.pdf")


    return smoothed_hist


def unroll(hist, rebin=1):
    if "TH1" in hist.ClassName():
        hist = hist.Rebin(rebin)
        hist = removeNegativeBins(hist)
        return hist

    elif "TH2" in hist.ClassName():
        # Get the number of bins in X and Y
        n_bins_x = hist.GetNbinsX()
        n_bins_y = hist.GetNbinsY()

        # Create a 1D histogram to hold the unrolled data
        n_bins_1d = n_bins_x * n_bins_y
        h1 = ROOT.TH1D("h1", "1D Unrolled Histogram", n_bins_1d, 0.5, n_bins_1d + 0.5)

        # Loop over all bins in the 2D histogram
        for bin_x in range(1, n_bins_x + 1):  # Bin indexing starts at 1
            for bin_y in range(1, n_bins_y + 1):
                # Calculate the global bin number for the 1D histogram
                bin_1d = (bin_y - 1) * n_bins_x + bin_x

                # Get the content and error from the 2D histogram
                content = hist.GetBinContent(bin_x, bin_y)
                error = hist.GetBinError(bin_x, bin_y)

                # Set the content and error in the 1D histogram
                h1.SetBinContent(bin_1d, content)
                h1.SetBinError(bin_1d, error)
                #print(hist.GetXaxis().GetBinCenter(bin_x), hist.GetYaxis().GetBinCenter(bin_y), content, error)
        return h1


    elif "TH3" in hist.ClassName():
        # Get binning information
        nbinsX = hist.GetNbinsX()
        nbinsY = hist.GetNbinsY()
        nbinsZ = hist.GetNbinsZ()
        nbins1D = nbinsX * nbinsY * nbinsZ

        # Create a 1D histogram with the correct number of bins
        h1 = ROOT.TH1D("h1_unrolled", "Unrolled 3D Histogram", nbins1D, 0, nbins1D)

        # Fill the 1D histogram by unrolling the 3D histogram
        bin1D = 1  # ROOT bins are 1-based
        for x in range(1, nbinsX + 1):
            for y in range(1, nbinsY + 1):
                for z in range(1, nbinsZ + 1):
                    content = hist.GetBinContent(x, y, z)
                    error = hist.GetBinError(x, y, z)  # Retrieve bin error
                    
                    if content < 0:
                        print("WARNING NEGATIVE CONTENT", content, hist.GetName())
                        content = 0
                        error = 0

                    h1.SetBinContent(bin1D, content)
                    h1.SetBinError(bin1D, error)  # Set the bin error
                    bin1D += 1  # Increment the 1D bin index
                    
        return h1

    else:
        return hist

def range_hist(hist_original, x_min, x_max):

    # Get the bin numbers corresponding to the range
    bin_min = hist_original.FindBin(x_min)
    bin_max = hist_original.FindBin(x_max)
    

    hist_selected = ROOT.TH1D(hist_original.GetName()+"new", "", bin_max - bin_min, x_min, x_max)

    for bin in range(bin_min, bin_max + 1):
        new_bin = bin - bin_min + 1  # Adjust for new histogram bin indexing
        hist_selected.SetBinContent(new_bin, hist_original.GetBinContent(bin))
        hist_selected.SetBinError(new_bin, hist_original.GetBinError(bin))  # Preserve errors

    print(bin_min, bin_max, hist_original.Integral(), hist_selected.Integral())
    return hist_selected


def make_pseudodata(procs, target="bb", variation=1.0):

    yields = {}
    xsec_tot = 0 # total cross-section
    xsec_target = 0  # nominal cross-section of the target process
    xsec_rest = 0 # cross-section of the rest

    sigProcs = procs_cfg[procs[0]] # get all signal processes
    for h_decay in h_decays:
        xsec = 0.
        for z_decay in z_decays:
            proc = f'wzp6_ee_{z_decay}H_H{h_decay}_ecm{ecm}'
            if not proc in sigProcs:
                continue
            xsec += getMetaInfo(proc)
        xsec_tot += xsec
        if h_decay != target:
            xsec_rest += xsec
        else:
            xsec_target += xsec

    xsec_new = variation*xsec_tot
    xsec_delta = xsec_new - xsec_tot # difference in cross-section
    print(xsec_target)
    scale_target = (xsec_target + xsec_delta)/xsec_target
    scale_rest = (xsec_rest - xsec_delta)/xsec_rest

    print("xsec_tot=", xsec_tot)
    print("xsec_new=", xsec_new)
    print("xsec_delta=", xsec_delta)
    print("scale_target=", scale_target)
    print("scale_rest=", scale_rest)

    hist_pseudo = None # start with all backgrounds
    for proc in procs[1:]:
        h = getHists(hName, proc)
        if hist_pseudo == None:
            hist_pseudo = h
        else:
            hist_pseudo.Add(h)

    xsec_tot_new = 0
    for h_decay in h_decays:
        xsec = 0.
        hist = None
        for z_decay in z_decays:
            proc = f'wzp6_ee_{z_decay}H_H{h_decay}_ecm{ecm}'
            if not proc in sigProcs:
                continue
            xsec += getMetaInfo(proc)
            h = getSingleHist(hName, proc)
            if hist == None:
                hist = h
            else:
                hist.Add(h)
        #xsec = sum([getMetaInfo(f'wzp6_ee_{z_decay}H_H{h_decay}_ecm240') for z_decay in z_decays])
        #hist = getSingleHist(hName, proc)
        #hist = getHist(hName, [f'wzp6_ee_{z_decay}H_H{h_decay}_ecm240' for z_decay in z_decays], rebin=rebin, scale=lumi)
        if h_decay == target:
            scale = (xsec+xsec_delta)/xsec
            hist.Scale(scale_target)
            xsec_tot_new += xsec*scale_target
        else:
            #hist.Scale(scale_rest)
            xsec_tot_new += xsec*scale_rest
        hist_pseudo.Add(hist)
    
    print(xsec_tot, xsec_tot_new) ## must be equal
    return hist_pseudo

if __name__ == "__main__":

    ecm = args.ecm
    outDir = f"output/h_zh/combine/ecm{ecm}/"
    plot_dir = "/home/submit/jaeyserm/public_html/fccee/h_zh/combine_smoothing/"
    sigma = 1
    bkg_unc = 1.01

    z_decays = ["qq", "bb", "cc", "ss", "ee", "mumu", "nunu" , "tautau"]
    h_decays = ["bb", "cc", "gg", "ss", "mumu", "tautau", "ZZ", "WW", "Za", "aa", "inv"]

    bbb = " --binByBinStat" if args.bbb else ""

    procs_cfg_240 = {
        "ZH"        : [f'wzp6_ee_{x}H_H{y}_ecm240' for x in z_decays for y in h_decays],
        "ZqqH"      : [f'wzp6_ee_{x}H_H{y}_ecm240' for x in ["qq", "bb", "cc", "ss"] for y in h_decays],
        "ZmumuH"    : [f'wzp6_ee_{x}H_H{y}_ecm240' for x in ["mumu"] for y in h_decays],
        "ZeeH"      : [f'wzp6_ee_{x}H_H{y}_ecm240' for x in ["ee"] for y in h_decays],
        "WW"        : ['p8_ee_WW_ecm240', 'p8_ee_WW_mumu_ecm240', 'p8_ee_WW_ee_ecm240'], # p8_ee_WW_mumu_ecm240 p8_ee_WW_ecm240
        "ZZ"        : ['p8_ee_ZZ_ecm240'],
        "Zgamma"    : ['wz3p6_ee_tautau_ecm240', 'wz3p6_ee_mumu_ecm240', 'wz3p6_ee_ee_Mee_30_150_ecm240', 'wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240', 'wz3p6_ee_nunu_ecm240'],
        "Zqqgamma"  : ['wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240'],
        "Rare"      : ['wzp6_egamma_eZ_Zmumu_ecm240', 'wzp6_gammae_eZ_Zmumu_ecm240', 'wzp6_gaga_mumu_60_ecm240', 'wzp6_egamma_eZ_Zee_ecm240', 'wzp6_gammae_eZ_Zee_ecm240', 'wzp6_gaga_ee_60_ecm240', 'wzp6_gaga_tautau_60_ecm240', 'wzp6_ee_nuenueZ_ecm240'],
    }

    procs_cfg_365 = {
        "ZH"        : [f'wzp6_ee_{x}H_H{y}_ecm365' for x in z_decays for y in h_decays],
        "ZqqH"      : [f'wzp6_ee_{x}H_H{y}_ecm365' for x in ["qq", "bb", "cc", "ss"] for y in h_decays],
        "ZmumuH"    : [f'wzp6_ee_{x}H_H{y}_ecm365' for x in ["mumu"] for y in h_decays],
        "ZeeH"      : [f'wzp6_ee_{x}H_H{y}_ecm365' for x in ["ee"] for y in h_decays],
        "WW"        : ['p8_ee_WW_ecm365', 'p8_ee_WW_mumu_ecm365', 'p8_ee_WW_ee_ecm365'], # p8_ee_WW_mumu_ecm365 p8_ee_WW_ecm365
        "ZZ"        : ['p8_ee_ZZ_ecm365'],
        "Zgamma"    : ['wz3p6_ee_tautau_ecm365', 'wz3p6_ee_mumu_ecm365', 'wz3p6_ee_ee_Mee_30_150_ecm365', 'wz3p6_ee_uu_ecm365', 'wz3p6_ee_dd_ecm365', 'wz3p6_ee_cc_ecm365', 'wz3p6_ee_ss_ecm365', 'wz3p6_ee_bb_ecm365', 'wz3p6_ee_nunu_ecm365'],
        "Zqqgamma"  : ['wz3p6_ee_uu_ecm365', 'wz3p6_ee_dd_ecm365', 'wz3p6_ee_cc_ecm365', 'wz3p6_ee_ss_ecm365', 'wz3p6_ee_bb_ecm365'],
        "Rare"      : ['wzp6_egamma_eZ_Zmumu_ecm365', 'wzp6_gammae_eZ_Zmumu_ecm365', 'wzp6_gaga_mumu_60_ecm365', 'wzp6_egamma_eZ_Zee_ecm365', 'wzp6_gammae_eZ_Zee_ecm365', 'wzp6_gaga_ee_60_ecm365', 'wzp6_gaga_tautau_60_ecm365', 'wzp6_ee_nuenueZ_ecm365'],
    }

    if ecm == 240:
        procs_cfg = procs_cfg_240
    if ecm == 365:
        procs_cfg = procs_cfg_365

    procs_scales_250 = {
        "ZH": 1.048,
        "WW": 0.971,
        "ZZ": 0.939,
        "Zgamma": 0.919,
    }
    
    procs_scales_polL = {
        "ZH": 1.554,
        "WW": 2.166,
        "ZZ": 1.330,
        "Zgamma": 1.263,
    }
    
    procs_scales_polR = {
        "ZH": 1.047,
        "WW": 0.219,
        "ZZ": 1.011,
        "Zgamma": 1.018,
    }
    
    proc_scales = procs_scales_250 ## change fit to ASIMOV -t -1 !!!
    proc_scales = {}



    cats = args.cats.split('-')
    p = -1 if args.floatBackgrounds else 1
    for cat in cats:
        hists = []
        if cat == "qq":
            inputDir = f"output/h_zh_hadronic/histmaker/ecm{ecm}/"
            hName = 'zqq_recoil_m_mqq_mva' # zqq_recoil_m zqq_recoil_m_mqq mva_score
            #hName = 'zqq_recoil_m_mqq' # zqq_recoil_m_mqq zqq_recoil_m
            procs = ["ZH", "WW", "ZZ", "Zgamma", "Rare"] # first must be signal
            proc_idx = [0, p*1, p*2, p*3, p*4]

            for proc in procs:
                h = getHists(hName, proc)
                h = unroll(h)
                h.SetName(f"{cat}_{proc}")
                hists.append(h)

            hist_pseudo = make_pseudodata(procs, target=args.target, variation=args.pert)
            hist_pseudo = unroll(hist_pseudo)
            hist_pseudo.SetName(f"{cat}_data")
            hists.append(hist_pseudo)

        if cat == "ee" or cat == "mumu":
            #procs_cfg["ZH"] = [f'wzp6_ee_{x}H_H{y}_ecm240' for x in [cat] for y in h_decays]

            # WITH COS THETA MISS
            #inputDir = "output/h_zh_leptonic/histmaker/ecm240/WithCosThetaMiss/"
            #hName, rebin = f'{cat}_zll_recoil', 50 # recoil 0.5 GeV bins

            ## FIT MVA SCORE (100 bins)
            #inputDir = f"output/h_zh_leptonic/histmaker/ecm{ecm}/"
            #hName, rebin = f'{cat}_mva_score', 10 # MVA score 100 bins

            # FIT RECOIL
            #inputDir = f"output/h_zh_leptonic/histmaker/ecm{ecm}/"
            #hName, rebin = f'{cat}_zll_recoil', 50 # recoil 0.5 GeV bins

            ## NOMINAL 2D MVA
            inputDir = f"output/h_zh_leptonic/histmaker/ecm{ecm}/"
            hName, rebin = f'{cat}_recoil_m_mva', 1

            procs = ["ZH", "WW", "ZZ", "Zgamma", "Rare"] # first must be signal
            proc_idx = [0, p*1, p*2, p*3, p*4]

            for proc in procs:
                h = getHists(hName, proc)
                h = unroll(h, rebin=rebin)
                h.SetName(f"{cat}_{proc}")
                hists.append(h)

            hist_pseudo = make_pseudodata(procs, target=args.target, variation=args.pert)
            hist_pseudo = unroll(hist_pseudo, rebin=rebin)
            hist_pseudo.SetName(f"{cat}_data")
            hists.append(hist_pseudo)


        fOut = ROOT.TFile(f"{outDir}/datacard_{cat}.root", "RECREATE")
        for hist in hists:
            hist.Write()
        fOut.Close()


        dc = ""
        dc += "imax *\n"
        dc += "jmax *\n"
        dc += "kmax *\n"
        dc += "####################\n"
        dc += f"shapes *        * datacard_{cat}.root $CHANNEL_$PROCESS\n"
        dc += f"shapes data_obs * datacard_{cat}.root $CHANNEL_data\n"
        dc += "####################\n"
        dc += f"bin          {cat}\n"
        dc += "observation  -1\n"
        dc += "####################\n"
        dc += "bin          {}\n".format('\t'.join([cat]*len(procs)))
        dc += "process      {}\n".format('\t'.join(procs))
        dc += "process      {}\n".format('\t'.join([str(idx) for idx in proc_idx]))
        dc += "rate         {}\n".format('\t'.join(['-1']*len(procs)))
        dc += "####################\n"
        if not args.freezeBackgrounds and not args.floatBackgrounds:
            for i,proc in enumerate(procs):
                if i == 0: continue # no signal
                str_before = '\t'.join(['-']*(i-1))
                str_after = '\t'.join(['-']*(len(procs)-i-1))
                dc += f"norm_{proc} lnN - {str_before}\t {bkg_unc} \t {str_after}\n"
        else:
            dc += "dummy lnN    1.000000005 {}\n".format('\t'.join(['-']*(len(procs)-1)))
        #dc += "* autoMCStats 0 0 1\n"

        f = open(f"{outDir}/datacard_{cat}.txt", 'w')
        f.write(dc)
        f.close()

        print(dc)

        if args.run and len(cats) == 1:
            cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; text2hdf5.py --X-allow-no-background datacard_{cat}.txt -o ws_{cat}.hdf5; combinetf.py ws_{cat}.hdf5 -o fit_output_{cat}.root -t 0  --expectSignal=1 {bbb}'"
            os.system(cmd)

    if len(cats) > 1 and args.run:
        cards = ' '.join([f'datacard_{x}.txt' for x in cats])
        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; combineCards.py {cards} > datacard_combined.txt'"
        #print(cmd)
        os.system(cmd)
        
        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; text2hdf5.py --X-allow-no-background datacard_combined.txt -o ws_combined.hdf5; combinetf.py ws_combined.hdf5 -o fit_output_combined.root -t 0  --expectSignal=1  {bbb}'" # --binByBinStat
        os.system(cmd)

