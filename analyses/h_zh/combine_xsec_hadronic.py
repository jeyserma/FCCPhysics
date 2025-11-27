import sys,array,ROOT,math,os,copy
import argparse
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

parser = argparse.ArgumentParser()
parser.add_argument("--cat", type=str, help="Category (mumu, ee)", default="mumu")
parser.add_argument("--run", help="Run combine", action='store_true')
parser.add_argument("--plot", help="Plot", action='store_true')
args = parser.parse_args()

def getMetaInfo(proc):
    if "Hinv" in proc:
        proc = proc.replace("wzp6", "wz3p6")
    fIn = ROOT.TFile(f"{inputDir}/{proc}.root")
    xsec = fIn.Get("crossSection").GetVal()
    
    #if "HZZ" in proc:
    #    xsec_inv = getMetaInfo(proc.replace("wzp6", "wz3p6").replace("HZZ", "Hinv"))
    #    print("REMOVE INV FROM ZZ XSEC", proc, xsec, xsec-xsec_inv)
    #    xsec = xsec - xsec_inv
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

def getHist(hName, procs, rebin=1, scale=1, scales=[]):
    hist = None
    for k,proc in enumerate(procs):
        if "Hinv" in proc:
            proc = proc.replace("wzp6", "wz3p6")
        fIn = ROOT.TFile(f"{inputDir}/{proc}.root")
        h = fIn.Get(hName)
        if len(scales) > 0:
            h.Scale(scales[k])
        h.SetDirectory(0)
        #if "HZZ" in proc:
        #    print("REMOVE INV FROM ZZ", proc)
        #    h_inv = getHist(hName, [proc.replace("wzp6", "wz3p6").replace("HZZ", "Hinv")], scales=scales)
        #    h.Add(h_inv, -1)
        h = removeNegativeBins(h) ## need to remove potential negative bins here at single process level (otherwise issues with pseudodata etc.)
        if hist == None:
            hist = h
        else:
            hist.Add(h)
        fIn.Close()
    #hist.Rebin(rebin)
    hist.Scale(scale)
    
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


def unroll_2d(hist):
    if "TH1" in hist.ClassName():
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


def make_pseudodata(target="bb", variation=1.0):

    xsec_dict = {}
    xsec_dict['bb'] = 0.03107+0.00394+0.003932+0.01745+0.0269+0.01745+0.01359+0.004171
    xsec_dict['cc'] = 0.001542
    xsec_dict['ss'] = 1.28e-05
    xsec_dict['gg'] = 0.004367
    xsec_dict['mumu'] = 1.161e-05
    xsec_dict['tautau'] = 0.003346
    xsec_dict['WW'] = 0.01148
    xsec_dict['ZZ'] = 0.001409
    xsec_dict['Za'] = 8.177e-05
    xsec_dict['aa'] = 0.0001211
    xsec_dict['inv'] = 5.8684e-05 # is in ZZ
    
    yields = {}
    xsec_tot = 0 # total cross-section
    xsec_target = 0  # nominal cross-section of the target process
    xsec_rest = 0 # cross-section of the rest
    
    for h_decay in higgs_decays:
        xsec = sum([getMetaInfo(f'wzp6_ee_{z_decay}H_H{h_decay}_ecm240') for z_decay in z_decays])
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
    
    #yield_tot_orig = sum(yields.values())
    #yield_tot_new = yield_tot_orig - yields[target] + variation*yields[target]
    #yield_others = yield_tot_orig - yields[target]
    #delta_yield = yield_tot_new - yield_tot_orig
    #yield_others_new = yield_others - delta_yield
    
    #variation_others = yield_others_new / yield_others

    
    
    #print("yield_tot_orig", yield_tot_orig)
    #print("yield_tot_new", yield_tot_new)
    
    #print("variation_others", variation_others, "ddd")

    #xsec_norm = sum(xsec_dict.values()) - xsec_dict[target]
    #delta = variation*xsec_dict[target] - xsec_dict[target] # difference in cross-section w.r.t nominal (signed)
    #print(delta, xsec_norm)
    hist_pseudo = getHist(hName, bkgs, rebin=rebin, scale=lumi)
    #hist_pseudo.Reset("ACE")
    #ss = 0
    #sss = 0

    xsec_tot_new = 0
    for h_decay in higgs_decays:
        xsec = sum([getMetaInfo(f'wzp6_ee_{z_decay}H_H{h_decay}_ecm240') for z_decay in z_decays])
        hist = getHist(hName, [f'wzp6_ee_{z_decay}H_H{h_decay}_ecm240' for z_decay in z_decays], rebin=rebin, scale=lumi)

        if h_decay == target:
            scale = (xsec+xsec_delta)/xsec
            hist.Scale(scale_target)
            xsec_tot_new += xsec*scale_target
        else:
            #hist.Scale(scale_rest)
            xsec_tot_new += xsec*scale_rest
        hist_pseudo.Add(hist)
    
    print(xsec_tot, xsec_tot_new) ## must be equal
    
    '''
    for sig in sigs:
        hist = getHist(hName, [sig], rebin=rebin, scale=lumi)
        z_decay =sig.split("_")[2][:-1]
        h_decay = sig.split("_")[3][1:]
        
        if z_decay in z_decays:
            #print(sig, z_decay, h_decay)
            xsec = getMetaInfo(sig)
            if h_decay == target:
                #hist.Scale(variation)
                var = (xsec+xsec_delta)/xsec
                hist.Scale(var)
                print("Scale", sig, var)
            else:
                pass
                #delta_s = xsec_dict[h_decay] / xsec_norm * delta
                #ss += delta_s
                #sss += xsec_dict[h_decay]
                #print(delta_s)
                #v = 1/variation * xsec_dict[h_decay] / xsec_norm
                ####hist.Scale(variation_others)
                #print(f"SCALE {sig} with {v}")
                hist.Scale(scale_rest)
        hist_pseudo.Add(hist)
    #print(hist_pseudo.Integral())
    #print(ss, sss)
    #quit()
    '''
    return hist_pseudo

if __name__ == "__main__":

    inputDir = "output/h_zh_hadronic/histmaker/" # histmaker_newcostheta_ilcsel histmaker_newcostheta
    outDir = "output/h_zh_hadronic/combine/"
    plot_dir = "/home/submit/jaeyserm/public_html/fccee/h_zh/combine_smoothing/"
    rebin = 1
    sigma = 1
    lumi = 1

    z_decays = ["qq", "bb", "cc", "ss", "ee", "mumu", "tautau", "nunu"]
    z_decays = ["qq", "bb", "cc", "ss"]
    higgs_decays = ["bb", "cc", "gg", "ss", "mumu", "tautau", "ZZ", "WW", "Za", "aa", "inv"]

    sigs = [f'wzp6_ee_{x}H_H{y}_ecm240' for x in z_decays for y in higgs_decays]
    bkgs = ['p8_ee_WW_ecm240', 'p8_ee_ZZ_ecm240', 'p8_ee_Zqq_ecm240']
    #sigs_scales, bkgs_scales = [1.554], [2.166, 1.330, 1.263, 1.263] # 250 + polL
    #sigs_scales, bkgs_scales = [1.047], [0.219, 1.011, 1.018, 1.018] # 250 + polR
    #sigs_scales, bkgs_scales = [1.048], [0.971, 0.939, 0.919, 0.919] # 250
    sigs_scales, bkgs_scales = [1.0]*len(sigs), [1.0, 1.0, 1.0, 1.0] # 240

    hName = 'zqq_recoil_m_mqq_mva' # zqq_recoil_m zqq_recoil_m_mqq mva_score
    cat = args.cat

    h_sig = getHist(hName, sigs, rebin=rebin, scale=lumi, scales=sigs_scales)
    h_bkg = getHist(hName, bkgs, rebin=rebin, scale=lumi, scales=bkgs_scales)
    
    h_ww = getHist(hName, [bkgs[0]], rebin=rebin, scale=lumi, scales=bkgs_scales)
    h_zz = getHist(hName, [bkgs[1]], rebin=rebin, scale=lumi, scales=bkgs_scales)
    h_zg = getHist(hName, [bkgs[2]], rebin=rebin, scale=lumi, scales=bkgs_scales)

    #h_sig = range_hist(h_sig, 80, 180) # zqq_recoil_m_best_nosel
    #h_bkg = range_hist(h_bkg, 80, 180)

    h_sig = unroll_2d(h_sig)
    h_bkg = unroll_2d(h_bkg)
    h_ww = unroll_2d(h_ww)
    h_zz = unroll_2d(h_zz)
    h_zg = unroll_2d(h_zg)
    #h_sig = smooth_histogram_gaussian_kernel(h_sig, sigma, f"{cat}_sig")
    #h_bkg = smooth_histogram_gaussian_kernel(h_bkg, sigma, f"{cat}_bkg")

    h_sig.SetName(f"{cat}_sig")
    h_bkg.SetName(f"{cat}_bkg")
    h_ww.SetName(f"{cat}_ww")
    h_zz.SetName(f"{cat}_zz")
    h_zg.SetName(f"{cat}_zg")
    
    hist_pseudo = make_pseudodata(target="inv", variation=1.05)
    hist_pseudo = unroll_2d(hist_pseudo)
    hist_pseudo.SetName(f"{cat}_data")

    fOut = ROOT.TFile(f"{outDir}/datacard_{cat}.root", "RECREATE")
    h_sig.Write()
    h_bkg.Write()
    hist_pseudo.Write()
    #h_data = h_sig.Clone(f"{cat}_data")
    #h_data.Add(h_bkg)
    #h_data.Write()
    h_ww.Write()
    h_zz.Write()
    h_zg.Write()

    fOut.Close()

    '''
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
    dc += f"bin          {cat}      {cat}\n"
    dc += "process      sig         bkg\n"
    dc += "process      0           1\n"
    dc += "rate         -1          -1\n"
    dc += "####################\n"
    dc += "dummy lnN    -     1.000000005\n"
    #dc += "bkg lnN      -           1.01\n"
    #dc += "* autoMCStats 0 0 1\n"
    '''
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
    dc += f"bin          {cat}      {cat}   {cat}   {cat}\n"
    dc += "process      sig         ww      zz      zg\n"
    dc += "process      0           -1       -2       -3 \n"
    dc += "rate         -1          -1      -1      -1\n"
    dc += "####################\n"
    dc += "dummy lnN    -     1.000000005 - -\n"
    #dc += "bkg lnN      -           1.0001 - -\n"
    #dc += "bkg lnN      -           - 1.0001 -\n"
    #dc += "bkg lnN      -           - - 1.0001\n"
    #dc += "bkg lnN      -           1.01\n"
    #dc += "* autoMCStats 0 0 1\n"

    f = open(f"{outDir}/datacard_{cat}.txt", 'w')
    f.write(dc)
    f.close()

    print(dc)

    if args.run:
        #cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd {outDir}; text2workspace.py datacard_{cat}.txt -o ws_{cat}.root; combine -M MultiDimFit -v 10 --rMin 0.5 --rMax 1.5 --setParameters r=1 ws_{cat}.root'"
        #cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd {outDir}; text2workspace.py datacard_{cat}.txt -o ws_{cat}.root; combine -v 10 -M FitDiagnostics -t 0 --setParameters r=1 ws_{cat}.root -n {cat} --rMin 0.5 --rMax 1.5 --cminDefaultMinimizerStrategy 0'"
        
        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; which python; text2hdf5.py --X-allow-no-background datacard_{cat}.txt -o ws_{cat}.hdf5; combinetf.py ws_{cat}.hdf5 -o fit_output_{cat}.root -t 0  --expectSignal=1'" # --binByBinStat
        
        os.system(cmd)
        
        '''
        fIn = ROOT.TFile(f"{outDir}/fitDiagnostics{cat}.root")
        t = fIn.Get("tree_fit_sb")
        t.GetEntry(0)
        err = t.rErr*100.
        status = t.fit_status
        print(f"{cat}\t{err:.3f} {status}")
        '''
        
    '''
    if args.run:
        cards = ' '.join([f'datacard_{x}.txt' for x in hists])
        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd {outDir}; combineCards.py {cards} > datacard_combined.txt'"
        os.system(cmd)

        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd {outDir}; text2workspace.py datacard_combined.txt -o ws_combined.root; combine -M FitDiagnostics -t -1 --setParameters r=1 ws_combined.root -n combined --cminDefaultMinimizerStrategy 0'"
        os.system(cmd)

        ## extract results
        for hist in hists + ['combined']:
            fIn = ROOT.TFile(f"{outDir}/fitDiagnostics{hist}.root")
            t = fIn.Get("tree_fit_sb")
            t.GetEntry(0)
            err = t.rErr
            status = t.fit_status
            print(f"{hist}\t{err:.3f} {status}")

    '''


