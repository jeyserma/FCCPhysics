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
    return hist




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



if __name__ == "__main__":

    ecm = args.ecm
    outDir = "output/h_zz_ww/combine/"
    bkg_unc = 1.01

    bbb = " --binByBinStat" if args.bbb else ""

    procs_cfg_240 = {
        "ZHWW"        : ['wzp6_ee_ccH_HWW_ecm240', 'wzp6_ee_bbH_HWW_ecm240', 'wzp6_ee_qqH_HWW_ecm240', 'wzp6_ee_ssH_HWW_ecm240', 'wzp6_ee_eeH_HWW_ecm240', 'wzp6_ee_tautauH_HWW_ecm240', 'wzp6_ee_mumuH_HWW_ecm240'],
        "ZHZZ"        : ['wzp6_ee_ccH_HZZ_ecm240', 'wzp6_ee_bbH_HZZ_ecm240', 'wzp6_ee_qqH_HZZ_ecm240', 'wzp6_ee_ssH_HZZ_ecm240', 'wzp6_ee_nunuH_HZZ_ecm240', 'wzp6_ee_eeH_HZZ_ecm240', 'wzp6_ee_tautauH_HZZ_ecm240', 'wzp6_ee_mumuH_HZZ_ecm240'],
        "ZqqHWW"        : ['wzp6_ee_ccH_HWW_ecm240', 'wzp6_ee_bbH_HWW_ecm240', 'wzp6_ee_qqH_HWW_ecm240', 'wzp6_ee_ssH_HWW_ecm240'],
        "ZqqHZZ"        : ['wzp6_ee_ccH_HZZ_ecm240', 'wzp6_ee_bbH_HZZ_ecm240', 'wzp6_ee_qqH_HZZ_ecm240', 'wzp6_ee_ssH_HZZ_ecm240'],

        "WW"        : ['p8_ee_WW_ecm240'],
        "ZZ"        : ['p8_ee_ZZ_ecm240'],
        "Zgamma"    : ['wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240'],
    }


    if ecm == 240:
        procs_cfg = procs_cfg_240
    if ecm == 365:
        procs_cfg = procs_cfg_365


    inputDir = f"output/h_zz_ww/histmaker/ecm{ecm}/stage2/"
    procs = ["ZHWW", "ZHZZ", "WW", "ZZ", "Zgamma"]
    #procs = ["ZHZZ", "ZHWW", "WW", "ZZ", "Zgamma"]
    proc_idx = [0, -1, 2, 3, 4]
    hNames = ["H_WW_m_Z_WW_m"] # H_WW_m_Z_WW_m  H_WW_m   WW_H_WW_m_Z_WW_m  ZZ_H_ZZ_m_Z_ZZ_m
    hNames = ["WW_H_WW_m_Z_WW_m", "ZZ_H_ZZ_m_Z_ZZ_m"]
    hNames = ["mva_score_ww", "mva_score_zz"]
    #hNames = ["final_mva_ww", "final_mva_zz"]
    #hNames = ["mva_score"] # ww inclusive
    #hNames = ["mva_score_1d"] # ww inclusive
    for hName in hNames:
        hists = []
        for proc in procs:
            h = getHists(hName, proc)
            h = unroll(h, rebin=10)
            h.SetName(f"{hName}_{proc}")
            hists.append(h)

        h_data = hists[0].Clone(f"{hName}_data")
        hists.append(h_data)

        fOut = ROOT.TFile(f"{outDir}/datacard_{hName}.root", "RECREATE")
        for hist in hists:
            hist.Write()
        fOut.Close()


        dc = ""
        dc += "imax *\n"
        dc += "jmax *\n"
        dc += "kmax *\n"
        dc += "####################\n"
        dc += f"shapes *        * datacard_{hName}.root $CHANNEL_$PROCESS\n"
        dc += f"shapes data_obs * datacard_{hName}.root $CHANNEL_data\n"
        dc += "####################\n"
        dc += f"bin          {hName}\n"
        dc += "observation  -1\n"
        dc += "####################\n"
        dc += "bin          {}\n".format('\t'.join([hName]*len(procs)))
        dc += "process      {}\n".format('\t'.join(procs))
        dc += "process      {}\n".format('\t'.join([str(idx) for idx in proc_idx]))
        dc += "rate         {}\n".format('\t'.join(['-1']*len(procs)))
        dc += "####################\n"
        if not args.freezeBackgrounds and not args.floatBackgrounds:
            for i,proc in enumerate(procs):
                if proc_idx[i] <= 0: continue # no signal
                str_before = '\t'.join(['-']*(i-1))
                str_after = '\t'.join(['-']*(len(procs)-i-1))
                dc += f"norm_{proc} lnN - {str_before}\t {bkg_unc} \t {str_after}\n"
        else:
            dc += "dummy lnN    1.000000005 {}\n".format('\t'.join(['-']*(len(procs)-1)))
        #dc += "* autoMCStats 0 0 1\n"

        f = open(f"{outDir}/datacard_{hName}.txt", 'w')
        f.write(dc)
        f.close()

        print(dc)

        if args.run and len(hNames) == 1:
            cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; text2hdf5.py --X-allow-no-background datacard_{hName}.txt -o ws_{hName}.hdf5; combinetf.py ws_{hName}.hdf5 -o fit_output_{hName}.root -t -1  --expectSignal=1 {bbb}'"
            os.system(cmd)



    if len(hNames) > 1 and args.run:
        cards = ' '.join([f'datacard_{x}.txt' for x in hNames])
        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; combineCards.py {cards} > datacard_combined.txt'"
        #print(cmd)
        os.system(cmd)
        
        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; text2hdf5.py --X-allow-no-background datacard_combined.txt -o ws_combined.hdf5; combinetf.py ws_combined.hdf5 -o fit_output_combined.root -t -1  --expectSignal=1  {bbb}'" # --binByBinStat
        os.system(cmd)

