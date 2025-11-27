
import ROOT
import os

def run_fit(cats, target, pert, extraArgs):
    cmd = f"python FCCPhysics/analyses/h_zh/combine.py --cats {cats} --run --ecm {ecm}  --target {target} --pert {pert} {extraArgs}"
    os.system(cmd)

    if '-' in cats:
        fIn = ROOT.TFile(f"output/h_zh/combine/ecm{ecm}/fit_output_combined.root") # combination
    else:
        fIn = ROOT.TFile(f"output/h_zh/combine/ecm{ecm}/fit_output_{cats}.root") # single category
    tree = fIn.Get("fitresults")
    tree.GetEntry(0)
    status = tree.status
    errstatus = tree.errstatus
    ZH_mu = tree.ZH_mu
    ZH_mu_err = tree.ZH_mu_err

    return [status, errstatus, ZH_mu, ZH_mu_err]

if __name__ == "__main__":

    ecm = 240 # 240 365
    extraArgs = "--freezeBackgrounds"
    #extraArgs = "--floatBackgrounds"
    extraArgs = ""
    pert = 1.01
    cats = "qq-mumu-ee"
    h_decays = ["bb", "cc", "ss", "gg", "mumu", "tautau", "ZZ", "WW", "Za", "aa", "inv"]

    res = []
    for h_decay in h_decays:
        r = run_fit(cats, h_decay, pert, extraArgs)
        res.append(r)


    for i,h_decay in enumerate(h_decays):
        status = res[i][0]
        errstatus = res[i][0] = res[i][1]
        ZH_mu = res[i][2]
        ZH_mu_err = res[i][3]
        #print(status, errstatus, ZH_mu, ZH_mu_err, pert)
        bias = 100.*(ZH_mu - pert)

        if status != 0 and errstatus != 0:
            print(f"{h_decay}\tERROR: status={status} errstatus={errstatus}")
        else:
            print(f"{h_decay}\t{bias:.2f}")

