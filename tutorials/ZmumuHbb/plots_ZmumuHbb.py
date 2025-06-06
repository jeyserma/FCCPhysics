import ROOT

# global parameters
intLumi        = 1. # assume histograms are scaled in previous step
intLumiLabel   = "L = 10.8 ab^{-1}"
ana_tex        = 'e^{+}e^{-} #rightarrow ZH #rightarrow #mu^{+}#mu^{-} + X'
energy         = 240.0
collider       = 'FCC-ee'
formats        = ['png','pdf']

outdir         = "output/ZmumuHbb/plots/"
inputDir       = "output/ZmumuHbb/histmaker/"
# You can also use the pre-generated files
# /ceph/submit/data/group/fcc/ee/tutorials/FNAL2025/ZmumuHbb/histmaker/ MIT
# /eos/project/f/fccsw-web/www/tutorials/FNAL2025/ZmumuHbb/histmaker/ CERN

plotStatUnc    = False
leg_position = [0.8, 0.7, 0.96, 0.88]

procs = ["", "", ""]

procs = {}
procs['signal'] = {'ZH':['wzp6_ee_mumuH_Hbb_ecm240', 'wzp6_ee_mumuH_Hcc_ecm240', 'wzp6_ee_mumuH_Hss_ecm240', 'wzp6_ee_mumuH_Hgg_ecm240', 'wzp6_ee_mumuH_Haa_ecm240', 'wzp6_ee_mumuH_HZa_ecm240', 'wzp6_ee_mumuH_HWW_ecm240', 'wzp6_ee_mumuH_HZZ_ecm240', 'wzp6_ee_mumuH_Hmumu_ecm240', 'wzp6_ee_mumuH_Htautau_ecm240']}
procs['backgrounds'] =  {
    'Zgamma':['wzp6_ee_mumu_ecm240', 'wzp6_ee_tautau_ecm240'],
    'WW':['p8_ee_WW_mumu_ecm240'],
    'ZZ':['p8_ee_ZZ_ecm240'],
    }

stacksig = False

colors = {}
colors['ZH'] = ROOT.kRed
colors['WW'] = ROOT.kBlue+1
colors['ZZ'] = ROOT.kGreen+2
colors['Zgamma'] = ROOT.kRed+2

legend = {}
legend['ZH'] = 'ZH'
legend['WW'] = 'WW'
legend['ZZ'] = 'ZZ'
legend['Zgamma'] = 'Z/#gamma*'


hists = {}



hists["zmumu_m"] = {
    "output":   "zmumu_m",
    "logy":     True,
    "stack":    False,
    "rebin":    1,
    "xmin":     0,
    "xmax":     200,
    "ymin":     1e1,
    "ymax":     1e7,
    "xtitle":   "m(#mu^{#plus}#mu^{#minus}) (GeV)",
    "ytitle":   "Events / 1 GeV",
}

hists["zmumu_p"] = {
    "output":   "zmumu_p",
    "logy":     True,
    "stack":    False,
    "rebin":    1,
    "xmin":     0,
    "xmax":     110,
    "ymin":     1e1,
    "ymax":     1e7,
    "xtitle":   "p(#mu^{#plus}#mu^{#minus}) (GeV)",
    "ytitle":   "Events / 1 GeV",
}

hists["zmumu_recoil_m"] = {
    "output":   "zmumu_recoil_m",
    "logy":     True,
    "stack":    False,
    "rebin":    1,
    "xmin":     100,
    "xmax":     150,
    "ymin":     1e1,
    "ymax":     1e7,
    "xtitle":   "Recoil mass (GeV)",
    "ytitle":   "Events / 1 GeV",
}


hists["dijet_m"] = {
    "output":   "dijet_m",
    "logy":     True,
    "stack":    False,
    "rebin":    1,
    "xmin":     100,
    "xmax":     150,
    "ymin":     1e1,
    "ymax":     1e7,
    "xtitle":   "Dijet mass (GeV)",
    "ytitle":   "Events / 1 GeV",
}







hists["recojet_isB"] = {
    "output":   "recojet_isB",
    "logy":     True,
    "stack":    False,
    "rebin":    1,
    "xmin":     0,
    "xmax":     1,
    "ymin":     1e1,
    "ymax":     1e7,
    "xtitle":   "B-jet probability",
    "ytitle":   "Events / 0.5 GeV",
}


hists["recojet_isC"] = {
    "output":   "recojet_isC",
    "logy":     True,
    "stack":    False,
    "rebin":    1,
    "xmin":     0,
    "xmax":     1,
    "ymin":     1e1,
    "ymax":     1e7,
    "xtitle":   "C-jet probability",
    "ytitle":   "Events / 0.5 GeV",
}


hists["zmumu_recoil_m_final"] = {
    "output":   "zmumu_recoil_m_final",
    "logy":     False,
    "stack":    True,
    "rebin":    1,
    "xmin":     120,
    "xmax":     140,
    "ymin":     0,
    "ymax":     2e4,
    "xtitle":   "Recoil mass (GeV)",
    "ytitle":   "Events / 0.5 GeV",
}


hists["cutFlow"] = {
    "output":   "cutFlow",
    "logy":     True,
    "stack":    False,
    "xmin":     0,
    "xmax":     7,
    "ymin":     1e3,
    "ymax":     1e10,
    "xtitle":   ["All events", "2 #mu^{#pm} + OS", "86 < m_{#mu^{+}#mu^{#minus}} < 96", "20 < p_{#mu^{+}#mu^{#minus}} < 70", "|cos#theta_{miss}| < 0.98", "120 < m_{rec} < 140", "B-tag score > 0.5"],
    "ytitle":   "Events ",
    "scaleSig": 100, # only for plotting
    "dumpTable": True
}