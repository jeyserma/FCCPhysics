
import ROOT

fIn = "FCCPhysics/analyses/h_zh/trainings/xgb_bdt_ee.root"
fOutName = "FCCPhysics/analyses/h_zh/trainings/xgb_bdt_ee_converted.root"

variables = ["zll_leading_p", "zll_leading_theta", "zll_subleading_p", "zll_subleading_theta", "acolinearity", "acoplanarity", "zll_m", "zll_p", "zll_theta"]

# append the variables
variables_ = ROOT.TList()
for var in variables:
     variables_.Add(ROOT.TObjString(var))
fOut = ROOT.TFile(fOutName, "UPDATE")
fOut.WriteObject(variables_, "variables")
