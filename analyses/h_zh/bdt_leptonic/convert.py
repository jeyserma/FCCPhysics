
import ROOT
import os

tag = "xgb_bdt_mumu_365"

fIn = f"FCCPhysics/analyses/h_zh/bdt_leptonic/{tag}.root"
fOutName = f"FCCPhysics/analyses/h_zh/bdt_leptonic/{tag}_converted.root"

os.system(f"cp {fIn} {fOutName}") # first make a copy of the original file

variables = ["zll_leading_p", "zll_leading_theta", "zll_subleading_p", "zll_subleading_theta", "acolinearity", "acoplanarity", "zll_m", "zll_p", "zll_theta"]

# append the variables
variables_ = ROOT.TList()
for var in variables:
     variables_.Add(ROOT.TObjString(var))
fOut = ROOT.TFile(fOutName, "UPDATE")
fOut.WriteObject(variables_, "variables")
