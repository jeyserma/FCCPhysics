import uproot
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score
import ROOT
import pickle
import os


ROOT.gROOT.SetBatch(True)
# e.g. https://root.cern/doc/master/tmva101__Training_8py.html

def load_process(proc, variables, target=0):
    weight = 1
    if ":" in proc:
        tmp = proc.split(":")
        proc = tmp[0]
        weight = float(tmp[1])

    fIn = f"{inputDir}/{proc}.root"

    f = uproot.open(fIn)
    print(f.keys())
    tree = f["events"]
    xsec = f["crossSection"].value
    ev_proc = f["eventsProcessed"].value
    weight *= xsec*3e6/ev_proc
    weight = 1.0
    #weight = 0.2*10.8e6 / ev_proc # 200 fb-1 for each Higgs decay == equal

    
    print(f"Load {fIn.replace('.root', '')} with {tree.num_entries} events and weight {weight} and cross-section {xsec}")

    df = tree.arrays(variables, library="pd") # convert the signal and background data to pandas DataFrames
    #df = df.head(int(0.1*ev_proc))
    df['target'] = target # add a target column to indicate signal (1) and background (0)
    df['weight'] = weight
    return df





################################################################################

tag = "zz_ww"

inputDir = "/ceph/submit/data/group/fcc/ee/analyses/h_zz_ww/treemaker/ecm240/stage2/"


variables = ["W1_WW_m", "W2_WW_m", "Z_WW_m", "Z_WW_p", "Z_WW_theta", "H_WW_m", "H_WW_p", "H_WW_theta", "Z1_ZZ_m", "Z2_ZZ_m", "Z_ZZ_m", "Z_ZZ_p", "Z_ZZ_theta", "H_ZZ_m", "H_ZZ_p", "H_ZZ_theta"]

# "sum_score_B_V1", "sum_score_B_V2", "sum_score_C_V1", "sum_score_C_V2", "sum_score_S_V1", "sum_score_S_V2", "sum_score_U_V1", "sum_score_U_V2", "sum_score_D_V1", "sum_score_D_V2", "sum_score_TAU_V1", "sum_score_TAU_V2"]
for i in range(0, 6):
    variables.append(f"jets_px_{i}")
    variables.append(f"jets_py_{i}")
    variables.append(f"jets_pz_{i}")
    variables.append(f"jets_phi_{i}")
    variables.append(f"jets_theta_{i}")
    variables.append(f"recojet_isB_{i}")
    variables.append(f"recojet_isC_{i}")
    variables.append(f"recojet_isS_{i}")
    variables.append(f"recojet_isU_{i}")
    variables.append(f"recojet_isD_{i}")
    variables.append(f"recojet_isTAU_{i}")

vars_kinematics = ["W1_WW_m", "W2_WW_m", "Z_WW_m", "Z_WW_p", "Z_WW_theta", "H_WW_m", "H_WW_p", "H_WW_theta", "Z1_ZZ_m", "Z2_ZZ_m", "Z_ZZ_m", "Z_ZZ_p", "Z_ZZ_theta", "H_ZZ_m", "H_ZZ_p", "H_ZZ_theta"]

# remove variables that should be similar between WW and ZZ
vars_kinematics = ["W1_WW_m", "W2_WW_m", "Z_WW_m", "Z1_ZZ_m", "Z2_ZZ_m", "Z_ZZ_m"]
            
vars_jets = ["W1_WW_jet1_p", "W1_WW_jet2_p", "W2_WW_jet1_p", "W2_WW_jet2_p", "W1_WW_jet1_theta", "W1_WW_jet2_theta", "W2_WW_jet1_theta", "W2_WW_jet2_theta", "Z1_ZZ_jet1_p", "Z1_ZZ_jet2_p", "Z2_ZZ_jet1_p", "Z2_ZZ_jet2_p", "Z1_ZZ_jet1_theta", "Z1_ZZ_jet2_theta", "Z2_ZZ_jet1_theta", "Z2_ZZ_jet2_theta"]

vars_flavor = ["W1_WW_jet1_isB", "W1_WW_jet2_isB", "W1_WW_jet1_isC", "W1_WW_jet2_isC", "W1_WW_jet1_isS", "W1_WW_jet2_isS", "W1_WW_jet1_isU", "W1_WW_jet2_isU", "W1_WW_jet1_isD", "W1_WW_jet2_isD", "W1_WW_jet1_isTAU", "W1_WW_jet2_isTAU", "W2_WW_jet1_isB", "W2_WW_jet2_isB", "W2_WW_jet1_isC", "W2_WW_jet2_isC", "W2_WW_jet1_isS", "W2_WW_jet2_isS", "W2_WW_jet1_isU", "W2_WW_jet2_isU", "W2_WW_jet1_isD", "W2_WW_jet2_isD", "W2_WW_jet1_isTAU", "W2_WW_jet2_isTAU", "Z1_ZZ_jet1_isB", "Z1_ZZ_jet2_isB", "Z1_ZZ_jet1_isC", "Z1_ZZ_jet2_isC", "Z1_ZZ_jet1_isS", "Z1_ZZ_jet2_isS", "Z1_ZZ_jet1_isU", "Z1_ZZ_jet2_isU", "Z1_ZZ_jet1_isD", "Z1_ZZ_jet2_isD", "Z1_ZZ_jet1_isTAU", "Z1_ZZ_jet2_isTAU", "Z2_ZZ_jet1_isB", "Z2_ZZ_jet2_isB", "Z2_ZZ_jet1_isC", "Z2_ZZ_jet2_isC", "Z2_ZZ_jet1_isS", "Z2_ZZ_jet2_isS", "Z2_ZZ_jet1_isU", "Z2_ZZ_jet2_isU", "Z2_ZZ_jet1_isD", "Z2_ZZ_jet2_isD", "Z2_ZZ_jet1_isTAU", "Z2_ZZ_jet2_isTAU"]

variables = vars_kinematics + vars_jets + vars_flavor

sig_procs = ['wzp6_ee_ccH_HWW_ecm240', 'wzp6_ee_bbH_HWW_ecm240', 'wzp6_ee_qqH_HWW_ecm240', 'wzp6_ee_ssH_HWW_ecm240']
bkg_procs = ['wzp6_ee_ccH_HZZ_ecm240', 'wzp6_ee_bbH_HZZ_ecm240', 'wzp6_ee_qqH_HZZ_ecm240', 'wzp6_ee_ssH_HZZ_ecm240']



parms = {
    'objective': 'binary:logistic',
    'eval_metric': 'logloss',
    'n_estimators': 2000,
    'max_depth': 6,
    #'eval_metric': ["merror", "mlogloss"],
    'early_stopping_rounds': 5,
}


parms = {
    'objective': 'binary:logistic',
    'eval_metric': 'logloss',
    'n_estimators': 350,
    'max_depth': 6, # 10 8
    #'eval_metric': ["merror", "mlogloss"],
    'early_stopping_rounds': 2,
}



################################################################################


print("Parse inputs")

dfs = []





for sig in sig_procs:
    df = load_process(sig, variables, target=1)
    dfs.append(df)



for bkg in bkg_procs:
    df = load_process(bkg, variables)
    dfs.append(df)



# Concatenate the dataframes into a single dataframe
data = pd.concat(dfs, ignore_index=True)


# split data in train/test events
train_data, test_data, train_labels, test_labels, train_weights, test_weights  = train_test_split(
    data[variables], data['target'], data['weight'], test_size=0.3, random_state=42
)



# conversion to numpy needed to have default feature_names (fN), needed for conversion to TMVA
train_data = train_data.to_numpy()
test_data = test_data.to_numpy()
train_labels = train_labels.to_numpy()
test_labels = test_labels.to_numpy()
train_weights = train_weights.to_numpy()
test_weights = test_weights.to_numpy()





# train the XGBoost model
print("Start training")
#eval_set = [(train_data, train_labels), (test_data, test_labels)]
eval_set = [(train_data, train_labels), (test_data, test_labels)]
bdt = xgb.XGBClassifier(**parms)
#bdt.fit(train_data, train_labels, verbose=True, eval_set=eval_set, sample_weight=train_weights)
bdt.fit(train_data, train_labels, verbose=True, eval_set=eval_set, sample_weight_eval_set=[train_weights, test_weights])



# export model (to ROOT and pkl)
print("Export model")
fOutName = f"/ceph/submit/data/group/fcc/ee/analyses/h_zz_ww/training/{tag}.root"
ROOT.TMVA.Experimental.SaveXGBoost(bdt, tag, fOutName, num_inputs=len(variables))

# append the variables
variables_ = ROOT.TList()
for var in variables:
     variables_.Add(ROOT.TObjString(var))
fOut = ROOT.TFile(fOutName, "UPDATE")
fOut.WriteObject(variables_, "variables")


save = {}
save['model'] = bdt
save['train_data'] = train_data
save['test_data'] = test_data
save['train_labels'] = train_labels
save['test_labels'] = test_labels
save['variables'] = variables
pickle.dump(save, open(f"{fOutName.replace('.root', '.pkl')}", "wb"))

print(tag)
print(f"Save to {fOutName}")

os.system(f"python FCCPhysics/analyses/h_zz_ww/training/evaluate.py --tag {tag}")