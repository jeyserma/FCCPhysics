
base_dir = "/ceph/submit/data/group/fcc/ee/beam_backgrounds/guineapig/ipc/"
base_dir_studies = "/ceph/submit/data/group/fcc/ee/beam_backgrounds/guineapig/ipc_studies/"
base_dir_ipc_sensitivity_studies = "/ceph/submit/data/group/fcc/ee/beam_backgrounds/guineapig/ipc_sensitivity_studies/"

ipc = {}


ipc['FCCee_Z_GHC_V25p1'] = {}


ipc['FCCee_Z_GHC_V25p1']['Z256_2T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV",
    "dir"       : f"{base_dir}/FCCee_Z_GHC_V25p1/Z256_2T_grids8/merged",
    "card"      : "cards/fccee.dat"
}


ipc['FCCee_Z_LCC_V105'] = {}
ipc['FCCee_Z_LCC_V105']['Z256_2T_grids8'] = {
    "label"     :"LCC V105, 91.2 GeV",
    "dir"       : f"{base_dir}/FCCee_Z_LCC_V105/Z256_2T_grids8/merged",
    "card"      : "cards/fccee.dat",
    "n_bunches" : 12000, # cross-check
}

ipc['FCCee_Z_LCC_V105_v2_25ns'] = {}
ipc['FCCee_Z_LCC_V105_v2_25ns']['Z256_2T_grids8'] = {
    "label"     :"LCC, 25 ns, 91.2 GeV",
    "dir"       : f"{base_dir}/FCCee_Z_LCC_V105_v2_25ns/Z256_2T_grids8/merged",
    "card"      : "cards/fccee.dat",
    "n_bunches" : 11200,
}

ipc['FCCee_Z_LCC_V105_v2_50ns'] = {}
ipc['FCCee_Z_LCC_V105_v2_50ns']['Z256_2T_grids8'] = {
    "label"     :"LCC, 50 ns, 91.2 GeV",
    "dir"       : f"{base_dir}/FCCee_Z_LCC_V105_v2_50ns/Z256_2T_grids8/merged",
    "card"      : "cards/fccee.dat",
    "n_bunches" : 5600,
}

ipc['FCCee_Z_CDR'] = {}
ipc['FCCee_Z_CDR']['Z256_2T_grids8'] = {
    "label"     :"CDR, 91.2 GeV",
    "dir"       : f"{base_dir}/FCCee_Z_CDR/Z256_2T_grids8/merged",
    "card"      : "cards/fccee.dat"
}



ipc['FCCee_Z_GHC_V23'] = {}
ipc['FCCee_Z_GHC_V23']['Z256_2T_grids8'] = {
    "label"     :"GHC V23 (MTR), 91.2 GeV",
    "dir"       : f"{base_dir}/FCCee_Z_GHC_V23/Z256_2T_grids8/merged",
    "card"      : "cards/fccee.dat"
}

ipc['FCCee_TOP_GHC_V25p1'] = {}
ipc['FCCee_TOP_GHC_V25p1']['Z256_2T_grids8'] = {
    "label"     :"GHC V25.1 (FSR), 365 GeV",
    "dir"       : f"{base_dir}/FCCee_TOP_GHC_V25p1/Z256_2T_grids8/merged",
    "card"      : "cards/fccee.dat"
}


ipc['FCCee_WW_GHC_V25p1'] = {}
ipc['FCCee_WW_GHC_V25p1']['Z256_2T_grids8'] = {
    "label"     :"GHC V25.1 (FSR), 160 GeV",
    "dir"       : f"{base_dir}/FCCee_WW_GHC_V25p1/Z256_2T_grids8/merged",
    "card"      : "cards/fccee.dat"
}


ipc['FCCee_ZH_GHC_V25p1'] = {}
ipc['FCCee_ZH_GHC_V25p1']['Z256_2T_grids8'] = {
    "label"     :"GHC V25.1 (FSR), 240 GeV",
    "dir"       : f"{base_dir}/FCCee_ZH_GHC_V25p1/Z256_2T_grids8/merged",
    "card"      : "cards/fccee.dat"
}


ipc['FCCee_Z_4IP_GHC_V24p4'] = {}
ipc['FCCee_Z_4IP_GHC_V24p4']['Z256_2T_grids8'] = {
    "label"     :"GHC V25.1 (FSR), 91.2 GeV",
    "dir"       : f"{base_dir}/FCCee_Z_4IP_GHC_V24p4/Z256_2T_grids8/merged",
    "card"      : "cards/fccee.dat"
}


ipc['FCCee_Z_GHC_V25p3_4'] = {}
ipc['FCCee_Z_GHC_V25p3_4']['Z256_2T_grids8'] = {
    "label"     :"GHC V25.3-4, 91.2 GeV",
    "dir"       : f"{base_dir}/FCCee_Z_GHC_V25p3_4/Z256_2T_grids8/merged",
    "card"      : "cards/fccee.dat"
}


############ STUDIES
ipc_studies = {}
ipc_studies['FCCee_Z_GHC_V25p1'] = {}

ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids1'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:1:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids1/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids2'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:2:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids2/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids3'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:3:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids3/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids4'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:4:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids4/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids5'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:5:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids5/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids6'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:6:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids6/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids7'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:7:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids7/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}

ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids1_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:1:0T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids1_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids2_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:2:0T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids2_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids3_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:3:0T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids3_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids4_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:4:0T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids4_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids5_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:5:0T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids5_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids6_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:6:0T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids6_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids7_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:7:0T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids7_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids8_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:0T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids8_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}

ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids1'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:1:2T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids1/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids2'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:2:2T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids2/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids3'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:3:2T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids3/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids4'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:4:2T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids4/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids5'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:5:2T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids5/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids6'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:6:2T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids6/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids7'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:7:2T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids7/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:2T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}

ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids1_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:1:2T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids1_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids2_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:2:2T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids2_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids3_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:3:2T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids3_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids4_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:4:2T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids4_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids5_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:5:2T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids5_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids6_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:6:2T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids6_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids7_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:7:2T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids7_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids8_nbf'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:2T:nbf)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids8_nbf/merged",
    "card"      : "cards/fccee_studies.dat"
}



ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0p1T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:0.1T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0p1T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0p2T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:0.2T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0p2T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0p3T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:0.3T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0p3T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0p4T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:0.4T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0p4T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0p5T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:0.5T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0p5T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0p6T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:0.6T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0p6T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0p7T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:0.7T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0p7T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0p75T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:0.75T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0p75T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0p8T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:0.8T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0p8T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0p9T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:0.9T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0p9T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_1T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:1T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_1T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_1p5T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:1.5T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_1p5T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_3T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:3T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_3T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z128_4T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (128:8:4T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_4T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}



ipc_studies['FCCee_Z_GHC_V25p1']['Z32_2T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (32:1:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z32_2T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z64_2T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (64:1:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z64_2T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z256_2T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (256:1:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z256_2T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z320_2T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (320:1:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z320_2T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}
ipc_studies['FCCee_Z_GHC_V25p1']['Z512_2T_grids8'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (320:1:0T)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z512_2T_grids8/merged",
    "card"      : "cards/fccee_studies.dat"
}


ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids8_theta0'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (320:1:0T:theta0)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids8_theta0/merged",
    "card"      : "cards/fccee_studies.dat"
}


ipc_studies['FCCee_Z_GHC_V25p1']['Z128_2T_grids8_noBoundary'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (320:1:2T:noBoundary)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_2T_grids8_noBoundary/merged",
    "card"      : "cards/fccee_studies.dat"
}



ipc_studies['FCCee_Z_GHC_V25p1']['Z128_0T_grids8_noBoundary'] = {
    "label"     : "GHC V25.1 (FSR) 91.2 GeV (320:1:0T:noBoundary)",
    "dir"       : f"{base_dir_studies}/FCCee_Z_GHC_V25p1/Z128_0T_grids8_noBoundary",
    "card"      : "cards/fccee_studies.dat"
}





# CAMPAIGN: sensitivity studies

grid = "Z256_2T_grids8"
card = "cards/fccee_GHC_V25p1_variations.dat"

ipc_sensitivity_studies = {
    "FCCee_Z_GHC_V25p1": {
        grid: {
            "label": "Nominal",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1/{grid}/merged",
            "card": card,
        }
    },

    "FCCee_Z_GHC_V25p1_sigmaX_10pctUp": {
        grid: {
            "label": "#sigma_{x} 10% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_sigmaX_10pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_sigmaX_10pctDown": {
        grid: {
            "label": "#sigma_{x} 10% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_sigmaX_10pctDown/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_sigmaX_25pctUp": {
        grid: {
            "label": "#sigma_{x} 25% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_sigmaX_25pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_sigmaX_25pctDown": {
        grid: {
            "label": "#sigma_{x} 25% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_sigmaX_25pctDown/{grid}/merged",
            "card": card,
        }
    },

    "FCCee_Z_GHC_V25p1_sigmaY_10pctUp": {
        grid: {
            "label": "#sigma_{y} 10% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_sigmaY_10pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_sigmaY_10pctDown": {
        grid: {
            "label": "#sigma_{y} 10% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_sigmaY_10pctDown/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_sigmaY_25pctUp": {
        grid: {
            "label": "#sigma_{y} 25% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_sigmaY_25pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_sigmaY_25pctDown": {
        grid: {
            "label": "#sigma_{y} 25% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_sigmaY_25pctDown/{grid}/merged",
            "card": card,
        }
    },

    "FCCee_Z_GHC_V25p1_sigmaZ_10pctUp": {
        grid: {
            "label": "#sigma_{z} 10% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_sigmaZ_10pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_sigmaZ_10pctDown": {
        grid: {
            "label": "#sigma_{z} 10% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_sigmaZ_10pctDown/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_sigmaZ_25pctUp": {
        grid: {
            "label": "#sigma_{z} 25% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_sigmaZ_25pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_sigmaZ_25pctDown": {
        grid: {
            "label": "#sigma_{z} 25% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_sigmaZ_25pctDown/{grid}/merged",
            "card": card,
        }
    },

    "FCCee_Z_GHC_V25p1_betaX_10pctUp": {
        grid: {
            "label": "#beta_{x} 10% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_betaX_10pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_betaX_10pctDown": {
        grid: {
            "label": "#beta_{x} 10% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_betaX_10pctDown/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_betaX_25pctUp": {
        grid: {
            "label": "#beta_{x} 25% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_betaX_25pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_betaX_25pctDown": {
        grid: {
            "label": "#beta_{x} 25% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_betaX_25pctDown/{grid}/merged",
            "card": card,
        }
    },

    "FCCee_Z_GHC_V25p1_betaY_10pctUp": {
        grid: {
            "label": "#beta_{y} 10% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_betaY_10pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_betaY_10pctDown": {
        grid: {
            "label": "#beta_{y} 10% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_betaY_10pctDown/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_betaY_25pctUp": {
        grid: {
            "label": "#beta_{y} 25% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_betaY_25pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_betaY_25pctDown": {
        grid: {
            "label": "#beta_{y} 25% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_betaY_25pctDown/{grid}/merged",
            "card": card,
        }
    },

    "FCCee_Z_GHC_V25p1_xing_10pctUp": {
        grid: {
            "label": "#alpha 10% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_xing_10pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_xing_10pctDown": {
        grid: {
            "label": "#alpha 10% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_xing_10pctDown/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_xing_5pctUp": {
        grid: {
            "label": "#alpha 5% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_xing_5pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_xing_5pctDown": {
        grid: {
            "label": "#alpha 5% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_xing_5pctDown/{grid}/merged",
            "card": card,
        }
    },

    "FCCee_Z_GHC_V25p1_np_10pctUp": {
        grid: {
            "label": "Bunch intensity 10% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_np_10pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_np_10pctDown": {
        grid: {
            "label": "Bunch intensity 10% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_np_10pctDown/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_np_25pctUp": {
        grid: {
            "label": "Bunch intensity 25% up",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_np_25pctUp/{grid}/merged",
            "card": card,
        }
    },
    "FCCee_Z_GHC_V25p1_np_25pctDown": {
        grid: {
            "label": "Bunch intensity 25% down",
            "dir": f"{base_dir_ipc_sensitivity_studies}/FCCee_Z_GHC_V25p1_np_25pctDown/{grid}/merged",
            "card": card,
        }
    },
}








