import ROOT
import os

# Input directory
input_dir = "output/h_zh_leptonic/histmaker/ecm240/"
input_dir = "output/h_zh_hadronic/histmaker/ecm240/"

# Processes to include
XX = ["qq", "bb", "ss", "cc"]

YY = ["mumu", "tautau", "Za", "aa", "bb", "cc", "ss", "gg", "WW", "ZZ"]

XX = ["qq", "bb", "ss", "cc"]

XX = ["ee"]

target = "cutFlow"

# Sum bin contents
n_initial = 0.0
n_selected = 0.0

for xx in XX:
    for yy in YY:

        filename = os.path.join(
            input_dir,
            f"wzp6_ee_{xx}H_H{yy}_ecm240.root"
        )

        if not os.path.exists(filename):
            print(f"File not found: {filename}")
            continue

        f = ROOT.TFile.Open(filename)

        if not f or f.IsZombie():
            print(f"Cannot open: {filename}")
            continue

        h = f.Get(target)

        if not h:
            print(f"Missing cutflow: {filename}")
            f.Close()
            continue

        n_initial += h.GetBinContent(1)
        n_selected += h.GetBinContent(7)

        f.Close()

# Overall efficiency
efficiency = n_selected / n_initial if n_initial > 0 else 0.0

print(f"Initial events:       {n_initial:.0f}")
print(f"Selected events:      {n_selected:.0f}")
print(f"Selection efficiency: {efficiency:.6f} ({efficiency * 100:.4f}%)")