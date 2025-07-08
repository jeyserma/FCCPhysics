# Material budget in Fast and FullSim

## FastSim
For Delphes cards it only works for the tracker definition with the track covariance module.

It takes an input Delphes card as argument. The configuration of the materials and grouping are defined in the python script.

	python delphesScanPlot.py -i delphes_card_IDEA.tcl --ymax 15

It produces the material budget plot for the defined in the Python file and a ROOT file with the histograms.

## FullSim

Need to run on XML files. First we run the scan for different angles (cos(theta)), which produces a ROOT file with the material budgets of all materials:

	k4run fullSimScan.py 

Then we can plot it:

	python fullSimPlot.py -f out_material_scan.root -i Air --ymax 15

It also produces the material budged plot for all the materials and a grouped definition as the defined in Python file. A ROOT file is also written containing the histograms.

## Compare

To compare a given detector group (e.g. beam pipe or vertex), you can use the following script:

	python compare.py

The relevant configurations for the plotting are written in the Python file.