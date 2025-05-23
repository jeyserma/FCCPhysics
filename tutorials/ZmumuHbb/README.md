


# Analyzing events for FCC-ee
In this tutorial, we're going to analyze events from the FCC-ee and measure cross-sections for selected processes. Throughout the tutorial, you'll learn how to perform a basic analysis, fit histograms, apply jet clustering and flavor tagging:

- **Part I:** Basic analysis and event selection, plotting, and cross-section measurement
- **Part II:** Fitting of histograms
- **Part III:** Jet clustering and flavor tagging



# Prerequisites

#### Computing environment

Instructions for account creation, logging in into the SubMIT cluster can be found [here](https://github.com/jeyserma/FCCPhysics/blob/main/tutorials/ZmumuHbb/SUBMIT.md).

An introduction to Linux and Git can be found [here](https://github.com/jeyserma/FCCPhysics/blob/main/tutorials/ZmumuHbb/LINUX.md).

#### Basic analysis tools (optional)

The following pages are not strictly required to complete the tutorial, but they provide helpful background on how analyses are typically performed in high energy physics. They cover useful tools, concepts, and workflows such as looping over events, making histograms, and plotting results.

- [Basics of High Energy Physics Computing and Analysis](https://github.com/jeyserma/FCCPhysics/blob/main/tutorials/ZmumuHbb/BASICS.md)
- [ROOT DataFrames](https://github.com/jeyserma/FCCPhysics/blob/main/tutorials/ZmumuHbb/RDATAFRAMES.md): ROOT’s high-level interface for processing events efficiently and in parallel across many files




#### FCCAnalysis framework
We will be working with the FCCAnalyses framework, available on [GitHub](https://github.com/HEP-FCC/FCCAnalyses). This is a common analysis framework developed for FCC-related studies, based on ROOT Dataframes. It lets you run full analyses over existing simulated samples, apply event selections, and produce plots and histograms. 

If you are running the tutorial on the MIT computing infrastructure, there is a pre-installed version of the analysis software available. You can set it up by running (this command must be executed every time you log into a new terminal session):

    source /work/submit/jaeyserm/software/FCCAnalyses/setup.sh


If you're running at CERN or prefer to install the analysis framework yourself, follow the instructions below (adapted from the [FCCAnalyses GitHub repository](https://github.com/HEP-FCC/FCCAnalyses)):

    cd go/to/my/directory
    source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10
    git clone --branch pre-edm4hep1 git@github.com:HEP-FCC/FCCAnalyses.git
    cd FCCAnalyses
    source ./setup.sh
    fccanalysis build -j 8

These setup steps need to be executed only once during the initial installation.

For subsequent sessions, simply run:

    cd go/to/my/directory/FCCAnalyses
    source setup.sh





#### Event samples 

For each physics process—such as signal, background, and specific Higgs decays—a sufficient number of events has been generated and is shared across all FCC analyses. These events are produced using various event generators (e.g., *Pythia*, *Whizard*) and processed with the [*Delphes* fast simulation framework](https://arxiv.org/abs/1307.6346) for detector simulation and reconstruction. The resulting events are stored as ROOT files in the [*edm4hep* data format](https://github.com/key4hep/EDM4hep), which is part of the [Key4HEP](https://key4hep.github.io/) software stack for future collider studies.

All samples in this tutorial are simulated using the [IDEA detector concept](https://arxiv.org/abs/2005.00022), which serves as the reference detector for FCC-ee physics studies.

You can find the full list of available samples [here](https://fcc-physics-events.web.cern.ch/fcc-ee/delphes/winter2023/idea/). 

In this tutorial, we will work with a selected subset of these samples, as described in a later section.

####  Getting the Tutorial Files

We’ve prepared a few Python files to guide you through this tutorial. You can find them in a separate repository  
[here](https://github.com/jeyserma/FCCPhysics/tree/main/tutorials/ZmumuHbb) [to be defined where we will store the tutorial files]

You can either download the files individually from the repository or clone the entire repository using:

    git clone https://github.com/jeyserma/FCCPhysics.git
    cd FCCPhysics/tutorials/ZmumuHbb

#### CMS Combination Tool (Combine)
We'll use the CMS Combination Tool for performing the statistical fits. This tool is built on top of the RooFit framework and is widely used within CMS for signal extraction, limit setting, and uncertainty estimation.

To simplify setup, we'll use a pre-compiled standalone version that can be run inside a Singularity image — an isolated Linux environment where we can install and run software without affecting the host system. This helps ensure the analysis runs the same way on different machines.

You can find the pre-built image at the following locations:

	/work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif # MIT
	/eos/project/f/fccsw-web/www/analysis/auxiliary/combine-standalone_v9.2.1.sif # CERN  

Instructions on how to run it with Singularity are provided below in the tutorial.

More information on how to compile the package locally and use all its features can be found in the [official documentation](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#oustide-of-cmssw-recommended-for-non-cms-users).




# Part I: Measuring a Cross-Section
In this first part of the tutorial, we will perform a basic event selection at the Higgs threshold. Our target process is:

```math
	\rm e⁺+e⁻ → ZH
```

with the Z boson decaying into a pair of muons (Z → μ⁺μ⁻), and the Higgs boson allowed to decay into any final state.

We aim to extract the total Higgs production cross-section using the recoil method. To isolate the signal, we perform a cutflow analysis to reduce the main backgrounds, which include:

- WW production  
- ZZ production  
- Z/γ* (Drell-Yan-like processes)

Finally, we construct a recoil mass histogram, which is used as input to the [CMS Combine tool](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) to estimate the uncertainty on the cross-section via a binned maximum likelihood fit.


### Step 1: Running the Event Selection

To generate the event selection and histograms for all signal and background samples (see below), run the following script from the main `FCCAnalyses` directory:

	fccanalysis run analysis_ZmumuH.py

The script typically takes a few minutes to run. It includes a main function, `build_graph(df, dataset)`, which is responsible for constructing the analysis graph. The `dataset` argument is a string that identifies the current dataset or process (see below for the list of processes used). The `df` argument is a dataframe object that you can modify within the function by applying selection cuts using `df.Filter()` or defining new variables on top of the existing ones using `df.Define()`. Histograms can be created with `df.Histo1D()` and should be added to the `hists` list. This list must be returned by the `build_graph()` function for the analysis to proceed correctly.

> *Exercise:* Inspect the `build_graph()` in the `analysis_ZmumuH.py` file. How many cuts are applied? What do they represent?

The processes that are included in during the analysis are specified in the processList dictionary:

    fraction = 0.2
    processList = {
        'wzp6_ee_mumuH_HXX_ecm240':      {'fraction': 1}, # signal
        'p8_ee_ZZ_ecm240':               {'fraction': fraction}, # background
        'p8_ee_WW_mumu_ecm240':          {'fraction': fraction}, # background
        'wzp6_ee_mumu_ecm240':           {'fraction': fraction}, # background
    }

Each process should be generated beforehand centrally (see section above for the Event Samples). For the Higgs signal sample `wzp6_ee_mumuH_HXX_ecm240`, we have several samples available for each individual Higgs decay.

The `fraction` parameter represents the fraction of the events you would like to run on. If you'd like to test the functionality of the script without processing the full datasets, you can reduce the number of events by adjusting the `fraction` values per process. A fraction of `1` means all events are processed. You can reduce this value (e.g., to `0.1`) to speed up the creation of histograms for testing purposes. By default, we run over the full signal statistics, but only 20% of the background statistics. With this configuration, in total you run over 100 million of events in just a few minutes!


Normalization of the histograms to a given integrated luminosity can be automated during the analysis. To enable this, modify the following parameters in `analysis_ZmumuH.py`:

    doScale = True
    intLumi = 10.8e6  # Integrated luminosity at the Higgs threshold (10.8 ab-1)

The cross-sections for each process are automatically taken into account in the background (the cross-sections are stored as meta information per sample).

The FCCAnalysis framework provides built-in functions—such as `FCCAnalyses::ReconstructedParticle::get_p()`—to access particles and their properties in the event samples, along with utilities for computing important observables in electron-positron collider physics. You can also define your own analysis-specific functions in C++, place them in a header file, and include them in your analysis script using the `includePaths` option, for example: `includePaths = ["utils.h"]`. This option accepts a list of header files that contain your custom functions and code snippets. In this tutorial, we use a `utils.h` file that defines several helper functions to compute the missing energy vector, its polar angle, and the acolinearity. Take a look at the contents of this file to understand how these quantities are calculated.

After running the script (which typically takes about 5 minutes depending on the load of the machine you're on), a ROOT file is created for each process in the output directory (`output/ZmumuH/histmaker/`). These files contain all the histograms needed for the next steps of the analysis. If you run into issues executing the setup or want to skip the processing step, we’ve pre-generated the files with full statistics. You can access them here:

	/ceph/submit/data/group/fcc/ee/tutorials/FNAL2025/ZmumuH/histmaker/ # MIT
	/eos/project/f/fccsw-web/www/tutorials/FNAL2025/ZmumuH/histmaker/ # CERN




### Step 2: Plotting Histograms
The FCCAnalysis framework provides a built-in plotting tool that you can use to visualize histograms.  
However, you are not limited to this tool—other plotting libraries such as `matplotlib` can also be used, especially for more customized or publication-ready plots (see the *Basic Analysis Tools* section for details).

To visualize the histograms and the cutflow plot, run the following command:

    fccanalysis plots plots_ZmumuH.py

Make sure to specify the input (by default `output/ZmumuH/histmaker/`) and output directories inside the `plots_ZmumuH.py` script (e.g. if you want to plot and run the pre-generated files as explained above, you should point to that directory).

The cutflow plot provides a clear illustration of how background events are progressively reduced by each selection cut, while ideally retaining the signal. To quantify the efficiency of this selection, we calculate the significance:

```math
\rm significance = \frac{S}{\sqrt{S + B}},
```

where `S` and `B` are the number of signal and background events, respectively, after each cut. As background is reduced and signal is preserved, the significance increases. This value reflects how confidently we can measure the signal in the presence of background. The uncertainty, in %, on the signal is calculated as 1/significance. In the special case where background is negligible, this expression simplifies to:

```math
\rm significance = \sqrt{S},
```

and the uncertainty on the signal yield becomes:

```math
\rm uncertainty =\frac{1}{\sqrt{S}}.
```

Maximizing the significance and mimnimizing the uncertainty is a key goal in any physics analysis — it's our job as physicists! More advanced methods, such as Machine Learning (covered later), can help you push this even further.

The significance after each cut has already been calculated for you and is listed in the file `cutFlow.txt`. Do you observe a final significance of 120? What is the corresponding uncertainty on the signal process?

> *Exercise:*  In the `build_graph` function, the **acolinearity** is computed between the two selected muons. Acolinearity measures how far the muons deviate from being perfectly back-to-back, which is characteristic of signal events from the **Z/γ\*** process.
> Update the plotting script to include this variable and produce a histogram of the acolinearity distribution.
> You can also explore whether applying a selection cut (a filter) on this variable helps suppress background events while retaining most of the signal. Try different cut values and compare the resulting distributions to evaluate the effectiveness of such a cut.

> *Exercise:*  In the `build_graph` function, the **missing energy** vector (`missingEnergy`) is defined and used to compute the **cosine of the polar angle** of the missing energy.  You can use this vector to calculate the **missing momentum** (i.e., its magnitude) with the built-in function `FCCAnalyses::ReconstructedParticle::get_p()` as we use for the muons, create a histogram, and add it to the plotting script.
> Once plotted, analyze the distribution of the missing momentum. Consider whether applying a cut on this variable can help suppress specific background processes. Which background can be effectively reduced by such a cut, and why? 


# Part II: Statistical Analysis
From the definition of significance above, we can also evaluate it bin-by-bin, treating each histogram bin as a separate measurement. Instead of tightening the recoil mass cut, we can instead combine information from all bins to extract a single uncertainty on the signal yield.

This is achieved using a likelihood fit, which does this combination optimally. The advantage of likelihood fits is that they also allow you to assign systematic uncertainties on both signal and background components — although that’s beyond the scope of this tutorial.

## Preparing the Datacards
To extract the uncertainty on the signal cross-section, we first need to prepare the Combine-compatible datacards. Use the following command to generate the required text and ROOT files containing the final histograms:

    fccanalysis combine combine_ZmumuH.py

The output files are saved in:

    output/ZmumuH/combine

This directory contains:
- A **text datacard**, also printed to screen  
- A **ROOT file** with the input histograms

The datacard is the key input to Combine — it tells the tool how to interpret the histograms in the context of a likelihood fit.


## Running the Fit and Output
To perform the likelihood fit using the CMS Combine tool inside a Singularity container, run (MIT and CERN respectively):

    singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd output/ZmumuH/combine; text2workspace.py datacard.txt -o ws.root; combine -M MultiDimFit -v 10 --rMin 0.9 --rMax 1.1 --setParameters r=1 ws.root'
    singularity exec /eos/project/f/fccsw-web/www/analysis/auxiliary/combine-standalone_v9.2.1.sif bash -c 'cd output/ZmumuH/combine; text2workspace.py datacard.txt -o ws.root; combine -M MultiDimFit -v 10 --rMin 0.9 --rMax 1.1 --setParameters r=1 ws.root'



After the fit runs, you’ll see output similar to this:

    Minuit2Minimizer : Valid minimum - status = 0
    FVAL  = -2.77185177588457066e-10
    Edm   = 3.21698353423223593e-13
    Nfcn  = 15
    r         = 1    +/-  0.00607216        (limited)
    Minimization finished with status=0
    Minimization success! status=0
    Minimized in 0.043089 seconds (0.040000 CPU time)
    FINAL NLL - NLL0 VALUE = -2.771851776e-10

From this output, you can extract:
- The best-fit value of `r` (signal strength), which is `1.0`
- The uncertainty on `r`, which is approximately `±0.006`, or 0.6%

> *Exercise:*  How does 0.6% uncertainty compares to the uncertainty from the significance? Is it better or worse? Think about:
> - Whether combining bins gives you more statistical power  
> - How much information you lost by applying a single tight cut  
> - The advantages of fitting a shape vs. counting events












# Part III: Higgs to a pair of b-quarks

## Extract the uncertainty for H → bb̄ 
The Higgs boson predominantly decays to a pair of b-quarks, with a branching ratio of approximately 58%. One of the key objectives of FCC-ee is to precisely measure the cross section for H → bb̄, which directly relates to the Higgs coupling to b-quarks.

So far, we have considered all Higgs decay modes (H → bb̄, cc̄, τ⁺τ⁻, WW, ZZ, etc.) as signal. To focus and measure specifically on H → bb̄, we must redefine the signal and background processes in the `combine_ZmumuH.py`

    sig_procs = {'sig': ['wzp6_ee_mumuH_Hbb_ecm240']}
    bkg_procs = {
        'ZHnoBB': [
            'wzp6_ee_mumuH_Hcc_ecm240',
            'wzp6_ee_mumuH_Hss_ecm240',
            'wzp6_ee_mumuH_Hgg_ecm240',
            'wzp6_ee_mumuH_Haa_ecm240',
            'wzp6_ee_mumuH_HZa_ecm240',
            'wzp6_ee_mumuH_HWW_ecm240',
            'wzp6_ee_mumuH_HZZ_ecm240',
            'wzp6_ee_mumuH_Hmumu_ecm240',
            'wzp6_ee_mumuH_Htautau_ecm240',
        ],
        'bkg': [
            'wzp6_ee_mumu_ecm240',
            'wzp6_ee_tautau_ecm240',
            'p8_ee_WW_mumu_ecm240',
            'p8_ee_ZZ_ecm240'
        ]
    }

This setup treats only H → bb̄ as signal. All other Higgs decays (e.g. H → cc̄, WW, gluons, taus, etc.) are treated as an additional background (ZHnoBB), along with non-Higgs processes. 

Let's run the Combine fit using this new signal/background definition: change the process definitions in ```combine_ZmumuH.py``` as defined above and re-run Combine. Extract the uncertainty on the cross section for H → bb̄:

- What value do you obtain?
- Compare it to the previous inclusive Higgs result (0.6% uncertainty).
- Is the result compatible with the expected statistical worsening by a factor of 1/√0.58?


While we might expect a statistical degradation of the uncertainty by a factor of 1/√0.58 ~ 1.31, that would naively lead to an 0.6%*1.31 ≈ 0.8% uncertainty), in practice, the result is worse. Why?

Because the non-bb Higgs decays form a significant background to H → bb̄ in our current selection, and they all exhibit the same shape of the recoil distribution. In order to overcome this, we will refine our selection and focus exclusively on the bb decays.



## Using the Flavor tagger
Both b-quarks from the Higgs decay manifest as jets in the detector — collimated sprays of particles resulting from hadronization. These final-state particles must be clustered into jets using a jet clustering algorithm. We use the [FastJet](https://indico.cern.ch/event/264054/contributions/592237/attachments/467910/648313/fastjet-doc-3.0.3.pdf) library, which provides a variety of jet algorithms commonly used in collider physics. For this analysis, we apply the Durham k<sub>T</sub> algorithm to cluster each event into exactly two jets, corresponding to the two b-quarks from the Higgs decay.

To distinguish b-jets from other types of jets (e.g. from gluons, c-quarks, or taus), we rely on jet flavor tagging. This process uses the properties of the jet — such as displaced vertices, track multiplicity, and invariant mass — to assign a probability that a given jet originated from a b-quark. This is done using a machine learning–based flavor tagger, specifically optimized for the FCC-ee environment. You can read more about the FCC-ee flavor tagger here: [FCC-ee Flavor Tagger – arXiv:2202.03285](https://arxiv.org/abs/2202.03285).


Both the jet clustering and flavor tagging can be easily performed within the FCCAnalyzer framework. To build on the previous `ZmumuH` analysis, we use a modified script:

    FCCPhysics run analysis_ZmumuHbb.py

Here are some important lines to inspect in comparison with the previous file:

    # Remove the muons from the event before clustering
    df = df.Define("rps_no_muons", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, muons)")

    # Cluster the remaining particles into exactly 2 jets
    jetClusteringHelper = ExclusiveJetClusteringHelper("rps_no_muons", 2, "")
    df = jetClusteringHelper.define(df)  # Performs clustering and creates jet collections

    # Setup and run the flavor tagger
    df = jetFlavourHelper.define(df)
    df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)  # Run inference

The jet clustering produces a collection of reconstructed jets, which can be converted to Lorentz vectors for physics analysis (e.g. dijet invariant mass). The flavor tagger returns per-jet probabilities for the jet being a `b`, `c`, `s`, `g`, or `τ` jet.

To select b-jets, we apply the following cut:

    df = df.Filter("recojet_isB[0] > 0.5 && recojet_isB[1] > 0.5")

This keeps only events where both jets have a b-tagging probability greater than 0.5.

Execute the full analysis pipeline:

    fccanalysis run analysis_ZmumuHbb.py
    fccanalysis plots plots_ZmumuHbb.py
    fccanalysis combine combine_ZmumuHbb.py
    singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd output/ZmumuHbb/combine; text2workspace.py datacard.txt -o ws.root; combine -M MultiDimFit -v 10 --rMin 0.9 --rMax 1.1 --setParameters r=1 ws.root'
    singularity exec /eos/project/f/fccsw-web/www/analysis/auxiliary/combine-standalone_v9.2.1.sif bash -c 'cd output/ZmumuHbb/combine; text2workspace.py datacard.txt -o ws.root; combine -M MultiDimFit -v 10 --rMin 0.9 --rMax 1.1 --setParameters r=1 ws.root'

Make sure to adapt the input/output directory in the scripts if needed. 

As you will notice, the first command might take a while to run — this is because the tagger inference step is relatively slow. Normally we parallelize this step and submit it to a batch system, which is outside of the scope of this tutorial. We’ve pre-generated the files with full statistics, and you can access them here:

	/ceph/submit/data/group/fcc/ee/tutorials/FNAL2025/ZmumuHbb/histmaker/ # MIT
	/eos/project/f/fccsw-web/www/tutorials/FNAL2025/ZmumuHbb/histmaker/ # CERN

If you use these samples, make sure to update the input directory in the plotting and combine scripts (or you can copy them to the input directory `output/ZmumuHbb/histmaker`).

> *Exercise:*  After running the full chain:
>- How much does the significance improve after applying the b-tagging probability cut?
>- How does the fit result (uncertainty on the H → bb̄ cross section) compare to the previous result without tagging?
>- Is it closer to the expected statistical limit?
>- Is the background from non-bb Higgs decays better suppressed? What about the backgrounds, in particular WW and Z/γ*?
>
> Use these comparisons to understand how flavor tagging enhances the precision of the measurement.











