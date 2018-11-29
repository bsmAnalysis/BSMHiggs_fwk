# Installation for 94X (2017) 
```bash
export SCRAM_ARCH=slc6_amd64_gcc630
#or
setenv SCRAM_ARCH slc6_amd64_gcc630
cmsrel CMSSW_9_4_9
cd CMSSW_9_4_9/src
cmsenv
git cms-init

#Packages to rerun electrons energy correction
#https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#2017%20MiniAOD%20V2
git cms-merge-topic cms-egamma:EgammaPostRecoTools_940
scram b -j 4

#Calculating corrections and uncertainties on MET
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription#Instructions%20for%20%209_4_X,%20X%20%3E=0%20f
git cms-merge-topic cms-met:METFixEE2017_949_v2
scram b -j 4

#Add Fall17_94X_V2 ID for electrons
#https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
git cms-merge-topic guitargeek:EgammaID_9_4_X
scram b -j 4

#Add Fall17_94X_V2 ID for photons
#https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP
scram b -j 10
# Add the area containing the MVA weights (from cms-data, to appear in "external").
# Note: the "external" area appears after "scram build" is run at least once, as above
cd $CMSSW_BASE/external/slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
cd data/RecoEgamma/PhotonIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/external/slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/src

git clone https://github.com/bsmAnalysis/BSMHiggs_fwk.git UserCode/bsmhiggs_fwk
cd $CMSSW_BASE/src/UserCode/bsmhiggs_fwk
git checkout 2017_Analysis
cd $CMSSW_BASE/src

#And compile
scram b -j 4
```

# Installation for 80X (2016) 

```bash
export SCRAM_ARCH=slc6_amd64_gcc530
#or:
setenv SCRAM_ARCH slc6_amd64_gcc530
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src
cmsenv
 
#Checkout Some Packages from Egamma 
# ( https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMRegression#Consistent_EGMSmearer )
git cms-init 
git cms-merge-topic cms-egamma:EGM_gain_v1
cd $CMSSW_BASE/src/EgammaAnalysis/ElectronTools/data
git clone -b Moriond17_gainSwitch_unc https://github.com/ECALELFS/ScalesSmearings.git
cd $CMSSW_BASE/src
scram b -j 8

#Adding DeepCSV
git cms-merge-topic -u mverzett:DeepFlavour-from-CMSSW_8_0_21 
mkdir $CMSSW_BASE/src/RecoBTag/DeepFlavour/data/
cd $CMSSW_BASE/src/RecoBTag/DeepFlavour/data/  
wget http://home.fnal.gov/~verzetti//DeepFlavour/training/DeepFlavourNoSL.json
cd $CMSSW_BASE/src/ 
scram b -j 8

# Check out packages for Hbb tagging (https://twiki.cern.ch/twiki/bin/viewauth/CMS/Hbbtagging#8_0_X)
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
#or
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch --tags btv-cmssw
git cms-merge-topic -u cms-btv-pog:BoostedDoubleSVTaggerV4-WithWeightFiles-v1_from-CMSSW_8_0_21
scram b -j 8

#git clone -b svFit_2015Apr03 https://github.com/veelken/SVfit_standalone.git TauAnalysis/SVfitStandalone

#https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription#Instructions_for_8_0_X_X_26_patc
git cms-merge-topic cms-met:METRecipe_8020 -u
git cms-merge-topic cms-met:METRecipe_80X_part2 -u


git clone https://github.com/bsmAnalysis/BSMHiggs_fwk.git UserCode/bsmhiggs_fwk
cd $CMSSW_BASE/src/UserCode/bsmhiggs_fwk
git checkout -b modified #copy the branch to a new one to host future modifications (ease pull request and code merging)
cd $CMSSW_BASE/src

## HiggsCombine in 80X 
cd $CMSSW_BASE/src
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
##Update to a reccomended tag - currently the reccomended tag is v7.0.6
git checkout v7.0.6
scramv1 b clean; scramv1 b # always make a clean build
cd ../../

#JetToolBox: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetToolbox
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_80X_V3

#And compile
scramv1 b -j 16 
```

# For developers

We have decided to use pull-request mode for the master development.

- Fork the code with your personal github ID. See [details](https://help.github.com/articles/fork-a-repo/)
- Make a clean git clone in the UserCode directory
```
cd $CMSSW_BASE/src/UserCode 
git clone git@github.com:yourgithubid/BSMHiggs_fwk.git bsmhiggs_fwk
cd bsmhiggs_fwk
git remote add upstream git@github.com:bsmAnalysis/BSMHiggs_fwk.git

- Update your repository and start your changes:
git remote update
git merge upstream/master
```
- Make your own change and commit
```
git commit -a -m "Added feature A, B, C"
git push
```
- Make a pull request against the bsmAnalysis. See [details](https://help.github.com/articles/using-pull-requests/)


# For creating private samples


# Analysing
See example files in UserCode/bsmhiggs_fwk/test/haa4b/ Executing the
file 'submit.sh' will run the analysis on all samples defined in
'samples.json'.

Please see the full usage of 'submit.sh' with

> ./submit.sh

# Analysis specific
If you need set your stuff under test/analysis_subdir

In general to run an executable or script you can do:
```
runLocalAnalysisOverSamples.py -e my_exe -j data/my_samples.json -d my_input_dir -o my_output_dir -c test/runAnalysis_cfg.py.templ -p "@par1=val1 @par2=val2" -s queue
```

 The easiest is however to create a submit.sh script in the directory
 of your analysis.  The executable can also be added in the
 UserCode/bsmhiggs_fwk/test/haa4b/bin directory (also add it to the
 Buildfile there)

