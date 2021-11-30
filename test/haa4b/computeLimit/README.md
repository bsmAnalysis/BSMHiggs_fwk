## Instructions for running limits and making plots and tables

Updated 2021-11-30

### Set up your release directory, check out code, and compile

Here's a sketch of what I do.
```
cd <your scratch area>
mkdir limits-combine-v8.1.0
cd limits-combine-v8.1.0
setenv SCRAM_ARCH slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src/
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.1.0
scramv1 b | & tee build-try1a.log
cd ../..
git clone https://github.com/bsmAnalysis/BSMHiggs_fwk.git UserCode/bsmhiggs_fwk
cd $CMSSW_BASE/src/UserCode/bsmhiggs_fwk
cd test/haa4b
#-- type 1 when prompted by converter.sh below
sh converter.sh
cd $CMSSW_BASE/src
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
scram b -j 4 | & tee build-try2a.log
```
Copy the plotter root files to the ```src/UserCode/bsmhiggs_fwk/test/haa4b``` directory and run hadd on the Wh files (needed for the DDQCD).
```
cp <where-your-plotter-files-are>/plotter*.root $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b
cd $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b
hadd all_plotter_forLimits.root plotter_WH*.root
```

### Run the Combine limit jobs in condor

The first argument is the year (2016, 2017, 2018, or all) and the second argument is which type of limit to run.

The most significant digit is:
 * 4 = Run limits for an individual year specified by the first argument (2016, 2017, or 2018)
 * 5 = Run limits for all of Run II combined (first argument is all)
 * 6 = Run limits for Wh and Zh combined for an individual year.
 * 7 = Run limits for Wh and Zh combined for all of Run II.

The least significant digit is:
 * 0 = Run limit for Wh channel
 * 1 = Run limit for Zh channel
 * 2 = Run limit for Wh+Zh channels combined

```
  cd $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit

  python optimize_haa.py 2016 4.0
  python optimize_haa.py 2016 4.1
  python optimize_haa.py 2016 4.2

  python optimize_haa.py 2017 4.0
  python optimize_haa.py 2017 4.1
  python optimize_haa.py 2017 4.2

  python optimize_haa.py 2018 4.0
  python optimize_haa.py 2018 4.1
  python optimize_haa.py 2018 4.2

#-- All of the above jobs need to finish before running the "all" jobs below.
#   This is because the script needs to extract the scale factors from the completed jobs above.

  python optimize_haa.py all  5.0
  python optimize_haa.py all  5.1
  python optimize_haa.py all  5.2
```

### Make the prefit plots

The first argument is the channel (wh or zh) and the second argument is the path to the 60 GeV mass point limit output directory.
The optional ```--do_liny``` and ```--do_linzoom``` arguments run plots with a linear y scale and a linear y scale zoomed in to show the region where the signal is (low numbers).
```
python prefitShapes.py zh cards_SB13TeV_SM_Zh_2016_noSoftb/0060/
python prefitShapes.py zh cards_SB13TeV_SM_Zh_2016_noSoftb/0060/ --do_liny
python prefitShapes.py zh cards_SB13TeV_SM_Zh_2016_noSoftb/0060/ --do_linzoom

python prefitShapes.py zh cards_SB13TeV_SM_Zh_2017_noSoftb/0060/
python prefitShapes.py zh cards_SB13TeV_SM_Zh_2017_noSoftb/0060/ --do_liny
python prefitShapes.py zh cards_SB13TeV_SM_Zh_2017_noSoftb/0060/ --do_linzoom

python prefitShapes.py zh cards_SB13TeV_SM_Zh_2018_noSoftb/0060/
python prefitShapes.py zh cards_SB13TeV_SM_Zh_2018_noSoftb/0060/ --do_liny
python prefitShapes.py zh cards_SB13TeV_SM_Zh_2018_noSoftb/0060/ --do_linzoom



python prefitShapes.py wh cards_SB13TeV_SM_Wh_2016_noSoftb/0060/
python prefitShapes.py wh cards_SB13TeV_SM_Wh_2016_noSoftb/0060/ --do_liny
python prefitShapes.py wh cards_SB13TeV_SM_Wh_2016_noSoftb/0060/ --do_linzoom

python prefitShapes.py wh cards_SB13TeV_SM_Wh_2017_noSoftb/0060/
python prefitShapes.py wh cards_SB13TeV_SM_Wh_2017_noSoftb/0060/ --do_liny
python prefitShapes.py wh cards_SB13TeV_SM_Wh_2017_noSoftb/0060/ --do_linzoom

python prefitShapes.py wh cards_SB13TeV_SM_Wh_2018_noSoftb/0060/
python prefitShapes.py wh cards_SB13TeV_SM_Wh_2018_noSoftb/0060/ --do_liny
python prefitShapes.py wh cards_SB13TeV_SM_Wh_2018_noSoftb/0060/ --do_linzoom
```

### Make the postfit plots

```
python postfitShapes2.py zh cards_SB13TeV_SM_Zh_2016_noSoftb/0060/
python postfitShapes2.py zh cards_SB13TeV_SM_Zh_2016_noSoftb/0060/ --do_liny
python postfitShapes2.py zh cards_SB13TeV_SM_Zh_2016_noSoftb/0060/ --do_linzoom

python postfitShapes2.py zh cards_SB13TeV_SM_Zh_2017_noSoftb/0060/
python postfitShapes2.py zh cards_SB13TeV_SM_Zh_2017_noSoftb/0060/ --do_liny
python postfitShapes2.py zh cards_SB13TeV_SM_Zh_2017_noSoftb/0060/ --do_linzoom

python postfitShapes2.py zh cards_SB13TeV_SM_Zh_2018_noSoftb/0060/
python postfitShapes2.py zh cards_SB13TeV_SM_Zh_2018_noSoftb/0060/ --do_liny
python postfitShapes2.py zh cards_SB13TeV_SM_Zh_2018_noSoftb/0060/ --do_linzoom



python postfitShapes2.py wh cards_SB13TeV_SM_Wh_2016_noSoftb/0060/
python postfitShapes2.py wh cards_SB13TeV_SM_Wh_2016_noSoftb/0060/ --do_liny
python postfitShapes2.py wh cards_SB13TeV_SM_Wh_2016_noSoftb/0060/ --do_linzoom

python postfitShapes2.py wh cards_SB13TeV_SM_Wh_2017_noSoftb/0060/
python postfitShapes2.py wh cards_SB13TeV_SM_Wh_2017_noSoftb/0060/ --do_liny
python postfitShapes2.py wh cards_SB13TeV_SM_Wh_2017_noSoftb/0060/ --do_linzoom

python postfitShapes2.py wh cards_SB13TeV_SM_Wh_2018_noSoftb/0060/
python postfitShapes2.py wh cards_SB13TeV_SM_Wh_2018_noSoftb/0060/ --do_liny
python postfitShapes2.py wh cards_SB13TeV_SM_Wh_2018_noSoftb/0060/ --do_linzoom

```

### Make the Brazilian flag limit plots

```
python optimize_haa.py 2016 6.0
python optimize_haa.py 2016 6.1
python optimize_haa.py 2016 6.2

python optimize_haa.py 2017 6.0
python optimize_haa.py 2017 6.1
python optimize_haa.py 2017 6.2

python optimize_haa.py 2018 6.0
python optimize_haa.py 2018 6.1
python optimize_haa.py 2018 6.2

python optimize_haa.py all 7.0
python optimize_haa.py all 7.1
python optimize_haa.py all 7.2
```

### Run the impacts in the condor batch system

The first argument is the Combine job output directory for the 60 GeV mass point.  The second argument is the step, where
 * 1 = do initial fit
 * 2 = submit batch jobs to condor
 * 3 = make the plots pdf file (run after all condor jobs finish from step 2)

Note:  If your fits take longer than the time limit for the default queue in condor, you may need to set the queue to workday.  To do that, edit this file
     ```CombineHarvester/CombineTools/python/combine/CombineToolBase.py``` and add the following line
```
+JobFlavour = "workday"
```
in the CONDOR_TEMPLATE section right after the line that sets log.  If step3 fails for you and you don't get the pdf file, this is probably the reason.

```
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2016_noSoftb/0060/ 1
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2017_noSoftb/0060/ 1
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2018_noSoftb/0060/ 1

sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2016_noSoftb/0060/ 1
sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2017_noSoftb/0060/ 1
sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2018_noSoftb/0060/ 1


sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2016_noSoftb/0060/ 2
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2017_noSoftb/0060/ 2
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2018_noSoftb/0060/ 2

sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2016_noSoftb/0060/ 2
sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2017_noSoftb/0060/ 2
sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2018_noSoftb/0060/ 2

#-- wait until all condor jobs finish before running step 3

sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2016_noSoftb/0060/ 3
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2017_noSoftb/0060/ 3
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2018_noSoftb/0060/ 3

sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2016_noSoftb/0060/ 3
sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2017_noSoftb/0060/ 3
sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2018_noSoftb/0060/ 3
```










