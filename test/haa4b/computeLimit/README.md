## Instructions for running limits and making plots and tables

Updated 2021-11-30, Owen Long

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
cp -p <where-your-plotter-files-are>/plotter*.root $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b
cd $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b
hadd all_plotter_forLimits.root plotter_WH*.root
```

### Run the Combine limit jobs in condor

The first argument is the year (2016, 2017, 2018, or all) and the second argument is which type of limit to run.
Note that the setup above uses Combine v8.1.0, which gives good results for us.  In that version, the FitDiagnostics
output file is named **fitDiagnostics.root**.  In more recent versions of Combine, it may be named **fitDiagnosticsTest.root**.
If you update to a more recent version of Combine, you will need to make that change inside **optimize_haa.py**.

The most significant digit is:
 * 4 = Run limits jobs for an individual year specified by the first argument (2016, 2017, or 2018)
 * 5 = Run limits jobs for all of Run II combined (first argument is all)
 * 6 = Make Brazilian flag limit plots for Wh and Zh combined for an individual year.
 * 7 = Make Brazilian flag limit plots for Wh and Zh combined for all of Run II.

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
The output files will go into a directory that's made by the python script (e.g. **prefit-plots-wh-2016**).
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

The arguments are similar to the prefit plot script, as described above.
The output files will go into a directory that's made by the python script (e.g. **postfit-plots-wh-2016**).
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


### Make the summary tables

This is done with the **yieldTableBDTcut3.py** python script.  Here's an example.
```
python yieldTableBDTcut3.py cards_SB13TeV_SM_Wh_2016_noSoftb
python yieldTableBDTcut3.py cards_SB13TeV_SM_Wh_2017_noSoftb
python yieldTableBDTcut3.py cards_SB13TeV_SM_Wh_2018_noSoftb

python yieldTableBDTcut3.py cards_SB13TeV_SM_Zh_2016_noSoftb
python yieldTableBDTcut3.py cards_SB13TeV_SM_Zh_2017_noSoftb
python yieldTableBDTcut3.py cards_SB13TeV_SM_Zh_2018_noSoftb
```
The **make-all-tables.sh** shell script does this and also renames the output files to have unique names.  To use that, do this
```
sh make-all-tables.sh
```
The output files will be in the cards directories.  To see them all, do this ```ls -l cards*/yield-table*.tex```

To look at the tables, try doing ```pdflatex yield-tables.tex```




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
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2016_noSoftb 1
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2017_noSoftb 1
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2018_noSoftb 1

sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2016_noSoftb 1
sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2017_noSoftb 1
sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2018_noSoftb 1


sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2016_noSoftb 2
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2017_noSoftb 2
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2018_noSoftb 2

sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2016_noSoftb 2
sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2017_noSoftb 2
sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2018_noSoftb 2

#-- wait until all condor jobs finish before running step 3

sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2016_noSoftb 3
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2017_noSoftb 3
sh run-impacts-batch.sh cards_SB13TeV_SM_Wh_2018_noSoftb 3

sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2016_noSoftb 3
sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2017_noSoftb 3
sh run-impacts-batch.sh cards_SB13TeV_SM_Zh_2018_noSoftb 3
```

### Running the impacts and GOF tests

Under the directory Unblinding/Stage3 you may find the scripts:

sh run-impacts.sh
sh run-goftest.sh , run-goftest-2.sh

### Redoing the BDT binning optimization

Merge the plotter root files for the data taking years into one root file each for WH and ZH.

```
  hadd all-WH.root plotter_WH_2016_2020_06_19_forLimits.root plotter_WH_2017_2020_02_05_forLimits.root plotter_WH_2018_2020_02_05_forLimits.root
  hadd all-ZH.root plotter_ZH_2016_2020_06_19_forLimits.root plotter_ZH_2017_2020_02_05_forLimits.root plotter_ZH_2018_2020_02_05_forLimits.root
```

For each WH/ZH + nb + a mass combination you want to look at, create a simple root file with two histograms (h_sig and h_bg_sum) that holds the sum over lepton flavors and the sum over background components (for the background).

```
  root
  .L draw_plotter_input.c
  draw_plotter_input( 3, "all-WH.root", 60 )
  draw_plotter_input( 4, "all-WH.root", 60 )
  draw_plotter_input( 3, "all-ZH.root", 60 )
  draw_plotter_input( 4, "all-ZH.root", 60 )
  .q
```

Run the optimization

```
  root
  .L optimization_from_plotter_hists.c
  optimization_from_plotter_hists( "binning-optimization-input-WH-3b-ma60.root", 0.1, 2 )
  optimization_from_plotter_hists( "binning-optimization-input-WH-4b-ma60.root", 0.1, 2 )
  optimization_from_plotter_hists( "binning-optimization-input-ZH-3b-ma60.root", 0.1, 2 )
  optimization_from_plotter_hists( "binning-optimization-input-ZH-4b-ma60.root", 0.1, 2 )
  .q
```





