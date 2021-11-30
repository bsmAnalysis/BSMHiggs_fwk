## Instructions for running limits and making plots and tables

Updated 2021-11-30

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










