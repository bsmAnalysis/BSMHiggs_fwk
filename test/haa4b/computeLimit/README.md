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

### Make the pre and post fit plots.

The first argument is the channel (wh or zh) and the second argument is the path to the 60 GeV mass point limit output directory.
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

