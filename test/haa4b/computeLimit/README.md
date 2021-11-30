## Test limit for one A mass point
> source TestLimit.sh

where :
* m = A mass point
* histo = 2D histo in shape variable vs cut index
* index = cut index for BDT cut value (cuts on 2D bdt_shapes histo)
* bins = 3b, 4b etc bins to consider

## Run the limit scan in all A mass points
1. configure the input variable jsonPath, inUrl_wh and inUrl_zh in the oprimize_haa.py
2. submit jobs to the batch mode:
> python optimize_haa.py -p 4.0
  , where '4.0' is to run Wh channel, use 4.1 to run Zh channel and 4.2 to run Wh+Zh combined
3. when jobs finished, produce limit plots:
> python optimize_haa.py -p 6.0
  , where '6.0' is to produce plots for Wh channel, use 6.1 for Zh channel and 6.2 for Wh+Zh combined 

## Run the limit scan in all A mass points using full RunII data (2016+2017+2018)
1. configure the input variable jsonPaths, inUrl_whs and inUrl_zhs in the oprimize_haa.py
2. submit jobs to the batch mode:
> python optimize_haa.py -p 5.0
  , where '5.0' is to run Wh channel, use 5.1 to run Zh channel and 5.2 to run Wh+Zh combined
3. when jobs finished, produce limit plots:
> python optimize_haa.py -p 7.0
  , where '7.0' is to produce plots for Wh channel, use 7.1 for Zh channel and 7.2 for Wh+Zh combined 

## Run the BDT cuts optimization scan
> python optimize_haa.py -p 1
> python optimize_haa.py -p 2

## Run the postfit plots
1. configure the variables limit_dir, dir1 and wz in the file: postfitShapes.py.
2. and then: python postfitShapes.py

## produce the postfit-b event yields
under the Wh limits folder, launch:
> python $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/dumpPostfits.py -u fitDiagnostics.root

under the Zh limits folder, launch:
> python $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/dumpPostfits.py  -u fitDiagnostics_e.root (or fitDiagnostics_mu.root for mumu channel)
