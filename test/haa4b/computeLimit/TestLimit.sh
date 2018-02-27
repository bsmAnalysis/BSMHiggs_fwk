mkdir -p TEST;
cd TEST;
computeLimit --m 60 --histo bdt_shapes --in $PWD/../../plotter_2017_09_21_forLimits.root --syst --shape --index 10 --json $PWD/../../samples2016.json  --shapeMin -9999 --shapeMax 9999 --bins 3b --systpostfix _13TeV --dropBckgBelow 0.0 ;
sh combineCards.sh;

text2workspace.py card_combined.dat -o workspace.root --PO verbose  
combine -M FitDiagnostics workspace.root --plots --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --rMin=-2 --rMax=2 --robustFit 1 --expectSignal 1

python ../mlfitNormsToText.py -u fitDiagnostics.root

combine -M Asymptotic -m 60 workspace.root --run blind -v 3 > COMB.log;

cd ..
