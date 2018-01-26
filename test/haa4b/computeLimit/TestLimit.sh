mkdir -p TEST;
cd TEST;
computeLimit --m 60 --histo bdt_shapes --in $PWD/../../plotter_2017_09_21.root  --index 5 --json $PWD/../../samples2016.json  --shapeMin -9999 --shapeMax 9999  --bins 3b --systpostfix _13TeV --rebin 8 --dropBckgBelow 0.0 ;
sh combineCards.sh;
combine -M Asymptotic -m 60 --run expected card_combined.dat > COMB.log;
cd ..
