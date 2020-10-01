rm -rf TEST_*

mkdir -p TEST_WH_2016;
cd TEST_WH_2016;

computeLimit --m 60 --histo bdt_shapes --lumi 35.9 --in $PWD/../../plotter_2016_WH_Sys_noSoftb_forLimits-replaced-signal-mc.root --syst --shape --index 1 --json $PWD/../../samples2016_legacy.json  --shapeMin -9999 --shapeMax 9999 --bins 3b,4b --systpostfix _13TeV --dropBckgBelow 0.015 ;
#sh combineCards.sh;
cd ../
mkdir -p TEST_ZH_2016;
cd TEST_ZH_2016;

computeLimit --m 60 --runZh --histo bdt_shapes --lumi 35.9 --in $PWD/../../plotter_2016_ZH_Sys_noSoftb_forLimits-replaced-signal-mc.root --syst --shape --index 1 --json $PWD/../../samples2016_legacy.json  --shapeMin -9999 --shapeMax 9999 --bins 3b,4b --systpostfix _13TeV --dropBckgBelow 0.015 ;

cd ../
mkdir TEST_WH_1718
cd TEST_WH_1718

computeLimit --m 60 --histo bdt_shapes --lumi 101.2 --in $PWD/../../plotter_1718_WH_Sys_noSoftb_forLimits-replaced-signal-mc.root --syst --shape --index 1 --json $PWD/../../samples2017.json  --shapeMin -9999 --shapeMax 9999 --bins 3b,4b --systpostfix _13TeV --dropBckgBelow 0.015 ;

cd ../
mkdir TEST_ZH_1718
cd TEST_ZH_1718 

computeLimit --m 60 --runZh --histo bdt_shapes --lumi 101.2 --in $PWD/../../plotter_1718_ZH_Sys_noSoftb_forLimits-replaced-signal-mc.root --syst --shape --index 1 --json $PWD/../../samples2017.json  --shapeMin -9999 --shapeMax 9999 --bins 3b,4b --systpostfix _13TeV --dropBckgBelow 0.015 ;

cd ..

