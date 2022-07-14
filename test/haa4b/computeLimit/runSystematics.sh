rm -rf TEST_*

mkdir -p TEST_WH_2016;
cd TEST_WH_2016;

computeLimit --m 60 --histo bdt_shapes --lumi 36.3 --in $PWD/../../plotter_WH_2016_2020_06_19_forLimits.root --syst --simfit --shape --index 1 --json $PWD/../../samples2016_legacy.json  --shapeMin -9999 --shapeMax 9999 --bins 3b,4b --systpostfix _13TeV --dropBckgBelow 0.015 ;
#sh combineCards.sh;
cd ../
mkdir -p TEST_ZH_2016;
cd TEST_ZH_2016;

computeLimit --m 60 --runZh --histo bdt_shapes --lumi 36.3 --in $PWD/../../plotter_ZH_2016_2020_06_19_forLimits.root --syst --simfit --shape --index 1 --json $PWD/../../samples2016_legacy.json  --shapeMin -9999 --shapeMax 9999 --bins 3b,4b --systpostfix _13TeV --dropBckgBelow 0.015 ;

#cd ../
#mkdir TEST_WH_1718
#cd TEST_WH_1718

#computeLimit --m 60 --histo bdt_shapes --lumi 101.2 --in $PWD/../../plotter_WH_201718_2020_02_05_forLimits.root --syst --shape --index 1 --json $PWD/../../samples2017.json  --shapeMin -9999 --shapeMax 9999 --bins 3b,4b --systpostfix _13TeV --dropBckgBelow 0.015 ;

#cd ../
#mkdir TEST_ZH_1718
#cd TEST_ZH_1718 

#computeLimit --m 60 --runZh --histo bdt_shapes --lumi 101.2 --in $PWD/../../plotter_ZH_201718_2020_02_05_forLimits.root --syst --shape --index 1 --json $PWD/../../samples2017.json  --shapeMin -9999 --shapeMax 9999 --bins 3b,4b --systpostfix _13TeV --dropBckgBelow 0.015 ;

cd ..

