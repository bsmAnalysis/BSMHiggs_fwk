#rm -rf TEST_*

mkdir -p TEST_WH_all;
cd TEST_WH_all;

computeLimit --m 20 --histo bdt_shapes --lumi 138.0 --in $PWD/../../plotter_WH_all.root --syst --shape --index 1 --json $PWD/../../samples2017.json  --shapeMin -9999 --shapeMax 9999 --bins 3b,4b --systpostfix _13TeV --dropBckgBelow 0.015 ;

cd ../
mkdir -p TEST_ZH_all;
cd TEST_ZH_all;

computeLimit --m 20 --runZh --histo bdt_shapes --lumi 138.0 --in $PWD/../../plotter_ZH_all.root --syst --shape --index 1 --json $PWD/../../samples2017.json  --shapeMin -9999 --shapeMax 9999 --bins 3b,4b --systpostfix _13TeV --dropBckgBelow 0.015 ;

cd ..

