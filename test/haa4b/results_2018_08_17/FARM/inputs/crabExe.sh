#!/bin/sh

echo "Here there are all the input arguments"
echo $@

# Copy files to lib and bin dir
cp ./lib/$SCRAM_ARCH/* $CMSSW_BASE/lib/$SCRAM_ARCH
cp -rd ./src/* $CMSSW_BASE/src/.
cp -rd ./python/* $CMSSW_BASE/python/.
cp  $CMSSW_BASE/bin/$SCRAM_ARCH
cp x509_proxy $CMSSW_BASE/
export X509_USER_PROXY=$CMSSW_BASE/x509_proxy

#If you are curious, you can have a look at the tweaked PSet.
echo "================= PSet.py file =================="
cat PSet.py
# This is what you need if you want to look at the tweaked parameter set!!
echo "================= Dumping PSet ===================="
python -c "import PSet; print PSet.process.dumpPython()"

#just needed to create the JobReport
cmsRun -j FrameworkJobReport.xml -p PSet.py

echo "================= Dumping Input files ===================="
python -c "import PSet; print 
.join(list(PSet.process.source.fileNames))"

#Actually run the script
 PSet.py
