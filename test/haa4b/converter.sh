#!/usr/bin/env bash
# Note: this script is designed to be run under the folder: test/haa4b

echo "Please select what you want to do:"
echo "1. Switch 2016 to 2017   2. Switch 2017 to 2016"
read select

if [ ${select} -eq 1 ]
then
    echo "Switching 2016 to 2017..."
    com='/#define YEAR_2017/{s/.*YEAR_2017/#define YEAR_2017/g}'
elif [ ${select} -eq 2 ]
then
    echo "Switching 2017 to 2016..."
    com='s/#define YEAR_2017/\/\/#define YEAR_2017/g'
else
    echo "Input invalid, exsting..."
    exit 0
fi

 path=`pwd`
 Files=("/../../interface/DataEvtSummaryHandler.h" "/../../src/DataEvtSummaryHandler.cc" "/../../plugins/mainNtuplizer.cc" "/../../interface/BSMPhysicsEvent.h" "/../../src/BSMPhysicsEvent.cc" "/../../bin/haa4b/runhaaAnalysis.cc" "/../../src/PatUtils.cc")
 for file in "${Files[@]}"
 do
    file="${path}${file}"
    if [ -f ${file} ]
    then
        sed -i "${com}" ${file}
        echo "sed -i \"${com}\" ${file}"
#        echo "Chaning file: "${file}
    else
        echo "Error! Cannot find the file: "${file}
        echo "Please make sure you run this script under the folder: test/haa4b"
    fi
 done

 echo "Done! Please compile again. If error happens, please fix it before compiling!"
