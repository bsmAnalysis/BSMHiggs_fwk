#!/bin/bash

year=$1 
#("2016" "2017" "2018" "all")

echo -e "|-------------------------------------------------------------------------------------------------------------|"  
format="%15s%20s%20s%20s\n"
printf "$format" "    p-value (mA) " "|  p-value (WH) " "|  p-value (ZH)   "  "|   p -value (WH+ZH)   "    
#printf "$format" "    p-value (mA) " "|  p-value (WH) " "|  p-value (ZH-e) " "| p-value (ZH-mu) " "|  P-value (WH+ZH)   "
echo -e "|-------------------------------------------------------------------------------------------------------------|" 

m_dirs=("15" "20" "25" "30" "40" "50" "60")

for mA in ${m_dirs[@]}; do                                                                                                                              

    file3=pvalues__${year}_m${mA}.txt ; file1=pvalues_Wh_${year}_m${mA}.txt ; file2=pvalues_Zh_${year}_m${mA}.txt 

    p1=`sed -n '1p' $file1` ; p1=`printf "%.4f" $p1`
    p2=`sed -n '1p' $file2` ;p2=`printf "%.4f" $p2`;
#    p2_e=`sed -n '1p' $file2` ; p2_mu=`sed -n '2p' $file2`
#    p2_e=`printf "%.4f" $p2_e`; p2_mu=`printf "%.4f" $p2_mu`
    p3=`sed -n '1p' $file3` ; p3=`printf "%.4f" $p3`;
#    p3_e=`sed -n '1p' $file3` ; p3_mu=`sed -n '2p' $file3`   
#    p3_e=`printf "%.4f" $p3_e`; p3_mu=`printf "%.4f" $p3_mu`   

    if [[ $year == *"all"* ]]; then
	printf "$format" $mA $p1 $p2 $p3 # $p2_e  $p3 # $p3_e     
    else 
	printf "$format" $mA $p1 $p2 $p3 # $p2_mu $p3 # $p3_mu
    fi
done
echo -e "|-------------------------------------------------------------------------------------------------------------|"   
