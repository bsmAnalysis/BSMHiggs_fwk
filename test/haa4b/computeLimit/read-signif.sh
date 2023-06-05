dirs=("cards_SB13TeV_SM_2016_noSoftb" "cards_SB13TeV_SM_2017_noSoftb" "cards_SB13TeV_SM_2018_noSoftb" "cards_SB13TeV_SM_Wh_2016_noSoftb" "cards_SB13TeV_SM_Wh_2017_noSoftb" "cards_SB13TeV_SM_Wh_2018_noSoftb" "cards_SB13TeV_SM_Zh_2016_noSoftb" "cards_SB13TeV_SM_Zh_2017_noSoftb" "cards_SB13TeV_SM_Zh_2018_noSoftb" "cards_SB13TeV_SM_Wh_all_noSoftb_Combined" "cards_SB13TeV_SM_Zh_all_noSoftb_Combined" "cards_SB13TeV_SM_all_noSoftb_Combined")      
#dirs=("cards_SB13TeV_SM_2018_noSoftb" "cards_SB13TeV_SM_Wh_2018_noSoftb" "cards_SB13TeV_SM_Zh_2018_noSoftb")   

m_dirs=("15" "20" "25" "30" "40" "50" "60") 
#m_dir=("0015" "0020" "0030" "0040" "0050" "0060")

rm pvalues_*.txt

for dir in ${dirs[@]}; do
    
    for subdir in ${m_dirs[@]}; do
	
	cd $dir/00$subdir

	if [[ $dir == *"2016"* ]]; then
	    year="2016"
	elif [[ $dir == *"2017"* ]]; then 
	    year="2017"
	elif [[ $dir == *"2018"* ]]; then  
	    year="2018"
	else 
	    year="all"
	fi
	
	if [[ $dir == *"all"* ]]; then
	    file=COMB-signif.log   
	    if [[ $dir == *"Wh"* ]]; then 
		ch="Wh"
	    elif [[ $dir == *"Zh"* ]]; then  
		ch="Zh"
	    else 
		ch=""
	    fi
	    outfile=pvalues_${ch}_${year}_m${subdir}.txt; rm $outfile  
	    more $file |grep "p-value" | awk '{print $4}' >>$outfile  
	    mv $outfile ../../   
	    
	elif [[ $dir == *"Wh"* ]]; then
	    ch="Wh"
	    file=COMB-signif-second.log
	    outfile=pvalues_${ch}_${year}_m${subdir}.txt; rm $outfile
	    
	    more $file |grep "p-value" | awk '{print $4}' >>$outfile    
	    mv $outfile ../../
	    
	elif [[ $dir == *"Zh"* ]]; then
	    ch="Zh"
	    outfile=pvalues_${ch}_${year}_m${subdir}.txt;   rm $outfile  
	    
	    file=COMB-signif.log
	     more $file |grep "p-value" | awk '{print $4}' >>$outfile      
#	    e_file=COMB-signif-e.log
#	    more $e_file |grep "p-value" | awk '{print $4}' >>$outfile       
#	    mu_file=COMB-signif-mu.log
#	    more $mu_file |grep "p-value" | awk '{print $4}' >>$outfile
	    mv $outfile ../../         
	    
	else
	    ch=""
	    outfile=pvalues_${ch}_${year}_m${subdir}.txt; rm $outfile  

	    file=COMB-signif.log
	    more $file |grep "p-value" | awk '{print $4}' >>$outfile 
#	    e_file=COMB_e.log
#	    more $e_file |grep "p-value" | awk '{print $4}' >>$outfile   
#	    mu_file=COMB_mu.log
#	    more $mu_file |grep "p-value" | awk '{print $4}' >>$outfile 
	    mv $outfile ../../         
	    
	fi
	cd ../../
    done
done
