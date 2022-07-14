

  printf "\n\n Number of arguments: %d\n\n" $#

  if [ $# -ne 2 ]
  then
     printf "\n\n Usage:\n"
     printf "  sh run-impacts-batch.sh <cards-dir> <step>\n"
     printf "\n"
     printf "    step:\n"
     printf "     1 : do initial fit\n"
     printf "     2 : submit batch jobs\n"
     printf "     3 : make the plots pdf file\n"
     printf "\n"
     printf "\n"
     printf "   Note: If your fits take longer than the time limit for the default queue in condor\n"
     printf "         you may need to set the queue to workday.  To do that, edit this file\n"
     printf "         CombineHarvester/CombineTools/python/combine/CombineToolBase.py\n"
     printf "         and add the following line\n\n"
     printf "+JobFlavour = \"workday\" \n\n"
     printf "         in the CONDOR_TEMPLATE section right after the line that sets log.\n\n"
     printf "         If step3 fails for you and you don't get the pdf file, this is probably the reason.\n\n\n\n\n"
     exit
  fi


  cards_dir=`echo $1 | sed s#\/##`

  if [ ! -d "$cards_dir" ]
  then
     printf "\n\n *** cards dir %s does not exist.\n\n" $cards_dir
     exit -1
  else
     printf "\n\n  cards dir = %s\n\n" $cards_dir
  fi
      
  step=$2

  printf "\n step = %d\n\n" $step



  impacts_dir=`echo $cards_dir | sed s/cards/impacts/`

  printf "\n\n impacts output directory: %s\n\n" $impacts_dir



# if [ -d "$impacts_dir" ]
# then
#    datestring=`date +%0G-%0m-%0d`
#    mv_dir=`printf "old-%s-%s" $impacts_dir $datestring`
#    printf "  impacts dir exists.  moving it to %s \n\n" $mv_dir
#    mv $impacts_dir $mv_dir
# fi


  if [[ $cards_dir == *"Wh_all"* ]]; then 
     printf "\n Found Wh_all\n\n"
     rangeString="tt_norm_e_2016=0.1,4.0:tt_norm_mu_2016=0.1,4.0:w_norm_e_2016=0.1,4:w_norm_mu_2016=0.1,4:tt_norm_e_2017=0.1,4.0:tt_norm_mu_2017=0.1,4.0:w_norm_e_2017=0.1,4:w_norm_mu_2017=0.1,4:tt_norm_e_2018=0.1,4.0:tt_norm_mu_2018=0.1,4.0:w_norm_e_2018=0.1,4:w_norm_mu_2018=0.1,4"
  elif [[ $cards_dir == *"Zh_all"* ]]; then
     printf "\n Found Zh_all\n\n"
     rangeString="z_norm_3b_e_2016=0.1,4.0:z_norm_3b_mu_2016=0.1,4.0:z_norm_4b_e_2016=0.1,4.0:z_norm_4b_mu_2016=0.1,4.0:tt_norm_e_2016=0.1,4.0:tt_norm_mu_2016=0.1,4.0:z_norm_3b_e_2017=0.1,4.0:z_norm_3b_mu_2017=0.1,4.0:z_norm_4b_e_2017=0.1,4.0:z_norm_4b_mu_2017=0.1,4.0:tt_norm_e_2017=0.1,4.0:tt_norm_mu_2017=0.1,4.0:z_norm_3b_e_2018=0.1,4.0:z_norm_3b_mu_2018=0.1,4.0:z_norm_4b_e_2018=0.1,4.0:z_norm_4b_mu_2018=0.1,4.0:tt_norm_e_2018=0.1,4.0:tt_norm_mu_2018=0.1,4.0"
  elif [[ $cards_dir == *"Wh_"* ]]; then
     printf "\n Found Wh_\n\n"
     rangeString="tt_norm_e=0.1,4.0:tt_norm_mu=0.1,4.0:w_norm_e=0.1,4:w_norm_mu=0.1,4"
  elif [[ $cards_dir == *"Zh_"* ]]; then
     printf "\n Found Zh_\n\n"
     rangeString="z_norm_3b_e=0.1,4.0:z_norm_3b_mu=0.1,4.0:z_norm_4b_e=0.1,4.0:z_norm_4b_mu=0.1,4.0:tt_norm_e=0.1,4.0:tt_norm_mu=0.1,4.0"
  elif [[ $cards_dir == *"SM_2"* ]]; then
     printf "\n Found SM_2\n\n"
     rangeString="tt_norm_e=0.1,4.0:tt_norm_mu=0.1,4.0:w_norm_e=0.1,4:w_norm_mu=0.1,4:z_norm_3b_e=0.1,4.0:z_norm_3b_mu=0.1,4.0:z_norm_4b_e=0.1,4.0:z_norm_4b_mu=0.1,4.0"
  elif [[ $cards_dir == *"SM_all"* ]]; then
     printf "\n Found SM_all\n\n"
     rangeString="tt_norm_e_2016=0.1,4.0:tt_norm_mu_2016=0.1,4.0:w_norm_e_2016=0.1,4:w_norm_mu_2016=0.1,4:z_norm_3b_e_2016=0.1,4.0:z_norm_3b_mu_2016=0.1,4.0:z_norm_4b_e_2016=0.1,4.0:z_norm_4b_mu_2016=0.1,4.0:tt_norm_e_2017=0.1,4.0:tt_norm_mu_2017=0.1,4.0:w_norm_e_2017=0.1,4:w_norm_mu_2017=0.1,4:z_norm_3b_e_2017=0.1,4.0:z_norm_3b_mu_2017=0.1,4.0:z_norm_4b_e_2017=0.1,4.0:z_norm_4b_mu_2017=0.1,4.0:tt_norm_e_2018=0.1,4.0:tt_norm_mu_2018=0.1,4.0:w_norm_e_2018=0.1,4:w_norm_mu_2018=0.1,4:z_norm_3b_e_2018=0.1,4.0:z_norm_3b_mu_2018=0.1,4.0:z_norm_4b_e_2018=0.1,4.0:z_norm_4b_mu_2018=0.1,4.0"
  else
     printf "\n\n *** can't figure out what kind of cards directory this is:  %s\n\n" $cards_dir
     exit -1
  fi

  printf "\n\n rangeString = %s\n\n" $rangeString
   




  mkdir -p $impacts_dir

  cd $impacts_dir

  this_dir=`pwd`
  printf "\n\n this directory: %s\n\n" $this_dir

  f_name=$impacts_dir


  if [ "$step" = "1" ]
  then
     printf "\n will run initial fit\n\n"
     combineTool.py -M Impacts -d ../$cards_dir/0060/workspace.root -m 60 --doInitialFit --robustFit 1 --rMin -20.0  --setParameterRanges $rangeString    |& tee step1.log
  fi


  if [ "$step" = "2" ]
  then
     printf "\n will submit condor jobs\n\n"
     combineTool.py -M Impacts -d ../$cards_dir/0060/workspace.root -m 60 --doFits  --robustFit 1 --rMin -20.0  --setParameterRanges $rangeString  --job-mode condor |& tee step2.log
  fi


  if [ "$step" = "3" ]
  then
     printf "\n will make the plots pdf file\n\n"
     combineTool.py -M Impacts -d ../$cards_dir/0060/workspace.root -m 60 -o impacts.json --rMin -20.0 -t -1  |& tee step3.log
     plotImpacts.py -i impacts.json -o $f_name |& tee step4.log
  fi

  printf "\n\n Done.\n\n"

  cd -

