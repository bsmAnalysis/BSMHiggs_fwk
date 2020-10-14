

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




  mkdir -p $impacts_dir

  cd $impacts_dir

  f_name=$impacts_dir


  if [ "$step" = "1" ]
  then
     printf "\n will run initial fit\n\n"
     combineTool.py -M Impacts -d ../$cards_dir/0060/workspace.root -m 60 --doInitialFit --robustFit 1 --rMin -1.0  |& tee step1.log
  fi


  if [ "$step" = "2" ]
  then
     printf "\n will submit condor jobs\n\n"
#    combineTool.py -M Impacts -d ../$cards_dir/0060/workspace.root -m 60 --robustFit 1 --doFits --rMin -1.0 --job-mode condor --dry-run 
     combineTool.py -M Impacts -d ../$cards_dir/0060/workspace.root -m 60 --robustFit 1 --doFits --rMin -1.0 --job-mode condor |& tee step2.log
  fi


  if [ "$step" = "3" ]
  then
     printf "\n will make the plots pdf file\n\n"
     combineTool.py -M Impacts -d ../$cards_dir/0060/workspace.root -m 60 -o impacts.json --rMin -1.0  |& tee step3.log
     plotImpacts.py -i impacts.json -o $f_name |& tee step4.log
  fi

  printf "\n\n Done.\n\n"

  cd -

