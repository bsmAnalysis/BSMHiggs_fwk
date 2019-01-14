#!/usr/bin/env bash

DIRS=(`find ./ -maxdepth 1 -name "results*" -type d`)
echo ${DIRS[@]}

if [ ${#DIRS[@]} -eq 0 ];  
then
  echo "Cannot find results* directories, nothing to do. Exiting..."
  exit 0
fi

if [ ${#DIRS[@]} -eq 1 ]
then
  BASE=${DIRS[0]}
else
  read -p "Find multiple results directories, which one you want to kill? Please Indicate with number: " choice
  while true
  do
    if [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -le ${#DIRS[@]} ]
    then
      break
    else
      read -p "Your input is illegal, either contains non-numeric characters or number is larger than the total number of directories listed above. Please try again: " choice
    fi
  done
  BASE=${DIRS[${choice}-1]}
fi

DIR=""$BASE"/FARM/inputs"

echo "Will kill all submitted jobs under: "$DIR""
echo "Please confirm [N/y]"
read answer

if [[ $answer == "y" ]];
then
  echo "Killing jobs..."
  for job in "$DIR"/crab*
  do
    if [[ -d $job ]];then
      echo " crab kill --dir="$job" "
      crab kill --dir=${job}
#      crab purge --dir=${job}
    fi
  done
fi

read -p "Please check above no errors happened. If not, do you want to delete the whole directory: "$BASE"? [N/y]: " answer
if [[ $answer == "y" ]];
then
  echo "Deleting the direcory..."
  `rm -rf "$BASE"`
fi

echo "All Done, bye"
