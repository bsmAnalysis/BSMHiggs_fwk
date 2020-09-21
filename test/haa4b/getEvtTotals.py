# The usage is: python computeNevents.py [dtag]

import ROOT as r
import glob
import sys
import commands
import os

GREEN    = '\033[92m'
FAIL     = '\033[91m'
END      = '\033[0m'

#dtags = [
#"MC13TeV_Wh_amass12_2016",
#"MC13TeV_Wh_amass15_2016",
#"MC13TeV_Wh_amass20_2016",
#"MC13TeV_Wh_amass25_2016",
#"MC13TeV_Wh_amass30_2016",
#"MC13TeV_Wh_amass40_2016",
#"MC13TeV_Wh_amass50_2016",
#"MC13TeV_Wh_amass60_2016",
#"MC13TeV_Zh_amass12_2016",
#"MC13TeV_Zh_amass15_2016",
#"MC13TeV_Zh_amass20_2016",
#"MC13TeV_Zh_amass25_2016",
#"MC13TeV_Zh_amass30_2016",
#"MC13TeV_Zh_amass40_2016",
#"MC13TeV_Zh_amass50_2016",
#"MC13TeV_Zh_amass60_2016",
#]

adtfile = open("mc-dtags-2018-2.txt","r")
dtags = adtfile.readlines()

printFiles = True
listF = []

os.system("mkdir -p evt-total-files2")

who = commands.getstatusoutput('whoami')[1]
num=1
for tag in dtags:
  tag = tag.rstrip("\n")
  print(str(num)+". Processing the tag: "+tag)
 #--- some of 2018 ntuples
  ntplpath = "/eos/cms/store/user/yuanc/results_2020_02_05/*/crab_"+tag+'*/*/*/'
 #--- other 2018 ntuples
  #ntplpath = "/eos/user/y/yuanc/backup_2018Analysis/*/crab_"+tag+'*/*/*/'
  nFiles = 0
  nTot = 0
  nPos = 0
  nNeg = 0
  outfile = "evt-total-files2/" + tag + ".txt"
  listF.append(outfile)
  with open(outfile,"w+") as _f:
    for filename in glob.glob(ntplpath+'*.root'):
      f = r.TFile(filename)
      if (f.IsZombie() or (not f.IsOpen())):
        print FAIL + "Error: cannot open " + filename + " or the file is not valid,please check if filename is valid!" + END
        continue
      this_file_nTot = f.Get("mainNtuplizer/nevents").GetBinContent(1)
      this_file_nPos = f.Get("mainNtuplizer/n_posevents").GetBinContent(1)
      this_file_nNeg = f.Get("mainNtuplizer/n_negevents").GetBinContent(1)
      nTot += this_file_nTot
      nPos += this_file_nPos
      nNeg += this_file_nNeg
      _f.write( "{:200}  {:15}  {:15}  {:15}\n".format(filename, this_file_nTot, this_file_nPos, this_file_nNeg) )
      f.Close()
    _f.write("Totals {:15}  {:15}  {:15}\n".format(nTot, nPos, nNeg))

  print("nTot: {0:15}, nPos: {1:15}, nNeg: {2:15}, nPos-nNeg: {3:15}\n".format(nTot,nPos,nNeg,nPos-nNeg))
  num += 1

  
