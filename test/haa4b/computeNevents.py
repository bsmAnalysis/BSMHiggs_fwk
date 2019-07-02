# The usage is: python computeNevents.py [dtag]

import ROOT as r
import glob
import sys
import commands

GREEN    = '\033[92m'
FAIL     = '\033[91m'
END      = '\033[0m'

inputdir = "results_2019_06_04"
#inputdir = "results_2019_06_04"
#dtag = sys.argv[1]
#dtags = ["MC13TeV_WJets_2017","MC13TeV_WJets_2017_ext1","MC13TeV_W1Jets_2017","MC13TeV_W2Jets_2017","MC13TeV_W3Jets_2017","MC13TeV_W4Jets_2017"]
#dtags = ["MC13TeV_DYJetsToLL_10to50_2017","MC13TeV_DYJetsToLL_M50_2017","MC13TeV_DYJetsToLL_M50_2017_ext1","MC13TeV_DY1JetsToLL_M50_2017","MC13TeV_DY1JetsToLL_M50_2017_ext1","MC13TeV_DY2JetsToLL_M50_2017_ext1","MC13TeV_DY2JetsToLL_M50_2017","MC13TeV_DY3JetsToLL_M50_2017","MC13TeV_DY3JetsToLL_M50_2017_ext1","MC13TeV_DY4JetsToLL_M50_2017","MC13TeV_WJets_2017","MC13TeV_WJets_2017_ext1","MC13TeV_W1Jets_2017","MC13TeV_W2Jets_2017","MC13TeV_W3Jets_2017","MC13TeV_W4Jets_2017"]
dtags = ["MC13TeV_DYJetsToLL_10to50_2016","MC13TeV_DY1JetsToLL_10to50_2016","MC13TeV_DY2JetsToLL_10to50_2016","MC13TeV_DY3JetsToLL_10to50_2016","MC13TeV_DY4JetsToLL_10to50_2016","MC13TeV_DYJetsToLL_50toInf_ext1_2016","MC13TeV_DYJetsToLL_50toInf_ext2_2016","MC13TeV_DY1JetsToLL_50toInf_2016","MC13TeV_DY2JetsToLL_50toInf_2016","MC13TeV_DY3JetsToLL_50toInf_2016","MC13TeV_DY4JetsToLL_50toInf_2016","MC13TeV_WJets_2016","MC13TeV_WJets_ext2_2016","MC13TeV_W1Jets_2016","MC13TeV_W2Jets_2016","MC13TeV_W2Jets_ext1_2016","MC13TeV_W3Jets_2016","MC13TeV_W3Jets_ext1_2016","MC13TeV_W4Jets_2016","MC13TeV_W4Jets_ext1_2016","MC13TeV_W4Jets_ext2_2016"]
#dtags = ["MC13TeV_DYJetsToLL_10to50_2017"]

who = commands.getstatusoutput('whoami')[1]
num=1
for tag in dtags:
  print(str(num)+". Processing the tag: "+tag)
  ntplpath = "/eos/cms/store/user/yuanc/"+inputdir+'/*/crab_'+tag+'*/*/*/'
  nFiles = 0
  nTot = 0
  nPos = 0
  nNeg = 0
  for file in glob.glob(ntplpath+'*.root'):
#    print((file))
    f = r.TFile(file)
    if (f.IsZombie() or (not f.IsOpen())):
      print FAIL + "Error: cannot open " + file + " or the file is not valid,please check if filename is valid!" + END
      continue
    nTot += f.Get("mainNtuplizer/nevents").GetBinContent(1)
    nPos += f.Get("mainNtuplizer/n_posevents").GetBinContent(1)
    nNeg += f.Get("mainNtuplizer/n_negevents").GetBinContent(1)
    f.Close()

  print("nTot: {0}, nPos: {1}, nNeg: {2}, nPos-nNeg: {3}\n".format(nTot,nPos,nNeg,nPos-nNeg))
  num += 1
