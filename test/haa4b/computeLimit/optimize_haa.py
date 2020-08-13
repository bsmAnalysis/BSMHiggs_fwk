#!/usr/bin/env python
import math
import os,sys
#import json
import getopt
import commands
#import ROOT
#from ROOT import TFile, TGraph, TCanvas, TF1, TH1
sys.path.append('../../../scripts/')
import LaunchOnCondor

import argparse

CMSSW_BASE=os.environ.get('CMSSW_BASE')

def help() :
   print '\n\n\n'
   print ' -p phase (no default value is assigned)'
   print '\t 1 --> compute limit for all possible selection cuts  (in view of optimization)'
   print '\t 2 --> wrap-up and order the results either by best significance or best limit'
   print '\t 3 --> choose the best cuts for the selection and produce limits and control plots for this selection'
   print '\t 4.0 --> run combine jobs for final limit plot in WH channel (brazilian flag plot)'
   print '\t 4.1 --> run combine jobs for final limit plot in ZH channel (brazilian flag plot)'
   print '\t 4.2 --> run combine jobs for final limit plot in WH and ZH combined channel (brazilian flag plot)'
   print '\t 5.0 --> run combine jobs for final limit plot in WH channel combined Run II data'
   print '\t 5.1 --> run combine jobs for final limit plot in ZH channel combined Run II data'
   print '\t 5.2 --> run combine jobs for final limit plot in WH and ZH combined channel with Run II data'
   print '\t 6.0, 6.1, 6.2 to make limit plots for WH channel, ZH Channel, combined WH+ZH channels respectively'
   print '\t 7.0, 7.1, 7.2 to make limit plots for WH channel, ZH Channel, combined WH+ZH channels respectively, with combined Run II data'
   print '\nNote: CMSSW_BASE must be set when launching optimize.py (current values is: ' + CMSSW_BASE + ')\n' 
   print '\n\n\n'
   
parser = argparse.ArgumentParser()

parser.add_argument( "year_to_run", type=str, help="Year to run jobs for: 2016, 2017, 2018, all") ;
parser.add_argument( "phase", type=float, help="Walue for phase.  See above for how to set it." )
parser.add_argument( "--noSubmit", dest='no_submit', help="Create condor scripts but do not execute condor_submit", action='store_true' )
parser.add_argument( "--usage", help="Print usage", action='help' )
parser.add_argument( "-i", type=str, help="input plotter.root file" )
parser.add_argument( "-o", type=str, help="set CWD" )
parser.add_argument( "-j", type=str, help=" set jsonUrl" )

help()

args = parser.parse_args()

no_submit = False
if args.no_submit:
   no_submit = True
   print("\n no_submit set to True.  Will not call condor_submit.\n\n")

year_to_run = args.year_to_run
if ( not (year_to_run=="2016" or year_to_run=="2017" or year_to_run=="2018" or year_to_run=="all")):
   print "\n\n *** invalid year_to_run : " + year_to_run + "\n\n"
   help()
   sys.exit(-1)

phase = args.phase

import ROOT
from ROOT import TFile, TGraph, TCanvas, TF1, TH1

#default value

signalSuffixVec = []
OUTName         = []
LandSArgOptions = []
BIN             = []
MODEL           = []

LaunchOnCondor.Jobs_Queue='cmscaf1nd'
FarmDirectory  = "FARM"
JobName        = "computeLimits"
CWD=os.getcwd()
#phase=-1

autoMCstats = False
# https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination/1425/1.html
#thredMCstat = 0.001
thredMCstat = 0
BackExtrapol = " --BackExtrapol"
###################################################
##   VALUES TO BE EDITED BY THE USE ARE BELLOW   ##
###################################################

MODELS=["SM"] 
based_key="haa_mcbased" #mcbased_" #to run limits on MC use: haa_mcbased_, to use data driven obj use: haa_datadriven_

jsonPath=''
inUrl_wh=''
inUrl_zh=''

if ( year_to_run == "2016" ):
   jsonPath='$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/samples2016_legacy.json'
   inUrl_wh='$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_2016_Wh_Sys_noSoftb_forLimits.root'
   inUrl_zh='$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_2016_Zh_Sys_noSoftb_forLimits.root'

if ( year_to_run == "2017" ):
   jsonPath='$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/samples2017.json'
   inUrl_wh='$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_2017_WH_Sys_noSoftb_forLimits.root'
   inUrl_zh='$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_2017_ZH_Sys_noSoftb_forLimits.root'

if ( year_to_run == "2018" ):
   jsonPath='$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/samples2018.json'
   inUrl_wh='$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_2018_WH_Sys_noSoftb_forLimits.root'
   inUrl_zh='$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_2018_ZH_Sys_noSoftb_forLimits.root'

# configure your  forLimits files and json files below in order to have combined RunII data limits
jsonPaths=[
  '$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/samples2016_legacy.json',
  '$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/samples2017.json',
  '$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/samples2018.json'
]
inUrl_whs=[
  '$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_2016_Wh_Sys_noSoftb_forLimits.root',
  '$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_2017_WH_Sys_noSoftb_forLimits.root',
  '$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_2018_WH_Sys_noSoftb_forLimits.root'
]
inUrl_zhs=[
  '$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_2016_Zh_Sys_noSoftb_forLimits.root',
  '$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_2017_ZH_Sys_noSoftb_forLimits.root',
  '$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_2018_ZH_Sys_noSoftb_forLimits.root'
]
years = ["2016", "2017", "2018"]

BESTDISCOVERYOPTIM=True #Set to True for best discovery optimization, Set to False for best limit optimization
ASYMTOTICLIMIT=True #Set to True to compute asymptotic limits (faster) instead of toy based hybrid-new limits
#BINS = ["3b","4b","3b,4b"] # list individual analysis bins to consider as well as combined bins (separated with a coma but without space)
BINS = ["3b,4b"] # list individual analysis bins to consider as well as combined bins (separated with a coma but without space)
#BINS = ["3b+4b"] # list individual analysis bins to consider as well as combined bins (separated with a coma but without space)

MASS = [20, 25, 30, 40, 50, 60]
SUBMASS = [20, 25, 30, 40, 50, 60]

#MASS = [12, 15, 20, 25, 30, 40, 50, 60]
#SUBMASS = [12, 15, 20, 25, 30, 40, 50, 60]

#-- owen, july 24: replace --statBinByBin with --autoMCStats to use Combine implementation of bin-by-bin MC stat errs.
LandSArgCommonOptions=" --dropBckgBelow 0.01  --autoMCStats  "
LandSArgCommonOptions_2016wh=" --dropBckgBelow 0.015  --autoMCStats   "

#--- these were the last ones from Yuan
#LandSArgCommonOptions=" --dropBckgBelow 0.01 --statBinByBin 0.001 " #--BackExtrapol " #--statBinByBin 0.00001 "
#LandSArgCommonOptions_2016wh=" --dropBckgBelow 0.015 --statBinByBin 0.001 " #--BackExtrapol " #--statBinByBin 0.00001 "

#LandSArgCommonOptions=" --dropBckgBelow 0.01 " # --statBinByBin 0.001 " #--BackExtrapol " #--statBinByBin 0.00001 "
#LandSArgCommonOptions="  --BackExtrapol --statBinByBin 0.00001 --dropBckgBelow 0.00001 --blind"

for model in MODELS:
   for shape in ["bdt_shapes"]: #here run all the shapes you want to test.
      for bin in BINS:
         if(model=="SM"):
            suffix = "" 
            signalSuffixVec += [ suffix ]
            #OUTName         += ["SB13TeV_SM_Wh_2016_noSoftb"]
            #OUTName         += ["SB13TeV_SM_Wh_2017_noSoftb"]
            #OUTName         += ["SB13TeV_SM_Wh_2018_noSoftb"]
            #OUTName         += ["SB13TeV_SM_Wh_AllRun2_noSoftb"]
            OUTName += ["SB13TeV_SM_Wh_"+year_to_run+"_noSoftb"]
            LandSArgOptions += [" --histo " + shape + "  --systpostfix _13TeV --shape "]
            BIN             += [bin]
            MODEL           += [model]

###################################################
##   MAIN  CODE : MODIFY AT YOUR OWN RISK        ##
###################################################

#parse the options
#try:
#   # retrive command line options
#   shortopts  = "p:f:m:i:s:j:o:h:t?"
#   opts, args = getopt.getopt( sys.argv[1:], shortopts )
#except getopt.GetoptError:
#   # print help information and exit:
#   print "ERROR: unknown options in argument %s" % sys.argv[1:]
#   help()
#   sys.exit(1)
#
#for o,a in opts:
#   if o in("-?", "-h"):
#      help()
#      sys.exit(1)
#   elif o in('-i'): inUrl = a
#   elif o in('-p'): phase = float(a)
#   elif o in('-o'): CWD=a
#   elif o in('-j'): jsonUrl=a
#   elif o in('-t'): autoMCstats=True
      
if(phase<0 or len(CMSSW_BASE)==0):
   help()
   sys.exit(1)


#auxiliary function
def findCutIndex(cutsH, Gcut, m):
   for i in range(1, cutsH.GetXaxis().GetNbins()+1):
      passAllCuts=True
      for y in range(1, cutsH.GetYaxis().GetNbins()+1):
         if( (cutsH.GetYaxis().GetBinLabel(y).find("<")>=0 and cutsH.GetBinContent(i,y)>(Gcut[y-1].Eval(m,0,"")+0.000001)) or (cutsH.GetYaxis().GetBinLabel(y).find(">")>=0 and cutsH.GetBinContent(i,y)<(Gcut[y-1].Eval(m,0,"")-0.000001)) or  (cutsH.GetYaxis().GetBinLabel(y).find("<")==-1 and cutsH.GetYaxis().GetBinLabel(y).find(">")==-1 and cutsH.GetBinContent(i,y)<(Gcut[y-1].Eval(m,0,"")-0.000001)) ):
            passAllCuts=False;
            break;
      if(passAllCuts==True): return i;   
   return cutsH.GetXaxis().GetNbins()+1;

def findSideMassPoint(mass):
   global MASS
   LMass=0
   RMass=9999
   for m in MASS:
      if(m<=mass and m>=LMass):LMass=m
      if(m>=mass and m<=RMass):RMass=m
   return [LMass,RMass]

#######################
print("autoMCstats = {}".format(autoMCstats))
print("BackExtrapol = {}".format(BackExtrapol))
#Loop over all configurations
iConf = -1
cp = 0
for signalSuffix in signalSuffixVec : 
   iConf+=1;
   if(phase<=3 and ',' in BIN[iConf]):continue #only need individual bin for these phases

   jsonUrl = jsonPath + " --key " + based_key #+ MODEL[iConf]
   LandSArg = LandSArgCommonOptions + ' ' + LandSArgOptions[iConf];
   LandSArg_2016wh = LandSArgCommonOptions_2016wh + ' ' + LandSArgOptions[iConf];
   if(signalSuffix != ""):LandSArg+=' --signalSufix \"' + signalSuffix +'\" '
   binSuffix = ""
   if(',' not in BIN[iConf]):binSuffix="_"+ BIN[iConf]   

   ##########if(phase == 4.0): # owen
   if(phase == 4.0 or phase == 6.0):
      inUrl = inUrl_wh
   if(phase == 4.1 or phase == 6.1):
      OUTName[iConf] = OUTName[iConf].replace('_Wh','_Zh')
      inUrl = inUrl_zh
   if(phase == 4.2 or phase == 6.2):
      OUTName[iConf] = OUTName[iConf].replace('_Wh','')
      inUrl = inUrl_wh
   if(phase == 5.0 or phase == 7.0):
      OUTName[iConf] = OUTName[iConf] + "_Combined"
      inUrl = inUrl_whs[0]
   if(phase == 5.1 or phase == 7.1):
      OUTName[iConf] = OUTName[iConf].replace('_Wh','_Zh') + "_Combined"
      inUrl = inUrl_zhs[0]
   if(phase == 5.2 or phase == 7.2):
      OUTName[iConf] = OUTName[iConf].replace('_Wh','') + "_Combined"
      inUrl = inUrl_whs[0]
   DataCardsDir='cards_'+OUTName[iConf]+signalSuffix+binSuffix

   #prepare the output
   OUT = CWD+'/JOBS/'+OUTName[iConf]+signalSuffix+binSuffix+'/'
   os.system('mkdir -p ' + OUT)

   #get the cuts
   file = ROOT.TFile(inUrl)
   cutsH = file.Get('Wh (20)/all_optim_cut') 

   #Cp
   if "100.0" in signalSuffix:
   	cp = 100.0
   elif "10.0" in signalSuffix:
	cp = 10.0
   elif "5.0" in signalSuffix:
        cp = 5.0     
 
   ###################################################
   ##   OPTIMIZATION LOOP                           ##
   ###################################################

   if( phase == 1 ):
      print '# RUN LIMITS FOR ALL POSSIBLE CUTS  for ' + DataCardsDir + '#\n'
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+binSuffix+OUTName[iConf])

      FILE = open(OUT+"/LIST.txt","w")
      i = 1 # cut index 8 corresponds to BDT cut at ~0.1 (please check!)
      while (i<cutsH.GetXaxis().GetNbins()+1):                   
          shapeCutMin_ = -9999
          shapeCutMax_ = 9999
          SCRIPT = open(OUT+'script_'+str(i)+'_'+str(shapeCutMin_)+'_'+str(shapeCutMax_)+'.sh',"w")
          SCRIPT.writelines('cd ' + CMSSW_BASE + '/src;\n')
          SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc6_amd64_gcc491")+";\n")
          SCRIPT.writelines("eval `scram r -sh`;\n")
          SCRIPT.writelines('cd -;\n')     
          for j in range(0, 1): #always run 1points per jobs
             for m in MASS:
                cardsdir = 'H'+ str(m) + '_' + OUTName[iConf] + '_' + str(i);
                SCRIPT.writelines('mkdir -p ' + cardsdir+';\ncd ' + cardsdir+';\n')
                #SCRIPT.writelines("computeLimit --m " + str(m) + " --in " + inUrl + " --index " + str(i)     + " --json " + jsonUrl + " --shapeMin " + str(shapeCutMin_) + " --shapeMax " + str(shapeCutMax_) + " " + LandSArg + " --bins " + BIN[iConf] + " ;\n")
                SCRIPT.writelines("computeLimit --replaceHighSensitivityBinsWithBG  --m " + str(m) + " --in " + inUrl + " --index " + str(i)     + " --json " + jsonUrl + " --shapeMin " + str(shapeCutMin_) + " --shapeMax " + str(shapeCutMax_) + " " + LandSArg + " --bins " + BIN[iConf] + " ;\n")
#                SCRIPT.writelines("computeLimit --m " + str(m) + " --in " + inUrl + " --syst " + " --index " + str(i)     + " --json " + jsonUrl + " --shapeMin " + str(shapeCutMin_) + " --shapeMax " + str(shapeCutMax_) + " " + LandSArg + " --bins " + BIN[iConf] + " ;\n")
                SCRIPT.writelines("sh combineCards.sh;\n")
                #SCRIPT.writelines("combine -M AsymptoticLimits -m " +  str(m) + " --run expected card_combined.dat > LIMIT.log;\n") #limit computation
                SCRIPT.writelines("combine -M AsymptoticLimits  -v 2  -m " +  str(m) + " --run expected card_combined.dat > LIMIT.log;\n") #limit computation
                SCRIPT.writelines("combine -M Significance  -m " +  str(m) + " --significance -t -1 --expectSignal=1 card_combined.dat  > SIGN.log;\n") #apriori significance computation
                SCRIPT.writelines('tail -n 100 LIMIT.log > ' +OUT+str(m)+'_'+str(i)+'_'+str(shapeCutMin_)+'_'+str(shapeCutMax_)+'.log;\n')
                SCRIPT.writelines('tail -n 100 SIGN.log >> ' +OUT+str(m)+'_'+str(i)+'_'+str(shapeCutMin_)+'_'+str(shapeCutMax_)+'.log;\n')
                SCRIPT.writelines('cat LIMIT.log; cat SIGN.log\n')
                SCRIPT.writelines('cd ..;\n')
                SCRIPT.writelines('mv ' + cardsdir + ' ' + OUT + '/.\n')
                SCRIPT.writelines('rm -rd ' + cardsdir+';\n')            
          SCRIPT.close()
          LaunchOnCondor.SendCluster_Push(["BASH", 'sh ' + OUT+'script_'+str(i)+'_'+str(shapeCutMin_)+'_'+str(shapeCutMax_)+'.sh'])
          i = i+1#increment the cut index
      FILE.close()
      if ( not no_submit ):
         LaunchOnCondor.SendCluster_Submit()
      
   ###################################################
   ##   WRAPPING UP RESULTS                         ##
   ###################################################
   elif(phase == 2):
      print '# SCANNING ALL SETS OF CUTS  for ' + DataCardsDir + ' (you may want to go for a coffee...)#\n'
      fileName = OUT + "/OPTIM"+signalSuffix
      FILE = open(fileName+".txt","w")
      for m in MASS:
         print 'Starting mass ' + str(m)
         FILE.writelines("------------------------------------------------------------------------------------\n")
         BestLimit = []
         fileList = commands.getstatusoutput("find " + OUT +" -name " + str(m)+"_*.log")[1].split();           
         for f in fileList:
            try:
               value = -1.0;
               if(BESTDISCOVERYOPTIM==True):
                  exp = commands.getstatusoutput("cat " + f + " | grep \"Significance:\"")[1]; 
                  if(len(exp)<=0):continue
                  value = exp.split()[1]
               else:
                  exp = commands.getstatusoutput("cat " + f + " | grep \"Expected 50.0%\"")[1];  
                  if(len(exp)<=0):continue
                  value = exp.split()[4]

               if(value=='matches'):continue             
               if(float(value)<=0.0):continue
               if(BESTDISCOVERYOPTIM and float(value)>1000.0):continue #too high for a significance --> something is going wrong here
               f = f.replace(".log","")
               fields = f.split('_')
               N = len(fields)
               index = fields[N-3] 
               Cuts = ''
               for c in range(1, cutsH.GetYaxis().GetNbins()+1): 
                  Cuts += str(cutsH.GetBinContent(int(index),c)).rjust(7) + " ("+str(cutsH.GetYaxis().GetBinLabel(c))+")   "
               BestLimit.append("mH= "+str(m)+ " --> Limit= " + ('%10.3f' % float(value)) + "  Index: " + str(index)   + "  Cuts: " + Cuts + "   CutsOnShape: " + str(fields[N-2]).rjust(5) + " " + str(fields[N-1]).rjust(5))
            except:
               print "File %s does not contain a valid limit" % f

         #sort the limits for this mass
         BestLimit.sort(reverse=BESTDISCOVERYOPTIM)
         for s in BestLimit:
            FILE.writelines(s+"\n")

      #all done
      FILE.close()
      print("file "+fileName+".txt is written: it contains all selection points ordered by exp limit")

      os.system("root -l -b -q plotOptim.C+'(\"\",\""+fileName+".txt\",\"\", false, true, 13 , 35914.143 )'")  

   ###################################################
   ##   CHOSE BEST SELECTION CUTS and SAVE IT       ##
   ###################################################

   elif(phase == 3 ):
      LaunchOnCondor.Jobs_RunHere        = 1
      print '# CHOOSE BEST SELECTION CUTS  for ' + DataCardsDir + '#\n'
      Gcut  = []
      for c in range(1, cutsH.GetYaxis().GetNbins()+1):
         Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #add a graph for each cut
         Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMin
         Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMax

      fileName = OUT+"/OPTIM"+signalSuffix
      fileName+=".txt"
         
      mi=0
      for m in MASS:
         #if you want to display more than 3 options edit -m3 field
         cut_lines=commands.getstatusoutput("cat " + fileName + " | grep 'mH="+str(m)+"' -m20")[1].split('\n')
         if(len(cut_lines)<1):continue  #make sure that we get some lines
         if(len(cut_lines[0])<1):continue # make sure the line is not empty
         print 'mH='+str(m)+'\tOption \tR \tLimits and Cuts' 
         ictr=1
         for c in cut_lines: 
            cplit = c.split()
            print '\t #'+ str(ictr) + '\t' + cplit[2] + '\tIndex=' +  cplit[4] + '\t' + str(cplit[6:len(cplit)])
            ictr+=1
         print "Which option you want to keep?"
         #opt = int(raw_input(">"))-1   # remove prompting
         opt = 7 #always take the best cut

         #save cut chosen
         cutString = ''
         for c in range(1, cutsH.GetYaxis().GetNbins()+1):
            cutString += str(cutsH.GetYaxis().GetBinLabel(c)) + cut_lines[opt].split()[6+(c-1)*2] + '\t'
            Gcut[c-1].SetPoint(mi, m, float(cut_lines[opt].split()[6+(c-1)*2]) );
         cutString += cut_lines[opt].split()[6+(cutsH.GetYaxis().GetNbins()-1)*2 + 3 ] + '<shape<' + cut_lines[opt].split()[6+(cutsH.GetYaxis().GetNbins()-1)*2 + 4] 
         print cutString
         Gcut[cutsH.GetYaxis().GetNbins()+0].SetPoint(mi, m, float(cut_lines[opt].split()[6+(cutsH.GetYaxis().GetNbins()-1)*2 + 3 ]) );
         Gcut[cutsH.GetYaxis().GetNbins()+1].SetPoint(mi, m, float(cut_lines[opt].split()[6+(cutsH.GetYaxis().GetNbins()-1)*2 + 4 ]) );
         mi+=1

      #run limits for the cuts chosen (for intermediate masses use spline interpolation)
      print("carefully check that the cuts you have choosen are well listed below:\n")
      for m in SUBMASS:
           index = findCutIndex(cutsH, Gcut, m);
           Cuts = ''
           for c in range(1, cutsH.GetYaxis().GetNbins()+1):
              Cuts += str(cutsH.GetBinContent(int(index),c)).rjust(7) + " ("+str(cutsH.GetYaxis().GetBinLabel(c))+")   "

           print "M=%04i : Index=% 5i --> Cuts: %s"  % (m, index, Cuts)

      while True:
           #ans = raw_input('Use this fit and compute final limits? (y or n)\n')
           ans='y'
           if(ans=='y' or ans == 'Y'): break;
           else:			       sys.exit(0);           
      print 'YES'

      listcuts = open(OUT+'cuts.txt',"w")
      #LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+binSuffix+OUTName[iConf])
      for m in SUBMASS:
           index = findCutIndex(cutsH, Gcut, m);
           listcuts.writelines(str(m)+' '+str(index)+' ');
           for c in range(1, cutsH.GetYaxis().GetNbins()+3):
              listcuts.writelines(str(Gcut[c-1].Eval(m,0,""))+' ');
           listcuts.writelines('\n');
      listcuts.close();

  #------------------------------------------------------------------------------------------------

   elif(phase == 4.0 ):
      LaunchOnCondor.Jobs_RunHere        = 0
      print '# FINAL COMBINED LIMITS  for ' + DataCardsDir + '#\n'
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+binSuffix+OUTName[iConf])
      for m in SUBMASS:

           SideMasses = findSideMassPoint(m)
           indexString = ' '
           indexLString = ' '
           indexRString = ' '
           for bin in BIN[iConf].split(',') :
               Gcut  = []
               for c in range(1, cutsH.GetYaxis().GetNbins()+1):
                 Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #add a graph for each cut
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMin
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMax

               INbinSuffix = "_" + bin 
               IN = CWD+'/JOBS/'+OUTName[iConf]+signalSuffix+INbinSuffix+'/'
               try:
                  listcuts = open(IN+'cuts.txt',"r")
                  mi=0
                  for line in listcuts :
                     vals=line.split(' ')
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
   #                     Gcut[c-1].SetPoint(mi, float(vals[0]), float(vals[c+1]));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);
                  listcuts.close();          
               except:
                  mi=0
                  for mtmp in SUBMASS:
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);

               #add comma to index string if it is not empty
               if(indexString!=' '):
                  indexString+=','
                  if(not (SideMasses[0]==SideMasses[1])):
                     indexLString+=','
                     indexRString+=','

               #find the cut index for the current mass point
               indexString += str(findCutIndex(cutsH, Gcut, m));
               if(not (SideMasses[0]==SideMasses[1])):
                  indexLString = str(findCutIndex(cutsH, Gcut, SideMasses[0]));
                  indexRString = str(findCutIndex(cutsH, Gcut, SideMasses[1]));

           #print indexString

           cutStr = " "
           SideMassesArgs = ""
           if(not (SideMasses[0]==SideMasses[1])):
               SideMassesArgs += "--mL " + str(SideMasses[0]) + " --mR " + str(SideMasses[1]) + " --indexL " + indexLString +  " --indexR " + indexRString + " "

           SCRIPT = open(OUT+'/script_mass_'+str(m)+'.sh',"w")
           SCRIPT.writelines('cd ' + CMSSW_BASE + ';\n')
           SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc6_amd64_gcc491")+";\n")
           SCRIPT.writelines("eval `scram r -sh`;\n")
           SCRIPT.writelines('cd -;\n')

           cardsdir=DataCardsDir+"/"+('%04.0f' % float(m));
           SCRIPT.writelines('mkdir -p out_{};\ncd out_{};\n'.format(m,m)) #--blind instead of --simfit

          #--- first pass

           SCRIPT.writelines("computeLimit  --noCorrelatedStatUnc  --verbose  --replaceHighSensitivityBinsWithBG  --m " + str(m) + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " --modeDD --shape --subFake " + LandSArg + cutStr  +" --sumInputFile $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/all_plotter_Sys_noSoftb_forLimits.root >& cl-first.log\n")
           SCRIPT.writelines("sh combineCards_wh.sh;\n"); 


	   if autoMCstats:
	      SCRIPT.writelines("\nsed -i '$a*	   autoMCStats	   {}' card_combined_wh.dat".format(thredMCstat))
           SCRIPT.writelines("\ntext2workspace.py card_combined_wh.dat -o workspace.root --PO verbose --channel-masks  --PO \'ishaa\' --PO m=\'" + str(m) + "\' >& t2w-first.log \n")  
           #compute pvalue
           SCRIPT.writelines("combine -M Significance --signif --pvalue -m " +  str(m) + " workspace.root >& COMB-signif-first.log;\n")
           SCRIPT.writelines("combine -M FitDiagnostics workspace.root -m " +  str(m) + " -v 3     --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL  --rMin=0 --rMax=20 --stepSize=0.05 --robustFit 1  >& log-first.txt \n") 

           SCRIPT.writelines("mv fitDiagnostics.root fitDiagnostics-first.root\n\n\n")
           SCRIPT.writelines("mkdir datacards-wh-first-pass\n") ;
           SCRIPT.writelines("cp *.dat datacards-wh-first-pass\n\n") ;


          #--- second pass

           SCRIPT.writelines("computeLimit  --noCorrelatedStatUnc  --verbose  --replaceHighSensitivityBinsWithBG  --m " + str(m) + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " --modeDD --shape --subFake " + LandSArg + cutStr  +" --sumInputFile $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/all_plotter_Sys_noSoftb_forLimits.root  --fitDiagnosticsInputFile fitDiagnostics-first.root >& cl-second.log\n")
           SCRIPT.writelines("sh combineCards_wh.sh;\n"); 


	   if autoMCstats:
	      SCRIPT.writelines("\nsed -i '$a*	   autoMCStats	   {}' card_combined_wh.dat".format(thredMCstat))
           SCRIPT.writelines("\ntext2workspace.py card_combined_wh.dat -o workspace.root --PO verbose --channel-masks  --PO \'ishaa\' --PO m=\'" + str(m) + "\' >& t2w-second.log \n")  
           #compute pvalue
           SCRIPT.writelines("combine -M Significance --signif --pvalue -m " +  str(m) + " workspace.root >& COMB-signif-second.log;\n")
           SCRIPT.writelines("combine -M FitDiagnostics workspace.root -m " +  str(m) + " -v 3     --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL  --rMin=0 --rMax=20 --stepSize=0.05 --robustFit 1  >& log.txt \n") 



           # save likelihood fit info
           SCRIPT.writelines("python " + CMSSW_BASE + "/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/print.py -u fitDiagnostics.root > simfit_m"+ str(m)+".txt \n")
	   SCRIPT.writelines("python " + CMSSW_BASE + "/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -A -a fitDiagnostics.root -g Nuisance_CrossCheck.root >> simfit_m"+ str(m) +".txt\n")

           ### THIS IS FOR Asymptotic fit
           if(ASYMTOTICLIMIT==True):
              SCRIPT.writelines("tt_e=`cat simfit_m"+ str(m) +".txt | grep 'tt_norm_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("tt_mu=`cat simfit_m"+ str(m) +".txt | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
	      SCRIPT.writelines("v_e=`cat simfit_m"+ str(m) +".txt | grep 'w_norm_e' | awk '{print $4;}'`;\n")
	      SCRIPT.writelines("v_mu=`cat simfit_m"+ str(m) +".txt | grep 'w_norm_mu' | awk '{print $4;}'`;\n")
	      SCRIPT.writelines("combine -M AsymptoticLimits   -v 2  -m " +  str(m) + " workspace.root -t -1 --setParameters tt_norm_e=$tt_e,w_norm_e=$v_e,tt_norm_mu=$tt_mu,w_norm_mu=$v_mu >& COMB.log;\n")


           ### THIS is for toy (hybridNew) fit
           else:
              SCRIPT.writelines("combine -M AsymptoticLimits  -v 2  -m " +  str(m) + " workspace.root > COMB.log;\n") #first run assymptotic limit to get quickly the range of interest
              SCRIPT.writelines("rm higgsCombineTest.AsymptoticLimits*.root;\n")
              SCRIPT.writelines("RMIN=`cat COMB.log | grep 'Expected  2.5%' | awk '{print $5;}'`;\n") #get the low edge 2sigma band from the assymptotic --> will be used to know where to put points
              SCRIPT.writelines("RMAX=`cat COMB.log | grep 'Expected 97.5%' | awk '{print $5;}'`;\n") #get the high edge 2sigma band from the assymptotic --> will be used to know where to put points
              SCRIPT.writelines('echo "expected limit from the assymptotic in the 2sigma range [$RMIN, $RMAX]";\n')
              SCRIPT.writelines('RMIN=$(echo "$RMIN*0.5" | bc);\n'); #DIVIDE RMIN   BY 2 to make sure we are considering large space enough
              SCRIPT.writelines('RMAX=$(echo "$RMAX*3.0" | bc);\n'); #MULTIPLY RMAX BY 3 to make sure we are considering large space enough
              SCRIPT.writelines('echo "for the hybridNew, consider r to be in the range [$RMIN, $RMAX]";\n')
              SCRIPT.writelines("makeGridUsingCrab.py card_combined.dat $RMIN $RMAX -n 40 -m "+str(m)+" -o grid ;\n")
              SCRIPT.writelines("rm grid.root;\n")
              SCRIPT.writelines("sh grid.sh 1 16 &> /dev/null;\n")
              SCRIPT.writelines("rm higgsCombinegrid.HybridNew.*;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+";\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.025;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.160;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.500;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.840;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.975;\n")
              SCRIPT.writelines("hadd -f higgsCombineTest.HybridNewMerged.mH"+str(m)+".root  higgsCombineTest.HybridNew.mH"+str(m)+"*.root;\n")
              SCRIPT.writelines("rm higgsCombineTest.HybridNew.mH"+str(m)+"*.root;\n")

           SCRIPT.writelines('mkdir -p ' + CWD+'/'+cardsdir+';\n')
           SCRIPT.writelines('mv * ' + CWD+'/'+cardsdir+'/.;\n')
           SCRIPT.writelines('cd ..;\n\n') 
           SCRIPT.close()
           #os.system('sh ' + OUT+'script_mass_'+str(m)+'.sh ')  #uncomment this line to launch interactively (this may take a lot of time)
           LaunchOnCondor.SendCluster_Push(["BASH", 'sh ' + OUT+'script_mass_'+str(m)+'.sh'])
      if ( not no_submit ):
         LaunchOnCondor.SendCluster_Submit()

  #------------------------------------------------------------------------------------------------

   elif(phase == 4.1 ):
      LaunchOnCondor.Jobs_RunHere        = 0
      print '# FINAL COMBINED LIMITS  for ' + DataCardsDir + '#\n'
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+binSuffix+OUTName[iConf])
      for m in SUBMASS:
           SideMasses = findSideMassPoint(m)
           indexString = ' '
           indexLString = ' '
           indexRString = ' '
           for bin in BIN[iConf].split(',') :
               Gcut  = []
               for c in range(1, cutsH.GetYaxis().GetNbins()+1):
                 Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #add a graph for each cut
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMin
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMax

               INbinSuffix = "_" + bin 
               IN = CWD+'/JOBS/'+OUTName[iConf]+signalSuffix+INbinSuffix+'/'
               try:
                  listcuts = open(IN+'cuts.txt',"r")
                  mi=0
                  for line in listcuts :
                     vals=line.split(' ')
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
   #                     Gcut[c-1].SetPoint(mi, float(vals[0]), float(vals[c+1]));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);
                  listcuts.close();          
               except:
                  mi=0
                  for mtmp in SUBMASS:
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);

               #add comma to index string if it is not empty
               if(indexString!=' '):
                  indexString+=','
                  if(not (SideMasses[0]==SideMasses[1])):
                     indexLString+=','
                     indexRString+=','

               #find the cut index for the current mass point
               indexString += str(findCutIndex(cutsH, Gcut, m));
               if(not (SideMasses[0]==SideMasses[1])):
                  indexLString = str(findCutIndex(cutsH, Gcut, SideMasses[0]));
                  indexRString = str(findCutIndex(cutsH, Gcut, SideMasses[1]));

           #print indexString

           cutStr = " "
           SideMassesArgs = ""
           if(not (SideMasses[0]==SideMasses[1])):
               SideMassesArgs += "--mL " + str(SideMasses[0]) + " --mR " + str(SideMasses[1]) + " --indexL " + indexLString +  " --indexR " + indexRString + " "

           SCRIPT = open(OUT+'/script_mass_'+str(m)+'.sh',"w")
           SCRIPT.writelines('cd ' + CMSSW_BASE + ';\n')
           SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc6_amd64_gcc491")+";\n")
           SCRIPT.writelines("eval `scram r -sh`;\n")
           SCRIPT.writelines('cd -;\n')

           cardsdir=DataCardsDir+"/"+('%04.0f' % float(m));
           SCRIPT.writelines('mkdir -p out_{};\ncd out_{};\n'.format(m,m)) #--blind instead of --simfit
           SCRIPT.writelines("computeLimit   --noCorrelatedStatUnc  --verbose  --replaceHighSensitivityBinsWithBG  --runZh --m " + str(m) + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --shape --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " " + LandSArg + cutStr  +"  >& cl.log\n")
           SCRIPT.writelines("sh combineCards_zh.sh;\n"); 

	   if autoMCstats:
	      SCRIPT.writelines("\nsed -i '$a*	   autoMCStats	   {}' card_e_zh.dat".format(thredMCstat))
           SCRIPT.writelines("\ntext2workspace.py card_e_zh.dat -o workspace_e.root --PO verbose --channel-masks  --PO \'ishaa\' --PO m=\'" + str(m) + "\' >& t2w-e.log \n")  
           SCRIPT.writelines("combine -M Significance --signif --pvalue -m " +  str(m) + " workspace_e.root >& COMB-signif-e.log;\n")
           SCRIPT.writelines("combine -M FitDiagnostics workspace_e.root -m " +  str(m) + " -v 3     --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL  --rMin=0 --rMax=10 --stepSize=0.05 --robustFit 1  >& log_e.txt \n") 
           SCRIPT.writelines("python " + CMSSW_BASE + "/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/print.py -u fitDiagnostics.root > simfit_m"+ str(m)+"_e.txt \n")
	   SCRIPT.writelines("python " + CMSSW_BASE + "/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -A -a fitDiagnostics.root -g Nuisance_CrossCheck.root >> simfit_m"+ str(m) +"_e.txt\n")
	   SCRIPT.writelines("mv fitDiagnostics.root fitDiagnostics_e.root\n")

	   if autoMCstats:
	      SCRIPT.writelines("\nsed -i '$a*	   autoMCStats	   {}' card_mu_zh.dat".format(thredMCstat))
           SCRIPT.writelines("\ntext2workspace.py card_mu_zh.dat -o workspace_mu.root --PO verbose --channel-masks  --PO \'ishaa\' --PO m=\'" + str(m) + "\' >& t2w-mu.log \n")  
           SCRIPT.writelines("combine -M Significance --signif --pvalue -m " +  str(m) + " workspace_mu.root > COMB-signif-mu.log;\n")
           SCRIPT.writelines("combine -M FitDiagnostics workspace_mu.root -m " +  str(m) + " -v 3     --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL  --rMin=0 --rMax=10 --stepSize=0.05 --robustFit 1  >& log_mu.txt \n") 
           SCRIPT.writelines("python " + CMSSW_BASE + "/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/print.py -u fitDiagnostics.root > simfit_m"+ str(m)+"_mu.txt \n")
	   SCRIPT.writelines("python " + CMSSW_BASE + "/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -A -a fitDiagnostics.root -g Nuisance_CrossCheck.root >> simfit_m"+ str(m) +"_mu.txt\n")
	   SCRIPT.writelines("mv fitDiagnostics.root fitDiagnostics_mu.root\n")


	   if autoMCstats:
	      SCRIPT.writelines("\nsed -i '$a*	   autoMCStats	   {}' card_combined_zh.dat".format(thredMCstat))
           SCRIPT.writelines("\ntext2workspace.py card_combined_zh.dat -o workspace.root --PO verbose --channel-masks  --PO \'ishaa\' --PO m=\'" + str(m) + "\' >& t2w.log \n")  
           ### THIS IS FOR Asymptotic fit
           if(ASYMTOTICLIMIT==True):
           ### THIS is for toy (hybridNew) fit
              SCRIPT.writelines("tt_e=`cat simfit_m"+ str(m) +"_e.txt | grep 'tt_norm_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("tt_mu=`cat simfit_m"+ str(m) +"_mu.txt | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("v_3b_e=`cat simfit_m"+ str(m) +"_e.txt | grep 'z_norm_3b_e' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("v_4b_e=`cat simfit_m"+ str(m) +"_e.txt | grep 'z_norm_4b_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("v_3b_mu=`cat simfit_m"+ str(m) +"_mu.txt | grep 'z_norm_3b_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("v_4b_mu=`cat simfit_m"+ str(m) +"_mu.txt | grep 'z_norm_4b_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("combine -M AsymptoticLimits   -v 2  -m " +  str(m) + " workspace.root -t -1 --freezeParameters tt_norm_e,tt_norm_mu --setParameters tt_norm_e=$tt_e,z_norm_3b_e=$v_3b_e,z_norm_4b_e=$v_4b_e,tt_norm_mu=$tt_mu,z_norm_3b_mu=$v_3b_mu,z_norm_4b_mu=$v_4b_mu >& COMB.log;\n")  

           else:
              SCRIPT.writelines("combine -M AsymptoticLimits   -v 2  -m " +  str(m) + " workspace.root >& COMB.log;\n") #first run assymptotic limit to get quickly the range of interest
              SCRIPT.writelines("rm higgsCombineTest.AsymptoticLimits*.root;\n")
              SCRIPT.writelines("RMIN=`cat COMB.log | grep 'Expected  2.5%' | awk '{print $5;}'`;\n") #get the low edge 2sigma band from the assymptotic --> will be used to know where to put points
              SCRIPT.writelines("RMAX=`cat COMB.log | grep 'Expected 97.5%' | awk '{print $5;}'`;\n") #get the high edge 2sigma band from the assymptotic --> will be used to know where to put points
              SCRIPT.writelines('echo "expected limit from the assymptotic in the 2sigma range [$RMIN, $RMAX]";\n')
              SCRIPT.writelines('RMIN=$(echo "$RMIN*0.5" | bc);\n'); #DIVIDE RMIN   BY 2 to make sure we are considering large space enough
              SCRIPT.writelines('RMAX=$(echo "$RMAX*3.0" | bc);\n'); #MULTIPLY RMAX BY 3 to make sure we are considering large space enough
              SCRIPT.writelines('echo "for the hybridNew, consider r to be in the range [$RMIN, $RMAX]";\n')
              SCRIPT.writelines("makeGridUsingCrab.py card_combined_zh.dat $RMIN $RMAX -n 40 -m "+str(m)+" -o grid ;\n")
              SCRIPT.writelines("rm grid.root;\n")
              SCRIPT.writelines("sh grid.sh 1 16 &> /dev/null;\n")
              SCRIPT.writelines("rm higgsCombinegrid.HybridNew.*;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+";\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.025;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.160;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.500;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.840;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.975;\n")
              SCRIPT.writelines("hadd -f higgsCombineTest.HybridNewMerged.mH"+str(m)+".root  higgsCombineTest.HybridNew.mH"+str(m)+"*.root;\n")
              SCRIPT.writelines("rm higgsCombineTest.HybridNew.mH"+str(m)+"*.root;\n")

           SCRIPT.writelines('mkdir -p ' + CWD+'/'+cardsdir+';\n')
           SCRIPT.writelines('mv * ' + CWD+'/'+cardsdir+'/.;\n')
           SCRIPT.writelines('cd ..;\n\n') 
           SCRIPT.close()
           #os.system('sh ' + OUT+'script_mass_'+str(m)+'.sh ')  #uncomment this line to launch interactively (this may take a lot of time)
           LaunchOnCondor.SendCluster_Push(["BASH", 'sh ' + OUT+'script_mass_'+str(m)+'.sh'])
      if ( not no_submit ):
         LaunchOnCondor.SendCluster_Submit()

   #########################ZH + WH combined####################################

   elif(phase == 4.2 ):
      LaunchOnCondor.Jobs_RunHere        = 0
      print '# FINAL COMBINED LIMITS  for ' + DataCardsDir + '#\n'
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+binSuffix+OUTName[iConf])
      for m in SUBMASS:
           SideMasses = findSideMassPoint(m)
           indexString = ' '
           indexLString = ' '
           indexRString = ' '
           for bin in BIN[iConf].split(',') :
               Gcut  = []
               for c in range(1, cutsH.GetYaxis().GetNbins()+1):
                 Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #add a graph for each cut
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMin
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMax

               INbinSuffix = "_" + bin 
               IN = CWD+'/JOBS/'+OUTName[iConf]+signalSuffix+INbinSuffix+'/'
               try:
                  listcuts = open(IN+'cuts.txt',"r")
                  mi=0
                  for line in listcuts :
                     vals=line.split(' ')
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
   #                     Gcut[c-1].SetPoint(mi, float(vals[0]), float(vals[c+1]));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);
                  listcuts.close();          
               except:
                  mi=0
                  for mtmp in SUBMASS:
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);

               #add comma to index string if it is not empty
               if(indexString!=' '):
                  indexString+=','
                  if(not (SideMasses[0]==SideMasses[1])):
                     indexLString+=','
                     indexRString+=','

               #find the cut index for the current mass point
               indexString += str(findCutIndex(cutsH, Gcut, m));
               if(not (SideMasses[0]==SideMasses[1])):
                  indexLString = str(findCutIndex(cutsH, Gcut, SideMasses[0]));
                  indexRString = str(findCutIndex(cutsH, Gcut, SideMasses[1]));

           #print indexString

           cutStr = " "
           SideMassesArgs = ""
           if(not (SideMasses[0]==SideMasses[1])):
               SideMassesArgs += "--mL " + str(SideMasses[0]) + " --mR " + str(SideMasses[1]) + " --indexL " + indexLString +  " --indexR " + indexRString + " "

           SCRIPT = open(OUT+'/script_mass_'+str(m)+'.sh',"w")
           SCRIPT.writelines('cd ' + CMSSW_BASE + ';\n')
           SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc6_amd64_gcc491")+";\n")
           SCRIPT.writelines("eval `scram r -sh`;\n")
           SCRIPT.writelines('cd -;\n')

           cardsdir=DataCardsDir+"/"+('%04.0f' % float(m));
           SCRIPT.writelines('mkdir -p out_{};\ncd out_{};\n'.format(m,m)) #--blind instead of --simfit

          #--- first pass for Wh

           SCRIPT.writelines("computeLimit  --noCorrelatedStatUnc  --verbose  --replaceHighSensitivityBinsWithBG  --m " + str(m) + BackExtrapol + " --in " + inUrl_wh + " " + "--syst --simfit --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " --modeDD --shape --subFake " + LandSArg + cutStr  +  " --sumInputFile $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/all_plotter_Sys_noSoftb_forLimits.root >& cl-wh-first.log\n")
           SCRIPT.writelines("sh combineCards_wh.sh;\n"); 
           SCRIPT.writelines("text2workspace.py card_combined_wh.dat -o workspace-wh-first.root --PO verbose --channel-masks  --PO \'ishaa\' --PO m=\'" + str(m) + "\' >& t2w-wh-first.log \n")  
           SCRIPT.writelines("combine -M FitDiagnostics workspace-wh-first.root -m " +  str(m) + " -v 3     --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL  --rMin=0 --rMax=20 --stepSize=0.05 --robustFit 1  >& log-wh-first.txt \n") 
           SCRIPT.writelines("mv fitDiagnostics.root fitDiagnostics-wh-first.root\n\n\n")
           SCRIPT.writelines("mkdir datacards-wh-first-pass\n") ;
           SCRIPT.writelines("cp *.dat datacards-wh-first-pass\n\n") ;

          #--- second pass for Wh
           SCRIPT.writelines("computeLimit  --noCorrelatedStatUnc  --verbose  --replaceHighSensitivityBinsWithBG  --m " + str(m) + BackExtrapol + " --in " + inUrl_wh + " " + "--syst --simfit --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " --modeDD --shape --subFake " + LandSArg + cutStr  +" --sumInputFile $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/all_plotter_Sys_noSoftb_forLimits.root  --fitDiagnosticsInputFile fitDiagnostics-wh-first.root >& cl-wh-second.log\n")
           SCRIPT.writelines("sh combineCards_wh.sh;\n"); 



           SCRIPT.writelines("computeLimit  --noCorrelatedStatUnc  --verbose  --replaceHighSensitivityBinsWithBG  --runZh --m " + str(m) + BackExtrapol + " --in " + inUrl_zh + " " + "--syst --simfit --shape --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " " + LandSArg + cutStr  +" >& cl-zh.log\n")

           SCRIPT.writelines("sh combineCards_zh.sh;\n"); 
           SCRIPT.writelines("combineCards.py card_e_wh.dat card_e_zh.dat > card_e.dat\n"); 
           SCRIPT.writelines("combineCards.py card_mu_wh.dat card_mu_zh.dat > card_mu.dat\n"); 
           SCRIPT.writelines("combineCards.py card_combined_wh.dat card_combined_zh.dat > card_combined.dat\n"); 
           SCRIPT.writelines("sed -i '/rateParam ch2_emu_/d' card_combined.dat\n\n"); 

	   if autoMCstats:
	      SCRIPT.writelines("\nsed -i '$a*	   autoMCStats	   {}' card_e.dat".format(thredMCstat))
           SCRIPT.writelines("text2workspace.py card_e.dat -o workspace_e.root --PO verbose --channel-masks  --PO \'ishaa\' --PO m=\'" + str(m) + "\' >& t2w-e.log  \n")  
           SCRIPT.writelines("combine -M Significance --signif --pvalue -m " +  str(m) + " workspace_e.root >& COMB_e.log;\n")
           SCRIPT.writelines("combine -M FitDiagnostics workspace_e.root -m " +  str(m) + " -v 3      --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL  --rMin=0 --rMax=20 --stepSize=0.05 --robustFit 1  >& log_e.txt \n") 
           SCRIPT.writelines("python " + CMSSW_BASE + "/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/print.py -u fitDiagnostics.root > simfit_m"+ str(m)+"_e.txt \n")
	   SCRIPT.writelines("python " + CMSSW_BASE + "/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -A -a fitDiagnostics.root -g Nuisance_CrossCheck.root >> simfit_m"+ str(m) +"_e.txt\n")
	   SCRIPT.writelines("mv fitDiagnostics.root fitDiagnostics_e.root\n\n")

	   if autoMCstats:
	      SCRIPT.writelines("\nsed -i '$a*	   autoMCStats	   {}' card_mu.dat".format(thredMCstat))
           SCRIPT.writelines("text2workspace.py card_mu.dat -o workspace_mu.root --PO verbose --channel-masks  --PO \'ishaa\' --PO m=\'" + str(m) + "\' >& t2w-mu.log  \n")  
           SCRIPT.writelines("combine -M Significance --signif --pvalue -m " +  str(m) + " workspace_mu.root > COMB_mu.log;\n")
           SCRIPT.writelines("combine -M FitDiagnostics workspace_mu.root -m " +  str(m) + " -v 3     --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL  --rMin=0 --rMax=20 --stepSize=0.05 --robustFit 1  >& log_mu.txt \n") 
           SCRIPT.writelines("python " + CMSSW_BASE + "/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/print.py -u fitDiagnostics.root > simfit_m"+ str(m)+"_mu.txt \n")
	   SCRIPT.writelines("python " + CMSSW_BASE + "/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -A -a fitDiagnostics.root -g Nuisance_CrossCheck.root >> simfit_m"+ str(m) +"_mu.txt\n")
	   SCRIPT.writelines("mv fitDiagnostics.root fitDiagnostics_mu.root\n\n")

	   if autoMCstats:
	      SCRIPT.writelines("\nsed -i '$a*	   autoMCStats	   {}' card_combined.dat".format(thredMCstat))
           SCRIPT.writelines("text2workspace.py card_combined.dat -o workspace.root --PO verbose --channel-masks  --PO \'ishaa\' --PO m=\'" + str(m) + "\' >& t2w.log  \n")
           ### THIS IS FOR Asymptotic fit
           if(ASYMTOTICLIMIT==True):
           ### THIS is for toy (hybridNew) fit
              SCRIPT.writelines("tt_e=`cat simfit_m"+ str(m) +"_e.txt | grep 'tt_norm_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("z_3b_e=`cat simfit_m"+ str(m) +"_e.txt | grep 'z_norm_3b_e' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("z_4b_e=`cat simfit_m"+ str(m) +"_e.txt | grep 'z_norm_4b_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("w_e=`cat simfit_m"+ str(m) +"_e.txt | grep 'w_norm_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("tt_mu=`cat simfit_m"+ str(m) +"_mu.txt | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("z_3b_mu=`cat simfit_m"+ str(m) +"_mu.txt | grep 'z_norm_3b_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("z_4b_mu=`cat simfit_m"+ str(m) +"_mu.txt | grep 'z_norm_4b_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("w_mu=`cat simfit_m"+ str(m) +"_mu.txt | grep 'w_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("combine -M AsymptoticLimits  -v 2  -m " +  str(m) + " workspace.root -t -1 --freezeParameters tt_norm_e,tt_norm_mu --setParameters tt_norm_e=$tt_e,z_norm_3b_e=$z_3b_e,z_norm_4b_e=$z_4b_e,w_norm_e=$w_e,tt_norm_mu=$tt_mu,z_norm_3b_mu=$z_3b_mu,z_norm_4b_mu=$z_4b_mu,w_norm_mu=$w_mu >& COMB.log;\n")  

           else:
	      print("Do not support this mode!!!")
	      exit(-1)

           SCRIPT.writelines('mkdir -p ' + CWD+'/'+cardsdir+';\n')
           SCRIPT.writelines('mv * ' + CWD+'/'+cardsdir+'/.;\n')
           SCRIPT.writelines('cd ..;\n\n') 
           SCRIPT.close()
           #os.system('sh ' + OUT+'script_mass_'+str(m)+'.sh ')  #uncomment this line to launch interactively (this may take a lot of time)
           LaunchOnCondor.SendCluster_Push(["BASH", 'sh ' + OUT+'script_mass_'+str(m)+'.sh'])
      if ( not no_submit ):
         LaunchOnCondor.SendCluster_Submit()


   ##------------------------------------------------------------------------
   
   elif(phase == 5.0 ):
      LaunchOnCondor.Jobs_RunHere        = 0
      print '# FINAL COMBINED LIMITS  for ' + DataCardsDir + '#\n'
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+binSuffix+OUTName[iConf])
      for m in SUBMASS:
           SideMasses = findSideMassPoint(m)
           indexString = ' '
           indexLString = ' '
           indexRString = ' '
           for bin in BIN[iConf].split(',') :
               Gcut  = []
               for c in range(1, cutsH.GetYaxis().GetNbins()+1):
                 Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #add a graph for each cut
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMin
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMax

               INbinSuffix = "_" + bin 
               IN = CWD+'/JOBS/'+OUTName[iConf]+signalSuffix+INbinSuffix+'/'
               try:
                  listcuts = open(IN+'cuts.txt',"r")
                  mi=0
                  for line in listcuts :
                     vals=line.split(' ')
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
   #                     Gcut[c-1].SetPoint(mi, float(vals[0]), float(vals[c+1]));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);
                  listcuts.close();          
               except:
                  mi=0
                  for mtmp in SUBMASS:
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);

               #add comma to index string if it is not empty
               if(indexString!=' '):
                  indexString+=','
                  if(not (SideMasses[0]==SideMasses[1])):
                     indexLString+=','
                     indexRString+=','

               #find the cut index for the current mass point
               indexString += str(findCutIndex(cutsH, Gcut, m));
               if(not (SideMasses[0]==SideMasses[1])):
                  indexLString = str(findCutIndex(cutsH, Gcut, SideMasses[0]));
                  indexRString = str(findCutIndex(cutsH, Gcut, SideMasses[1]));

           #print indexString

	   if not len(inUrl_whs) == len(jsonPaths):
	      print("The number of json files is not equal to the number of plotter files!!")
	      exit(-1)
           cutStr = " "
           SideMassesArgs = ""
           if(not (SideMasses[0]==SideMasses[1])):
               SideMassesArgs += "--mL " + str(SideMasses[0]) + " --mR " + str(SideMasses[1]) + " --indexL " + indexLString +  " --indexR " + indexRString + " "

           SCRIPT = open(OUT+'/script_mass_'+str(m)+'.sh',"w")
           SCRIPT.writelines('cd ' + CMSSW_BASE + ';\n')
           SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc6_amd64_gcc491")+";\n")
           SCRIPT.writelines("eval `scram r -sh`;\n")
           SCRIPT.writelines('cd -;\n')

           cardsdir=DataCardsDir+"/"+('%04.0f' % float(m));
           SCRIPT.writelines('mkdir -p out_{};\ncd out_{};\n'.format(m,m)) #--blind instead of --simfit
	   
	   for i in range(0, len(jsonPaths)):
	      inUrl = inUrl_whs[i]
	      jsonUrl = jsonPaths[i] + " --key " + based_key
	      year = years[i]
	      datacard_e = "card_e_{}_wh.dat".format(year)
	      workspace_e = "workspace_{}_e.root".format(year)
	      datacard_mu = "card_mu_{}_wh.dat".format(year)
	      workspace_mu = "workspace_{}_mu.root".format(year)
	      datacard = "card_{}_wh.dat".format(year)
         
	      SCRIPT.writelines("\n#****************** {} *****************\n".format(year)) 
	      if "2016" in inUrl:
	          #SCRIPT.writelines("\ncomputeLimit --m " + str(m) + " --year " + year + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " --modeDD --shape --subFake " + LandSArg_2016wh + cutStr  +" ;\n")  
	          SCRIPT.writelines("\ncomputeLimit   --noCorrelatedStatUnc  --verbose  --replaceHighSensitivityBinsWithBG  --m " + str(m) + " --year " + year + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " --modeDD --shape --subFake " + LandSArg_2016wh + cutStr  +"   --sumInputFile $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/all_plotter_Sys_noSoftb_forLimits.root   >& cl-"+ year + ".log\n")  
	      else:
	          #SCRIPT.writelines("\ncomputeLimit --m " + str(m) + " --year " + year + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " --modeDD --shape --subFake " + LandSArg + cutStr  +" ;\n")  
	          SCRIPT.writelines("\ncomputeLimit   --noCorrelatedStatUnc  --verbose  --replaceHighSensitivityBinsWithBG  --m " + str(m) + " --year " + year + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " --modeDD --shape --subFake " + LandSArg + cutStr  +"   --sumInputFile $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/all_plotter_Sys_noSoftb_forLimits.root   >& cl-" + year + ".log\n")  
	      SCRIPT.writelines("sh combineCards_"+year+"_wh.sh;\n"); 

	   SCRIPT.writelines("\ncombineCards.py card_combined_2016_wh.dat card_combined_2017_wh.dat card_combined_2018_wh.dat > card_combined_wh.dat")
	   SCRIPT.writelines("\ntext2workspace.py card_combined_wh.dat -o workspace.root --PO verbose --channel-masks  --PO \'ishaa\' --PO m=\'" + str(m) + "\' >& t2w.log \n")  
	   SCRIPT.writelines("combine -M Significance --signif --pvalue -m " +  str(m) + " workspace.root > COMB-signif.log;\n\n\n")

	    
	   if(ASYMTOTICLIMIT==True):
           ### THIS is for toy (hybridNew) fit

              wh2016simfit  = CWD + "/cards_SB13TeV_SM_Wh_2016_noSoftb/00" + str(m) + "/simfit_m" + str(m) + ".txt"
              wh2017simfit  = CWD + "/cards_SB13TeV_SM_Wh_2017_noSoftb/00" + str(m) + "/simfit_m" + str(m) + ".txt"
              wh2018simfit  = CWD + "/cards_SB13TeV_SM_Wh_2018_noSoftb/00" + str(m) + "/simfit_m" + str(m) + ".txt"

              SCRIPT.writelines("if [ -f " + wh2016simfit + " ]; then\n")
              SCRIPT.writelines("   tt_e_2016=`cat "  + wh2016simfit + " | grep 'tt_norm_e'  | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   tt_mu_2016=`cat " + wh2016simfit + " | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   v_e_2016=`cat "   + wh2016simfit + " | grep 'w_norm_e'   | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   v_mu_2016=`cat "  + wh2016simfit + " | grep 'w_norm_mu'  | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find wh2016simfit file " + wh2016simfit + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2016,   tt_e = %s , tt_mu = %s , w_e = %s , w_mu = %s\\n\" $tt_e_2016, $tt_mu_2016, $v_e_2016, $v_mu_2016\n\n")

              SCRIPT.writelines("if [ -f " + wh2017simfit + " ]; then\n")
              SCRIPT.writelines("   tt_e_2017=`cat "  + wh2017simfit + " | grep 'tt_norm_e'  | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   tt_mu_2017=`cat " + wh2017simfit + " | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   v_e_2017=`cat "   + wh2017simfit + " | grep 'w_norm_e'   | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   v_mu_2017=`cat "  + wh2017simfit + " | grep 'w_norm_mu'  | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find wh2017simfit file " + wh2017simfit + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2017,   tt_e = %s , tt_mu = %s , w_e = %s , w_mu = %s\\n\" $tt_e_2017, $tt_mu_2017, $v_e_2017, $v_mu_2017\n\n")

              SCRIPT.writelines("if [ -f " + wh2018simfit + " ]; then\n")
              SCRIPT.writelines("   tt_e_2018=`cat "  + wh2018simfit + " | grep 'tt_norm_e'  | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   tt_mu_2018=`cat " + wh2018simfit + " | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   v_e_2018=`cat "   + wh2018simfit + " | grep 'w_norm_e'   | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   v_mu_2018=`cat "  + wh2018simfit + " | grep 'w_norm_mu'  | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find wh2018simfit file " + wh2018simfit + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2018,   tt_e = %s , tt_mu = %s , w_e = %s , w_mu = %s\\n\" $tt_e_2018, $tt_mu_2018, $v_e_2018, $v_mu_2018\n\n")

	      SCRIPT.writelines("combine -M AsymptoticLimits  -v 2  -m " +  str(m) + " workspace.root -t -1 --setParameters tt_norm_e_2016=$tt_e_2016,w_norm_e_2016=$v_e_2016,tt_norm_mu_2016=$tt_mu_2016,w_norm_mu_2016=$v_mu_2016,tt_norm_e_2017=$tt_e_2017,w_norm_e_2017=$v_e_2017,tt_norm_mu_2017=$tt_mu_2017,w_norm_mu_2017=$v_mu_2017,tt_norm_e_2018=$tt_e_2018,w_norm_e_2018=$v_e_2018,tt_norm_mu_2018=$tt_mu_2018,w_norm_mu_2018=$v_mu_2018 >& COMB.log;\n")

           else:
	      print("Do not support this mode!!!")
	      exit(-1)

           SCRIPT.writelines('\nmkdir -p ' + CWD+'/'+cardsdir+';\n')
           SCRIPT.writelines('mv * ' + CWD+'/'+cardsdir+'/.;\n')
           SCRIPT.writelines('cd ..;\n\n') 
           SCRIPT.close()
           #os.system('sh ' + OUT+'script_mass_'+str(m)+'.sh ')  #uncomment this line to launch interactively (this may take a lot of time)
           LaunchOnCondor.SendCluster_Push(["BASH", 'sh ' + OUT+'script_mass_'+str(m)+'.sh'])
      if ( not no_submit ):
         LaunchOnCondor.SendCluster_Submit()
   

  #------------------------------------------------------------------------------------------------

   elif(phase == 5.1 ):
      LaunchOnCondor.Jobs_RunHere        = 0
      print '# FINAL COMBINED LIMITS  for ' + DataCardsDir + '#\n'
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+binSuffix+OUTName[iConf])
      for m in SUBMASS:
           SideMasses = findSideMassPoint(m)
           indexString = ' '
           indexLString = ' '
           indexRString = ' '
           for bin in BIN[iConf].split(',') :
               Gcut  = []
               for c in range(1, cutsH.GetYaxis().GetNbins()+1):
                 Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #add a graph for each cut
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMin
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMax

               INbinSuffix = "_" + bin 
               IN = CWD+'/JOBS/'+OUTName[iConf]+signalSuffix+INbinSuffix+'/'
               try:
                  listcuts = open(IN+'cuts.txt',"r")
                  mi=0
                  for line in listcuts :
                     vals=line.split(' ')
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
   #                     Gcut[c-1].SetPoint(mi, float(vals[0]), float(vals[c+1]));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);
                  listcuts.close();          
               except:
                  mi=0
                  for mtmp in SUBMASS:
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);

               #add comma to index string if it is not empty
               if(indexString!=' '):
                  indexString+=','
                  if(not (SideMasses[0]==SideMasses[1])):
                     indexLString+=','
                     indexRString+=','

               #find the cut index for the current mass point
               indexString += str(findCutIndex(cutsH, Gcut, m));
               if(not (SideMasses[0]==SideMasses[1])):
                  indexLString = str(findCutIndex(cutsH, Gcut, SideMasses[0]));
                  indexRString = str(findCutIndex(cutsH, Gcut, SideMasses[1]));

           #print indexString

	   if not len(inUrl_whs) == len(jsonPaths):
	      print("The number of json files is not equal to the number of plotter files!!")
	      exit(-1)
           cutStr = " "
           SideMassesArgs = ""
           if(not (SideMasses[0]==SideMasses[1])):
               SideMassesArgs += "--mL " + str(SideMasses[0]) + " --mR " + str(SideMasses[1]) + " --indexL " + indexLString +  " --indexR " + indexRString + " "

           SCRIPT = open(OUT+'/script_mass_'+str(m)+'.sh',"w")
           SCRIPT.writelines('cd ' + CMSSW_BASE + ';\n')
           SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc6_amd64_gcc491")+";\n")
           SCRIPT.writelines("eval `scram r -sh`;\n")
           SCRIPT.writelines('cd -;\n')

           cardsdir=DataCardsDir+"/"+('%04.0f' % float(m));
           SCRIPT.writelines('mkdir -p out_{};\ncd out_{};\n'.format(m,m)) #--blind instead of --simfit
	   
	   for i in range(0, len(jsonPaths)):
	      inUrl = inUrl_zhs[i]
	      jsonUrl = jsonPaths[i] + " --key " + based_key
	      year = years[i]
	      datacard_e = "card_e_{}_zh.dat".format(year)
	      workspace_e = "workspace_{}_e.root".format(year)
	      datacard_mu = "card_mu_{}_zh.dat".format(year)
	      workspace_mu = "workspace_{}_mu.root".format(year)
	      datacard = "card_{}_zh.dat".format(year)
         
	      SCRIPT.writelines("\n#****************** {} *****************\n".format(year)) 
	      #SCRIPT.writelines("\ncomputeLimit --runZh --m " + str(m) + " --year " + year + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --shape --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " " + LandSArg + cutStr  +" ;\n")
	      SCRIPT.writelines("\ncomputeLimit  --noCorrelatedStatUnc  --replaceHighSensitivityBinsWithBG  --runZh --m " + str(m) + " --year " + year + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --shape --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " " + LandSArg + cutStr  +" >& cl-" + year + ".log \n")
	      SCRIPT.writelines("sh combineCards_{}_zh.sh;\n".format(year))
           
	   SCRIPT.writelines("\ncombineCards.py card_combined_2016_zh.dat card_combined_2017_zh.dat card_combined_2018_zh.dat > card_combined_zh.dat")
	   SCRIPT.writelines("\ntext2workspace.py card_combined_zh.dat -o workspace.root --PO verbose --channel-masks  --PO \'ishaa\' --PO m=\'" + str(m) + "\' >& t2w.log \n")  
	   SCRIPT.writelines("combine -M Significance --signif --pvalue -m " +  str(m) + " workspace.root >& COMB-signif.log;\n")

	    
	   if(ASYMTOTICLIMIT==True):
           ### THIS is for toy (hybridNew) fit

              zh2016simfit_e  = CWD + "/cards_SB13TeV_SM_Zh_2016_noSoftb/00" + str(m) + "/simfit_m" + str(m) + "_e.txt"
              zh2017simfit_e  = CWD + "/cards_SB13TeV_SM_Zh_2017_noSoftb/00" + str(m) + "/simfit_m" + str(m) + "_e.txt"
              zh2018simfit_e  = CWD + "/cards_SB13TeV_SM_Zh_2018_noSoftb/00" + str(m) + "/simfit_m" + str(m) + "_e.txt"
              zh2016simfit_mu = CWD + "/cards_SB13TeV_SM_Zh_2016_noSoftb/00" + str(m) + "/simfit_m" + str(m) + "_mu.txt"
              zh2017simfit_mu = CWD + "/cards_SB13TeV_SM_Zh_2017_noSoftb/00" + str(m) + "/simfit_m" + str(m) + "_mu.txt"
              zh2018simfit_mu = CWD + "/cards_SB13TeV_SM_Zh_2018_noSoftb/00" + str(m) + "/simfit_m" + str(m) + "_mu.txt"

              SCRIPT.writelines("if [ -f " + zh2016simfit_e + " ]; then\n")
              SCRIPT.writelines("   tt_e_2016=`cat " + zh2016simfit_e + " | grep 'tt_norm_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   v_3b_e_2016=`cat " + zh2016simfit_e + " | grep 'z_norm_3b_e' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   v_4b_e_2016=`cat " + zh2016simfit_e + " | grep 'z_norm_4b_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find zh2016simfit file " + zh2016simfit_e + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2016,  e:  tt = %s , v_3b = %s , v_4b = %s\\n\" $tt_e_2016, $v_3b_e_2016, $v_4b_e_2016\n\n")

              SCRIPT.writelines("if [ -f " + zh2016simfit_mu + " ]; then\n")
              SCRIPT.writelines("   tt_mu_2016=`cat " + zh2016simfit_mu + " | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   v_3b_mu_2016=`cat " + zh2016simfit_mu + " | grep 'z_norm_3b_mu' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   v_4b_mu_2016=`cat " + zh2016simfit_mu + " | grep 'z_norm_4b_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find zh2016simfit file " + zh2016simfit_mu + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2016, mu:  tt = %s , v_3b = %s , v_4b = %s\\n\" $tt_mu_2016, $v_3b_mu_2016, $v_4b_mu_2016\n\n")


              SCRIPT.writelines("if [ -f " + zh2017simfit_e + " ]; then\n")
              SCRIPT.writelines("   tt_e_2017=`cat " + zh2017simfit_e + " | grep 'tt_norm_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   v_3b_e_2017=`cat " + zh2017simfit_e + " | grep 'z_norm_3b_e' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   v_4b_e_2017=`cat " + zh2017simfit_e + " | grep 'z_norm_4b_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find zh2017simfit file " + zh2017simfit_e + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2017,  e:  tt = %s , v_3b = %s , v_4b = %s\\n\" $tt_e_2017, $v_3b_e_2017, $v_4b_e_2017\n\n")

              SCRIPT.writelines("if [ -f " + zh2017simfit_mu + " ]; then\n")
              SCRIPT.writelines("   tt_mu_2017=`cat " + zh2017simfit_mu + " | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   v_3b_mu_2017=`cat " + zh2017simfit_mu + " | grep 'z_norm_3b_mu' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   v_4b_mu_2017=`cat " + zh2017simfit_mu + " | grep 'z_norm_4b_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find zh2017simfit file " + zh2017simfit_mu + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2017, mu:  tt = %s , v_3b = %s , v_4b = %s\\n\" $tt_mu_2017, $v_3b_mu_2017, $v_4b_mu_2017\n\n")


              SCRIPT.writelines("if [ -f " + zh2018simfit_e + " ]; then\n")
              SCRIPT.writelines("   tt_e_2018=`cat " + zh2018simfit_e + " | grep 'tt_norm_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   v_3b_e_2018=`cat " + zh2018simfit_e + " | grep 'z_norm_3b_e' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   v_4b_e_2018=`cat " + zh2018simfit_e + " | grep 'z_norm_4b_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find zh2018simfit file " + zh2018simfit_e + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2018,  e:  tt = %s , v_3b = %s , v_4b = %s\\n\" $tt_e_2018, $v_3b_e_2018, $v_4b_e_2018\n\n")

              SCRIPT.writelines("if [ -f " + zh2018simfit_mu + " ]; then\n")
              SCRIPT.writelines("   tt_mu_2018=`cat " + zh2018simfit_mu + " | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   v_3b_mu_2018=`cat " + zh2018simfit_mu + " | grep 'z_norm_3b_mu' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   v_4b_mu_2018=`cat " + zh2018simfit_mu + " | grep 'z_norm_4b_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find zh2018simfit file " + zh2018simfit_mu + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2018, mu:  tt = %s , v_3b = %s , v_4b = %s\\n\" $tt_mu_2018, $v_3b_mu_2018, $v_4b_mu_2018\n\n")


	      SCRIPT.writelines("combine -M AsymptoticLimits  -v 2  -m " +  str(m) + " workspace.root -t -1 --freezeParameters tt_norm_e_2016,tt_norm_mu_2016,tt_norm_e_2017,tt_norm_mu_2017,tt_norm_e_2018,tt_norm_mu_2018 --setParameters tt_norm_e_2016=$tt_e_2016,z_norm_3b_e_2016=$v_3b_e_2016,z_norm_4b_e_2016=$v_4b_e_2016,tt_norm_mu_2016=$tt_mu_2016,z_norm_3b_mu_2016=$v_3b_mu_2016,z_norm_4b_mu_2016=$v_4b_mu_2016,tt_norm_e_2017=$tt_e_2017,z_norm_3b_e_2017=$v_3b_e_2017,z_norm_4b_e_2017=$v_4b_e_2017,tt_norm_mu_2017=$tt_mu_2017,z_norm_3b_mu_2017=$v_3b_mu_2017,z_norm_4b_mu_2017=$v_4b_mu_2017,tt_norm_e_2018=$tt_e_2018,z_norm_3b_e_2018=$v_3b_e_2018,z_norm_4b_e_2018=$v_4b_e_2018,tt_norm_mu_2018=$tt_mu_2018,z_norm_3b_mu_2018=$v_3b_mu_2018,z_norm_4b_mu_2018=$v_4b_mu_2018 >& COMB.log;\n") 

           else:
	      print("Do not support this mode!!!")
	      exit(-1)

           SCRIPT.writelines('\nmkdir -p ' + CWD+'/'+cardsdir+';\n')
           SCRIPT.writelines('mv * ' + CWD+'/'+cardsdir+'/.;\n')
           SCRIPT.writelines('cd ..;\n\n') 
           SCRIPT.close()
           #os.system('sh ' + OUT+'script_mass_'+str(m)+'.sh ')  #uncomment this line to launch interactively (this may take a lot of time)
           LaunchOnCondor.SendCluster_Push(["BASH", 'sh ' + OUT+'script_mass_'+str(m)+'.sh'])
      if ( not no_submit ):
         LaunchOnCondor.SendCluster_Submit()
   

  #------------------------------------------------------------------------------------------------

   elif(phase == 5.2 ):
      LaunchOnCondor.Jobs_RunHere        = 0
      print '# FINAL COMBINED LIMITS  for ' + DataCardsDir + '#\n'
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+binSuffix+OUTName[iConf])
      for m in SUBMASS:
           SideMasses = findSideMassPoint(m)
           indexString = ' '
           indexLString = ' '
           indexRString = ' '
           for bin in BIN[iConf].split(',') :
               Gcut  = []
               for c in range(1, cutsH.GetYaxis().GetNbins()+1):
                 Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #add a graph for each cut
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMin
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMax

               INbinSuffix = "_" + bin 
               IN = CWD+'/JOBS/'+OUTName[iConf]+signalSuffix+INbinSuffix+'/'
               try:
                  listcuts = open(IN+'cuts.txt',"r")
                  mi=0
                  for line in listcuts :
                     vals=line.split(' ')
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
   #                     Gcut[c-1].SetPoint(mi, float(vals[0]), float(vals[c+1]));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);
                  listcuts.close();          
               except:
                  mi=0
                  for mtmp in SUBMASS:
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);

               #add comma to index string if it is not empty
               if(indexString!=' '):
                  indexString+=','
                  if(not (SideMasses[0]==SideMasses[1])):
                     indexLString+=','
                     indexRString+=','

               #find the cut index for the current mass point
               indexString += str(findCutIndex(cutsH, Gcut, m));
               if(not (SideMasses[0]==SideMasses[1])):
                  indexLString = str(findCutIndex(cutsH, Gcut, SideMasses[0]));
                  indexRString = str(findCutIndex(cutsH, Gcut, SideMasses[1]));

           #print indexString

	   if not len(inUrl_whs) == len(jsonPaths):
	      print("The number of json files is not equal to the number of plotter files!!")
	      exit(-1)
           cutStr = " "
           SideMassesArgs = ""
           if(not (SideMasses[0]==SideMasses[1])):
               SideMassesArgs += "--mL " + str(SideMasses[0]) + " --mR " + str(SideMasses[1]) + " --indexL " + indexLString +  " --indexR " + indexRString + " "

           SCRIPT = open(OUT+'/script_mass_'+str(m)+'.sh',"w")
           SCRIPT.writelines('cd ' + CMSSW_BASE + ';\n')
           SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc6_amd64_gcc491")+";\n")
           SCRIPT.writelines("eval `scram r -sh`;\n")
           SCRIPT.writelines('cd -;\n')

           cardsdir=DataCardsDir+"/"+('%04.0f' % float(m));
           SCRIPT.writelines('mkdir -p out_{};\ncd out_{};\n'.format(m,m)) #--blind instead of --simfit
	   
	   for i in range(0, len(jsonPaths)):
	      inUrl = inUrl_whs[i]
	      jsonUrl = jsonPaths[i] + " --key " + based_key
	      year = years[i]
         
	      SCRIPT.writelines("\n#****************** {} Wh *****************\n".format(year)) 
	      if "2016" in inUrl:
	          #SCRIPT.writelines("\n#computeLimit --m " + str(m) + " --year " + year + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " --modeDD --shape --subFake " + LandSArg_2016wh + cutStr  +" ;\n")  
	          SCRIPT.writelines("\ncomputeLimit  --noCorrelatedStatUnc  --verbose  --replaceHighSensitivityBinsWithBG  --m " + str(m) + " --year " + year + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " --modeDD --shape --subFake " + LandSArg_2016wh + cutStr  +"  --sumInputFile $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/all_plotter_Sys_noSoftb_forLimits.root >& cl-wh-" + year + ".log \n")  
	      else:
	          #SCRIPT.writelines("\n#computeLimit --m " + str(m) + " --year " + year + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " --modeDD --shape --subFake " + LandSArg + cutStr  +" ;\n")  
	          SCRIPT.writelines("\ncomputeLimit  --noCorrelatedStatUnc  --verbose  --replaceHighSensitivityBinsWithBG  --m " + str(m) + " --year " + year + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " --modeDD --shape --subFake " + LandSArg + cutStr  +"  --sumInputFile $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/all_plotter_Sys_noSoftb_forLimits.root >& cl-wh-" + year + ".log \n")  
	      SCRIPT.writelines("sh combineCards_"+year+"_wh.sh;\n"); 
	      
	      inUrl = inUrl_zhs[i]
	      jsonUrl = jsonPaths[i] + " --key " + based_key
	      year = years[i]
         
	      SCRIPT.writelines("\n#****************** {} Zh *****************\n".format(year)) 
	      #SCRIPT.writelines("\n#computeLimit --runZh --m " + str(m) + " --year " + year + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --shape --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " " + LandSArg + cutStr  +" ;\n")
	      SCRIPT.writelines("\ncomputeLimit  --noCorrelatedStatUnc  --verbose  --replaceHighSensitivityBinsWithBG  --runZh --m " + str(m) + " --year " + year + BackExtrapol + " --in " + inUrl + " " + "--syst --simfit --shape --index 1 --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " " + LandSArg + cutStr  +" >& cl-zh-" + year + ".log ;\n")
	      SCRIPT.writelines("sh combineCards_{}_zh.sh;\n".format(year))

	      SCRIPT.writelines("combineCards.py card_combined_{}_wh.dat card_combined_{}_zh.dat > card_combined_{}.dat\n".format(year, year, year))

	   SCRIPT.writelines("\ncombineCards.py card_combined_2016.dat card_combined_2017.dat card_combined_2018.dat > card_combined.dat")
	   SCRIPT.writelines("\ntext2workspace.py card_combined.dat -o workspace.root --PO verbose --channel-masks  --PO \'ishaa\' --PO m=\'" + str(m) + "\' >& t2w.log \n")  
	   SCRIPT.writelines("combine -M Significance --signif --pvalue -m " +  str(m) + " workspace.root >& COMB-signif.log;\n")

	    
	   if(ASYMTOTICLIMIT==True):
           ### THIS is for toy (hybridNew) fit
	      # put normalizations below	      

              wh2016simfit  = CWD + "/cards_SB13TeV_SM_Wh_2016_noSoftb/00" + str(m) + "/simfit_m" + str(m) + ".txt"
              wh2017simfit  = CWD + "/cards_SB13TeV_SM_Wh_2017_noSoftb/00" + str(m) + "/simfit_m" + str(m) + ".txt"
              wh2018simfit  = CWD + "/cards_SB13TeV_SM_Wh_2018_noSoftb/00" + str(m) + "/simfit_m" + str(m) + ".txt"

              SCRIPT.writelines("if [ -f " + wh2016simfit + " ]; then\n")
              SCRIPT.writelines("   tt_e_2016=`cat "  + wh2016simfit + " | grep 'tt_norm_e'  | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   tt_mu_2016=`cat " + wh2016simfit + " | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   w_e_2016=`cat "   + wh2016simfit + " | grep 'w_norm_e'   | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   w_mu_2016=`cat "  + wh2016simfit + " | grep 'w_norm_mu'  | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find wh2016simfit file " + wh2016simfit + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2016,   tt_e = %s , tt_mu = %s , w_e = %s , w_mu = %s\\n\" $tt_e_2016, $tt_mu_2016, $w_e_2016, $w_mu_2016\n\n")

              SCRIPT.writelines("if [ -f " + wh2017simfit + " ]; then\n")
              SCRIPT.writelines("   tt_e_2017=`cat "  + wh2017simfit + " | grep 'tt_norm_e'  | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   tt_mu_2017=`cat " + wh2017simfit + " | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   w_e_2017=`cat "   + wh2017simfit + " | grep 'w_norm_e'   | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   w_mu_2017=`cat "  + wh2017simfit + " | grep 'w_norm_mu'  | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find wh2017simfit file " + wh2017simfit + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2017,   tt_e = %s , tt_mu = %s , w_e = %s , w_mu = %s\\n\" $tt_e_2017, $tt_mu_2017, $w_e_2017, $w_mu_2017\n\n")

              SCRIPT.writelines("if [ -f " + wh2018simfit + " ]; then\n")
              SCRIPT.writelines("   tt_e_2018=`cat "  + wh2018simfit + " | grep 'tt_norm_e'  | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   tt_mu_2018=`cat " + wh2018simfit + " | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   w_e_2018=`cat "   + wh2018simfit + " | grep 'w_norm_e'   | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   w_mu_2018=`cat "  + wh2018simfit + " | grep 'w_norm_mu'  | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find wh2018simfit file " + wh2018simfit + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2018,   tt_e = %s , tt_mu = %s , w_e = %s , w_mu = %s\\n\" $tt_e_2018, $tt_mu_2018, $w_e_2018, $w_mu_2018\n\n")


              zh2016simfit_e  = CWD + "/cards_SB13TeV_SM_Zh_2016_noSoftb/00" + str(m) + "/simfit_m" + str(m) + "_e.txt"
              zh2017simfit_e  = CWD + "/cards_SB13TeV_SM_Zh_2017_noSoftb/00" + str(m) + "/simfit_m" + str(m) + "_e.txt"
              zh2018simfit_e  = CWD + "/cards_SB13TeV_SM_Zh_2018_noSoftb/00" + str(m) + "/simfit_m" + str(m) + "_e.txt"
              zh2016simfit_mu = CWD + "/cards_SB13TeV_SM_Zh_2016_noSoftb/00" + str(m) + "/simfit_m" + str(m) + "_mu.txt"
              zh2017simfit_mu = CWD + "/cards_SB13TeV_SM_Zh_2017_noSoftb/00" + str(m) + "/simfit_m" + str(m) + "_mu.txt"
              zh2018simfit_mu = CWD + "/cards_SB13TeV_SM_Zh_2018_noSoftb/00" + str(m) + "/simfit_m" + str(m) + "_mu.txt"

              SCRIPT.writelines("if [ -f " + zh2016simfit_e + " ]; then\n")
              SCRIPT.writelines("   tt_e_2016=`cat " + zh2016simfit_e + " | grep 'tt_norm_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   z_3b_e_2016=`cat " + zh2016simfit_e + " | grep 'z_norm_3b_e' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   z_4b_e_2016=`cat " + zh2016simfit_e + " | grep 'z_norm_4b_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find zh2016simfit file " + zh2016simfit_e + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2016,  e:  tt = %s , z_3b = %s , z_4b = %s\\n\" $tt_e_2016, $z_3b_e_2016, $z_4b_e_2016\n\n")

              SCRIPT.writelines("if [ -f " + zh2016simfit_mu + " ]; then\n")
              SCRIPT.writelines("   tt_mu_2016=`cat " + zh2016simfit_mu + " | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   z_3b_mu_2016=`cat " + zh2016simfit_mu + " | grep 'z_norm_3b_mu' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   z_4b_mu_2016=`cat " + zh2016simfit_mu + " | grep 'z_norm_4b_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find zh2016simfit file " + zh2016simfit_mu + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2016, mu:  tt = %s , z_3b = %s , z_4b = %s\\n\" $tt_mu_2016, $z_3b_mu_2016, $z_4b_mu_2016\n\n")


              SCRIPT.writelines("if [ -f " + zh2017simfit_e + " ]; then\n")
              SCRIPT.writelines("   tt_e_2017=`cat " + zh2017simfit_e + " | grep 'tt_norm_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   z_3b_e_2017=`cat " + zh2017simfit_e + " | grep 'z_norm_3b_e' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   z_4b_e_2017=`cat " + zh2017simfit_e + " | grep 'z_norm_4b_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find zh2017simfit file " + zh2017simfit_e + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2017,  e:  tt = %s , z_3b = %s , z_4b = %s\\n\" $tt_e_2017, $z_3b_e_2017, $z_4b_e_2017\n\n")

              SCRIPT.writelines("if [ -f " + zh2017simfit_mu + " ]; then\n")
              SCRIPT.writelines("   tt_mu_2017=`cat " + zh2017simfit_mu + " | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   z_3b_mu_2017=`cat " + zh2017simfit_mu + " | grep 'z_norm_3b_mu' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   z_4b_mu_2017=`cat " + zh2017simfit_mu + " | grep 'z_norm_4b_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find zh2017simfit file " + zh2017simfit_mu + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2017, mu:  tt = %s , z_3b = %s , z_4b = %s\\n\" $tt_mu_2017, $z_3b_mu_2017, $z_4b_mu_2017\n\n")


              SCRIPT.writelines("if [ -f " + zh2018simfit_e + " ]; then\n")
              SCRIPT.writelines("   tt_e_2018=`cat " + zh2018simfit_e + " | grep 'tt_norm_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   z_3b_e_2018=`cat " + zh2018simfit_e + " | grep 'z_norm_3b_e' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   z_4b_e_2018=`cat " + zh2018simfit_e + " | grep 'z_norm_4b_e' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find zh2018simfit file " + zh2018simfit_e + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2018,  e:  tt = %s , z_3b = %s , z_4b = %s\\n\" $tt_e_2018, $z_3b_e_2018, $z_4b_e_2018\n\n")

              SCRIPT.writelines("if [ -f " + zh2018simfit_mu + " ]; then\n")
              SCRIPT.writelines("   tt_mu_2018=`cat " + zh2018simfit_mu + " | grep 'tt_norm_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("   z_3b_mu_2018=`cat " + zh2018simfit_mu + " | grep 'z_norm_3b_mu' | awk '{print $4;}'`;\n") 
              SCRIPT.writelines("   z_4b_mu_2018=`cat " + zh2018simfit_mu + " | grep 'z_norm_4b_mu' | awk '{print $4;}'`;\n")
              SCRIPT.writelines("else\n")
              SCRIPT.writelines("   printf \" cant find zh2018simfit file " + zh2018simfit_mu + ".  quitting\\n\\n \" \n")
              SCRIPT.writelines("   exit -1\n")
              SCRIPT.writelines("fi\n")
              SCRIPT.writelines("printf \" 2018, mu:  tt = %s , z_3b = %s , z_4b = %s\\n\" $tt_mu_2018, $z_3b_mu_2018, $z_4b_mu_2018\n\n")


	      SCRIPT.writelines("combine -M AsymptoticLimits  -v 2  -m " +  str(m) + " workspace.root -t -1 --freezeParameters tt_norm_e_2016,tt_norm_mu_2016,tt_norm_e_2017,tt_norm_mu_2017,tt_norm_e_2018,tt_norm_mu_2018 --setParameters tt_norm_e_2016=$tt_e_2016,z_norm_3b_e_2016=$z_3b_e_2016,z_norm_4b_e_2016=$z_4b_e_2016,w_norm_e_2016=$w_e_2016,tt_norm_mu_2016=$tt_mu_2016,z_norm_3b_mu_2016=$z_3b_mu_2016,z_norm_4b_mu_2016=$z_4b_mu_2016,w_norm_mu_2016=$w_mu_2016,tt_norm_e_2017=$tt_e_2017,z_norm_3b_e_2017=$z_3b_e_2017,z_norm_4b_e_2017=$z_4b_e_2017,w_norm_e_2017=$w_e_2017,tt_norm_mu_2017=$tt_mu_2017,z_norm_3b_mu_2017=$z_3b_mu_2017,z_norm_4b_mu_2017=$z_4b_mu_2017,w_norm_mu_2017=$w_mu_2017,tt_norm_e_2018=$tt_e_2018,z_norm_3b_e_2018=$z_3b_e_2018,z_norm_4b_e_2018=$z_4b_e_2018,w_norm_e_2018=$w_e_2018,tt_norm_mu_2018=$tt_mu_2018,z_norm_3b_mu_2018=$z_3b_mu_2018,z_norm_4b_mu_2018=$z_4b_mu_2018,w_norm_mu_2018=$w_mu_2018  >& COMB.log;\n") 

           else:
	      print("Do not support this mode!!!")
	      exit(-1)

           SCRIPT.writelines('\nmkdir -p ' + CWD+'/'+cardsdir+';\n')
           SCRIPT.writelines('mv * ' + CWD+'/'+cardsdir+'/.;\n')
           SCRIPT.writelines('cd ..;\n\n') 
           SCRIPT.close()
           #os.system('sh ' + OUT+'script_mass_'+str(m)+'.sh ')  #uncomment this line to launch interactively (this may take a lot of time)
           LaunchOnCondor.SendCluster_Push(["BASH", 'sh ' + OUT+'script_mass_'+str(m)+'.sh'])
      if ( not no_submit ):
         LaunchOnCondor.SendCluster_Submit()
   ######################################################################
   
   elif(phase == 6.0 or phase == 6.1 or phase == 6.2 or phase == 7.0 or phase == 7.1 or phase == 7.2):
      print '# FINAL PLOT for ' + DataCardsDir + '#\n'
#      os.system("hadd -f "+DataCardsDir+"/PValueTree.root "+DataCardsDir+"/*/higgsCombineTest.ProfileLikelihood.*.root > /dev/null")
      os.system("hadd -f "+DataCardsDir+"/PValueTree.root "+DataCardsDir+"/*/higgsCombineTest.Significance.*.root > /dev/null")

      #THIS IS FOR ASYMPTOTIC
      if(ASYMTOTICLIMIT==True):
         #os.system("hadd -f "+DataCardsDir+"/LimitTree.root "+DataCardsDir+"/*/higgsCombineTest.Asymptotic.*.root > /dev/null")
         os.system("hadd -f "+DataCardsDir+"/LimitTree.root "+DataCardsDir+"/*/higgsCombineTest.AsymptoticLimits.*.root > /dev/null")
      #THIS IS FOR HYBRIDNEW
      else:
         os.system("hadd -f "+DataCardsDir+"/LimitTree.root "+DataCardsDir+"/*/higgsCombineTest.HybridNewMerged.*.root > /dev/null")

      integrated_luminosity = 0
      if year_to_run == "2016":
         integrated_luminosity =  35914.143
      if year_to_run == "2017":
         integrated_luminosity =  41529.152
      if year_to_run == "2018":
         integrated_luminosity =  59740.565
      if year_to_run == "all":
         integrated_luminosity = 137183.86

      os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Strength_\",\""+DataCardsDir+"/LimitTree.root\",\"\", false, true, 13 , "+integrated_luminosity+" )'")
      os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Strength_\",\""+DataCardsDir+"/LimitTree.root\",\"\", false, true, 13 , "+integrated_luminosity+" , \"Wh channels\" ,false)'")


#      os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Strength_\",\""+DataCardsDir+"/LimitTree.root\",\"\", false, true, 13 , 35914.143 )'")
#      os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Strength_\",\""+DataCardsDir+"/LimitTree.root\",\"\", false, true, 13 , 35914.143 , \"Wh channels\" ,false)'")
#      os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Strength_\",\""+DataCardsDir+"/LimitTree.root\",\"\", false, true, 13 , 41529.152 )'")
#      os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Strength_\",\""+DataCardsDir+"/LimitTree.root\",\"\", false, true, 13 , 41529.152 , \"Wh channels\" ,false)'")
#      os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Strength_\",\""+DataCardsDir+"/LimitTree.root\",\"\", false, true, 13 , 59740.565 )'")
#      os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Strength_\",\""+DataCardsDir+"/LimitTree.root\",\"\", false, true, 13 , 59740.565 , \"Wh channels\" ,false)'")
#      os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Strength_\",\""+DataCardsDir+"/LimitTree.root\",\"\", false, true, 13 , 137183.86 )'")
#      os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Strength_\",\""+DataCardsDir+"/LimitTree.root\",\"\", false, true, 13 , 137183.86 , \"Wh channels\" ,false)'")

   ######################################################################

if(phase>8):
      help()

