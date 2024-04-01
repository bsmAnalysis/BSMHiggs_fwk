
import ROOT
ROOT.gROOT.SetBatch(True)

years = [ "2016", "2017", "2018" ]
amasses = [ "0020", "0025", "0030", "0040", "0050", "0060" ]
leps = [ "e", "mu" ]

outfile = ROOT.TFile.Open( "all-roo-fit-results.root", "recreate" )

#--- wh

for yr in years:

   fname = "cards_SB13TeV_SM_Wh_" + yr + "_noSoftb/0040/fitDiagnostics.root"
   print  "wh : yr = ", yr, " file = ", fname  

   infile = ROOT.TFile.Open( fname, "read" )

   fit_b = infile.Get("fit_b")

   if fit_b != None:

       outfile.cd()

       rfr_name = "fit_b_wh_" + yr
       fit_b.Write( rfr_name )

   for amass in amasses:

      fname = "cards_SB13TeV_SM_Wh_" + yr + "_noSoftb/" + amass + "/fitDiagnostics.root"
      print  "wh : yr = ", yr, " amass = ", amass, " file = ", fname  

      infile = ROOT.TFile.Open( fname, "read" )

      fit_s = infile.Get("fit_s")

      if fit_s != None:

          outfile.cd()

          rfr_name = "fit_s_wh_" + yr + "_" + amass
          fit_s.Write( rfr_name )

print "\n\n"



#--- zh

for lep in leps:

   for yr in years:

      fname = "cards_SB13TeV_SM_Zh_" + yr + "_noSoftb/0040/fitDiagnostics_" + lep + ".root"
      print  "zh " + lep + "  : yr = ", yr, " file = ", fname  

      infile = ROOT.TFile.Open( fname, "read" )

      fit_b = infile.Get("fit_b")

      if fit_b != None:

          outfile.cd()

          rfr_name = "fit_b_zh_" + lep + "_" + yr
          fit_b.Write( rfr_name )

      for amass in amasses:

         fname = "cards_SB13TeV_SM_Zh_" + yr + "_noSoftb/" + amass + "/fitDiagnostics_" + lep + ".root"
         print  "zh " + lep + "  : yr = ", yr, " amass = ", amass, " file = ", fname  

         infile = ROOT.TFile.Open( fname, "read" )

         fit_s = infile.Get("fit_s")

         if fit_s != None:

             outfile.cd()

             rfr_name = "fit_s_zh_" + lep + "_" + yr + "_" + amass
             fit_s.Write( rfr_name )



outfile.Close()


print "\n\n"


