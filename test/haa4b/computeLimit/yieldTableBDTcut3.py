import ROOT as rt
import CMS_lumi, tdrstyle
from array import array
import os
import sys
import ctypes
import math

import argparse

#verbose = True
verbose = False

from ROOT import gROOT, gBenchmark, gRandom, gSystem, gStyle

import ROOT

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()

parser.add_argument( "--prefit", dest='use_prefit', help="Use the prefit histograms instead of the fit_b histograms", action='store_true' )
parser.add_argument( "limit_dir", type=str, help="Input directory for the mass point.  Example = new-results-2020-10-01/cards_SB13TeV_SM_Wh_2017_noSoftb/")

args = parser.parse_args()

limit_dir = args.limit_dir

use_prefit = False
if args.use_prefit:
   use_prefit = True


rrv = rt.RooRealVar()  # Get Wouter's print out of the way first


if use_prefit:
   print("\n\n Will use prefit histograms\n\n")
   dir1 = "shapes_prefit"
else:
   print("\n\n Will use fit_b histograms\n\n")
   dir1="shapes_fit_b"

print("\n\n dir1 = ", dir1, "\n\n")

chans = [ "A_SR_3b", "A_SR_4b" ]
chan_label = { "A_SR_3b" : "Signal region, 3b",
               "A_SR_4b" : "Signal region, 4b"
             }


a_mass_list = [ 20, 60 ]

proc_list = [ "otherbkg", "ddqcd", "zll", "wlnu", "ttbarbba", "ttbarcba", "ttbarlig", "total_background", "data" ]
proc_tex = { "otherbkg" : "Other Bkgs",
             "ddqcd"    : "DD QCD",
             "zll"      : "$Z \\to \\ell\\ell$",
             "wlnu"     : "$W \\to \\ell\\nu$",
             "ttbarbba" : "$t\\bar{t} + b\\bar{b}$",
             "ttbarcba" : "$t\\bar{t} + c\\bar{c}$",
             "ttbarlig" : "$t\\bar{t} + $ light ",
             "total_background" : "Total Bkg",
             "data"     : "Data"
           }

bin_group_list = [ {"first":1, "last":3, "blind":False}, {"first":4, "last":4, "blind":True}, {"first":5, "last":5, "blind":True} ]

if "2016" in limit_dir: year = "2016"
if "2017" in limit_dir: year = "2017"
if "2018" in limit_dir: year = "2018"

if ( "Wh" in limit_dir ):
   wz = "wh"
   wz_label = "Wh"
   lf_list = [ "e", "mu" ]
   lf_tex = { "e":"$e$",  "mu":"$\\mu$" }
   fd_file_dict = { "e": "fitDiagnostics.root", "mu": "fitDiagnostics.root" }
   #fd_file_dict = { "e": "fitDiagnosticsTest.root", "mu": "fitDiagnosticsTest.root" }
elif ( "Zh" in limit_dir ):
   wz = "zh"
   wz_label = "Zh"
   lf_list = [ "ee", "mumu" ]
   lf_tex = { "ee":"$ee$",  "mumu":"$\\mu\\mu$" }
   fd_file_dict = { "ee": "fitDiagnostics_e.root", "mumu": "fitDiagnostics_mu.root" }
else:
   print("\n\n *** Did not find Wh or Zh in limit_dir: {}.  Quitting.\n\n".format(limit_dir) )



for chan in chans:

   tex_file_name = "{}/yield-table-in-bdt-bins-{}.tex".format(limit_dir, chan)
   print("\n\n Opening output tex file: {}\n".format(tex_file_name) )
   tex_file = open( tex_file_name, "w" )

   tex_file.write("\\begin{table}\n")
   tex_file.write("\\begin{center}\n")
   #tex_file.write("\\resizebox{\\textwidth}{!}{\n")
   tex_file.write("\\small\n")
   tex_file.write("\\begin{tabular}{||l||")
   for lf in lf_list:
      for bin_group in bin_group_list:
         tex_file.write(">{\\collectcell{\\num[group-minimum-digits=4, group-separator = {,}]}}r<{\\endcollectcell}\n")
         tex_file.write("@{${}\pm{}$}\n")
         tex_file.write(">{\\collectcell{\\num[group-minimum-digits=4, group-separator = {,}]}}r<{\\endcollectcell}|\n")


      tex_file.write("|")
   tex_file.write("}\n")

   tex_file.write("\\hline\n")
   tex_file.write( "\\multicolumn{{13}}{{||c||}}{{ \\bf {} : {}, {} }} \\\\ \n".format( year, wz_label, chan_label[chan] ) )
   tex_file.write("\\hline\n")


   if (verbose): print("\n\n ====== {}".format(chan) )

   print("\n\n ====== {}".format(chan) )

   print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
   tex_file.write("\\hline\n")
   print(" Process             : ", end="")
   tex_file.write(" Process ")
   for lf in lf_list:
      for bin_group in bin_group_list:
         print("             {:4}[{},{}]     ".format(lf,bin_group["first"],bin_group["last"]),end="")
         if bin_group_list.index(bin_group) == (len(bin_group_list)-1):
            tex_file.write(" & \\multicolumn{{2}}{{c||}}{{ {:4} [{},{}] }} ".format(lf_tex[lf],bin_group["first"],bin_group["last"]))
         else:
            tex_file.write(" & \\multicolumn{{2}}{{c|}}{{ {:4} [{},{}] }} ".format(lf_tex[lf],bin_group["first"],bin_group["last"]))

      print("|", end="")
   print()
   tex_file.write(" \\\\ \n")
   print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
   tex_file.write("\\hline\n")


  #------ backgrounds and data
   for proc in proc_list:

      if proc == "data" :
         print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
         tex_file.write("\\hline\n")

      printline = "{0: <19} : ".format(proc)
      texline   = "{0: <23}   ".format(proc_tex[proc])


      for lf in lf_list:

         indent = "   "

         chan_name = "{}_{}".format( lf, chan )
         if ( verbose ): print("{} {}, {}, {} chan_name = {}".format( indent, chan, proc, lf, chan_name ) )

         chan_dir = "{}/{}".format( dir1, chan_name )
         if ( verbose ): print("{} {}, {}, {} chan_dir = {}".format( indent, chan, proc, lf, chan_dir ) )

         file_name = "{}/0060/{}".format(limit_dir, fd_file_dict[lf])
         if ( verbose ): print("{} {}, {}, {} file_name = {}".format( indent, chan, proc, lf, file_name ) )

         root_file = rt.TFile( file_name,"READ")

         hist_name = "{}/{}".format( chan_dir, proc )
         if ( verbose ): print("{} {}, {}, {}  hist_name = {}".format( indent, chan, proc, lf, hist_name ), end="" )

         hist = root_file.Get( hist_name )
         if ( verbose ) :
            if hist:
               print(" found")
            else:
               print(" *NOT* found")

         found = False

         if hist:

            found = True


          #-------- data

            if proc == "data":


               for bin_group in bin_group_list:


                  if ( not bin_group["blind"] ) :

                     if ( verbose ) : print( "   unblind   ")

                     if bin_group["first"] == bin_group["last"] :
                        val = hist.GetPointY( bin_group["first"]-1 )
                     else :
                        val = 0
                        i = bin_group["first"]
                        while ( i <= bin_group["last"]) :
                           val += hist.GetPointY( i-1 )
                           i += 1

                     printline += "  {:10,.0f}         ".format( val )
                     if bin_group_list.index(bin_group) == (len(bin_group_list)-1):
                        texline  += " & \\multicolumn{{2}}{{c||}}{{ \\num[group-minimum-digits=4, group-separator = {{,}}]{{ {:10.0f} }} }}   ".format( val )
                     else:
                        texline  += " & \\multicolumn{{2}}{{c|}}{{ \\num[group-minimum-digits=4, group-separator = {{,}}]{{ {:10.0f} }} }}   ".format( val )

                  else:
                     if ( verbose ) : print( "   blind   ")
                     printline += "                blind         "
                     if bin_group_list.index(bin_group) == (len(bin_group_list)-1):
                        texline  += " & \\multicolumn{2}{c||}{ }   "
                     else:
                        texline  += " & \\multicolumn{2}{c|}{ }   "


            else:

          #-------- backgrounds

               for bin_group in bin_group_list:

                  if ( verbose ) : print( " bin group [{},{}] ".format(bin_group["first"],bin_group["last"]), end="")

                  if bin_group["first"] == bin_group["last"] :

                     if ( verbose ) : print( " first and last same ", end="" )
                     val = hist.GetBinContent( bin_group["first"] )
                     err = hist.GetBinError( bin_group["first"] )

                  else:

                     if ( verbose ) : print( " first and last NOT same ", end="" )

                     val = 0.
                     i = bin_group["first"]
                     while ( i <= bin_group["last"] ):
                        val += hist.GetBinContent( i )
                        i += 1

                     cov_mat_name = "{}/total_covar".format( chan_dir )
                     cov_mat = root_file.Get( cov_mat_name )
                     if not cov_mat:
                        print("\n\n\n *** can't find {} in {}.  I quit.\n\n".format( cov_mat_name, file_name ) )
                        exit()

                     err2 = 0.
                     i = bin_group["first"]
                     while ( i <= bin_group["last"] ):
                        j = bin_group["first"]
                        while ( j <= bin_group["last"] ):
                           erri = hist.GetBinError( i )
                           errj = hist.GetBinError( j )
                           covij = cov_mat.GetBinContent( i, j )
                           covii = cov_mat.GetBinContent( i, i )
                           covjj = cov_mat.GetBinContent( j, j )
                           rho = 0
                           if covii * covjj > 0 :
                              rho = covij / math.sqrt( covii * covjj )
                           err2 += erri * errj * rho
                           j += 1
                        i += 1
                     err = math.sqrt( err2 )

                  if ( verbose ) : print( "{:.1f} +/- {:.1f}".format( val, err ) )

                  printline += "    {:10,.1f} +/- {:6,.1f}  ".format( val, err )
                  texline   += " &  {:10.1f} & {:6.1f}  ".format( val, err )


         else:

           #printline += " ------------------------------------------------------------------------------- "
            printline += "                                                                                 "
            texline   += "     & \\multicolumn{6}{c||}{ }                                                           "

         root_file.Close()

         printline += "|"

      if found :
         print(" {} ".format( printline ) )
         tex_file.write( "{} \\\\ \n".format( texline ) )

   print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
   tex_file.write("\\hline\n")

  #------ signals
   for a_mass in a_mass_list:

      printline = "signal({})          : ".format(a_mass)
      texline   = "signal({})     ".format(a_mass)

      for lf in lf_list:

         indent = "   "

         chan_name = "{}_{}".format( lf, chan )

         sig_file_name = "{}/00{}/haa4b_{}_13TeV_{}.root".format( limit_dir, a_mass, a_mass, wz )
         if ( verbose ): print("{} {}, {}, {} sig_file_name = {}".format( indent, chan, a_mass, lf, sig_file_name ) )

         sig_dir = chan_name
         if ( verbose ): print("{} {}, {}, {} sig_dir = {}".format( indent, chan, a_mass, lf, sig_dir ) )

         root_file = rt.TFile( sig_file_name,"READ")

         sig_hist_name = "{}/wh".format( sig_dir )

         hist = root_file.Get( sig_hist_name )

         if not hist:
            printf("\n\n *** missing {} in {}.  I quit.\n\n".format( sig_hist_name, sig_file_name ) )
            exit()


         for bin_group in bin_group_list:

            if ( verbose ) : print( " bin group [{},{}] ".format(bin_group["first"],bin_group["last"]), end="")

            if bin_group["first"] == bin_group["last"] :

               if ( verbose ) : print( " first and last same ", end="" )
               val = hist.GetBinContent( bin_group["first"] )
               err = hist.GetBinError( bin_group["first"] )

            else:

               if ( verbose ) : print( " first and last NOT same ", end="" )

               val = 0.
               err2 = 0.
               i = bin_group["first"]
               while ( i <= bin_group["last"] ):
                  val += hist.GetBinContent( i )
                  err2 += hist.GetBinError( i )
                  i += 1
               err = math.sqrt( err2 )

            if ( verbose ) : print( "{:.1f} +/- {:.1f}".format( val, err ) )

            printline += "    {:10,.1f} +/- {:6,.1f}  ".format( val, err )
            texline   += " &  {:10.1f} & {:6.1f}  ".format( val, err )

         root_file.Close()

         printline += "|"

      print( " {}".format( printline ) )
      tex_file.write( " {} \\\\ \n".format( texline ) )

   print("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
   tex_file.write("\\hline\n")

   tex_file.write("\\end{tabular}\n")
   #tex_file.write("}\n")
   tex_file.write("\\end{center}\n")
   tex_file.write("\\end{table}\n")


exit()


