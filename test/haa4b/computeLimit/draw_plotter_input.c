
#include "utils.c"


   void draw_plotter_input( 
           int nb = 3,
           const char* infile = "../limits-work/plotter-files-2022-03-08/all-ZH.root",
           int sig_mass = 60
         ) {

      gDirectory -> Delete("h*") ;

      char hist_name[100] ;
      char bg_hist_name[100] ;
      char sig_hist_name[100] ;



      bool is_wh(false) ;
      bool is_zh(false) ;
      TString tsif( infile ) ;
      if ( tsif.Contains("WH") ) is_wh = true ;
      if ( tsif.Contains("ZH") ) is_zh = true ;

      if ( !is_wh && !is_zh ) { printf("\n\n *** can't tell from input file name if this is WH or ZH.\n\n") ; return ; }

      if ( is_wh && is_zh ) { printf("\n\n *** can't tell from input file name if this is WH or ZH.\n\n") ; return ; }

      if ( is_wh ) printf("\n\n Input file is WH : %s\n\n", infile ) ;
      if ( is_zh ) printf("\n\n Input file is ZH : %s\n\n", infile ) ;

      char wz_string[10] ;
      if ( is_wh ) sprintf( wz_string, "WH" ) ;
      if ( is_zh ) sprintf( wz_string, "ZH" ) ;

      char lepton_flav[2][10] ;
      char bg_components[10][100] ;
      char bg_shortname[10][20] ;

      int nbg ;

      if ( is_zh ) {

         sprintf( lepton_flav[0], "ee" ) ;
         sprintf( lepton_flav[1], "mumu" ) ;

         sprintf( bg_components[0], "Other Bkgds" ) ;
         sprintf( bg_components[1], "t#bar{t} + b#bar{b}_filt5" ) ;
         sprintf( bg_components[2], "t#bar{t} + c#bar{c}_filt4" ) ;
         sprintf( bg_components[3], "t#bar{t} + light_filt1" ) ;
         sprintf( bg_components[4], "Z#rightarrow ll" ) ;

         sprintf( bg_shortname[0], "other" ) ;
         sprintf( bg_shortname[1], "ttbb" ) ;
         sprintf( bg_shortname[2], "ttcc" ) ;
         sprintf( bg_shortname[3], "ttlf" ) ;
         sprintf( bg_shortname[4], "zll" ) ;

         nbg = 5 ;

      }

      if ( is_wh ) {

         sprintf( lepton_flav[0], "e" ) ;
         sprintf( lepton_flav[1], "mu" ) ;

         sprintf( bg_components[0], "Other Bkgds" ) ;
         sprintf( bg_components[1], "t#bar{t} + b#bar{b}_filt5" ) ;
         sprintf( bg_components[2], "t#bar{t} + c#bar{c}_filt4" ) ;
         sprintf( bg_components[3], "t#bar{t} + light_filt1" ) ;
         sprintf( bg_components[4], "Z#rightarrow ll" ) ;
         sprintf( bg_components[5], "W#rightarrow l#nu" ) ;

         sprintf( bg_shortname[0], "other" ) ;
         sprintf( bg_shortname[1], "ttbb" ) ;
         sprintf( bg_shortname[2], "ttcc" ) ;
         sprintf( bg_shortname[3], "ttlf" ) ;
         sprintf( bg_shortname[4], "zll" ) ;
         sprintf( bg_shortname[5], "wlnu" ) ;

         nbg = 6 ;

      }



      TFile* tf_in = new TFile( infile ,"READ" ) ;
      if ( tf_in == 0x0 ) { printf("\n\n *** bad input file : %s\n\n", infile ) ; return ; }
      if ( !(tf_in -> IsOpen()) ) { printf("\n\n *** problem opening input file : %s\n\n", infile ) ; return ; }

      tf_in -> pwd() ;
      gDirectory -> pwd() ;

      TCanvas* can1 = get_canvas( "can1", "", 50, 50, 900, 900 ) ;
      can1 -> Clear() ;
      can1 -> Divide(3,3) ;
      int ci = 1 ;

      TH1F* h_bg_sum = 0x0 ;

      for ( int bgi=0; bgi<nbg; bgi++ ) {

         tf_in -> cd("") ;
         printf("\n\n Changing to directory %s\n", bg_components[bgi] ) ;
         tf_in -> cd( bg_components[bgi] ) ;

         TH1F* h_comp_lf_sum = 0x0 ;

         for ( int lfi=0; lfi<2; lfi++ ) {

            sprintf( hist_name, "%s_A_SR_%db_bdt", lepton_flav[lfi], nb ) ;
            TH1F* h_bg = (TH1F*) gDirectory -> Get( hist_name ) ;
            if ( h_bg == 0x0 ) {
               printf("\n\n *** %s not found.\n", hist_name ) ;
               continue ;
            } else {
               printf("  %s : %p\n", hist_name, h_bg ) ;
            }
            if ( h_comp_lf_sum == 0x0 ) {
               sprintf( bg_hist_name, "h_A_SR_%db_bdt_%s", nb, bg_shortname[bgi] ) ;
               h_comp_lf_sum = (TH1F*) h_bg -> Clone( bg_hist_name ) ;
            } else {
               h_comp_lf_sum -> Add( h_bg ) ;
            }
            if ( h_bg_sum == 0x0 ) {
               h_bg_sum = (TH1F*) h_bg -> Clone( "h_bg_sum" ) ;
            } else {
               h_bg_sum -> Add( h_bg ) ;
            }

         } // lfi

         can1 -> cd( ci++ ) ;
         if ( h_comp_lf_sum == 0x0 ) continue ;
         h_comp_lf_sum -> Draw() ;

      } // bgi

      can1 -> cd( ci++ ) ;
      h_bg_sum -> Draw() ;

      tf_in -> cd("") ;
      char sig_dir[100] ;
      sprintf( sig_dir, "Wh (%d)", sig_mass ) ;
      tf_in -> cd( sig_dir ) ;
      TH1F* h_sig_lf_sum = 0x0 ;
      for ( int lfi=0; lfi<2; lfi++ ) {
         TH1F* h_sig = get_hist( hist_name ) ;
         if ( h_sig_lf_sum == 0x0 ) {
            sprintf( sig_hist_name, "h_sig" ) ;
            h_sig_lf_sum = (TH1F*) h_sig -> Clone( sig_hist_name ) ;
         } else {
            h_sig_lf_sum -> Add( h_sig ) ;
         }
      }

      h_sig_lf_sum -> Draw("same") ;

      can1 -> cd( ci++ ) ;
      h_sig_lf_sum -> Draw() ;
      h_bg_sum -> Draw("same") ;



      char out_filename[1000] ;
      sprintf( out_filename, "binning-optimization-input-%s-%db-ma%d.root", wz_string, nb, sig_mass ) ;
      TFile* tf_out = new TFile( out_filename, "RECREATE" ) ;
      h_sig_lf_sum -> Write() ;
      h_bg_sum -> Write() ;

      tf_out -> Close() ;




   }


