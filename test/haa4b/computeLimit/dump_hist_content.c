

#include "utils.c"


   void dump_hist_content( const char* infile = "binning-optimization-input-WH-3b-ma60.root" ) {

      gDirectory -> Delete( "h*") ;

      TFile* tf_in = new TFile( infile, "READ" ) ;
      if ( tf_in == 0x0 ) { printf("\n\n *** bad input file : %s\n\n", infile ) ; return ; }
      if ( !(tf_in -> IsOpen()) ) { printf("\n\n *** problem opening input file : %s\n\n", infile ) ; return ; }

      tf_in -> ls() ;

      TH1F* h_sig = get_hist( "h_sig" ) ;
      TH1F* h_bg_sum = get_hist( "h_bg_sum" ) ;

      int nbins = h_sig -> GetNbinsX() ;
      printf("\n\n Histogram has %d bins.\n\n", nbins ) ;

      for ( int bi=1; bi<=nbins; bi++ ) {

         float bin_low_edge = h_sig -> GetXaxis() -> GetBinLowEdge( bi ) ;
         float bin_up_edge = h_sig -> GetXaxis() -> GetBinUpEdge( bi ) ;

         float sig_val = h_sig -> GetBinContent( bi ) ;
         float sig_err = h_sig -> GetBinError( bi ) ;

         float bg_val = h_bg_sum -> GetBinContent( bi ) ;
         float bg_err = h_bg_sum -> GetBinError( bi ) ;

         printf("  %3d : [%5.2f,%5.2f]   BG = %9.2f +/- %6.2f    Sig = %9.2f +/- %6.2f\n",
            bi, bin_low_edge, bin_up_edge, bg_val, bg_err,  sig_val, sig_err ) ;

      } // bi

   }





