

#include "utils.c"


//=====================================================================

   void optimization_from_plotter_hists(
           const char* infile = "binning-optimization-input-WH-3b-ma60.root",
           double bg_syst_mult = 0.1, // fraction of background
           double bg_syst_add = 2.0 // bg events
           ) {

      gStyle -> SetPaintTextFormat( "5.2f" ) ;
      gStyle -> SetMarkerSize( 1.4 ) ;
      gStyle -> SetOptStat(0) ;

      TString tsif( infile ) ;

      bool is_wh(false) ;
      bool is_zh(false) ;
      if ( tsif.Contains("WH") ) is_wh = true ;
      if ( tsif.Contains("ZH") ) is_zh = true ;

      if ( !is_wh && !is_zh ) { printf("\n\n *** can't tell from input file name if this is WH or ZH.\n\n") ; return ; }

      if ( is_wh && is_zh ) { printf("\n\n *** can't tell from input file name if this is WH or ZH.\n\n") ; return ; }

      if ( is_wh ) printf("\n\n Input file is WH : %s\n\n", infile ) ;
      if ( is_zh ) printf("\n\n Input file is ZH : %s\n\n", infile ) ;

      int nb(0) ;
      if ( tsif.Contains("3b") ) nb = 3 ;
      if ( tsif.Contains("4b") ) nb = 4 ;
      if ( nb < 3 ) { printf("\n\n *** problem getting nb from input file name.\n\n") ; return ; }

      char answ ;
      char outfile[1000] ;
      char label[1000] ;

      TFile* tf_in = new TFile( infile, "READ" ) ;
      if ( tf_in == 0x0 ) { printf("\n\n *** bad input file : %s\n\n", infile ) ; return ; }
      if ( !(tf_in -> IsOpen()) ) { printf("\n\n *** problem opening input file : %s\n\n", infile ) ; return ; }



      gStyle -> SetOptTitle(0) ;
      gStyle -> SetPadBottomMargin( 0.13 ) ;

      TText* label_tt = new TText() ;
      TText* xtitle_tt = new TText() ;
      TText* ytitle_tt = new TText() ;
      ytitle_tt -> SetTextAngle( 90 ) ;
      ytitle_tt -> SetTextAlign( 31 ) ;
      xtitle_tt -> SetTextAlign( 31 ) ;
      label_tt -> SetTextSize( 0.04 ) ;
      xtitle_tt -> SetTextSize( 0.04 ) ;
      ytitle_tt -> SetTextSize( 0.04 ) ;

      gSystem -> Exec( "mkdir -p outputfiles" ) ;


      gDirectory -> Delete( "*" ) ;


      char hname[1000] ;

      TH1F* h_sig = get_hist( "h_sig" ) ;

      TH1F* h_bg = get_hist( "h_bg_sum" ) ;

      int nvals_highbinedge(0) ;
      int nvals_lowbinedge(0) ;
      double highbinedge_vals[100] ;
      double lowbinedge_vals[100] ;

         int i(0) ;
    //----------
         highbinedge_vals[i] = 0.07 ; i++ ;
         highbinedge_vals[i] = 0.09 ; i++ ;
         highbinedge_vals[i] = 0.11 ; i++ ;
         highbinedge_vals[i] = 0.13 ; i++ ;
         highbinedge_vals[i] = 0.15 ; i++ ;
         highbinedge_vals[i] = 0.17 ; i++ ;
         highbinedge_vals[i] = 0.19 ; i++ ;
         highbinedge_vals[i] = 0.21 ; i++ ;
         highbinedge_vals[i] = 0.23 ; i++ ;
         highbinedge_vals[i] = 0.25 ; i++ ;
         highbinedge_vals[i] = 0.27 ; i++ ;
         highbinedge_vals[i] = 0.29 ; i++ ;
         highbinedge_vals[i] = 0.31 ; i++ ;
         highbinedge_vals[i] = 0.33 ; i++ ;
         nvals_highbinedge = i ;
         i = 0 ;
    //----------
         lowbinedge_vals[i] = -0.09 ; i++ ;
         lowbinedge_vals[i] = -0.07 ; i++ ;
         lowbinedge_vals[i] = -0.05 ; i++ ;
         lowbinedge_vals[i] = -0.03 ; i++ ;
         lowbinedge_vals[i] = -0.01 ; i++ ;
         lowbinedge_vals[i] =  0.01 ; i++ ;
         lowbinedge_vals[i] =  0.03 ; i++ ;
         lowbinedge_vals[i] =  0.05 ; i++ ;
         lowbinedge_vals[i] =  0.07 ; i++ ;
         lowbinedge_vals[i] =  0.09 ; i++ ;
         lowbinedge_vals[i] =  0.11 ; i++ ;
         lowbinedge_vals[i] =  0.13 ; i++ ;
         lowbinedge_vals[i] =  0.15 ; i++ ;
         lowbinedge_vals[i] =  0.17 ; i++ ;
         lowbinedge_vals[i] =  0.19 ; i++ ;
         lowbinedge_vals[i] =  0.21 ; i++ ;
         lowbinedge_vals[i] =  0.23 ; i++ ;
         lowbinedge_vals[i] =  0.25 ; i++ ;
         lowbinedge_vals[i] =  0.27 ; i++ ;
         lowbinedge_vals[i] =  0.29 ; i++ ;
         nvals_lowbinedge = i ;

      TH2F* h_s_over_sigmab_quadsum = new TH2F( "h_s_over_sigmab_quadsum", "h_s_over_sigmab_quadsum",
                 nvals_highbinedge, 0.5, nvals_highbinedge+0.5,
                 nvals_lowbinedge, 0.5, nvals_lowbinedge+0.5 ) ;

      for ( int i=1; i<=nvals_highbinedge; i++ ) {
         char label[100] ;
         sprintf( label, "%.2f", highbinedge_vals[i-1] ) ;
         h_s_over_sigmab_quadsum -> GetXaxis() -> SetBinLabel( i, label ) ;
         h_s_over_sigmab_quadsum -> SetXTitle( "low edge of highest BDT bin" ) ;
      }
      for ( int i=1; i<=nvals_lowbinedge; i++ ) {
         char label[100] ;
         sprintf( label, "%.2f", lowbinedge_vals[i-1] ) ;
         h_s_over_sigmab_quadsum -> GetYaxis() -> SetBinLabel( i, label ) ;
         h_s_over_sigmab_quadsum -> SetYTitle( "low edge of 2nd highest BDT bin" ) ;
      }

      TH2F* h_s_over_sigmab_highbin = (TH2F*) h_s_over_sigmab_quadsum -> Clone( "h_s_over_sigmab_highbin" ) ;
      TH2F* h_s_over_sigmab_lowbin = (TH2F*) h_s_over_sigmab_quadsum -> Clone( "h_s_over_sigmab_lowbin" ) ;

      TH2F* h_s_yield_highbin = (TH2F*) h_s_over_sigmab_quadsum -> Clone( "h_s_yield_highbin" ) ;
      TH2F* h_bg_yield_highbin = (TH2F*) h_s_over_sigmab_quadsum -> Clone( "h_bg_yield_highbin" ) ;

      TH2F* h_s_yield_lowbin = (TH2F*) h_s_over_sigmab_quadsum -> Clone( "h_s_yield_lowbin" ) ;
      TH2F* h_bg_yield_lowbin = (TH2F*) h_s_over_sigmab_quadsum -> Clone( "h_bg_yield_lowbin" ) ;

      int nbins_1d_hist = h_sig -> GetNbinsX() ;

      int best_hbi(-1) ;
      int best_lbi(-1) ;

      float best_quad_sum_s_over_sigmab(0.) ;

      for ( int hbi=1; hbi<=nvals_highbinedge; hbi++ ) {

         int high_bdt_bin_low_edge_1d_hist_bin = h_sig -> FindBin( highbinedge_vals[hbi-1] + 0.001 ) ;
         printf("   hbi = %2d : high_bdt_bin_low_edge_1d_hist_bin = %2d, edge = %.2f\n", hbi, high_bdt_bin_low_edge_1d_hist_bin, highbinedge_vals[hbi-1] ) ;

         for ( int lbi=1; lbi<=nvals_lowbinedge; lbi++ ) {

            if ( lowbinedge_vals[lbi-1] >= highbinedge_vals[hbi-1] ) continue ;

            int low_bdt_bin_low_edge_1d_hist_bin = h_sig -> FindBin( lowbinedge_vals[lbi-1] + 0.001 ) ;
            printf("   lbi = %2d :  low_bdt_bin_low_edge_1d_hist_bin = %2d, edge = %.2f\n", lbi, low_bdt_bin_low_edge_1d_hist_bin, lowbinedge_vals[lbi-1] ) ;

            double high_bdt_bin_sig_err_mc_stat(0.) ;
            double  low_bdt_bin_sig_err_mc_stat(0.) ;

            double high_bdt_bin_bg_err_mc_stat(0.) ;
            double  low_bdt_bin_bg_err_mc_stat(0.) ;

            float high_bdt_bin_sig_val = h_sig -> IntegralAndError( high_bdt_bin_low_edge_1d_hist_bin, nbins_1d_hist, high_bdt_bin_sig_err_mc_stat ) ;
            float  low_bdt_bin_sig_val = h_sig -> IntegralAndError( low_bdt_bin_low_edge_1d_hist_bin, high_bdt_bin_low_edge_1d_hist_bin-1, low_bdt_bin_sig_err_mc_stat ) ;

            float high_bdt_bin_bg_val = h_bg -> IntegralAndError( high_bdt_bin_low_edge_1d_hist_bin, nbins_1d_hist, high_bdt_bin_bg_err_mc_stat ) ;
            float  low_bdt_bin_bg_val = h_bg -> IntegralAndError( low_bdt_bin_low_edge_1d_hist_bin, high_bdt_bin_low_edge_1d_hist_bin-1, low_bdt_bin_bg_err_mc_stat ) ;



            if ( high_bdt_bin_bg_val <= 0 ) high_bdt_bin_bg_val = 1. ; // always have at least one bg event.
            ///////float high_bdt_bin_bg_err = sqrt( high_bdt_bin_bg_val + pow( high_bdt_bin_bg_val * bg_syst_mult, 2 ) + pow( bg_syst_add, 2. ) ) ;
            float high_bdt_bin_bg_err = sqrt( pow( high_bdt_bin_bg_err_mc_stat, 2. ) + pow( high_bdt_bin_bg_val * bg_syst_mult, 2 ) + pow( bg_syst_add, 2. ) ) ;

            if ( low_bdt_bin_bg_val <= 0 ) low_bdt_bin_bg_val = 1. ; // always have at least one bg event.
            ///////float low_bdt_bin_bg_err = sqrt( low_bdt_bin_bg_val + pow( low_bdt_bin_bg_val * bg_syst_mult, 2 ) + pow( bg_syst_add, 2. ) ) ;
            float low_bdt_bin_bg_err = sqrt( pow( low_bdt_bin_bg_err_mc_stat, 2. ) + pow( low_bdt_bin_bg_val * bg_syst_mult, 2 ) + pow( bg_syst_add, 2. ) ) ;

            printf("     High bin, S = %8.1f, B = %8.1f +/- %7.1f\n", high_bdt_bin_sig_val, high_bdt_bin_bg_val, high_bdt_bin_bg_err ) ;
            printf("      Low bin, S = %8.1f, B = %8.1f +/- %7.1f\n", low_bdt_bin_sig_val, low_bdt_bin_bg_val, low_bdt_bin_bg_err ) ;

            float high_bdt_bin_s_over_sigmab = high_bdt_bin_sig_val / high_bdt_bin_bg_err ;
            float low_bdt_bin_s_over_sigmab = low_bdt_bin_sig_val / low_bdt_bin_bg_err ;
            float quad_sum_s_over_sigmab = sqrt( high_bdt_bin_s_over_sigmab * high_bdt_bin_s_over_sigmab + low_bdt_bin_s_over_sigmab * low_bdt_bin_s_over_sigmab ) ;

            if ( quad_sum_s_over_sigmab > best_quad_sum_s_over_sigmab ) {
               best_quad_sum_s_over_sigmab = quad_sum_s_over_sigmab ;
               best_hbi = hbi ;
               best_lbi = lbi ;
            }

            printf("     S / sigma B :   high bin = %7.2f, low bin = %7.2f , quad sum = %7.2f\n\n",
               high_bdt_bin_s_over_sigmab, low_bdt_bin_s_over_sigmab, quad_sum_s_over_sigmab ) ;

            h_s_over_sigmab_quadsum -> SetBinContent( hbi, lbi, quad_sum_s_over_sigmab ) ;
            h_s_over_sigmab_highbin -> SetBinContent( hbi, lbi, high_bdt_bin_s_over_sigmab ) ;
            h_s_over_sigmab_lowbin -> SetBinContent( hbi, lbi, low_bdt_bin_s_over_sigmab ) ;


            h_s_yield_highbin -> SetBinContent( hbi, lbi, high_bdt_bin_sig_val ) ;
            h_s_yield_lowbin  -> SetBinContent( hbi, lbi,  low_bdt_bin_sig_val ) ;

            h_s_yield_highbin -> SetBinError( hbi, lbi, high_bdt_bin_sig_err_mc_stat ) ;
            h_s_yield_lowbin  -> SetBinError( hbi, lbi,  low_bdt_bin_sig_err_mc_stat ) ;


            h_bg_yield_highbin -> SetBinContent( hbi, lbi, high_bdt_bin_bg_val ) ;
            h_bg_yield_lowbin  -> SetBinContent( hbi, lbi,  low_bdt_bin_bg_val ) ;

            h_bg_yield_highbin -> SetBinError( hbi, lbi, high_bdt_bin_bg_err ) ;
            h_bg_yield_lowbin  -> SetBinError( hbi, lbi,  low_bdt_bin_bg_err ) ;


         } // lbi

      } // hbi

    //--------

      double best_high_bin_low_edge = highbinedge_vals[best_hbi-1] ;
      double best_low_bin_low_edge = lowbinedge_vals[best_lbi-1] ;

      float best_high_bin_sig_val = h_s_yield_highbin->GetBinContent( best_hbi, best_lbi ) ;
      float best_high_bin_sig_err = h_s_yield_highbin->GetBinError( best_hbi, best_lbi ) ;

      float best_low_bin_sig_val = h_s_yield_lowbin->GetBinContent( best_hbi, best_lbi ) ;
      float best_low_bin_sig_err = h_s_yield_lowbin->GetBinError( best_hbi, best_lbi ) ;

      float best_high_bin_bg_val = h_bg_yield_highbin->GetBinContent( best_hbi, best_lbi ) ;
      float best_high_bin_bg_err = h_bg_yield_highbin->GetBinError( best_hbi, best_lbi ) ;

      float best_low_bin_bg_val = h_bg_yield_lowbin->GetBinContent( best_hbi, best_lbi ) ;
      float best_low_bin_bg_err = h_bg_yield_lowbin->GetBinError( best_hbi, best_lbi ) ;



      double var_bins[6] ;

      if ( is_wh && nb == 3 ) {
         var_bins[0] = -0.31 ;
         var_bins[1] = -0.13 ;
         var_bins[2] = -0.01 ;
         var_bins[3] = best_low_bin_low_edge ;
         var_bins[4] = best_high_bin_low_edge ;
         var_bins[5] = 0.35 ;
      }
      if ( is_wh && nb == 4 ) {
         var_bins[0] = -0.31 ;
         var_bins[1] = -0.13 ;
         var_bins[2] = -0.07 ;
         var_bins[3] = best_low_bin_low_edge ;
         var_bins[4] = best_high_bin_low_edge ;
         var_bins[5] = 0.35 ;
      }

      if ( is_zh && nb == 3 ) {
         var_bins[0] = -0.31 ;
         var_bins[1] = -0.09 ;
         var_bins[2] =  0.01 ;
         var_bins[3] = best_low_bin_low_edge ;
         var_bins[4] = best_high_bin_low_edge ;
         var_bins[5] = 0.35 ;
      }
      if ( is_zh && nb == 4 ) {
         var_bins[0] = -0.31 ;
         var_bins[1] = -0.15 ;
         var_bins[2] = -0.07 ;
         var_bins[3] = best_low_bin_low_edge ;
         var_bins[4] = best_high_bin_low_edge ;
         var_bins[5] = 0.35 ;
      }

      TH1* h_sig_5bins = h_sig -> Rebin( 5, "h_sig_5bins", var_bins ) ;
      TH1* h_bg_5bins = h_bg -> Rebin( 5, "h_bg_5bins", var_bins ) ;

      printf("\n\n ============================================================================== \n\n" ) ;
      printf("  Results for optimal S / sigmaB\n\n") ;
      printf("    Lower edge of high bin:  %.2f\n", best_high_bin_low_edge ) ;
      printf("    Lower edge of low  bin:  %.2f\n", best_low_bin_low_edge ) ;
      printf("\n") ;
      printf("    High bin content:   Signal = %6.1f +/- %4.1f        BG = %6.1f +/- %4.1f\n",
          best_high_bin_sig_val, best_high_bin_sig_err,
          best_high_bin_bg_val, best_high_bin_bg_err ) ;
      printf("    Low bin content:    Signal = %6.1f +/- %4.1f        BG = %6.1f +/- %4.1f\n",
          best_low_bin_sig_val, best_low_bin_sig_err,
          best_low_bin_bg_val, best_low_bin_bg_err ) ;
      printf("\n") ;
      printf("       Signal error is MC stat only.\n") ;
      printf("       BG error is quadrature sum of MC stats, %.0f %% multiplicative syst, %.1f events additive syst.\n", 100*bg_syst_mult, bg_syst_add ) ;


      printf("\n  Bin edges: ") ;
      for ( int i=0; i<6; i++ ) printf(" %5.2f ", var_bins[i] ) ;

      printf("\n\n ============================================================================== \n\n" ) ;




    //--------

      gDirectory -> cd( "Rint:/") ;

      TEllipse* te = new TEllipse() ;
      te -> SetLineColor(2) ;
      te -> SetLineWidth(3) ;

      int ci ;

      TCanvas* can2 = get_canvas( "can2", "", 150, 50, 2200, 1000 ) ;
      can2 -> Clear() ;
      can2 -> Divide(2,2) ;
      can2 -> cd(1) ;

      ci = 1 ;

      can2 -> cd( ci ) ; ci ++ ;
      h_s_yield_highbin -> Draw("colz") ;
      h_s_yield_highbin -> Draw("text same") ;
      label_tt -> DrawTextNDC( 0.1, 0.93, "High bin, signal yield" ) ;

      can2 -> cd( ci ) ; ci ++ ;
      h_bg_yield_highbin -> Draw("colz") ;
      h_bg_yield_highbin -> Draw("text same") ;
      label_tt -> DrawTextNDC( 0.1, 0.93, "High bin, BG yield" ) ;

      can2 -> cd( ci ) ; ci ++ ;
      h_s_yield_lowbin -> Draw("colz") ;
      h_s_yield_lowbin -> Draw("text same") ;
      label_tt -> DrawTextNDC( 0.1, 0.93, "Low bin, signal yield" ) ;

      can2 -> cd( ci ) ; ci ++ ;
      h_bg_yield_lowbin -> Draw("colz") ;
      h_bg_yield_lowbin -> Draw("text same") ;
      label_tt -> DrawTextNDC( 0.1, 0.93, "Low bin, BG yield" ) ;

    //--------

      TCanvas* can1 = get_canvas( "can1", "", 50, 50, 1500, 1200 ) ;
      can1 -> Clear() ;
      can1 -> Divide(2,2) ;
      can1 -> cd(1) ;

      can1 -> Clear() ;
      can1 -> Divide(2,2) ;

      ci = 1 ;

      can1 -> cd( ci ) ; ci ++ ;
      h_s_over_sigmab_quadsum -> Draw("colz2") ;
      h_s_over_sigmab_quadsum -> Draw("text same") ;
      //sprintf( label, "%db : mh = %2d : Bsyst = %4.2f : S / sigma B, sum in quad.", nb, sigmass, bg_syst_mult ) ;
      sprintf( label, "Bsyst = %4.2f : S / sigma B, sum in quad.", bg_syst_mult ) ;
      label_tt -> DrawTextNDC( 0.1, 0.93, label ) ;

      can1 -> cd( ci ) ; ci ++ ;
      h_s_over_sigmab_highbin -> Draw("colz2") ;
      h_s_over_sigmab_highbin -> Draw("text same") ;
      //sprintf( label, "%db : mh = %2d : Bsyst = %4.2f : S / sigma B, high bin", nb, sigmass, bg_syst_mult ) ;
      sprintf( label, "Bsyst = %4.2f : S / sigma B, high bin", bg_syst_mult ) ;
      label_tt -> DrawTextNDC( 0.1, 0.93, label ) ;

      can1 -> cd( ci ) ; ci ++ ;
      h_s_over_sigmab_lowbin -> Draw("colz2") ;
      h_s_over_sigmab_lowbin -> Draw("text same") ;
      //sprintf( label, "%db : mh = %2d : Bsyst = %4.2f : S / sigma B, low bin", nb, sigmass, bg_syst_mult ) ;
      sprintf( label, "Bsyst = %4.2f : S / sigma B, low bin", bg_syst_mult ) ;
      label_tt -> DrawTextNDC( 0.1, 0.93, label ) ;

    //--------


      TLine* tl = new TLine() ;
      tl -> SetLineColor(2) ;
      tl -> SetLineWidth(2) ;

      TCanvas* can3 = get_canvas( "can3", "", 150, 150, 1500, 1200 ) ;
      can3 -> Clear() ;
      can3 -> Divide(2,3) ;



      TH1F* h_sig_zoom = (TH1F*) h_sig -> Clone( "h_sig_zoom" ) ;
      h_sig_zoom -> SetMaximum( 3.0 * h_sig->GetMaximum() ) ;

      TH1* h_sig_5bins_zoom = (TH1*) h_sig_5bins -> Clone( "h_sig_5bins_zoom" ) ;
      h_sig_5bins_zoom -> SetMaximum( 3.0 * h_sig_5bins ->GetMaximum() ) ;

      h_bg -> SetMinimum(0.5) ;
      h_bg_5bins -> SetMinimum(0.5) ;

      ci = 1 ;

      can3 -> cd( ci ) ; ci ++ ;
      h_bg -> Draw() ;
      h_sig -> Draw( "same" ) ;
      tl -> DrawLine( best_low_bin_low_edge, 0., best_low_bin_low_edge, h_bg -> GetMaximum() ) ;
      tl -> DrawLine( best_high_bin_low_edge, 0., best_high_bin_low_edge, h_bg -> GetMaximum() ) ;

      can3 -> cd( ci ) ; ci ++ ;
      h_bg -> Draw() ;
      h_sig -> Draw( "same" ) ;
      gPad -> SetLogy(1) ;
      tl -> DrawLine( best_low_bin_low_edge, 0., best_low_bin_low_edge, h_bg -> GetMaximum() ) ;
      tl -> DrawLine( best_high_bin_low_edge, 0., best_high_bin_low_edge, h_bg -> GetMaximum() ) ;

      can3 -> cd( ci ) ; ci ++ ;
      h_sig_zoom -> Draw() ;
      h_bg -> Draw("same") ;
      tl -> DrawLine( best_low_bin_low_edge, 0., best_low_bin_low_edge, h_sig_zoom -> GetMaximum() ) ;
      tl -> DrawLine( best_high_bin_low_edge, 0., best_high_bin_low_edge, h_sig_zoom -> GetMaximum() ) ;

      can3 -> cd( ci ) ; ci ++ ;
      h_bg_5bins -> Draw() ;
      h_sig_5bins -> Draw( "same" ) ;
      gPad -> SetLogy(1) ;

      can3 -> cd( ci ) ; ci ++ ;
      h_bg_5bins -> Draw() ;
      h_sig_5bins -> Draw( "same" ) ;

      can3 -> cd( ci ) ; ci ++ ;
      h_sig_5bins_zoom -> Draw() ;
      h_bg_5bins -> Draw( "same" ) ;


   } // draw_bdt

//=====================================================================

