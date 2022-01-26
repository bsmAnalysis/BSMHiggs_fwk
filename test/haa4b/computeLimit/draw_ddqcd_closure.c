
   TH1* get_hist( const char* hname, TFile& tf, bool require_found = false ) ;

   //
   //  The jobs_dir argument is a directory that holds all of the job output directories.
   //  That is, the jobs_dir directory should have these directories in it
   //
   //     cards_SB13TeV_SM_Wh_2016_noSoftb
   //     cards_SB13TeV_SM_Wh_2017_noSoftb
   //     cards_SB13TeV_SM_Wh_2018_noSoftb
   //

   void draw_ddqcd_closure( bool pause = false, const char* jobs_dir = "" ) {

      gStyle -> SetPadBottomMargin(0.25) ;
      gStyle -> SetPadTopMargin(0.05) ;
      gStyle -> SetOptTitle(0) ;
      gStyle -> SetPadLeftMargin(0.20) ;
      //      gStyle -> SetStats(0) ;

      vector<string> sel ;
      vector<string> lep ;
      vector<string> nb ;
      vector<string> year ;

      sel.emplace_back( "CR" ) ;
      //sel.emplace_back( "SR" ) ;

      lep.emplace_back( "e" ) ;
      lep.emplace_back( "mu" ) ;

      nb.emplace_back( "3b" ) ;
      nb.emplace_back( "4b" ) ;

      year.emplace_back( "2016" ) ;
      year.emplace_back( "2017" ) ;
      year.emplace_back( "2018" ) ;


      float hist_max =  50000 ;
      float hist_min = -10000 ;



      //int nbins = ( sel.size() * lep.size() + 1 ) * nb.size() * year.size() ;
      int nbins = ( year.size() * sel.size() * lep.size() + 1 ) * nb.size() ;
      printf("\n nbins = %d\n", nbins ) ;

      TH1F* h_coverd_vs_channel = new TH1F( "h_coverd_vs_channel", "C/D ratio vs channel", nbins+2, -0.5, nbins+1.5 ) ;

      TH1F* h_pred_over_obs_vs_channel =  new TH1F( "h_pred_over_obs_vs_channel", "DDQCD over (data - nonQCD) in A vs channel", nbins+2, -0.5, nbins+1.5 ) ;


      TH1F* h_pred_vs_channel =  new TH1F( "h_pred_vs_channel", "DDQCD in A vs channel", nbins+1, -0.5, nbins+0.5 ) ;
      TH1F* h_pred_vs_channel_total_error =  new TH1F( "h_pred_vs_channel_total_error", "DDQCD in A vs channel", nbins+1, -0.5, nbins+0.5 ) ;

      TH1F* h_data_minus_nonqcd_vs_channel =  new TH1F( "h_data_minus_nonqcd_vs_channel", "(data - nonQCD) in A vs channel", nbins+1, -0.5, nbins+0.5 ) ;
      TH1F* h_data_minus_nonqcd_vs_channel_total_error =  new TH1F( "h_data_minus_nonqcd_vs_channel_total_error", "(data - nonQCD) in A vs channel", nbins+1, -0.5, nbins+0.5 ) ;


      vector<string> non_qcd_procs ;
      vector<string> scale_factor_name ;

      non_qcd_procs.emplace_back("wlnu")     ; scale_factor_name.emplace_back("w_norm") ;
      non_qcd_procs.emplace_back("otherbkg") ; scale_factor_name.emplace_back("") ;
      non_qcd_procs.emplace_back("zll")      ; scale_factor_name.emplace_back("") ;
      non_qcd_procs.emplace_back("ttbarbba") ; scale_factor_name.emplace_back("tt_norm") ;
      non_qcd_procs.emplace_back("ttbarcba") ; scale_factor_name.emplace_back("tt_norm") ;
      non_qcd_procs.emplace_back("ttbarlig") ; scale_factor_name.emplace_back("") ;

      {
         int bini = 2 ;
         for ( int bi=0; bi<nb.size(); bi++ ) {

            for ( int si=0; si<sel.size(); si++ ) {

               for ( int li=0; li<lep.size(); li++ ) {

                  for ( int yi=0; yi<year.size(); yi++ ) {

                      char infile[10000] ;

                      sprintf( infile, "%scards_SB13TeV_SM_Wh_%s_noSoftb/0040/haa4b_40_13TeV_wh.root", jobs_dir, year[yi].c_str() ) ;
                      TFile tf( infile, "READ" ) ;
                      if ( !tf.IsOpen() ) { printf("\n\n *** problem opening %s\n\n", infile ) ; gSystem -> Exit(-1) ; }

                      sprintf( infile, "%scards_SB13TeV_SM_Wh_%s_noSoftb/0040/fitDiagnostics-first.root", jobs_dir, year[yi].c_str() ) ;
                      TFile tf_fd( infile, "READ" ) ;
                      if ( !tf_fd.IsOpen() ) { printf("\n\n *** problem opening %s\n\n", infile ) ; gSystem -> Exit(-1) ; }

                      RooFitResult* rfr = (RooFitResult*) tf_fd.Get( "fit_b" ) ;
                      if ( rfr == 0x0 ) { printf("\n\n *** can't find RooFitResult for fit_b in %s\n\n", infile ) ; gSystem -> Exit(-1); }

                      char hname[1000] ;
                      char htitle[1000] ;

                      char label[100] ;
                      sprintf( label, "%s, %s, %s, %s", sel[si].c_str(), nb[bi].c_str(), lep[li].c_str(), year[yi].c_str() ) ;

                     //-- in boxes B, C, and C, data_obs already has the non-QCD stuff subtracted.

                      sprintf( hname, "%s_C_%s_%s/data_obs", lep[li].c_str(), sel[si].c_str(), nb[bi].c_str() ) ;
                      TH1* h_data_minus_nonqcd_c = get_hist( hname, tf, 1 ) ;

                      sprintf( hname, "%s_D_%s_%s/data_obs", lep[li].c_str(), sel[si].c_str(), nb[bi].c_str() ) ;
                      TH1* h_data_minus_nonqcd_d = get_hist( hname, tf, 1 ) ;

                      int integral_last_bin = h_data_minus_nonqcd_c -> GetNbinsX() ;
                      if ( strcmp( sel[si].c_str(), "SR" ) == 0 ) { integral_last_bin = 3 ; }
                      double ierr ;


                      sprintf( hname, "h_all_non_qcd_%s_A_%s_%s", lep[li].c_str(), sel[si].c_str(), nb[bi].c_str() ) ;
                      TH1* h_all_non_qcd = (TH1*) h_data_minus_nonqcd_d -> Clone( hname ) ;
                      h_all_non_qcd -> Reset() ;
                      sprintf( htitle, "All non-QCD BG, %s A %s %s", lep[li].c_str(), sel[si].c_str(), nb[bi].c_str() ) ;
                      h_all_non_qcd -> SetTitle( htitle ) ;

                      double sf_tt_val(0.) ;
                      double sf_tt_err(0.) ;
                      double sf_w_val(0.) ;
                      double sf_w_err(0.) ;
                      double tt_sum_val(0.) ;
                      double w_sum_val(0.) ;
                      for ( int pi=0; pi<non_qcd_procs.size(); pi++ ) {
                         sprintf( hname, "%s_A_%s_%s/%s", lep[li].c_str(), sel[si].c_str(), nb[bi].c_str(), non_qcd_procs[pi].c_str() ) ;
                         TH1* h_nqcdbg = get_hist( hname, tf, 0 ) ;
                         if ( h_nqcdbg != 0x0 ) {
                            float sf_val = 1. ;
                            float sf_err = 0. ;
                            if ( scale_factor_name[pi].length() > 0 ) {
                               char pname[100] ;
                               sprintf( pname, "%s_%s", scale_factor_name[pi].c_str(), lep[li].c_str() ) ;
                               RooRealVar* rrv = (RooRealVar*)(rfr -> floatParsFinal()).find( pname ) ;
                               sf_val = rrv -> getVal() ;
                               sf_err = rrv -> getError() ;
                               ///// printf("   hist  %s ,  sf  %s  =  %6.3f +/- %6.3f\n", hname, pname, sf_val, sf_err ) ;
                               if ( strcmp( scale_factor_name[pi].c_str(), "tt_norm" ) == 0 ) {
                                  tt_sum_val += h_nqcdbg -> IntegralAndError( 1, integral_last_bin, ierr ) ;
                                  sf_tt_val = sf_val ;
                                  sf_tt_err = sf_err ;
                               }
                               if ( strcmp( scale_factor_name[pi].c_str(), "w_norm" ) == 0 ) {
                                  w_sum_val += h_nqcdbg -> IntegralAndError( 1, integral_last_bin, ierr ) ;
                                  sf_w_val = sf_val ;
                                  sf_w_err = sf_err ;
                               }
                            }
                            h_all_non_qcd -> Add( h_nqcdbg, sf_val ) ;
                         }
                      } // pi

                      /////double all_non_qcd_sf_err = sqrt( pow( tt_sum_val * sf_tt_val * sf_tt_err, 2. ) + pow( w_sum_val * sf_w_val * sf_w_err, 2. ) ) ;
                      char pname_tt[100] ;
                      char pname_w[100] ;
                      sprintf( pname_tt, "tt_norm_%s", lep[li].c_str() ) ;
                      sprintf( pname_w, "w_norm_%s", lep[li].c_str() ) ;
                      double sf_rho = rfr -> correlation( pname_tt, pname_w ) ;
                      printf( "     rho for %s and %s = %7.3f\n", pname_tt, pname_w, sf_rho ) ;
                      double all_non_qcd_sf_err = sqrt( pow( tt_sum_val * sf_tt_val * sf_tt_err, 2. ) + pow( w_sum_val * sf_w_val * sf_w_err, 2. )  +  2. * sf_tt_val * tt_sum_val * sf_w_val * w_sum_val * sf_tt_err * sf_w_err * sf_rho ) ;


                      sprintf( hname, "%s_A_%s_%s/ddqcd", lep[li].c_str(), sel[si].c_str(), nb[bi].c_str() ) ;
                      TH1* h_ddqcd_pred = get_hist( hname, tf, 0 ) ;



                      double c_integral_val, c_integral_err ;
                      double d_integral_val, d_integral_err ;

                      c_integral_val = h_data_minus_nonqcd_c -> IntegralAndError( 1, h_data_minus_nonqcd_c ->GetNbinsX(), c_integral_err ) ;
                      d_integral_val = h_data_minus_nonqcd_d -> IntegralAndError( 1, h_data_minus_nonqcd_d ->GetNbinsX(), d_integral_err ) ;

                      double coverd_val = c_integral_val / d_integral_val ;
                      double coverd_err = coverd_val * sqrt( pow( c_integral_err/c_integral_val, 2. ) + pow( d_integral_err/d_integral_val, 2. ) ) ;






                      sprintf( hname, "%s_A_%s_%s/data_obs", lep[li].c_str(), sel[si].c_str(), nb[bi].c_str() ) ;
                      TH1* h_data_minus_nonqcd_a = get_hist( hname, tf, 1 ) ;
                      ///////printf("  Integral of %s : %9.1f\n", hname, h_data_minus_nonqcd_a -> Integral() ) ;
                      ///////printf("  Integral of %s : %9.1f\n", h_all_non_qcd->GetName(), h_all_non_qcd -> Integral() ) ;

                      double all_non_qcd_integral_val(1.), all_non_qcd_integral_stat_err(0.) ;
                      all_non_qcd_integral_val = h_all_non_qcd -> IntegralAndError( 1, integral_last_bin, all_non_qcd_integral_stat_err ) ;
                      double data_val = h_data_minus_nonqcd_a -> Integral() ;
                      double data_err = sqrt( data_val ) ;
                      double data_minus_nonqcd_val =  data_val - all_non_qcd_integral_val ;
                      double data_minus_nonqcd_stat_err = sqrt( pow( data_err, 2. ) + pow( all_non_qcd_integral_stat_err, 2. ) ) ;
                      double data_minus_nonqcd_sf_err = all_non_qcd_sf_err ;

                      h_data_minus_nonqcd_vs_channel -> SetBinContent( bini, data_minus_nonqcd_val ) ;
                      h_data_minus_nonqcd_vs_channel -> SetBinError( bini, data_minus_nonqcd_stat_err ) ;
                      h_data_minus_nonqcd_vs_channel_total_error -> SetBinContent( bini, data_minus_nonqcd_val ) ;
                      h_data_minus_nonqcd_vs_channel_total_error -> SetBinError( bini, sqrt( pow( data_minus_nonqcd_stat_err, 2. ) + pow( data_minus_nonqcd_sf_err, 2. ) ) ) ;

                      h_data_minus_nonqcd_vs_channel -> GetXaxis() -> SetBinLabel( bini, label ) ;
                      h_data_minus_nonqcd_vs_channel_total_error -> GetXaxis() -> SetBinLabel( bini, label ) ;

                      printf("\n  %20s :  data - nonQCD = (%9.1f +/- %6.1f) - (%9.1f +/- %6.1f +/- %6.1f sf )   =   %9.1f +/- %6.1f +/- %6.1f\n",
                         label,
                         data_val, data_err,
                         all_non_qcd_integral_val, all_non_qcd_integral_stat_err, all_non_qcd_sf_err,
                         data_minus_nonqcd_val, data_minus_nonqcd_stat_err, data_minus_nonqcd_sf_err
                         ) ;


                      h_data_minus_nonqcd_a -> Add( h_all_non_qcd, -1. ) ;
                      ///////printf("  Integral after subtraction : %9.1f\n", h_data_minus_nonqcd_a -> Integral() ) ;

                      double ddqcd_pred_val(0.), ddqcd_pred_err(1.), ddqcd_pred_syst(1.) ;
                      /////////double data_minus_nonqcd_val(1.), data_minus_nonqcd_err(1.) ;
                      double pred_over_obs_val(0.), pred_over_obs_err(0.), pred_over_obs_systsf(0.) ;


                      if ( h_ddqcd_pred != 0x0 ) {
                         ddqcd_pred_val = h_ddqcd_pred -> IntegralAndError( 1, integral_last_bin, ddqcd_pred_err ) ;
                         ddqcd_pred_syst = 0.5 * ddqcd_pred_val ;
                         ////////data_minus_nonqcd_val = h_data_minus_nonqcd_a -> IntegralAndError( 1, integral_last_bin, data_minus_nonqcd_stat_err ) ;
                         pred_over_obs_val = ddqcd_pred_val / data_minus_nonqcd_val ;
                         pred_over_obs_err = pred_over_obs_val * sqrt( pow( ddqcd_pred_err/ddqcd_pred_val, 2. ) + pow( data_minus_nonqcd_stat_err/data_minus_nonqcd_val, 2. ) ) ;
                         pred_over_obs_systsf = pred_over_obs_val * sqrt( pow( ddqcd_pred_syst/ddqcd_pred_val, 2. ) + pow( data_minus_nonqcd_sf_err/data_minus_nonqcd_val, 2. ) ) ;
                         h_pred_over_obs_vs_channel -> SetBinContent( bini, pred_over_obs_val ) ;
                         h_pred_over_obs_vs_channel -> SetBinError( bini, pred_over_obs_err ) ;
                      }
                      h_pred_over_obs_vs_channel -> GetXaxis() -> SetBinLabel( bini, label ) ;

                      h_pred_vs_channel -> SetBinContent( bini, ddqcd_pred_val ) ;
                      h_pred_vs_channel -> SetBinError( bini, ddqcd_pred_err ) ;
                      h_pred_vs_channel_total_error -> SetBinContent( bini, ddqcd_pred_val ) ;
                      h_pred_vs_channel_total_error -> SetBinError( bini, sqrt( pow( ddqcd_pred_err, 2. ) + pow( ddqcd_pred_syst, 2. ) ) ) ;

                      h_pred_vs_channel -> GetXaxis() -> SetBinLabel( bini, label ) ;
                      h_pred_vs_channel_total_error -> GetXaxis() -> SetBinLabel( bini, label ) ;



                      h_coverd_vs_channel -> GetXaxis() -> SetBinLabel( bini, label ) ;


                      printf("  %20s :  C/D = (%9.1f +/- %6.1f) / (%9.1f +.- %6.1f) = %6.3f +/- %6.3f   |",
                          label,
                          c_integral_val, c_integral_err,
                          d_integral_val, d_integral_err,
                          coverd_val, coverd_err ) ;
                      if ( h_ddqcd_pred != 0x0 ) {
                         printf("   pred / obs = (%9.1f +/- %6.1f +/- %6.1f) / (%9.1f +/- %6.1f +/- %6.1f) = %6.3f +/- %6.3f +/- %6.3f |",
                            ddqcd_pred_val, ddqcd_pred_err, ddqcd_pred_syst,
                            data_minus_nonqcd_val, data_minus_nonqcd_stat_err, data_minus_nonqcd_sf_err,
                            pred_over_obs_val, pred_over_obs_err, pred_over_obs_systsf ) ;
                      }
                      printf("\n") ;

                      h_coverd_vs_channel -> SetBinContent( bini, coverd_val ) ;
                      h_coverd_vs_channel -> SetBinError( bini, coverd_err ) ;


                      bini ++ ;

                  } // yi

               } // li

               bini ++ ;

            } // si

         } // bi

      } // scoping for bini


      h_pred_vs_channel_total_error->SetFillColor( kRed-9 ) ;
      h_pred_vs_channel_total_error->SetFillStyle( 3144 ) ;
      h_data_minus_nonqcd_vs_channel_total_error -> SetFillColor( kBlue-10 ) ;

      h_pred_vs_channel -> SetLineColor(2) ;

      h_data_minus_nonqcd_vs_channel -> SetLineWidth(2) ;

      h_data_minus_nonqcd_vs_channel -> SetMarkerStyle(20) ;
      h_pred_vs_channel -> SetMarkerStyle(22) ;
      h_pred_vs_channel -> SetMarkerColor(2) ;

      //h_pred_vs_channel -> LabelsDeflate( "X" ) ;
      h_pred_vs_channel -> LabelsOption( "v" ) ;
      //h_pred_vs_channel_total_error -> LabelsDeflate( "X" ) ;
      h_pred_vs_channel_total_error -> LabelsOption( "v" ) ;


      //h_data_minus_nonqcd_vs_channel -> LabelsDeflate( "X" ) ;
      h_data_minus_nonqcd_vs_channel -> LabelsOption( "v" ) ;
      //h_data_minus_nonqcd_vs_channel_total_error -> LabelsDeflate( "X" ) ;
      h_data_minus_nonqcd_vs_channel_total_error -> LabelsOption( "v" ) ;

      TLine blue_line ;
      blue_line.SetLineColor(4) ;

      TLine red_line ;
      red_line.SetLineColor(2) ;

      TLine green_line ;
      green_line.SetLineColor(8) ;



      TCanvas* can = new TCanvas("can", "C/D vs channel", 50, 50, 1000, 1200 ) ;

      h_data_minus_nonqcd_vs_channel_total_error -> SetMaximum( hist_max ) ;
      h_data_minus_nonqcd_vs_channel_total_error -> SetMinimum( hist_min ) ;
      h_data_minus_nonqcd_vs_channel_total_error -> SetYTitle( "Number of QCD background events" ) ;
      h_data_minus_nonqcd_vs_channel_total_error -> SetTitleOffset( 2.2, "y" ) ;

      TLegend* legend = new TLegend( 0.45, 0.80, 0.88, 0.93 ) ;
      legend -> AddEntry( h_data_minus_nonqcd_vs_channel_total_error, "Data - non QCD" ) ;
      legend -> AddEntry( h_pred_vs_channel_total_error, "DD QCD prediction" ) ;

      h_data_minus_nonqcd_vs_channel_total_error -> SetStats(0) ;

      h_data_minus_nonqcd_vs_channel_total_error -> Draw("e2" ) ;
      h_data_minus_nonqcd_vs_channel -> Draw("same" ) ;
      h_pred_vs_channel_total_error -> Draw("e2 same" ) ;
      h_pred_vs_channel -> Draw("same") ;
      h_data_minus_nonqcd_vs_channel -> Draw("same" ) ;
      legend -> Draw() ;

      //float ymax = 1.9 * h_data_minus_nonqcd_vs_channel_total_error -> GetMaximum() ;
      //float ymin = -0.2 * ymax ;
      float ymax = hist_max ;
      float ymin = hist_min ;

      //for ( int i=1; i<6; i++ ) { blue_line.DrawLine( 6*i+0.5, ymin, 6*i+0.5, ymax ) ; }
      //for ( int i=1; i<6; i++ ) { blue_line.DrawLine( 7*i+0.5, ymin, 7*i+0.5, ymax ) ; }

      ///////for ( int i=0; i<7; i++ ) { blue_line.DrawLine( 7*i, ymin, 7*i, ymax ) ; }



      can -> SaveAs( "ddqcd-closure.pdf" ) ;


 //   float ymin = -10. ;
 //   float ymax = 10. ;

 //   h_coverd_vs_channel -> LabelsDeflate( "X" ) ;
 //   h_coverd_vs_channel -> LabelsOption( "v" ) ;

 //   h_coverd_vs_channel -> SetMarkerStyle(20) ;

 //   h_coverd_vs_channel -> Draw() ;

 //   red_line.DrawLine( 0.5, 0., nbins+0.5, 0. ) ;


 //   gPad -> SetGridy(1) ;

 //  //--------

 //   h_pred_over_obs_vs_channel -> SetMarkerStyle(20) ;

 //   h_pred_over_obs_vs_channel -> LabelsDeflate( "X" ) ;
 //   h_pred_over_obs_vs_channel -> LabelsOption( "v" ) ;



///   TCanvas* can2 = new TCanvas("can2", "pred/obs in A vs channel", 150, 350, 1000, 600 ) ;

///   h_pred_over_obs_vs_channel -> Draw("e0") ;

///   green_line.DrawLine( 0.5, 1., nbins+0.5, 1. ) ;

///   h_pred_over_obs_vs_channel -> Draw("e0 same") ;



///   for ( int i=1; i<6; i++ ) { blue_line.DrawLine( 6*i+0.5, ymin, 6*i+0.5, ymax ) ; }

///   gPad -> SetGridy(1) ;

   }

  //--------------

   TH1* get_hist( const char* hname, TFile& tf, bool require_found = false ) {

      TH1* rp = (TH1*) tf.Get( hname ) ;
      if ( require_found && rp == 0x0 ) { printf("\n\n *** get_hist: hist %s required but not found.\n\n", hname ) ; gSystem -> Exit(-1) ; }
      return rp ;

   } // get_hist

