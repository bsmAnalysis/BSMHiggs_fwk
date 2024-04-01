#ifndef utils_c
#define utils_c


   TCanvas* get_canvas( const char* cname, const char* ctitle = "",
                        int px = 50, int py = 50,
                        int w = 900, int h = 900 ) {
       TCanvas* rp = (TCanvas*) gDirectory -> FindObject( cname ) ;
       if ( rp != 0 ) return rp ;
       return new TCanvas( cname, ctitle, px, py, w, h ) ;
   }

   TH1F* get_hist( const char* hname ) {
      TH1F* rp = (TH1F*) gDirectory -> Get( hname ) ;
      if ( rp == 0x0 ) { printf("\n\n *** Missing hist %s\n\n", hname ) ; gSystem -> Exit(-1) ; }
      return rp ;
   }

   TH2F* get_hist2d( const char* hname ) {
      TH2F* rp = (TH2F*) gDirectory -> Get( hname ) ;
      if ( rp == 0x0 ) { printf("\n\n *** Missing hist %s\n\n", hname ) ; gSystem -> Exit(-1) ; }
      return rp ;
   }

   TH1F* fget_hist( const char* hname, TFile* tfp ) {
      if ( tfp == 0x0 ) { printf("\n\n *** fget_hist : bad file pointer.\n\n") ; gSystem -> Exit(-1) ; }
      TH1F* rp = (TH1F*) tfp -> Get( hname ) ;
      if ( rp == 0x0 ) { printf("\n\n *** Missing hist %s\n\n", hname ) ; gSystem -> Exit(-1) ; }
      return rp ;
   }

   TH2F* fget_hist2d( const char* hname, TFile* tfp ) {
      if ( tfp == 0x0 ) { printf("\n\n *** fget_hist : bad file pointer.\n\n") ; gSystem -> Exit(-1) ; }
      TH2F* rp = (TH2F*) tfp -> Get( hname ) ;
      if ( rp == 0x0 ) { printf("\n\n *** Missing hist %s\n\n", hname ) ; gSystem -> Exit(-1) ; }
      return rp ;
   }

   TGraph* get_graph( const char* gname ) {
      TGraph* rp = (TGraph*) gDirectory -> Get( gname ) ;
      if ( rp == 0x0 ) { printf("\n\n *** Missing TGraph %s\n\n", gname ) ; gSystem -> Exit(-1) ; }
      return rp ;
   }

   void grid_on() {
      gPad -> SetGridx(1) ;
      gPad -> SetGridy(1) ;
   }

   void SetupCorrelationPalette() {

      static Bool_t initialized = kFALSE ;
      static Int_t colors[90] ;

      const Int_t Number = 3 ;
      Double_t Length[Number] = {0.,0.5, 1.} ;
      Double_t Red[Number] = {0.,1.,1.} ;
      Double_t Green[Number] = {0.,1.,0.} ;
      Double_t Blue[Number] = {1.,1.,0.} ;

      if ( !initialized ) {
         Int_t fi = TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,90);
         for ( int i=0; i<90; i++ ) colors[i] = fi+i ;
         initialized = kTRUE ;
         return ;
      }
      gStyle -> SetPalette( 90, colors ) ;

   }

   void Setup2DhistPalette() {
      gStyle -> SetPalette( kBird ) ;
   }




#endif

















