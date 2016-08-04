void compareCSVS(TString inputfile, TString tag){
 

    TString mkdir = ".! mkdir -p plots/plots_"+tag;
    gROOT->ProcessLine(mkdir);

    TString rocCurveDefault = "DQMData/Run 1/Btag/Run summary/CSVv2_GLOBAL/FlavEffVsBEff_DUSG_discr_CSVv2_GLOBAL";
    TFile* f = new TFile(inputfile);
    TH1F* hDefault;
    
    gDirectory->GetObject(rocCurveDefault,hDefault);

    Float_t xmin = 0.0;
    Float_t xmax = 0.2;
    Int_t nstep = 20;
    
    for (Int_t i = 0; i <= nstep; i++)
    {
      TString rocCurveCustom; 
      rocCurveCustom.Form("DQMData/Run 1/Btag/Run summary/CSVv2Custom%i_GLOBAL/FlavEffVsBEff_DUSG_discr_CSVv2Custom%i_GLOBAL",i, i); 
      //cout<<rocCurveCustom<<endl;

    // Check if rocCurveDefault and rocCurveCustom are in the file

      Float_t cutval = xmin + (xmax - xmin)/Float_t(nstep)*i;
      TString legString;
      legString.Form("cut value = %.3f",cutval);
     
      TString cutvalString;
      cutvalString.Form("%.3i",i);
      /*cutvalString.Form("%.2f",cutval);
      cutvalString.ReplaceAll(".","_");
 
      if(cutval >= 0) 
        cutvalString = "pos"+ cutvalString;
      else
        cutvalString.ReplaceAll("-","neg");
      */

      cout<<cutvalString<<endl;
      //gDirectory->pwd();
      TH1F* hCustom;
      gDirectory->GetObject(rocCurveCustom,hCustom);

      hCustom->SetLineColor(kRed);
      hCustom->SetMarkerColor(kRed);

      Double_t w = 1200;
      Double_t h = 800;
      TCanvas c1("c", "c", w, h);
      c1.SetLogy();
    
      hDefault->GetXaxis()->SetRangeUser(0.5,1.0);
      hDefault->Draw();
      hCustom->Draw("same");

      TLegend leg(0.1,0.7,0.48,0.9);
      leg.AddEntry(hDefault,"Default","l");
      leg.AddEntry(hCustom,legString,"l");
      leg.Draw();

      c1.Print("plots/plots_"+tag+"/"+tag+"_"+cutvalString+".png","png");
      
    }
    TString tarzip = ".! tar czvf plots/plots_"+tag+".tgz plots/plots_"+tag+"/*";
    gROOT->ProcessLine(tarzip);
         

}
