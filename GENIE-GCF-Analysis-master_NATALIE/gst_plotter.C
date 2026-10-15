#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include "TLatex.h"
#include "TVector3.h" 

#include <iostream>
#include <vector>

void PcmPlots(vector<TH1F*> x_hists, vector<TH1F*> y_hists, vector<TH1F*> z_hists,vector<TH1F*> data_hist, vector<TH1F*> mag_hists,TString varName, TLegend *legend,TString CutName, int order){

        double FontSize = 0.08;
        double label_size = FontSize;
        double title_size = FontSize*1.2;
     //   int FontStyle = 132;
        int FontStyle = 62;
   

        gStyle->SetOptStat(0);

        TString x_axis_title_x = "P_{c.m.}^{X}";
        TString x_axis_title_y = "P_{c.m.}^{Y}";
        TString x_axis_title_z = "P_{c.m.}^{Z}";
        TString x_axis_title_mag = "|P_{c.m.}|";

        TF1 COMGaus("COMGaus","gaus",-.5,.5);
        COMGaus.SetParameters(1,0,0.15); 
  
	TCanvas* canvas = new TCanvas("canvas_"+varName+CutName,"canvas_"+varName+CutName,205,34,3600,1600);				// ---------------------------------------------------------------------------------------------------------------------------				
        TPad* pad1 = nullptr;
	TPad* pad2 = nullptr;
	TPad* pad3 = nullptr;
	TPad* pad4 = nullptr;
	TPad* pad5 = nullptr;
	TPad* pad6 = nullptr;
	TPad* pad7 = nullptr;
	TPad* pad8 = nullptr;
        
        pad5 = new TPad(varName+CutName+"_5",varName+CutName+"_5",0.0,.04,0.27,.52,22);
	pad5->SetFillColor(kWhite); pad5->Draw(); 
        pad6 = new TPad(varName+CutName+"_6",varName+CutName+"_6",0.27,.04,.5,.52,22);
	pad6->SetFillColor(kWhite); pad6->Draw(); 
        pad7 = new TPad(varName+CutName+"_7",varName+CutName+"_7",0.5,.04,.73,.52,22);
	pad7->SetFillColor(kWhite); pad7->Draw(); 
        pad8 = new TPad(varName+CutName+"_8",varName+CutName+"_8",0.73,.04,1.,.52,22);
	pad8->SetFillColor(kWhite); pad8->Draw(); 

        pad1 = new TPad(varName+CutName+"_1",varName+CutName+"_1",0.0,.52,.27,1.,21);
	pad1->SetFillColor(kWhite); pad1->Draw();
	pad2 = new TPad(varName+CutName+"_2",varName+CutName+"_2",0.27,.52,0.5,1,22);
	pad2->SetFillColor(kWhite); pad2->Draw();
	pad3 = new TPad(varName+CutName+"_3",varName+CutName+"_3",0.5,.52,.73,1.,22);
	pad3->SetFillColor(kWhite); pad3->Draw(); 
        pad4 = new TPad(varName+CutName+"_4",varName+CutName+"_4",0.73,.52,1.,1.,22);
	pad4->SetFillColor(kWhite); pad4->Draw(); 					
	
    //    pad1->SetBottomMargin(0.1);
//	pad2->SetBottomMargin(0.1);
//	pad3->SetBottomMargin(0.1);

	pad7->SetBottomMargin(0.25);					
	pad8->SetBottomMargin(0.25);					
	pad5->SetBottomMargin(0.25);					
	pad6->SetBottomMargin(0.25);
        
        pad1->SetTopMargin(.18);
        pad2->SetTopMargin(.18);
        pad3->SetTopMargin(.18);
        pad4->SetTopMargin(.18);
        
        pad1->SetBottomMargin(0.);					
	pad2->SetBottomMargin(0.);					
	pad3->SetBottomMargin(0.);
	pad4->SetBottomMargin(0.);
        pad7->SetTopMargin(0.);				
        pad8->SetTopMargin(0.);				
        pad5->SetTopMargin(0.);				
        pad6->SetTopMargin(0.);				
	
        pad1->SetRightMargin(0.);
        pad2->SetRightMargin(0.);
        pad3->SetRightMargin(0.);
        pad5->SetRightMargin(0.);
        pad6->SetRightMargin(0.);
        pad7->SetRightMargin(0.);
	
        pad2->SetLeftMargin(0);
	pad3->SetLeftMargin(0);
	pad4->SetLeftMargin(0);
	pad6->SetLeftMargin(0);
	pad7->SetLeftMargin(0);
	pad8->SetLeftMargin(0);

        pad4->SetRightMargin(0.25);
        pad8->SetRightMargin(0.25);
	pad1->SetLeftMargin(0.25);
	pad5->SetLeftMargin(0.25);

        double bordersize = 2.;
        pad1->SetFrameLineWidth(bordersize);
        pad2->SetFrameLineWidth(bordersize);
        pad3->SetFrameLineWidth(bordersize);
        pad4->SetFrameLineWidth(bordersize);
        pad5->SetFrameLineWidth(bordersize);
        pad6->SetFrameLineWidth(bordersize);
        pad7->SetFrameLineWidth(bordersize);
        pad8->SetFrameLineWidth(bordersize);
   //     canvas->Divide(2,1);

        gStyle->SetOptTitle(0);
        pad1->cd();


        x_hists[0]->GetYaxis()->SetTitleFont(FontStyle);
        x_hists[0]->GetYaxis()->SetTitleSize(title_size);
        x_hists[0]->GetYaxis()->SetTitleOffset(1.);

        x_hists[0]->GetYaxis()->SetTitle("Counts");

      //  x_hists[0]->GetYaxis()->SetRangeUser(0,1.05*x_hists[0]->GetMaximum());
        double scaley;
        double scaley_acc = 1.2*mag_hists[2]->GetMaximum();
        double scaley_epp = 1.2*mag_hists[4]->GetMaximum();

        scaley = scaley_acc;
        x_hists[0]->GetYaxis()->SetRangeUser(0,scaley);

        

        x_hists[0]->GetYaxis()->CenterTitle();        

        x_hists[0]->GetXaxis()->SetTitle("");
        x_hists[0]->GetXaxis()->SetTitleFont(FontStyle);
        x_hists[0]->GetXaxis()->SetTitleSize(title_size);
        x_hists[0]->GetXaxis()->SetTitleOffset(.9);        
        x_hists[0]->GetXaxis()->CenterTitle(); 

        if(varName != "pRec")
       // x_hists[0]->GetXaxis()->ChangeLabel(-1,-1,0,-1,-1,-1,"");


        x_hists[0]->GetYaxis()->SetTitleFont(FontStyle);
        x_hists[0]->GetYaxis()->SetTitleSize(title_size);
        x_hists[0]->GetYaxis()->SetTitleOffset(1.2);
        x_hists[0]->GetYaxis()->SetTitle("Counts");
        x_hists[0]->SetLabelFont(FontStyle,"xy");
        x_hists[0]->SetLabelSize(label_size,"xy");


        x_hists[0]->GetYaxis()->SetNdivisions(6);
        x_hists[0]->GetXaxis()->SetNdivisions(6);

        x_hists[0]->Draw("SAME C HIST");
     //   x_hists[1]->Draw("SAME C HIST"); 
        x_hists[2]->Draw("SAME C HIST");
        if(order == 1){
     //           data_hist[3]->Draw("SAME PE X0");
        }

        legend->SetBorderSize(0);
        legend->SetTextSize(.08);
        legend->SetTextFont(FontStyle-20);
        



      
        pad2->cd();



      //  y_hists[0]->GetYaxis()->SetRangeUser(0,1.05*y_hists[0]->GetMaximum());
        y_hists[0]->GetYaxis()->SetRangeUser(0,scaley);



        y_hists[0]->GetYaxis()->SetTitleFont(FontStyle);
        y_hists[0]->GetYaxis()->SetTitleSize(title_size);
        y_hists[0]->GetYaxis()->SetTitleOffset(.95);


        y_hists[0]->GetYaxis()->CenterTitle();        

        y_hists[0]->GetXaxis()->SetTitle("");
        y_hists[0]->GetXaxis()->SetTitleFont(FontStyle);
        y_hists[0]->GetXaxis()->SetTitleSize(title_size);
        y_hists[0]->GetXaxis()->SetTitleOffset(.9);        
        y_hists[0]->GetXaxis()->CenterTitle();   
     
        if(varName != "pRec")
          //      y_hists[0]->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,"");

        
        
        y_hists[0]->SetLabelSize(label_size,"xy");
        y_hists[0]->SetLabelFont(FontStyle,"xy");

        std::cout << order << std::endl;

        y_hists[0]->GetYaxis()->SetNdivisions(6);
        y_hists[0]->GetXaxis()->SetNdivisions(6);

        y_hists[0]->Draw("SAME C HIST Y+"); 
   //     y_hists[1]->Draw("SAME C HIST Y+");
        y_hists[2]->Draw("SAME C HIST Y+");
     //   if(order == 1) data_hist[1]->Draw("SAME PE Y+ X0");
        legend->Draw("SAME");
 


        pad3->cd();

       // z_hists[0]->GetYaxis()->SetRangeUser(0,1.05*z_hists[0]->GetMaximum());
        z_hists[0]->GetYaxis()->SetRangeUser(0,scaley);



        z_hists[0]->GetYaxis()->SetTitleFont(FontStyle);
        z_hists[0]->GetYaxis()->SetTitleSize(title_size);
        z_hists[0]->GetYaxis()->SetTitleOffset(.95);


        z_hists[0]->GetYaxis()->CenterTitle();        

        z_hists[0]->GetXaxis()->SetTitle("");
        z_hists[0]->GetXaxis()->SetTitleFont(FontStyle);
        z_hists[0]->GetXaxis()->SetTitleSize(title_size);
        z_hists[0]->GetXaxis()->SetTitleOffset(.9);        
        z_hists[0]->GetXaxis()->CenterTitle();   
     
        if(varName != "pRec")
                z_hists[0]->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,"");
                x_hists[0]->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,"");
                mag_hists[0]->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,"");

        
        z_hists[0]->GetXaxis()->SetRangeUser(-.3,.7);
        
        z_hists[0]->SetLabelSize(label_size,"xy");
        z_hists[0]->SetLabelFont(FontStyle,"xy");

        std::cout << order << std::endl;

        z_hists[0]->GetYaxis()->SetNdivisions(6);
        z_hists[0]->GetXaxis()->SetNdivisions(6);


        z_hists[0]->Draw("SAME C HIST Y+"); 
    //    z_hists[1]->Draw("SAME C HIST Y+");
        z_hists[2]->Draw("SAME C HIST Y+");
      //  if(order == 1) data_hist[2]->Draw("SAME PE Y+ X0");

        TLatex* cut = new TLatex(0.1,0.1,"");
        cut->SetTextFont(FontStyle);
        cut->SetTextSize(0.08);
        cut->DrawLatexNDC(0.6,0.7,"^{12}C(e,e'pp)");

        pad5->cd();

        scaley = scaley_epp;

        x_hists[3]->GetYaxis()->SetTitleFont(FontStyle);
        x_hists[3]->GetYaxis()->SetTitleSize(title_size);
        x_hists[3]->GetYaxis()->SetTitleOffset(1.);

        x_hists[3]->GetYaxis()->SetTitle("Counts");

      //  x_hists[3]->GetYaxis()->SetRangeUser(0,1.05*x_hists[3]->GetMaximum());
        x_hists[3]->GetYaxis()->SetRangeUser(0,scaley);

        

        x_hists[3]->GetYaxis()->CenterTitle();        

        x_hists[3]->GetXaxis()->SetTitle(x_axis_title_x);
        x_hists[3]->GetXaxis()->SetTitleFont(FontStyle);
        x_hists[3]->GetXaxis()->SetTitleSize(title_size);
        x_hists[3]->GetXaxis()->SetTitleOffset(1.2);        
        x_hists[3]->GetXaxis()->CenterTitle(); 

        if(varName != "pRec")
       // x_hists[3]->GetXaxis()->ChangeLabel(-1,-1,0,-1,-1,-1,"");


        x_hists[3]->GetYaxis()->SetTitleFont(FontStyle);
        x_hists[3]->GetYaxis()->SetTitleSize(title_size);
        x_hists[3]->GetYaxis()->SetTitleOffset(1.2);
        x_hists[3]->GetYaxis()->SetTitle("");
        x_hists[3]->SetLabelFont(FontStyle,"xy");
        x_hists[3]->SetLabelSize(label_size,"xy");


        x_hists[3]->GetYaxis()->SetNdivisions(6);
        x_hists[3]->GetXaxis()->SetNdivisions(6);

        x_hists[3]->Draw("SAME C HIST");
        x_hists[4]->Draw("SAME C HIST"); 
        x_hists[5]->Draw("SAME C HIST");
        if(order == 1){
                data_hist[0]->Draw("SAME PE X0");
        }

        legend->SetBorderSize(0);
        legend->SetTextSize(.08);
        legend->SetTextFont(FontStyle-20);
        



      
        pad6->cd();



      //  y_hists[3]->GetYaxis()->SetRangeUser(0,1.05*y_hists[3]->GetMaximum());
        y_hists[3]->GetYaxis()->SetRangeUser(0,scaley);



        y_hists[3]->GetYaxis()->SetTitleFont(FontStyle);
        y_hists[3]->GetYaxis()->SetTitleSize(title_size);
        y_hists[3]->GetYaxis()->SetTitleOffset(.95);


        y_hists[3]->GetYaxis()->CenterTitle();        

        y_hists[3]->GetXaxis()->SetTitle(x_axis_title_y);
        y_hists[3]->GetXaxis()->SetTitleFont(FontStyle);
        y_hists[3]->GetXaxis()->SetTitleSize(title_size);
        y_hists[3]->GetXaxis()->SetTitleOffset(1.2);        
        y_hists[3]->GetXaxis()->CenterTitle();   
     
        if(varName != "pRec")
          //      y_hists[3]->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,"");

        
        
        y_hists[3]->SetLabelSize(label_size,"xy");
        y_hists[3]->SetLabelFont(FontStyle,"xy");

        std::cout << order << std::endl;

        y_hists[3]->GetYaxis()->SetNdivisions(6);
        y_hists[3]->GetXaxis()->SetNdivisions(6);

        y_hists[3]->Draw("SAME C HIST Y+"); 
        y_hists[4]->Draw("SAME C HIST Y+");
        y_hists[5]->Draw("SAME C HIST Y+");
        if(order == 1) data_hist[1]->Draw("SAME PE Y+ X0");
      //  legend->Draw("SAME");
 


        pad7->cd();

       // z_hists[3]->GetYaxis()->SetRangeUser(0,1.05*z_hists[3]->GetMaximum());
        z_hists[3]->GetYaxis()->SetRangeUser(0,scaley);



        z_hists[3]->GetYaxis()->SetTitleFont(FontStyle);
        z_hists[3]->GetYaxis()->SetTitleSize(title_size);
        z_hists[3]->GetYaxis()->SetTitleOffset(.95);


        z_hists[3]->GetYaxis()->CenterTitle();        

        z_hists[3]->GetXaxis()->SetTitle(x_axis_title_z);
        z_hists[3]->GetXaxis()->SetTitleFont(FontStyle);
        z_hists[3]->GetXaxis()->SetTitleSize(title_size);
        z_hists[3]->GetXaxis()->SetTitleOffset(1.2);        
        z_hists[3]->GetXaxis()->CenterTitle();   
     
        if(varName != "pRec")
                z_hists[3]->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,"");

        
        z_hists[3]->GetXaxis()->SetRangeUser(-.3,.7);
        
        z_hists[3]->SetLabelSize(label_size,"xy");
        z_hists[3]->SetLabelFont(FontStyle,"xy");

        std::cout << order << std::endl;

        z_hists[3]->GetYaxis()->SetNdivisions(6);
        z_hists[3]->GetXaxis()->SetNdivisions(6);

        z_hists[3]->Draw("SAME C HIST Y+"); 
        z_hists[4]->Draw("SAME C HIST Y+");
        z_hists[5]->Draw("SAME C HIST Y+");
        if(order == 1) data_hist[2]->Draw("SAME PE Y+ X0");

        pad4->cd();

        scaley = scaley_acc;


        mag_hists[0]->GetYaxis()->SetTitleFont(FontStyle);
        mag_hists[0]->GetYaxis()->SetTitleSize(title_size);
        mag_hists[0]->GetYaxis()->SetTitleOffset(1.);

        mag_hists[0]->GetYaxis()->SetTitle("");

      //  mag_hists[0]->GetYaxis()->SetRangeUser(0,1.05*mag_hists[0]->GetMaximum());
        mag_hists[0]->GetYaxis()->SetRangeUser(0,scaley);

        

        mag_hists[0]->GetYaxis()->CenterTitle();        

        mag_hists[0]->GetXaxis()->SetTitle("");
        mag_hists[0]->GetXaxis()->SetTitleFont(FontStyle);
        mag_hists[0]->GetXaxis()->SetTitleSize(title_size);
        mag_hists[0]->GetXaxis()->SetTitleOffset(.9);        
        mag_hists[0]->GetXaxis()->CenterTitle(); 

        if(varName != "pRec")
       // mag_hists[0]->GetYaxis()->ChangeLabel(-1,-1,0,-1,-1,-1,"");


        mag_hists[0]->GetYaxis()->SetTitleFont(FontStyle);
        mag_hists[0]->GetYaxis()->SetTitleSize(title_size);
        mag_hists[0]->GetYaxis()->SetTitleOffset(1.2);
        mag_hists[0]->GetYaxis()->SetTitle("");
        mag_hists[0]->SetLabelFont(FontStyle,"xy");
        mag_hists[0]->SetLabelSize(label_size,"xy");


        mag_hists[0]->GetYaxis()->SetNdivisions(6);
        mag_hists[0]->GetXaxis()->SetNdivisions(6);

        mag_hists[0]->Draw("SAME C HIST Y+");
     //   mag_hists[1]->Draw("SAME C HIST"); 
        mag_hists[2]->Draw("SAME C HIST");
        if(order == 1){
     //           data_hist[0]->Draw("SAME PE X0");
        }

        pad8->cd();

        scaley = scaley_epp;

        mag_hists[3]->GetYaxis()->SetTitleFont(FontStyle);
        mag_hists[3]->GetYaxis()->SetTitleSize(title_size);
        mag_hists[3]->GetYaxis()->SetTitleOffset(1.);

        mag_hists[3]->GetYaxis()->SetTitle("");

      //  mag_hists[3]->GetYaxis()->SetRangeUser(0,1.05*mag_hists[3]->GetMaximum());
        mag_hists[3]->GetYaxis()->SetRangeUser(0,scaley);

        

        mag_hists[3]->GetYaxis()->CenterTitle();        

        mag_hists[3]->GetXaxis()->SetTitle(x_axis_title_mag);
        mag_hists[3]->GetXaxis()->SetTitleFont(FontStyle);
        mag_hists[3]->GetXaxis()->SetTitleSize(title_size);
        mag_hists[3]->GetXaxis()->SetTitleOffset(1.2);        
        mag_hists[3]->GetXaxis()->CenterTitle(); 

        if(varName != "pRec")
         mag_hists[3]->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,"");



        mag_hists[3]->GetYaxis()->SetTitleFont(FontStyle);
        mag_hists[3]->GetYaxis()->SetTitleSize(title_size);
        mag_hists[3]->GetYaxis()->SetTitleOffset(1.2);
        mag_hists[3]->GetYaxis()->SetTitle("");
        mag_hists[3]->SetLabelFont(FontStyle,"xy");
        mag_hists[3]->SetLabelSize(label_size,"xy");


        mag_hists[3]->GetYaxis()->SetNdivisions(6);
        mag_hists[3]->GetXaxis()->SetNdivisions(6);

        mag_hists[3]->Draw("SAME C HIST Y+");
        mag_hists[4]->Draw("SAME C HIST"); 
        mag_hists[5]->Draw("SAME C HIST");
        if(order == 1){
                data_hist[3]->Draw("SAME PE X0");
        }
        



        canvas->Update();

        if(order == 0)
                canvas->SaveAs("Plots/pcom_components_acceptance_Overlay_FS_1D_EMQE_carbon.pdf");
        if(order == 1)
                canvas->SaveAs("Plots/pcom_components_epp_Overlay_FS_1D_EMQE_carbon.pdf");

     /*   cout << "PWIA X: " << endl;
        x_hists[3]->Fit("COMGaus");
        cout << "PWIA ACC X: " << endl;
        x_hists[0]->Fit("COMGaus");
        cout << "FULL ACC X: " << endl;
        x_hists[2]->Fit("COMGaus");
        cout << "FULL X: " << endl;
        x_hists[5]->Fit("COMGaus");
       /* cout << "PWIA Y: " << endl;
        y_hists[3]->Fit("COMGaus");
        cout << "FULL X: " << endl;
        x_hists[5]->Fit("COMGaus");
        */
}

void OverlayPlots(TH1F *hist_1, TH1F *hist_2, TH1F *hist_3,TH1F *hist_4,TH1F *hist_5,TH1F *hist_6, TH1F *hist_7,TString varName, TLegend *legend,TString CutName, int order){

        double FontSize = 0.08;
        double label_size = .8*FontSize;
        double title_size = FontSize;
     //   int FontStyle = 132;
        int FontStyle = 62;
        bool pn_ep = true;
        bool plead_ep = false;

        gStyle->SetOptStat(0);

     
        TString x_axis_title = "";
        if(varName == "pMiss"){
                x_axis_title = "P_{miss} [GeV/c]";
        } else if(varName == "EMiss"){
                x_axis_title = "E_{miss} [GeV]";
        } else if(varName == "pLead"){
                x_axis_title = "P_{Lead} [GeV/c]";
        } else if(varName == "pRec"){
                x_axis_title = "P_{Recoil} [GeV/c]";
        } else if(varName == "pcom_x"){
                x_axis_title = "P_{COM,x}";
               
        } else if(varName == "pcom_y"){
                x_axis_title = "P_{COM,y}";
              
        } else if(varName == "pcom_z"){
                x_axis_title = "P_{COM,z}";
             /*   hist_1->Rebin();
                hist_2->Rebin();
                hist_3->Rebin();*/
        } else if(varName == "prel"){
                x_axis_title = "P_{REL} [GeV/c]";
               // hist_1->GetXaxis()->SetRangeUser(0.,1.1);

        } else if(varName == "cth_prel"){
                x_axis_title = "P_{REL} Cos(#theta)";
                hist_1->Rebin();
                hist_2->Rebin();
                hist_3->Rebin();
             /*   hist_4->Rebin();
                hist_5->Rebin();
                */
        } else if(varName == "cth_pcom"){
                x_axis_title = "P_{com} Cos(#theta)";
                hist_1->Rebin();
                hist_2->Rebin();
                hist_3->Rebin();
            /*    hist_4->Rebin();
                hist_5->Rebin();
                */
        } else if(varName == "phi_pcom"){
                x_axis_title = "P_{com} #phi";
                hist_1->Rebin();
                hist_2->Rebin();
                hist_3->Rebin();
              /*  hist_4->Rebin();
                hist_5->Rebin();
                */
        } else if(varName == "phi_prel"){
                x_axis_title = "P_{rel} #phi";
                hist_1->Rebin();
                hist_2->Rebin();
                hist_3->Rebin();
              /*  hist_4->Rebin();
                hist_5->Rebin();
                */
        } else if(varName == "pq"){
                x_axis_title = "|p|/|q| [GeV/c]";
        }  else if(varName == "theta_pq"){
                x_axis_title = "#theta(P_{Lead},q) [Deg.]";
        }  else if(varName == "cthl"){
                x_axis_title = "Lepton Cos(#theta)";
                hist_1->Rebin();
                hist_2->Rebin();
                hist_3->Rebin();
             /*   hist_4->Rebin();
                hist_5->Rebin();
                */
              /*  hist_1->GetXaxis()->SetRangeUser(.985,1);
                hist_2->GetXaxis()->SetRangeUser(.985,1);
                hist_3->GetXaxis()->SetRangeUser(.985,1);
*/
        }  else if(varName == "th_pmiq"){
                x_axis_title = "#theta(P_{Miss},q) [Deg.]";
        }   else if(varName == "pn"){
                x_axis_title = "P_{Miss} [GeV/c]";
            /*   hist_1->GetXaxis()->SetRangeUser(0,1);
                hist_2->GetXaxis()->SetRangeUser(0,1);
                hist_3->GetXaxis()->SetRangeUser(0,1);
                */
        }    else if(varName == "En"){
                x_axis_title = "Struck Nucleon Energy [GeV/c]";
                hist_1->GetXaxis()->SetRangeUser(0,1);
                hist_2->GetXaxis()->SetRangeUser(0,1);
                hist_3->GetXaxis()->SetRangeUser(0,1);
        }     else if(varName == "pnx"){
                x_axis_title = "(Struck Nucleon Momentum)_{x} [GeV/c]";
               
        }     else if(varName == "pny"){
                x_axis_title = "(Struck Nucleon Momentum)_{y} [GeV/c]";
               
        }     else if(varName == "pnz"){
                x_axis_title = "(Struck Nucleon Momentum)_{z} [GeV/c]";
               
        }       else if(varName == "p4n"){
                x_axis_title = "Struck Nucleon Invariant Mass";
               hist_1->GetXaxis()->SetRangeUser(0,1);
                hist_2->GetXaxis()->SetRangeUser(0,1);
                hist_3->GetXaxis()->SetRangeUser(0,1);
        }       else if(varName == "qx"){
                x_axis_title = "q_{x} [GeV/c]";    
        }        else if(varName == "qy"){
                x_axis_title = "q_{y} [GeV/c]";    
        }         else if(varName == "qz"){
                x_axis_title = "q_{z} [GeV/c]";    
        }    else if(varName == "spectral_func"){
                x_axis_title = "Struck Nucleon/Cluster Momentum [GeV/c]";
        }   else if(varName == "pcom"){
                x_axis_title = "P_{c.m.} [GeV/c]";
        }    else if(varName == "pMiss"){
                x_axis_title = "|P_{Miss}| [GeV/c]";
        }    else if(varName == "pMiss_epp"){
                x_axis_title = "|P_{Miss}| [GeV/c]";
        }   else if(varName == "pcom_q"){
                x_axis_title = "P_{COM} in Direction of q [GeV/c]";
            

        } else if(varName == "pcom_qt"){
                x_axis_title = "P_{COM} in Direction of q_{T} [GeV/c]";   
        } else if(varName == "EmissMPmiss"){
                x_axis_title = "E_{Miss} - P_{Miss} in Direction of q";   
        }  else if(varName == "Emiss2MPmiss2"){
                x_axis_title = "E^{2}_{Miss} - P^{2}_{Miss} in Direction of q";   
        }  else if(varName == "pn_genie"){
                x_axis_title = "Struck Nucleon Momentum (From GENIE)";   
        }                    

  

        TCanvas *canvas = new TCanvas("canvas_"+varName+CutName+"_1"," ");
        gROOT->ForceStyle(true);
        TPad *pad1 = new TPad(varName+CutName+"_1",varName+CutName+"_1",0,0,0.5,1., 21);
	pad1->SetFillColor(kWhite); pad1->Draw();
	TPad *pad2 = new TPad(varName+CutName+"_2",varName+CutName+"_2",0.5,0,1.,1,22);
	pad2->SetFillColor(kWhite); pad2->Draw();
	pad1->SetBottomMargin(0.18);
	pad2->SetBottomMargin(0.18);
	
        pad1->SetRightMargin(0.);
	pad2->SetLeftMargin(0);

        pad2->SetRightMargin(0.25);
	pad1->SetLeftMargin(0.25);

        double bordersize = 2.;
        pad1->SetFrameLineWidth(bordersize);
        pad2->SetFrameLineWidth(bordersize);

   //     canvas->Divide(2,1);

        gStyle->SetOptTitle(0);
        pad1->cd();

        

        if(order == 2) hist_1 = hist_2;

        if(varName == "pLead"){
                if(plead_ep == false){
                        hist_1->GetXaxis()->SetRangeUser(.75,2.4);
                        hist_5->GetXaxis()->SetRangeUser(.75,2.4);
                } else {
                        hist_1->GetXaxis()->SetRangeUser(.75,2.4);
                        hist_5->GetXaxis()->SetRangeUser(.75,2.4);
                }
        }

        hist_1->GetYaxis()->SetTitleFont(FontStyle);
        hist_1->GetYaxis()->SetTitleSize(title_size);
        hist_1->GetYaxis()->SetTitleOffset(1.);

        hist_1->GetYaxis()->SetTitle("Counts");

        
        if(varName == "pn"){
            //    hist_1->GetYaxis()->SetTitle("(e,e'pp)/(e,e'p) Events");

                hist_1->GetYaxis()->SetRangeUser(hist_1->GetMinimum(),1.6*hist_1->GetMaximum());
                gPad->SetLogy();
        } else {

        hist_1->GetYaxis()->SetRangeUser(0,1.05*hist_1->GetMaximum());
        if(varName == "pLead" || varName == "pcom_x" || varName == "pcom_y" ||varName == "pcom_z")
                hist_1->GetYaxis()->SetRangeUser(0,1.2*hist_1->GetMaximum());

        }        

        
      //  if(varName == "pMiss")
        //        hist_1->GetYaxis()->SetRangeUser(0,6000);

        if(varName == "prel"){
                hist_1->GetXaxis()->SetRangeUser(.2,1.);
                hist_5->GetXaxis()->SetRangeUser(.2,1.);
        }

      //  if(varName == "pLead") hist_1->GetYaxis()->SetRangeUser(0,450);


        hist_1->GetYaxis()->CenterTitle();        

        hist_1->GetXaxis()->SetTitle(x_axis_title);
        hist_1->GetXaxis()->SetTitleFont(FontStyle);
        hist_1->GetXaxis()->SetTitleSize(title_size);
        hist_1->GetXaxis()->SetTitleOffset(.9);        
        hist_1->GetXaxis()->CenterTitle(); 

        if(varName == "pLead") hist_1->GetYaxis()->SetRangeUser(0,420);

        if( varName != "pLead")
        hist_1->GetXaxis()->ChangeLabel(-1,-1,0,-1,-1,-1,"");


        hist_1->GetYaxis()->SetTitleFont(FontStyle);
        hist_1->GetYaxis()->SetTitleSize(title_size);
        hist_1->GetYaxis()->SetTitleOffset(1.5);
        hist_1->GetYaxis()->SetTitle("Counts");
        hist_1->SetLabelFont(FontStyle,"xy");
        hist_1->SetLabelSize(label_size,"xy");


        hist_1->GetYaxis()->SetNdivisions(6);
        hist_1->GetXaxis()->SetNdivisions(6);

        if(order != 3){

                hist_1->Draw("SAME C HIST");
                hist_3->Draw("SAME C HIST");

        } else {

                hist_1->Draw("SAME C HIST");
                hist_2->Draw("SAME C HIST");
                
        }
       
        
        

        legend->SetBorderSize(0);
        legend->SetTextSize(.08);
        legend->SetTextFont(FontStyle-20);
        legend->Draw("SAME");


        pad2->cd();


        if(hist_4->GetMaximum() < hist_5->GetMaximum()) hist_5->GetYaxis()->SetRangeUser(0,1.05*hist_5->GetMaximum());
        else hist_5->GetYaxis()->SetRangeUser(0,1.05*hist_4->GetMaximum()+hist_4->GetBinError(hist_4->GetMaximumBin()));

        if(varName == "pn"){
                gPad->SetLogy();
                hist_5->GetYaxis()->SetRangeUser(.9,1.05*hist_5->GetMaximum());
                if(pn_ep == true)
                        hist_5->GetYaxis()->SetRangeUser(2,1.05*hist_5->GetMaximum());
        } 

        if(varName == "th_pmiq"){
               hist_5->GetYaxis()->SetRangeUser(0,1.05*hist_5->GetMaximum()+hist_4->GetBinError(hist_4->GetMaximumBin()));
        }


        hist_5->GetYaxis()->SetTitleFont(FontStyle);
        hist_5->GetYaxis()->SetTitleSize(title_size);
        hist_5->GetYaxis()->SetTitleOffset(.95);


        hist_5->GetYaxis()->CenterTitle();        

        hist_5->GetXaxis()->SetTitle(x_axis_title);
        hist_5->GetXaxis()->SetTitleFont(FontStyle);
        hist_5->GetXaxis()->SetTitleSize(title_size);
        hist_5->GetXaxis()->SetTitleOffset(.9);        
        hist_5->GetXaxis()->CenterTitle();   
     
        if( varName != "pRec") hist_5->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,"");

        
        
        hist_5->SetLabelSize(label_size,"xy");
        hist_5->SetLabelFont(FontStyle,"xy");


        hist_5->GetYaxis()->SetNdivisions(6);
        hist_5->GetXaxis()->SetNdivisions(6);

        hist_5->Draw("SAME C HIST Y+");
        hist_6->Draw("SAME C HIST Y+");        
        hist_7->Draw("SAME C HIST Y+");
        hist_4->Draw("SAME PE Y+ X0");
        
     


        TLatex* cut = new TLatex(0.1,0.1,"");
        cut->SetTextFont(FontStyle);
        cut->SetTextSize(0.08);
        if(varName == "pRec" || varName == "pcom" || varName == "pcom_x" || varName == "pcom_y" || varName == "pcom_z" || (varName == "pn" && pn_ep == false) || varName == "prel" || (varName == "pLead" && plead_ep == false))
                cut->DrawLatexNDC(0.35,0.82,"^{12}C(e,e'pp)");
        else
                cut->DrawLatexNDC(0.38,0.82,"^{12}C(e,e'p)");


        canvas->Update();

        if(varName == "th_pmiq") varName = "theta_pmiq";


        if(varName == "pn" && pn_ep ){ canvas->SaveAs("Plots/" + varName +CutName+"_ep_Overlay_FS_1D_EMQE_carbon.pdf");}
        else {if(varName == "pLead" && !plead_ep ) canvas->SaveAs("Plots/" + varName +CutName+"_epp_Overlay_FS_1D_EMQE_carbon.pdf");
        else canvas->SaveAs("Plots/" + varName +CutName+"_Overlay_FS_1D_EMQE_carbon.pdf");
        }
}



void Hist_draw_order(TFile *inFile_1, TFile *inFile_2, TFile *inFile_3,TFile *inFile_4, vector<double> norm_factors,/*TFile *inFile_5,*/TString VarName, TString cutName){


        TH1F *h_1;
        TH1F *h_2; 
        TH1F *h_3;
        TH1F *h_4;
        TH1F *h_5;
        TH1F *h_6;
        TH1F *h_7;


        //cutName is acceptance for paper    
        h_1 = (TH1F*)inFile_1->Get(VarName+cutName+"_hist");
        h_2 = (TH1F*)inFile_2->Get(VarName+cutName+"_hist");
        h_3 = (TH1F*)inFile_3->Get(VarName+cutName+"_hist");

        h_4 = (TH1F*)inFile_4->Get(VarName+"_ep"+"_hist");
        h_5 = (TH1F*)inFile_1->Get(VarName+"_ep"+"_hist");
        h_6 = (TH1F*)inFile_2->Get(VarName+"_ep"+"_hist");
        h_7 = (TH1F*)inFile_3->Get(VarName+"_ep"+"_hist");


        


        bool epp = false;
        if(VarName == "pRec" || VarName == "pLead" || VarName == "pcom" || VarName == "pcom_z" || VarName == "pcom_y"|| VarName == "pcom_x" || VarName == "prel"){
                epp = true;
                h_4 = (TH1F*)inFile_4->Get(VarName+"_epp"+"_hist");
                h_5 = (TH1F*)inFile_1->Get(VarName+"_epp"+"_hist");
                h_6 = (TH1F*)inFile_2->Get(VarName+"_epp"+"_hist");
                h_7 = (TH1F*)inFile_3->Get(VarName+"_epp"+"_hist");
        }
        if(epp){

                //Get Acceptance Hists with detected recoil

                h_1 = (TH1F*)inFile_1->Get(VarName+"_acc_epp"+"_hist");
                h_2 = (TH1F*)inFile_2->Get(VarName+"_acc_epp"+"_hist");
                h_3 = (TH1F*)inFile_3->Get(VarName+"_acc_epp"+"_hist");

                int rebin_num = 2;
                if(VarName == "pLead") rebin_num = 4;
                if(VarName == "pn") rebin_num = 4;

                h_1->Rebin(rebin_num);
                h_2->Rebin(rebin_num);
                h_3->Rebin(rebin_num);

                

        }

        


        if(VarName == "pcom" || VarName == "prel" || VarName == "pcom_x" || VarName == "pcom_y" || VarName == "pcom_z"){
                h_1->Rebin();
                h_2->Rebin();
                h_3->Rebin();
                h_5->Rebin();
                h_6->Rebin();
                h_7->Rebin();
                if(VarName == "pcom_x" || VarName == "pcom_y" || VarName == "pcom_z")
                        h_4->Rebin();
        }





        if(epp == false){

                double scale_fac = .5;
                if(VarName == "pLead"){
                        h_1->Rebin();
                        h_2->Rebin();
                        h_3->Rebin();
                }
                h_1->Rebin();
                h_2->Rebin();
                h_3->Rebin();
                h_5->Rebin();
                h_6->Rebin();
                h_7->Rebin();

                h_5->Scale(scale_fac*norm_factors[0]); // pwia e'p
                h_6->Scale(scale_fac*norm_factors[1]); 
                h_7->Scale(scale_fac*norm_factors[2]);

                if(VarName == "pLead") scale_fac/=2;
                                
                h_1->Scale(scale_fac*norm_factors[0]); // pwia acceptance 
                h_2->Scale(scale_fac*norm_factors[1]); 
                h_3->Scale(scale_fac*norm_factors[2]);


        } else {


                h_5->Rebin();
                h_6->Rebin();
                h_7->Rebin();

                double scale_fac = .5;

                if(VarName == "pLead" || VarName == "pn"){ 
                        scale_fac = .25;
                        h_5->Rebin();
                        h_6->Rebin();
                        h_7->Rebin();
                }

                if(VarName == "prel") scale_fac = .5;
                
                h_1->Scale(scale_fac*norm_factors[0]);
                h_2->Scale(scale_fac*norm_factors[1]);
                h_3->Scale(scale_fac*norm_factors[2]);

                if(VarName != "prel" && VarName != "pcom" && VarName != "pRec"){
                        
                        h_5->Scale(scale_fac*norm_factors[0]);      
                        h_6->Scale(scale_fac*norm_factors[1]);      
                        h_7->Scale(scale_fac*norm_factors[2]);      

                        
                } else{

                        h_5->Scale(scale_fac*norm_factors[0]);      
                        h_6->Scale(scale_fac*norm_factors[1]);      
                        h_7->Scale(scale_fac*norm_factors[2]);

                }

                
                
        
        }

        


        double hist_1_ymax = 0;
        double hist_2_ymax = 0;
        double hist_3_ymax = 0;
        double hist_4_ymax = 0;
        double hist_5_ymax = 0;
        double hist_6_ymax = 0;
        double hist_7_ymax = 0;



        if(h_1 != NULL){ 
                hist_1_ymax = h_1->GetMaximum();  
        }
        if(h_2 != NULL){ 
                hist_2_ymax = h_2->GetMaximum();   
        }
        if(h_3 != NULL){ 
                hist_3_ymax = h_3->GetMaximum();   
        }
        if(h_4 != NULL){ 
                hist_4_ymax = h_4->GetMaximum();   
        }
        if(h_5 != NULL){ 
                hist_5_ymax = h_5->GetMaximum();   
        }
        if(h_6 != NULL){ 
                hist_6_ymax = h_6->GetMaximum();   
        }
        if(h_7 != NULL){ 
                hist_7_ymax = h_7->GetMaximum();   
        }       

        int color = kAzure+2;

        h_1->SetLineStyle(kDotted);
        h_5->SetLineStyle(kDotted);

        h_2->SetLineStyle(kDashed);
        h_6->SetLineStyle(kDashed);

        h_4->SetMarkerColor(kBlack);
        h_4->SetMarkerStyle(20);
        h_4->SetMarkerSize(1.2);
        h_4->SetLineWidth(2);
        h_4->SetLineColor(kBlack);

        h_1->SetLineColor(color);
        h_2->SetLineColor(color);
        h_3->SetLineColor(color);
        h_5->SetLineColor(color);
        h_6->SetLineColor(color);
        h_7->SetLineColor(color);

        h_1->SetLineWidth(3);
        h_2->SetLineWidth(3);
        h_3->SetLineWidth(3);
        h_5->SetLineWidth(3);
        h_6->SetLineWidth(3);
        h_7->SetLineWidth(3);

        TLegend *leg = new TLegend(0.65,0.65,0.95,0.87);
        TLegend *leg_plead = new TLegend(0.28,0.67,0.63,0.88);
        if(VarName == "pLead")
                leg = leg_plead;

        leg->SetFillStyle(0);
        leg->AddEntry(h_1,"PWIA","l");
        leg->AddEntry(h_2,"T+SCX","l");
        leg->AddEntry(h_3,"FULL","l");
        leg->AddEntry(h_4,"Data","pe");
    

        std::vector<double> hist_ymax_vec;
        hist_ymax_vec.clear();

        hist_ymax_vec.push_back(hist_1_ymax);
        hist_ymax_vec.push_back(hist_2_ymax);
        hist_ymax_vec.push_back(hist_3_ymax);

        double h_ymax = 0;
        for(int j = 0; j < hist_ymax_vec.size(); j++){
                if(hist_ymax_vec[j] > h_ymax){
                        h_ymax = hist_ymax_vec[j];
                }
        }


        std::vector<double> hist_ymax_vec_2;
        hist_ymax_vec_2.clear();
        hist_ymax_vec_2.push_back(hist_5_ymax);
        hist_ymax_vec_2.push_back(hist_6_ymax);
        hist_ymax_vec_2.push_back(hist_7_ymax);
 
        double h_ymax_2 = 0;
        for(int j = 0; j < hist_ymax_vec_2.size(); j++){
                if(hist_ymax_vec_2[j] > h_ymax_2){
                        h_ymax_2 = hist_ymax_vec_2[j];
                }
        }


        if(hist_1_ymax == h_ymax){
                if(hist_5_ymax == h_ymax_2) OverlayPlots(h_1,h_2,h_3,h_4,h_5,h_6,h_7,VarName,leg,cutName,1);       
                if(hist_6_ymax == h_ymax_2) OverlayPlots(h_1,h_2,h_3,h_4,h_6,h_5,h_7,VarName,leg,cutName,1);       
                if(hist_7_ymax == h_ymax_2) OverlayPlots(h_1,h_2,h_3,h_4,h_7,h_5,h_6,VarName,leg,cutName,1);       
        }
        if(hist_2_ymax == h_ymax){
                if(hist_5_ymax == h_ymax_2) OverlayPlots(h_2,h_1,h_3,h_4,h_5,h_6,h_7,VarName,leg,cutName,2);
                if(hist_6_ymax == h_ymax_2) OverlayPlots(h_2,h_1,h_3,h_4,h_6,h_5,h_7,VarName,leg,cutName,2);
                if(hist_7_ymax == h_ymax_2) OverlayPlots(h_2,h_1,h_3,h_4,h_7,h_5,h_6,VarName,leg,cutName,2);
        } 
        if(hist_3_ymax == h_ymax){
                if(hist_5_ymax == h_ymax_2) OverlayPlots(h_3,h_1,h_2,h_4,h_5,h_6,h_7,VarName,leg,cutName,3);
                if(hist_6_ymax == h_ymax_2) OverlayPlots(h_3,h_1,h_2,h_4,h_6,h_5,h_7,VarName,leg,cutName,3);
                if(hist_7_ymax == h_ymax_2) OverlayPlots(h_3,h_1,h_2,h_4,h_7,h_5,h_6,VarName,leg,cutName,3);

        } 
    

        if(VarName == "prel"){
                TLegend *leg_pcom = new TLegend(0.05,0.45,0.4,0.75);
                leg_pcom->SetFillStyle(0);

                vector<TH1F*> data_hists;
                data_hists.push_back((TH1F*)inFile_4->Get("pcom_x_epp_hist"));
                data_hists.push_back((TH1F*)inFile_4->Get("pcom_y_epp_hist"));
                data_hists.push_back((TH1F*)inFile_4->Get("pcom_z_epp_hist"));
                data_hists.push_back((TH1F*)inFile_4->Get("pcom_epp_hist"));
                TString cut;
                vector<TFile*> inFiles;
                inFiles.push_back(inFile_1);
                inFiles.push_back(inFile_2);
                inFiles.push_back(inFile_3);
                vector<TH1F*> pcom_x,pcom_y,pcom_z;
                vector<TH1F*> pcom_x_epp,pcom_y_epp,pcom_z_epp;
                vector<TH1F*> mag_hists;

                

                for(int j = 0; j < data_hists.size(); ++j){
                        data_hists[j]->SetMarkerColor(kBlack);
                        data_hists[j]->SetMarkerStyle(20);
                        data_hists[j]->SetMarkerSize(3);
                        data_hists[j]->SetLineWidth(2);
                        data_hists[j]->SetLineColor(kBlack);
                        if(j != 3)
                                data_hists[j]->Rebin();
                }

                h_1->SetLineWidth(2);
                h_2->SetLineWidth(2);
                h_3->SetLineWidth(2);
                leg_pcom->SetNColumns(1);
                leg_pcom->AddEntry(h_1,"PWIA","l");
                leg_pcom->AddEntry(h_2,"T+SCX","l");
                leg_pcom->AddEntry(h_3,"FULL","l");
                        
                for(int order = 1; order < 2; ++order){
                        cut = "_epp";
                        if(order == 1) leg_pcom->AddEntry(data_hists[0],"Data","pe");

                        
                        if(pcom_x.size() != 0) pcom_x.clear();
                        if(pcom_y.size() != 0) pcom_y.clear();
                        if(pcom_z.size() != 0) pcom_z.clear();
                            

                        

                        for(int j = 0; j < inFiles.size(); ++j){
                                pcom_x.push_back((TH1F*)inFiles[j]->Get("pcom_x_acc_epp_hist"));
                                pcom_y.push_back((TH1F*)inFiles[j]->Get("pcom_y_acc_epp_hist"));
                                pcom_z.push_back((TH1F*)inFiles[j]->Get("pcom_z_acc_epp_hist"));
                                mag_hists.push_back((TH1F*)inFiles[j]->Get("pcom_acc_epp_hist"));
                        }

                        for(int j = 0; j < inFiles.size(); ++j){
                                pcom_x_epp.push_back((TH1F*)inFiles[j]->Get("pcom_x"+cut+"_hist"));
                                pcom_y_epp.push_back((TH1F*)inFiles[j]->Get("pcom_y"+cut+"_hist"));
                                pcom_z_epp.push_back((TH1F*)inFiles[j]->Get("pcom_z"+cut+"_hist"));
                                mag_hists.push_back((TH1F*)inFiles[j]->Get("pcom"+cut+"_hist"));

                        }

                       
        

                        for(int j = 0; j < inFiles.size(); ++j){
                            pcom_x[j]->Scale(norm_factors[j+6]);
                            pcom_y[j]->Scale(norm_factors[j+6]);
                            pcom_z[j]->Scale(norm_factors[j+6]);

                            pcom_x_epp[j]->Scale(norm_factors[j+6]);
                            pcom_y_epp[j]->Scale(norm_factors[j+6]);
                            pcom_z_epp[j]->Scale(norm_factors[j+6]);

                            mag_hists[j]->Scale(norm_factors[j+6]);
                            mag_hists[j+3]->Scale(norm_factors[j+6]);

                        }

                        int rebin_num = 4;

                        for(int j = 0; j < inFiles.size(); ++j){

                                
                                pcom_x[j]->Rebin(rebin_num);
                                pcom_y[j]->Rebin(rebin_num);
                                pcom_z[j]->Rebin(rebin_num);
                           
                                mag_hists[j]->Rebin(rebin_num/2);
                                mag_hists[j+3]->Rebin(rebin_num);

                                pcom_x[j]->SetLineColor(kAzure+2);
                                pcom_y[j]->SetLineColor(kAzure+2);
                                pcom_z[j]->SetLineColor(kAzure+2);

                                pcom_x[j]->SetLineWidth(2);
                                pcom_y[j]->SetLineWidth(2);
                                pcom_z[j]->SetLineWidth(2);

                                pcom_x_epp[j]->Rebin(rebin_num*2);
                                pcom_y_epp[j]->Rebin(rebin_num*2);
                                pcom_z_epp[j]->Rebin(rebin_num);
             

                                pcom_x_epp[j]->SetLineColor(kAzure+2);
                                pcom_y_epp[j]->SetLineColor(kAzure+2);
                                pcom_z_epp[j]->SetLineColor(kAzure+2);

                                pcom_x_epp[j]->SetLineWidth(2);
                                pcom_y_epp[j]->SetLineWidth(2);
                                pcom_z_epp[j]->SetLineWidth(2);

                                mag_hists[j]->SetLineWidth(2);
                                mag_hists[j+3]->SetLineWidth(2);

                                mag_hists[j]->SetLineColor(kAzure+2);
                                mag_hists[j+3]->SetLineColor(kAzure+2);

                                if(j == 0){
                                        pcom_x[j]->SetLineStyle(kDotted);
                                        pcom_y[j]->SetLineStyle(kDotted);
                                        pcom_z[j]->SetLineStyle(kDotted);

                                        pcom_x_epp[j]->SetLineStyle(kDotted);
                                        pcom_y_epp[j]->SetLineStyle(kDotted);
                                        pcom_z_epp[j]->SetLineStyle(kDotted);

                                        mag_hists[j]->SetLineStyle(kDotted);
                                        mag_hists[j+3]->SetLineStyle(kDotted);
                                }
                                if(j == 1){
                                        pcom_x[j]->SetLineStyle(kDashed);
                                        pcom_y[j]->SetLineStyle(kDashed);
                                        pcom_z[j]->SetLineStyle(kDashed);

                                        pcom_x_epp[j]->SetLineStyle(kDashed);
                                        pcom_y_epp[j]->SetLineStyle(kDashed);
                                        pcom_z_epp[j]->SetLineStyle(kDashed);

                                        mag_hists[j]->SetLineStyle(kDashed);
                                        mag_hists[j+3]->SetLineStyle(kDashed);

                                }
                        }


        

                         for(int j = 0; j < inFiles.size(); ++j){
                                double ratio = 1;
                                ratio = .25;
                                ratio *= 2;

                                pcom_x[j]->Scale(ratio);
                                pcom_y[j]->Scale(ratio);
                                pcom_z[j]->Scale(ratio);

                                ratio = .25;


                                pcom_x_epp[j]->Scale(ratio);
                                pcom_y_epp[j]->Scale(ratio);
                                ratio*=2;
                                
                                pcom_z_epp[j]->Scale(ratio);

                                ratio *= 2;

                                mag_hists[j]->Scale(ratio);
                                ratio/=2;
                                mag_hists[j+3]->Scale(ratio);

                        }


                        cout << "scale 6: " << norm_factors[0] << endl; 
                        cout << "scale 7: " << norm_factors[1] << endl;
                        cout << "scale 8: " << norm_factors[2] << endl;


                        for(int j = 0; j < pcom_x_epp.size(); ++j){
                                pcom_x.push_back(pcom_x_epp[j]);
                                pcom_y.push_back(pcom_y_epp[j]);
                                pcom_z.push_back(pcom_z_epp[j]);
                        }

                        
                        cout << "mag size: " <<  mag_hists.size() << endl;
                        
                        PcmPlots(pcom_x,pcom_y,pcom_z,data_hists,mag_hists,VarName,leg_pcom,cut,order);                        
                        
                }
        }
        
}

void gst_plotter_test(){


 
        TFile *f_input_1= TFile::Open("/pnfs/genie/scratch/users/nwright/devel_gaussian_q2/hists_emqe_gcf_G18_10a_eppstudy_501_1bil_weight.root","readonly");
        //TFile *f_input_2 = TFile::Open("/pnfs/genie/scratch/users/nwright/devel_gaussian_q2/hists_emqe_gcf_G18_10a_eppstudy_plead_500mil_cx_3.root","readonly");
        TFile *f_input_2 = TFile::Open("/pnfs/genie/scratch/users/nwright/devel_gaussian_q2/hists_emqe_gcf_G18_10a_eppstudy_501_500mil_cx.root","readonly");
        TFile *f_input_3 = TFile::Open("/pnfs/genie/scratch/users/nwright/devel_gaussian_q2/hists_emqe_gcf_G18_10a_eppstudy_501_1bil_reco_f.root","readonly");
        TFile *f_input_4 = TFile::Open("Macros/ep_data.root","readonly"); 

        vector<TString> cut;
  //      cut.push_back("");
        cut.push_back("_acceptance");
   /*     cut.push_back("_bkx");
        cut.push_back("_thpq");
        cut.push_back("_pmq");
        cut.push_back("_ep");
        cut.push_back("_epp");
*/


        std::vector<TString> varNames;
        varNames.clear();

    

        varNames.push_back("cth_pcom");
        varNames.push_back("cth_prel");
        varNames.push_back("th_pmiq");
        varNames.push_back("theta_pq");
      //  varNames.push_back("pcom");
        varNames.push_back("pLead");
        varNames.push_back("pn");
        varNames.push_back("prel");



        int varNames_size = varNames.size();
        int cut_size = cut.size();
       // nfiles = 1;
        vector<TFile*> inFiles;
        inFiles.push_back(f_input_1);
        inFiles.push_back(f_input_2);
        inFiles.push_back(f_input_3);
        vector<double> norm_factors;


        TH1F *h;
        TH1F *data;
        double transparency;
        vector<double> cuts;

        for(int j = 0; j < inFiles.size(); ++j){
                transparency = 1.;
                if(j == 1) transparency = (.44/.53);
                h = (TH1F*)inFiles[j]->Get("pn_ep_hist");
                data = (TH1F*)f_input_4->Get("pn_ep_hist");
                norm_factors.push_back(transparency*data->Integral()/h->Integral());
                // 0 = PWIA ep
                // 1 = T+SCX ep
                // 2 = FULL ep
        
        }
 
        
        for(int j = 0; j< cut_size; j++){
                
                for(int i = 0; i < varNames_size; i++){
                        if(varNames[i] == "spectral_func" && cut[j] != ""){continue;}
                        cout << varNames[i] << endl;
                        Hist_draw_order(f_input_1, f_input_2, f_input_3,  f_input_4,/* f_input_5,*/norm_factors, varNames[i], cut[j]);    

                }
                
        }


}