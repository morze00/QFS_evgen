#include "info/libs.hh"
#include "info/info.hh"

#define NP 61
#define nfiles 8

Double_t a[nfiles+1][NP];
Double_t b[nfiles+1][NP];
Double_t c[nfiles+1][NP];

Double_t slp[nfiles+1][NP]={0};
Double_t offs[nfiles+1][NP]={0};

Int_t nop[nfiles+1];
Double_t energy[nfiles+1];

Double_t funct(Double_t *x, Double_t *par){
	
	Double_t dxs = -1;
	Double_t dxssin = -1;
    Int_t pp = Int_t(par[0]);
    Int_t nn = nop[pp];

    //cout << par[0] << "\t" << pp << endl;
	
	if(x[0]>a[pp][nn])
		dxs = (c[pp][nn]+c[pp][nn-1])/2.;
	else if(x[0]<a[pp][0])
		dxs = c[pp][0];
	
	for(Int_t i=0; i<nn; i++){
		if(x[0]>a[pp][i] && x[0]<a[pp][i+1]){
			dxs = slp[pp][i]*x[0]+offs[pp][i];
		break;
		}
	}
	
	dxssin = (dxs*sin(x[0]/180.*PI));
	return dxssin;
}


Double_t fsin(Double_t *x, Double_t *par){

	return (sin(x[0]/180.*PI));
}

	TF1** fp;
	TF1* fsin0;

//Main function
void show_parametrization()
{

    Int_t inp_en[nfiles] = {300,325,350,375,385,400,425,450};
	FILE *parfile;

    fsin0 = new TF1("sinus",fsin,0.,180.);
    fsin0->SetLineColor(kGray);

    fsin0->GetYaxis()->SetRangeUser(0.,12.);
	fsin0->Draw();
	   
	fp = new TF1*[nfiles+1];
    TLegend * L = new TLegend(0.8,0.7,0.9,0.9);

    char line[64];
    Double_t buffer_a[128];
    Double_t buffer_b[128];
    Double_t buffer_c[128];
    Double_t buffer_e;
    Int_t tt;

    for(Int_t cc=0; cc<nfiles; cc++)
    {
   
        parfile = fopen(Form("data_filtered_calc_%dMeV.txt",inp_en[cc]),"r");
        tt=0;
     	while(fgets(line, 64, parfile) != NULL)
    	{
    	    sscanf (line, "%lf\t%lf\t%lf\t%lf", &buffer_c[tt], &buffer_b[tt], &buffer_e, &buffer_a[tt] );
            a[cc][tt] = buffer_a[tt];
            b[cc][tt] = buffer_b[tt];
            c[cc][tt] = buffer_c[tt];
    	    tt++;
    	}  
        energy[cc] = buffer_e;
        nop[cc] = tt;
    
        fclose(parfile);
    
        for(Int_t i=0; i<nop[cc]; i++){
    		cout << a[cc][i] << " " << b[cc][i] << " " << c[cc][i] << endl;
    		slp[cc][i]=(c[cc][i+1]-c[cc][i])/(a[cc][i+1]-a[cc][i]);
    		offs[cc][i]=c[cc][i]-slp[cc][i]*a[cc][i];
    	}
    
        fp[cc] = new TF1(Form("f%.0f",energy[cc]),funct,0.,180.,2);
        fp[cc]->SetParameters(Double_t(cc),0);
        fp[cc]->SetLineColor(cc+1);
        fp[cc]->SetLineWidth(3);
        //fp[cc]->SetLineStyle(cc+1);
    	//Cross section function
        fp[cc]->Draw("same");
        L->AddEntry(fp[cc], Form("%.0f MeV",energy[cc]), "f");
    }
  
    parfile = fopen("data_filtered.txt","r");
    tt=0;
 	while(fgets(line, 64, parfile) != NULL)
	{
	    sscanf (line, "%lf\t%lf\t%lf\t%lf", &buffer_c[tt], &buffer_b[tt], &buffer_e, &buffer_a[tt] );
        a[nfiles][tt] = buffer_a[tt];
        b[nfiles][tt] = buffer_b[tt];
        c[nfiles][tt] = buffer_c[tt];
	    tt++;
	}  
    energy[nfiles] = buffer_e;
    nop[nfiles] = tt;

    fclose(parfile);

    for(Int_t i=0; i<nop[nfiles]; i++){
		cout << a[nfiles][i] << " " << b[nfiles][i] << " " << c[nfiles][i] << endl;
		slp[nfiles][i]=(c[nfiles][i+1]-c[nfiles][i])/(a[nfiles][i+1]-a[nfiles][i]);
		offs[nfiles][i]=c[nfiles][i]-slp[nfiles][i]*a[nfiles][i];
	}

    fp[nfiles] = new TF1("fexp",funct,0.,180.,2);
    fp[nfiles]->SetParameters(nfiles,0);
    fp[nfiles]->SetLineColor(nfiles+1);
    fp[nfiles]->SetLineWidth(3);
    fp[nfiles]->SetLineStyle(nfiles+1);
	//Cross section function
    fp[nfiles]->Draw("same");
    L->AddEntry(fp[nfiles], "Exp. p,n", "f");
  
    L->AddEntry(fsin0, "sin(x)", "f");
    L->SetFillColor(kWhite);
    L->Draw();

}

