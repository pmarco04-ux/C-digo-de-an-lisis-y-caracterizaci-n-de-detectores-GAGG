// Empezamos metiendo las librerias necesarias 
#include <string>
#include <utility> 

#include <fstream>
#include <sstream>

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <filesystem> 

// ROOT libraries
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TSpectrum.h>
#include <TVirtualFitter.h>
#include <TSystem.h>
#include <TLine.h>
#include "TMatrixD.h"


#include "General.cpp"

// Si por ejemplo este código llamara a otros códigos/scripts (es común y lo harás en el futuro)
// Se escribe: 
// #include "ruta al codigo" (por ejemplo "src/Codigo2.cpp")

/////////////////////////////////////////
/////// Ir al final del Código /////////
////////////////////////////////////////

void Sort_GAGGs(TChain* chain, TChain* chainBKG, TFile* outputFile, 
    const std::string filename, const std::string isotope, const std::string Calib_filename ,const bool Substract_BKG, const std::string BKG_Calib_filename
    ){
    // En nuestro caso solo tenemos un MDPP16 conectado. Si te fijas en el fichero que te habrá dado Dani
    // Si hace root -l Nombre.root y new TBrowser veras las variables

    //Numero de modulos MDPP16 por si metemos alguno mas
    const int mdpp16_counter = 1;

    //Numero de canales
    const Int_t NumberOfChan_mdpp16 = 16; 

    //Maxima multiplicidad de los eventos (cuantos canales golpean) 
    const Int_t maxMul = 64;

    TRandom3 RandomR3;

    //Defino las variables del tree que nos interesan
    Int_t mul_MDPP16_Long[mdpp16_counter],mul_MDPP16_Short[mdpp16_counter];
    int detCh_MDPP16_Long[mdpp16_counter][maxMul],detCh_MDPP16_Short[mdpp16_counter][maxMul];
    Long64_t values_MDPP16_Long[mdpp16_counter][maxMul],values_MDPP16_Short[mdpp16_counter][maxMul],EventTimeStamp_MDPP16[mdpp16_counter],ExtTimeStamp_MDPP16[mdpp16_counter];

    std::cout << "Reading TTree for MDPP16" << std::endl;
    for (int i = 0; i < mdpp16_counter; i++) {
        // Hay que coger el nombre de cada una de las ramas 
        TString baseName = TString::Format("MDPP16_QDC_%d_", i);
        
        // Ponemos las ramas que aparecen en el TBrowser con el mismo nombre en el EventTree. Estos nombres deberían ser lo bueno, pero puede que Dani haya tocado algo de su código
        // que traduce de Mesytec a ROOT así que por si acaso comprueba que los nombres sean estos

        chain->SetBranchAddress(baseName + "ilong_Multiplicity", &mul_MDPP16_Long[i]); // Multiplicidad Integral Larga. A la variable con nombre baseName + "ilong_Multiplicity" le asigno una variable en mi código  &mul_MDPP16_Long[i]. Por el & no te preocupes son cosas de punteros.
        chain->SetBranchAddress(baseName + "ilong_Channels", detCh_MDPP16_Long[i]); // Canal (GAGG) que ha golpeado
        chain->SetBranchAddress(baseName + "ilong_Values", values_MDPP16_Long[i]); // Canal electrónico donde ha golpead
        chain->SetBranchAddress(baseName + "ishort_Multiplicity", &mul_MDPP16_Short[i]); //Lo mismo pero para el short
        chain->SetBranchAddress(baseName + "ishort_Channels", detCh_MDPP16_Short[i]); 
        chain->SetBranchAddress(baseName + "ishort_Values", values_MDPP16_Short[i]);
        chain->SetBranchAddress(baseName + "Event_timestamps", &EventTimeStamp_MDPP16[i]); //Tiempo 
        chain->SetBranchAddress(baseName + "Extended_timestamps", &ExtTimeStamp_MDPP16[i]);

        chainBKG->SetBranchAddress(baseName + "ilong_Multiplicity",  &mul_MDPP16_Long[i]);
        chainBKG->SetBranchAddress(baseName + "ilong_Channels",       detCh_MDPP16_Long[i]);
        chainBKG->SetBranchAddress(baseName + "ilong_Values",         values_MDPP16_Long[i]);
        chainBKG->SetBranchAddress(baseName + "ishort_Multiplicity", &mul_MDPP16_Short[i]);
        chainBKG->SetBranchAddress(baseName + "ishort_Channels",      detCh_MDPP16_Short[i]);
        chainBKG->SetBranchAddress(baseName + "ishort_Values",        values_MDPP16_Short[i]);
        chainBKG->SetBranchAddress(baseName + "Event_timestamps", &EventTimeStamp_MDPP16[i]); 
        chainBKG->SetBranchAddress(baseName + "Extended_timestamps", &ExtTimeStamp_MDPP16[i]);
    }


    //Cojo los valores de Energía
    std::vector<double> Calibration_Values = Select_Calibration_Values(isotope);
    std::vector<double> Ranges_values = Select_Fit_Ranges(isotope);
    const int nPeaks = Calibration_Values.size();

    //Defino 8 Histogramas uno por GAGG
    // Raw 
    char GAGGS_Raw_name_Total_Long[8][100]; 
    char GAGGS_Raw_name_histo_Total_Long[8][200]; 
    static TH1F* GAGGS_Raw_hist_Total_Long[8] = {nullptr}; 

    char GAGGS_Raw_name_Total_Short[8][100]; //Lo mismo pero para el short
    char GAGGS_Raw_name_histo_Total_Short[8][200]; 
    static TH1F* GAGGS_Raw_hist_Total_Short[8] = {nullptr};

    //Calib
    char GAGGS_Calib_name_Total_Long[8][100]; 
    char GAGGS_Calib_name_histo_Total_Long[8][200]; 
    static TH1F* GAGGS_Calib_hist_Total_Long[8] = {nullptr}; 

    char GAGGS_Calib_name_Total_Short[8][100]; 
    char GAGGS_Calib_name_histo_Total_Short[8][200]; 
    static TH1F* GAGGS_Calib_hist_Total_Short[8] = {nullptr};    

    //Residuo
    char GAGGS_Residuo_name_Long[8][100];//Histogramas para mostrar los residuos long
    char GAGGS_Residuo_name_histo_Long[8][200];
    static TH2F* GAGGS_Residuo_Long[8] = {nullptr};

    char GAGGS_Residuo_name_Short[8][100];//Histogramas para mostrar los residuos short
    char GAGGS_Residuo_name_histo_Short[8][200];
    static TH2F* GAGGS_Residuo_Short[8] = {nullptr};

    //FWHM
    char GAGGS_FWHM_name_Long[8][100];//Histogramas para mostrar FWHM
    char GAGGS_FWHM_name_histo_Long[8][200];
    static TH2F* GAGGS_FWHM_Long[8] = {nullptr};

    char GAGGS_FWHM_name_Short[8][100];//Histogramas para mostrar FWHM
    char GAGGS_FWHM_name_histo_Short[8][200];
    static TH2F* GAGGS_FWHM_Short[8] = {nullptr};

    //Coinci 
    char GAGGS_Coinci_Two_name[8][100];
    char GAGGS_Coinci_Two_name_histo[8][200];
    static TH1F* GAGGS_Coinci_Two_hist[8] = {nullptr};

    char GAGGS_Coinci_Two_Full_name[100];
    char GAGGS_Coinci_Two_Full_name_histo[200];
    static TH1F* GAGGS_Coinci_Two_Full_hist = {nullptr};

    //GAGGs Coincidence with 511
    char GAGGS_Coinci_with_511_name[8][100];
    char GAGGS_Coinci_with_511_name_histo[8][200];
    static TH1F* GAGGS_Coinci_with_511_hist[8] = {nullptr};

    char GAGGS_Coinci_with_511_Full_name[100];
    char GAGGS_Coinci_with_511_name_Full_histo[200];
    static TH1F* GAGGS_Coinci_with_511_Full_hist = {nullptr};

    //GAGGs Coincidence with 1 GAGG in 511 keV
    char GAGGS_Multi_2_name[100];
    char GAGGS_Multi_2_name_histo[200];
    static TH2F* GAGGS_Multi_2_hist = {nullptr};
    
    //GAGGs in Coincidence (Like P-P)
    char GAGGS_Coinci_2_name[100];
    char GAGGS_Coinci_2_name_histo[200];
    static TH2F* GAGGS_Coinci_2_hist = {nullptr};

    //Short vs Long
    char GAGGS_Short_vs_Long_name[8][100];
    char GAGGS_Short_vs_Long_name_histo[8][200];
    static TH2F* GAGGS_Short_vs_Long_hist[8] = {nullptr};
    
    //Check Calib
    char GAGGs_Check_Calib_name[100];
    char GAGGs_Check_Calib_name_histo[100];
    static TH2F* GAGGs_Check_Calib_hist = {nullptr};

    //Add-back
    char GAGGS_Add_back_name_Total_Long[8][100]; 
    char GAGGS_Add_back_name_histo_Total_Long[8][200]; 
    static TH1F* GAGGS_Add_back_hist_Total_Long[8] = {nullptr}; 

    //Numero de eventos en coincidencias
    static TH1F* GAGGS_n_Coinci_hist = {nullptr};

    //Comparación de 2 GAGGS con el tercero en coincidencias
    static TH2F* GAGGS_3hit_PairSum_vs_Third_hist = {nullptr};

    static TH2F* GAGGS_3hit_Emax_vs_2Emin = {nullptr};

    for (int i = 0; i < 8 ; i++){
        sprintf(GAGGS_Raw_name_Total_Long[i],"GAGG_%i Raw Long", i+1); 
        sprintf(GAGGS_Raw_name_histo_Total_Long[i],"GAGG_%i Raw Long; Channel ; Counts", i+1); 
        GAGGS_Raw_hist_Total_Long[i] = new TH1F(GAGGS_Raw_name_Total_Long[i],GAGGS_Raw_name_histo_Total_Long[i],65536.0/4.0,0,65536);    
    
        sprintf(GAGGS_Raw_name_Total_Short[i],"GAGG_%i Raw Short", i+1); //Lo mismo pero con short
        sprintf(GAGGS_Raw_name_histo_Total_Short[i],"GAGG_%i Raw Short; Channel ; Counts", i+1); 
        GAGGS_Raw_hist_Total_Short[i] = new TH1F(GAGGS_Raw_name_Total_Short[i],GAGGS_Raw_name_histo_Total_Short[i],65536.0/4.0,0,65536); 

        sprintf(GAGGS_Calib_name_Total_Long[i],"GAGG_%i Calib Long", i+1);
        sprintf(GAGGS_Calib_name_histo_Total_Long[i],"GAGG_%i Calib Long; Energia(KeV) ; Counts", i+1); 
        GAGGS_Calib_hist_Total_Long[i] = new TH1F(GAGGS_Calib_name_Total_Long[i],GAGGS_Calib_name_histo_Total_Long[i],5000.0/2.0,0,5000);

        sprintf(GAGGS_Calib_name_Total_Short[i],"GAGG_%i Calib Short", i+1); 
        sprintf(GAGGS_Calib_name_histo_Total_Short[i],"GAGG_%i Calib Short; Energia (KeV); Counts", i+1); 
        GAGGS_Calib_hist_Total_Short[i] = new TH1F(GAGGS_Calib_name_Total_Short[i],GAGGS_Calib_name_histo_Total_Short[i],5000.0/2.0,0,5000);

        sprintf(GAGGS_Add_back_name_Total_Long[i],"GAGG_%i Add-back", i+1);
        sprintf(GAGGS_Add_back_name_histo_Total_Long[i],"GAGG_%i Add-back; Energia(KeV) ; Counts", i+1); 
        GAGGS_Add_back_hist_Total_Long[i] = new TH1F(GAGGS_Add_back_name_Total_Long[i],GAGGS_Add_back_name_histo_Total_Long[i],5000.0/2.0,0,5000);


        for(int j = 0; j < nPeaks; j ++){
            Draw_Calib_lines(GAGGS_Calib_hist_Total_Long[i], Calibration_Values[j]);
            Draw_Calib_lines(GAGGS_Calib_hist_Total_Short[i], Calibration_Values[j]);
            Draw_Calib_lines(GAGGS_Add_back_hist_Total_Long[i], Calibration_Values[j]);
        }

        for(int j = 0; j < nPeaks; j ++){
            Draw_Compton_lines(GAGGS_Calib_hist_Total_Long[i], Calibration_Values[j]);
            Draw_Compton_lines(GAGGS_Calib_hist_Total_Short[i], Calibration_Values[j]);
            Draw_Compton_lines(GAGGS_Add_back_hist_Total_Long[i], Calibration_Values[j]);
        }

        for(int j = 0; j < nPeaks; j ++){
            Draw_BackScatter_lines(GAGGS_Calib_hist_Total_Long[i], Calibration_Values[j]);
            Draw_BackScatter_lines(GAGGS_Calib_hist_Total_Short[i], Calibration_Values[j]);
            Draw_BackScatter_lines(GAGGS_Add_back_hist_Total_Long[i], Calibration_Values[j]);
        }

        sprintf(GAGGS_Residuo_name_Long[i],"Residuo_GAGG_%i Long", i+1); 
        sprintf(GAGGS_Residuo_name_histo_Long[i],"Residuo_GAGG_%i Long; Energia(KeV);Residuo (KeV) ; Counts", i+1); 
        GAGGS_Residuo_Long[i] = new TH2F(GAGGS_Residuo_name_Long[i],GAGGS_Residuo_name_histo_Long[i],2000,0,2000,300,-10,10);
        add_zero_line(GAGGS_Residuo_Long[i]); 

        sprintf(GAGGS_Residuo_name_Short[i],"Residuo_GAGG_%i Short", i+1); 
        sprintf(GAGGS_Residuo_name_histo_Short[i],"Residuo_GAGG_%i Short; Energia(KeV);Residuo (KeV) ; Counts", i+1); 
        GAGGS_Residuo_Short[i] = new TH2F(GAGGS_Residuo_name_Short[i],GAGGS_Residuo_name_histo_Short[i],2000,0,2000,300,-10,10);
        add_zero_line(GAGGS_Residuo_Short[i]);

        sprintf(GAGGS_FWHM_name_Long[i],"FWHM_GAGG_%i Long", i+1); 
        sprintf(GAGGS_FWHM_name_histo_Long[i],"FWHM_GAGG_%i Long; Energia(KeV);FWHM (KeV) ; Counts", i+1); 
        GAGGS_FWHM_Long[i] = new TH2F(GAGGS_FWHM_name_Long[i],GAGGS_FWHM_name_histo_Long[i],2000,0,2000,300,0,100); 
        add_zero_line(GAGGS_FWHM_Long[i]);

        sprintf(GAGGS_FWHM_name_Short[i],"FWHM_GAGG_%i Short", i+1); 
        sprintf(GAGGS_FWHM_name_histo_Short[i],"FWHM_GAGG_%i Short; Energia(KeV);FWHM (KeV) ; Counts", i+1); 
        GAGGS_FWHM_Short[i] = new TH2F(GAGGS_FWHM_name_Short[i],GAGGS_FWHM_name_histo_Short[i],2000,0,2000,300,0,100);
        add_zero_line(GAGGS_FWHM_Short[i]);

        sprintf(GAGGS_Short_vs_Long_name[i],"GAGG_%i Short vs Long", i+1);
        sprintf(GAGGS_Short_vs_Long_name_histo[i],"GAGG_%i Short vs Long; E Long (keV) ; (E Short-E_Long)/E_Long ; Counts", i+1);
        GAGGS_Short_vs_Long_hist[i] = new TH2F(GAGGS_Short_vs_Long_name[i],GAGGS_Short_vs_Long_name_histo[i],500,0,5000,1000,-1,1);

        sprintf(GAGGS_Coinci_Two_name[i],"Two GAGGs in Coincidence GAGG_%i",i+1);
        sprintf(GAGGS_Coinci_Two_name_histo[i],"Two GAGGs in Coincidence GAGG_%i; Energy(keV) ; Counts",i+1);
        GAGGS_Coinci_Two_hist[i] = new TH1F(GAGGS_Coinci_Two_name[i],GAGGS_Coinci_Two_name_histo[i],3500,0,7000);

        sprintf(GAGGS_Coinci_with_511_name[i],"Two GAGGs in Coincidence one is 511 GAGG_%i",i+1);
        sprintf(GAGGS_Coinci_with_511_name_histo[i],"Two GAGGs in Coincidence  one is 511 GAGG_%i; Energy(keV) ; Counts",i+1);
        GAGGS_Coinci_with_511_hist[i] = new TH1F(GAGGS_Coinci_with_511_name[i],GAGGS_Coinci_with_511_name_histo[i],3500,0,7000);
    }

    //Coincidences
    sprintf(GAGGS_Coinci_Two_Full_name,"Two GAGGs in Coincidence");
    sprintf(GAGGS_Coinci_Two_Full_name_histo,"Two GAGGs in Coincidence; Energy(keV) ; Counts");
    GAGGS_Coinci_Two_Full_hist = new TH1F(GAGGS_Coinci_Two_Full_name,GAGGS_Coinci_Two_Full_name_histo,3500,0,7000);

    sprintf(GAGGS_Coinci_with_511_Full_name,"Two GAGGs in Coincidence one is 511");
    sprintf(GAGGS_Coinci_with_511_name_Full_histo,"Two GAGGs in Coincidence  one is 511; Energy(keV) ; Counts");
    GAGGS_Coinci_with_511_Full_hist = new TH1F(GAGGS_Coinci_with_511_Full_name,GAGGS_Coinci_with_511_name_Full_histo,3500,0,7000);

    //Single Double Escape
    sprintf(GAGGS_Multi_2_name,"2 GAGGs in coincidence Energy");
    sprintf(GAGGS_Multi_2_name_histo,"2 GAGGs in coincidence Energy; E GAGG_Max (keV) ; E GAGG_Min (keV); Counts");
    GAGGS_Multi_2_hist = new TH2F(GAGGS_Multi_2_name,GAGGS_Multi_2_name_histo,350,0,7000,350,0,7000);

    sprintf(GAGGS_Coinci_2_name,"2 GAGGs in coincidence Number");
    sprintf(GAGGS_Coinci_2_name_histo,"2 GAGGs in coincidence Number; GAGG_Number ;  GAGG_Number; Counts");
    GAGGS_Coinci_2_hist = new TH2F(GAGGS_Coinci_2_name,GAGGS_Coinci_2_name_histo, 8, 1, 9, 8, 1, 9);

    //Check Calib
    sprintf(GAGGs_Check_Calib_name,"Check Calib GAGG");
    sprintf(GAGGs_Check_Calib_name_histo,"Check Calib GAGG; Energy(keV) ; Counts");
    GAGGs_Check_Calib_hist = new TH2F(GAGGs_Check_Calib_name, GAGGs_Check_Calib_name_histo, 8, 1, 9, 2500, 0, 5000);

    // Total Histo
    TH1F* GAGGS_Calib_Total_Long = new TH1F("GAGG_Total_Calib_Long", "GAGG Total Calibrated Long; Energia (KeV); Counts", 5000,0,5000);
    TH1F* GAGGS_Calib_Total_Short = new TH1F("GAGG_Total_Calib_Short", "GAGG Total Calibrated Short; Energia (KeV); Counts", 5000,0,5000);

    GAGGS_n_Coinci_hist= new TH1F ("N_GAGGs_in_Coincidence", "N GAGGs in Coincidence; N GAGGs; Counts",8,1,9);

    GAGGS_3hit_PairSum_vs_Third_hist = new TH2F("3hit_PairSum_vs_Third","3 GAGGs: Pair Sum vs Third; E_{pair} (keV); E_{third} (keV); Counts",500, 0, 5000, 250, 0, 2500);

    GAGGS_3hit_Emax_vs_2Emin= new TH2F("3hit_Emax_vs_2Emin ","3 GAGGs: Emax vs 2Emin; E_{max} (keV); 2E_{min} (keV); Counts",500, 0, 5000, 250, 0, 2500);

    for(int j = 0; j < nPeaks; j ++){
        Draw_Check_Calib_lines(GAGGs_Check_Calib_hist, Calibration_Values[j]);
        Draw_Calib_lines(GAGGS_Calib_Total_Long, Calibration_Values[j]);
        Draw_Calib_lines(GAGGS_Calib_Total_Short, Calibration_Values[j]);
    }

    double a_MDPP16_Long[8], b_MDPP16_Long[8];
    double a_MDPP16_Short[8], b_MDPP16_Short[8];
    LoadMDPP16CalibParameters(Calib_filename,a_MDPP16_Long,b_MDPP16_Long,a_MDPP16_Short,b_MDPP16_Short);

    //Cogemos el número de entradas para hacernos una idea
    Long64_t nEntries = chain->GetEntries();
    cout << "Number of Entries: " << nEntries << endl;
    int printFrequency = nEntries / 100;


    int Loops_Time_Stamp = 0; //Por si no tengo Extended Timestamp
    int Last_Timestamp = 0;

    int n1=0;
    int n2=0;
    int n3=0;
    //Bucle sobre las entradas
    for (int i = 0; i < int(nEntries); i++){

        int hitCounter_GAGG_Long[8] ={0}; int hitCounter_GAGG_Short[8] = {0}; // Teóricamente esto debería coincidir siempre
        double Energy_GAGG_Long[8] = {0}; double Energy_GAGG_Short[8] = {0}; // En teoría los GAGGs solo pueden coger una señal por evento, en caso contrario deberíamos cambiar este double por un vector
        
        std::vector<int> GAGGs_hit; //Los GAGGs que han recibido una señal
        std::vector<int> GAGGs_hit_511; //Los GAGGs que han recibido una señal cerca de 511

        chain -> GetEvent(i);

        if(Last_Timestamp > EventTimeStamp_MDPP16[0]){
            Loops_Time_Stamp++;
        }
        Last_Timestamp = EventTimeStamp_MDPP16[0];

        for (int j = 0; j < mdpp16_counter; j++){

            for (Int_t k = 0; k < mul_MDPP16_Long[j]; k++){ 
                int GAGG_Number = detCh_MDPP16_Long[j][k];
                if (GAGG_Number<=7){
                    GAGGS_Raw_hist_Total_Long[GAGG_Number]->Fill(values_MDPP16_Long[j][k]);

                    double CalibratedEnergy_MDPP16_Long = (a_MDPP16_Long[GAGG_Number]*(values_MDPP16_Long[j][k]+RandomR3.Uniform()-0.5)+b_MDPP16_Long[GAGG_Number]);

                    if(CalibratedEnergy_MDPP16_Long > 100.0){ //Ponemos un Noise Cut 
                        GAGGs_hit.push_back(GAGG_Number);
                        hitCounter_GAGG_Long[GAGG_Number]++;
                        Energy_GAGG_Long[GAGG_Number] = CalibratedEnergy_MDPP16_Long;
                        GAGGS_Calib_hist_Total_Long[GAGG_Number]->Fill(CalibratedEnergy_MDPP16_Long);
                        GAGGS_Calib_Total_Long->Fill(CalibratedEnergy_MDPP16_Long);

                        GAGGs_Check_Calib_hist->Fill((double)GAGG_Number, CalibratedEnergy_MDPP16_Long);
                    }   
                }
            }
            for (Int_t k = 0; k < mul_MDPP16_Short[j]; k++){ 
                int GAGG_Number = detCh_MDPP16_Short[j][k];
                if (GAGG_Number<=7){
                    GAGGS_Raw_hist_Total_Short[GAGG_Number]->Fill(values_MDPP16_Short[j][k]);
                    
                    double CalibratedEnergy_MDPP16_Short = (a_MDPP16_Short[GAGG_Number]*(values_MDPP16_Short[j][k]+RandomR3.Uniform()-0.5)+b_MDPP16_Short[GAGG_Number]);
                    if(CalibratedEnergy_MDPP16_Short > 100.0){ //Ponemos un Noise Cut 
                        hitCounter_GAGG_Short[GAGG_Number]++;
                        Energy_GAGG_Short[GAGG_Number] = CalibratedEnergy_MDPP16_Short;
                        GAGGS_Calib_hist_Total_Short[GAGG_Number]->Fill(CalibratedEnergy_MDPP16_Short);
                        GAGGS_Calib_Total_Short->Fill(CalibratedEnergy_MDPP16_Short);
                    }
                }   
            }
        }

        // Ahora viene el análisis de verdad
        for(int j = 0; j < 8; j++){
            if(hitCounter_GAGG_Long[j] == 1 && hitCounter_GAGG_Short[j] == 1){
                GAGGS_Short_vs_Long_hist[j]->Fill(Energy_GAGG_Long[j],(Energy_GAGG_Long[j]-Energy_GAGG_Short[j])/Energy_GAGG_Long[j]);
            }        
        }

        GAGGS_n_Coinci_hist->Fill(GAGGs_hit.size());

        if(GAGGs_hit.size() == 1){
            double Energy_0 = Energy_GAGG_Long[GAGGs_hit[0]];
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0);
            n1+=1;
        }
        bool matched_sum_2hit=false; //Energy_0+Energy_1 cerca de un pico
        if(GAGGs_hit.size() == 2){
            double Energy_0 = Energy_GAGG_Long[GAGGs_hit[0]];
            double Energy_1 = Energy_GAGG_Long[GAGGs_hit[1]];
            GAGGS_Coinci_Two_hist[GAGGs_hit[0]]->Fill(Energy_0);
            GAGGS_Coinci_Two_hist[GAGGs_hit[1]]->Fill(Energy_1);
            GAGGS_Coinci_Two_Full_hist->Fill(Energy_0);
            GAGGS_Coinci_Two_Full_hist->Fill(Energy_1);

            GAGGS_Coinci_2_hist->Fill((double)(GAGGs_hit[0]+1),(double)(GAGGs_hit[1]+1));
            GAGGS_Coinci_2_hist->Fill((double)(GAGGs_hit[1]+1),(double)(GAGGs_hit[0]+1));

            GAGGS_Multi_2_hist->Fill(Energy_0, Energy_1);
            GAGGS_Multi_2_hist->Fill(Energy_1, Energy_0);

            for (int j=0;j<nPeaks && matched_sum_2hit==0 ;j++ ){
                if (std::abs(Energy_0+Energy_1-Calibration_Values[j])<(0.043*Calibration_Values[j]+23)){
                    if (Energy_0>Energy_1){
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0+Energy_1);
                    }
                    else {
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_0+Energy_1);
                    }
                    matched_sum_2hit=true;
                    n2+=1;
                }
            }
            if (matched_sum_2hit==0){
                GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0);
                GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_1);
                n1+=1;
            }

            /*
            if(GAGGs_hit[1]-GAGGs_hit[0] ==1 ||GAGGs_hit[1]-GAGGs_hit[0] == 7){
                GAGGS_Coinci_2_hist->Fill(Energy_0, Energy_1);
                GAGGS_Coinci_2_hist->Fill(Energy_1, Energy_0);
            }
            */


            if(Energy_0 > 480 && Energy_0 < 550){
                GAGGS_Coinci_with_511_hist[GAGGs_hit[0]]->Fill(Energy_0);
                GAGGS_Coinci_with_511_hist[GAGGs_hit[1]]->Fill(Energy_1);
                GAGGS_Coinci_with_511_Full_hist->Fill(Energy_0);
                GAGGS_Coinci_with_511_Full_hist->Fill(Energy_1);
            }
            if(Energy_1 > 480 && Energy_1 < 550){
                GAGGS_Coinci_with_511_hist[GAGGs_hit[0]]->Fill(Energy_0);
                GAGGS_Coinci_with_511_hist[GAGGs_hit[1]]->Fill(Energy_1);
                GAGGS_Coinci_with_511_Full_hist->Fill(Energy_0);
                GAGGS_Coinci_with_511_Full_hist->Fill(Energy_1);
            }
        }

        bool matched_sum_3hit=false; //Energy_0+Energy_1+Energy_2 cerca de un pico
        matched_sum_2hit=false;
        if(GAGGs_hit.size() == 3){
            double Energy_0 = Energy_GAGG_Long[GAGGs_hit[0]];
            double Energy_1 = Energy_GAGG_Long[GAGGs_hit[1]];
            double Energy_2 = Energy_GAGG_Long[GAGGs_hit[2]];
            // Las 3 combinaciones posibles de pares
            GAGGS_3hit_PairSum_vs_Third_hist->Fill(Energy_0 + Energy_1, Energy_2);
            if (Energy_0>Energy_1 && Energy_0>Energy_2){
                GAGGS_3hit_Emax_vs_2Emin->Fill(Energy_0,Energy_1+Energy_2);
            }
            else if (Energy_1>Energy_0 && Energy_1>Energy_2){
                GAGGS_3hit_Emax_vs_2Emin->Fill(Energy_1,Energy_0+Energy_2);                
            }
            else {
                GAGGS_3hit_Emax_vs_2Emin->Fill(Energy_2,Energy_1+Energy_0);
            }
            for (int j=0;j<nPeaks && matched_sum_3hit==0 && matched_sum_2hit==0 ;j++ ){

                if (std::abs(Energy_0+Energy_1+Energy_2-Calibration_Values[j])<(0.043*Calibration_Values[j]+23)){
                    if (Energy_0>Energy_1 && Energy_0>Energy_2){
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0+Energy_1+Energy_2);
                    }
                    else if (Energy_1>Energy_0 && Energy_1>Energy_2){
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_0+Energy_1+Energy_2);
                    }
                    else {
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[2]]->Fill(Energy_0+Energy_1+Energy_2);

                    }
                    matched_sum_3hit=true;
                    n3+=1;
                }
                else if (std::abs(Energy_0+Energy_1-Calibration_Values[j])<(0.043*Calibration_Values[j]+23)){
                    if (Energy_0>Energy_1){
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0+Energy_1);
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[2]]->Fill(Energy_2);
                    }
                    else {
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_0+Energy_1);
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[2]]->Fill(Energy_2);
                    }
                    matched_sum_2hit=true;
                    n2+=1;
                    n1+=1;
                }
                else if (std::abs(Energy_0+Energy_2-Calibration_Values[j])<(0.043*Calibration_Values[j]+23)){
                    if (Energy_0>Energy_2){
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0+Energy_2);
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_1);
                    }
                    else {
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[2]]->Fill(Energy_0+Energy_2);
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_1);
                    }
                    matched_sum_2hit=true;
                    n2+=1;
                    n1+=1;
                }
                else if (std::abs(Energy_1+Energy_2-Calibration_Values[j])<(0.043*Calibration_Values[j]+23)){
                    if (Energy_1>Energy_2){
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_1+Energy_2);
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0);
                    }
                    else {
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[2]]->Fill(Energy_1+Energy_2);
                        GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0);
                    }
                    matched_sum_2hit=true;
                    n2+=1;
                    n1+=1;
                }
            }

            if (matched_sum_3hit==0 && matched_sum_2hit==0){
                GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0);
                GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_1);
                GAGGS_Add_back_hist_Total_Long[GAGGs_hit[2]]->Fill(Energy_2);
                n1+=1;
            }
        }
        if(GAGGs_hit.size() == 4){
            double Energy_0 = Energy_GAGG_Long[GAGGs_hit[0]];
            double Energy_1 = Energy_GAGG_Long[GAGGs_hit[1]];
            double Energy_2 = Energy_GAGG_Long[GAGGs_hit[2]];
            double Energy_3 = Energy_GAGG_Long[GAGGs_hit[3]];

            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_1);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[2]]->Fill(Energy_2);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[3]]->Fill(Energy_3);
        }
        if(GAGGs_hit.size() == 5){
            double Energy_0 = Energy_GAGG_Long[GAGGs_hit[0]];
            double Energy_1 = Energy_GAGG_Long[GAGGs_hit[1]];
            double Energy_2 = Energy_GAGG_Long[GAGGs_hit[2]];
            double Energy_3 = Energy_GAGG_Long[GAGGs_hit[3]];
            double Energy_4 = Energy_GAGG_Long[GAGGs_hit[4]];


            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_1);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[2]]->Fill(Energy_2);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[3]]->Fill(Energy_3);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[4]]->Fill(Energy_4);
        }

        if(GAGGs_hit.size() == 6){
            double Energy_0 = Energy_GAGG_Long[GAGGs_hit[0]];
            double Energy_1 = Energy_GAGG_Long[GAGGs_hit[1]];
            double Energy_2 = Energy_GAGG_Long[GAGGs_hit[2]];
            double Energy_3 = Energy_GAGG_Long[GAGGs_hit[3]];
            double Energy_4 = Energy_GAGG_Long[GAGGs_hit[4]];
            double Energy_5 = Energy_GAGG_Long[GAGGs_hit[5]];


            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_1);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[2]]->Fill(Energy_2);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[3]]->Fill(Energy_3);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[4]]->Fill(Energy_4);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[5]]->Fill(Energy_5);
        }
        if(GAGGs_hit.size() == 7){
            double Energy_0 = Energy_GAGG_Long[GAGGs_hit[0]];
            double Energy_1 = Energy_GAGG_Long[GAGGs_hit[1]];
            double Energy_2 = Energy_GAGG_Long[GAGGs_hit[2]];
            double Energy_3 = Energy_GAGG_Long[GAGGs_hit[3]];
            double Energy_4 = Energy_GAGG_Long[GAGGs_hit[4]];
            double Energy_5 = Energy_GAGG_Long[GAGGs_hit[5]];
            double Energy_6 = Energy_GAGG_Long[GAGGs_hit[6]];

            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_1);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[2]]->Fill(Energy_2);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[3]]->Fill(Energy_3);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[4]]->Fill(Energy_4);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[5]]->Fill(Energy_5);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[6]]->Fill(Energy_6);
        }
        if(GAGGs_hit.size() == 8){
            double Energy_0 = Energy_GAGG_Long[GAGGs_hit[0]];
            double Energy_1 = Energy_GAGG_Long[GAGGs_hit[1]];
            double Energy_2 = Energy_GAGG_Long[GAGGs_hit[2]];
            double Energy_3 = Energy_GAGG_Long[GAGGs_hit[3]];
            double Energy_4 = Energy_GAGG_Long[GAGGs_hit[4]];
            double Energy_5 = Energy_GAGG_Long[GAGGs_hit[5]];
            double Energy_6 = Energy_GAGG_Long[GAGGs_hit[6]];
            double Energy_7 = Energy_GAGG_Long[GAGGs_hit[7]];

            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[0]]->Fill(Energy_0);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[1]]->Fill(Energy_1);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[2]]->Fill(Energy_2);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[3]]->Fill(Energy_3);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[4]]->Fill(Energy_4);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[5]]->Fill(Energy_5);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[6]]->Fill(Energy_6);
            GAGGS_Add_back_hist_Total_Long[GAGGs_hit[7]]->Fill(Energy_7);
        }

        if (i % printFrequency == 0 || i == nEntries - 1) {  // Imprimir cada 1% o en la última iteración
            double progress = (double(i) / nEntries) * 100.0;
            std::cout << " Building histos Raw: " << progress << "% c.\r";
            std::cout.flush();  // Asegurar que la salida se actualice en la consola
        }
    
    }
    std::cout<<"El número de eventos que han tocado una única vez es de: "<<n1<<std::endl;
    std::cout<<"El número de eventos que han tocado dos veces es de: "<<n2<<std::endl;
    std::cout<<"El número de eventos que han tocado tres veces es de: "<<n3<<std::endl;


    chain->GetEvent(nEntries);//Vamos a hacer que aparezca por pantalla el tiempo que ha estado midiendo.
    std::cout << "Time Stamp: " << EventTimeStamp_MDPP16[0] << "\n";
    std::cout << "Extended Time Stamp: " << ExtTimeStamp_MDPP16[0] << "\n";
    std::cout << "Loops: " << Loops_Time_Stamp << "\n";
    double Time;
    Time = (EventTimeStamp_MDPP16[0] + Loops_Time_Stamp*pow(2,30))/16.0; // En us
    double Time_hours = Time*pow(10,-6)/3600.0;
    std::cout<<"Ha tardado "<< Time_hours << " horas"<<"\n";

    //Vamos a llenar los histogramas de residuos y FWHM, También vamos a ajustar los picos a gaussianas+rectas
    // Abrir el fichero TXT antes del bucle principal
    std::ofstream outFile("Calibracion/Ajuste/Fit_Gauss_Lineal_Long" + filename + ".txt");
    outFile << "GAGG \t Peak_Energy \t Mean \t Sigma \t Slope \t Y_Intercept \t FWHM \t Area \n";
    outFile << std::fixed << std::setprecision(4);

    std::ofstream Areas_Long("Calibracion/Areas/Area_Long" + filename + ".txt");
    std::ofstream Err_Areas_Long("Calibracion/Areas/Area_error_Long" + filename + ".txt");
    // Cabecera


    for (int i = 0; i < 8; i++) {
        //TH1F* histo = GAGGS_Calib_hist_Total_Long[i];
        TH1F* histo = GAGGS_Add_back_hist_Total_Long[i];
        double residuo_acumulado_Long = 0;
        double residuo_acumulado_add_back = 0;

        for (size_t j = 0; j < nPeaks; j++) {
            double peakEnergy = Calibration_Values[j];

            double amplitud, mean, sigma, slope, y_inter,area, cov_01;
            double err_amp, err_mean, err_sigma, err_slope, err_y_inter, err_area_gaus,err_area=1;

            Ajuste_Gauss_Lineal(histo, peakEnergy, Ranges_values[j], Ranges_values[j], amplitud, mean, sigma, slope, y_inter,area, err_amp, err_mean, err_sigma, err_slope, err_y_inter, err_area_gaus, cov_01);

            double xmin = mean - Ranges_values[j];
            double xmax = mean + Ranges_values[j];

            TF1* fit_gaus = new TF1("fit_gaus", "gaus", xmin, xmax);
            fit_gaus->SetParameter(0, amplitud);
            fit_gaus->SetParameter(1, mean);
            fit_gaus->SetParameter(2, sigma);
            int color_id = TColor::GetColor(0xFF, 0x57, 0x33);
            fit_gaus->SetLineColor(color_id);
            fit_gaus->SetLineWidth(2);
            fit_gaus->SetLineStyle(2);

            TF1* fit_lin = new TF1("fit_lin", "pol1", xmin, xmax);
            fit_lin->SetParameter(0, y_inter);
            fit_lin->SetParameter(1, slope);
            fit_lin->SetLineColor(kGreen);
            fit_lin->SetLineWidth(2);
            fit_lin->SetLineStyle(2);

            //GAGGS_Calib_hist_Total_Long[i]->GetListOfFunctions()->Add(fit_gaus);
            //GAGGS_Calib_hist_Total_Long[i]->GetListOfFunctions()->Add(fit_lin);
            GAGGS_Add_back_hist_Total_Long[i]->GetListOfFunctions()->Add(fit_gaus);
            GAGGS_Add_back_hist_Total_Long[i]->GetListOfFunctions()->Add(fit_lin);

            double residuo = mean - peakEnergy;
            double fwhm = sigma * 2.355;
            residuo_acumulado_Long += residuo;

            int bin1 = histo->FindBin(xmin);
            int bin2 = histo->FindBin(xmax);
            double bin_width = histo->GetBinWidth(1);

            double area_hist = histo->Integral(bin1,bin2);
            double area_fondo = (xmax-xmin)*(y_inter+slope*(xmin+xmax)/2)/bin_width;
            double area_total = area_hist-area_fondo;

            double err_area_hist = TMath::Sqrt(area_hist);
            double err_area_fondo1 = (xmax-xmin)/bin_width*TMath::Sqrt(err_y_inter*err_y_inter+(err_slope/2)*(xmin+xmax)*(err_slope/2)*(xmin+xmax)+(xmax+xmin)*cov_01);
            double err_area_fondo2 = TMath::Sqrt(area_fondo);
            double err_area_total1 = TMath::Sqrt(err_area_hist*err_area_hist+err_area_fondo1*err_area_fondo1);
            double err_area_total2 = TMath::Sqrt(err_area_hist*err_area_hist+err_area_fondo2*err_area_fondo2);    

            GAGGS_Residuo_Long[i]->Fill(peakEnergy, residuo);
            GAGGS_FWHM_Long[i]->Fill(peakEnergy, fwhm);

            // Escribir fila en el TXT
            outFile << i+1        << "\t"
            << peakEnergy << "\t"
            << mean       << "\t" << err_mean  << "\t"
            << sigma      << "\t" << err_sigma << "\t"
            << slope      << "\t" << err_slope << "\t"
            << y_inter    << "\t" << err_y_inter<< "\t"
            << fwhm       << "\t" << err_sigma * 2.355<<"\t"
            <<area<<"\t"<<err_area_total1<<"\t"<<err_area_total2<<"\t"<<err_area_gaus<<"\n";

        }

        // Residuo acumulado por GAGG como fila especial
        outFile << "# GAGG " << i+1 << " Residuo_Acumulado\t" << residuo_acumulado_Long << "\n";
    }
        for (int i = 0; i < 8; i++) {
        TH1F* histo = GAGGS_Calib_hist_Total_Long[i];
        //TH1F* histo = GAGGS_Add_back_hist_Total_Long[i];
        double residuo_acumulado_Long = 0;
        double residuo_acumulado_add_back = 0;

        for (size_t j = 0; j < nPeaks; j++) {
            double peakEnergy = Calibration_Values[j];

            double amplitud, mean, sigma, slope, y_inter,area, cov_01;
            double err_amp, err_mean, err_sigma, err_slope, err_y_inter, err_area_gaus,err_area=1;

            Ajuste_Gauss_Lineal(histo, peakEnergy, Ranges_values[j], Ranges_values[j], amplitud, mean, sigma, slope, y_inter,area, err_amp, err_mean, err_sigma, err_slope, err_y_inter, err_area_gaus, cov_01);

            double xmin = mean - Ranges_values[j];
            double xmax = mean + Ranges_values[j];

            TF1* fit_gaus = new TF1("fit_gaus", "gaus", xmin, xmax);
            fit_gaus->SetParameter(0, amplitud);
            fit_gaus->SetParameter(1, mean);
            fit_gaus->SetParameter(2, sigma);
            int color_id = TColor::GetColor(0xFF, 0x57, 0x33);
            fit_gaus->SetLineColor(color_id);
            fit_gaus->SetLineWidth(2);
            fit_gaus->SetLineStyle(2);

            TF1* fit_lin = new TF1("fit_lin", "pol1", xmin, xmax);
            fit_lin->SetParameter(0, y_inter);
            fit_lin->SetParameter(1, slope);
            fit_lin->SetLineColor(kGreen);
            fit_lin->SetLineWidth(2);
            fit_lin->SetLineStyle(2);

            GAGGS_Calib_hist_Total_Long[i]->GetListOfFunctions()->Add(fit_gaus);
            GAGGS_Calib_hist_Total_Long[i]->GetListOfFunctions()->Add(fit_lin);
            //GAGGS_Add_back_hist_Total_Long[i]->GetListOfFunctions()->Add(fit_gaus);
            //GAGGS_Add_back_hist_Total_Long[i]->GetListOfFunctions()->Add(fit_lin);

            double residuo = mean - peakEnergy;
            double fwhm = sigma * 2.355;
            residuo_acumulado_Long += residuo;

            int bin1 = histo->FindBin(xmin);
            int bin2 = histo->FindBin(xmax);
            double bin_width = histo->GetBinWidth(1);

            double area_hist = histo->Integral(bin1,bin2);
            double area_fondo = (xmax-xmin)*(y_inter+slope*(xmin+xmax)/2)/bin_width;
            double area_total = area_hist-area_fondo;

            double err_area_hist = TMath::Sqrt(area_hist);
            double err_area_fondo1 = (xmax-xmin)/bin_width*TMath::Sqrt(err_y_inter*err_y_inter+(err_slope/2)*(xmin+xmax)*(err_slope/2)*(xmin+xmax)+(xmax+xmin)*cov_01);
            double err_area_fondo2 = TMath::Sqrt(area_fondo);
            double err_area_total1 = TMath::Sqrt(err_area_hist*err_area_hist+err_area_fondo1*err_area_fondo1);
            double err_area_total2 = TMath::Sqrt(err_area_hist*err_area_hist+err_area_fondo2*err_area_fondo2);    

            GAGGS_Residuo_Long[i]->Fill(peakEnergy, residuo);
            GAGGS_FWHM_Long[i]->Fill(peakEnergy, fwhm);

            // Escribir fila en el TXT
            outFile << i+1        << "\t"
            << peakEnergy << "\t"
            << mean       << "\t" << err_mean  << "\t"
            << sigma      << "\t" << err_sigma << "\t"
            << slope      << "\t" << err_slope << "\t"
            << y_inter    << "\t" << err_y_inter<< "\t"
            << fwhm       << "\t" << err_sigma * 2.355<<"\t"
            <<area<<"\t"<<err_area_total1<<"\t"<<err_area_total2<<"\t"<<err_area_gaus<<"\n";

        }

        // Residuo acumulado por GAGG como fila especial
        outFile << "# GAGG " << i+1 << " Residuo_Acumulado\t" << residuo_acumulado_Long << "\n";
    }
    outFile.close();
        
    //Caso short
    std::ofstream outFileShort("Calibracion/Ajuste/Fit_Gauss_Lineal_Short" + filename + ".txt");
    outFileShort << std::fixed << std::setprecision(4);
        // Cabecera
    outFileShort << "GAGG\tPeak_Energy\tMean\tSigma\tSlope\tY_Intercept\tFWHM\n";

    for (int i = 0;i < 8; i++){//Lo mismo para Short
        TH1F* histo = GAGGS_Calib_hist_Total_Short[i];
        double residuo_acumulado_Short = 0;
        for (int j = 0;j < nPeaks;j++){
            double peakEnergy = Calibration_Values[j];
            double amplitud,mean,sigma,slope,y_inter,area, cov_01;
            double err_amp, err_mean, err_sigma, err_slope, err_y_inter, err_area_gaus,err_area_hist,err_area=1;


            Ajuste_Gauss_Lineal(histo, peakEnergy, Ranges_values[j], Ranges_values[j], amplitud, mean, sigma, slope, y_inter,area, err_amp, err_mean, err_sigma, err_slope, err_y_inter, err_area_gaus, cov_01);
            double xmin = mean - Ranges_values[j];
            double xmax = mean + Ranges_values[j];

            TF1* fit_gaus = new TF1("fit_gaus", "gaus", xmin, xmax);
            fit_gaus->SetParameter(0, amplitud);
            fit_gaus->SetParameter(1, mean);
            fit_gaus->SetParameter(2, sigma);
            fit_gaus->SetLineColor(kPink);
            fit_gaus->SetLineWidth(2);
            fit_gaus->SetLineStyle(2);

            TF1* fit_lin = new TF1("fit_lin", "pol1", xmin, xmax);
            fit_lin->SetParameter(0, y_inter);
            fit_lin->SetParameter(1, slope);
            fit_lin->SetLineColor(kGreen);
            fit_lin->SetLineWidth(2);
            fit_lin->SetLineStyle(2);

            GAGGS_Calib_hist_Total_Short[i]->GetListOfFunctions()->Add(fit_gaus);
            GAGGS_Calib_hist_Total_Short[i]->GetListOfFunctions()->Add(fit_lin);

            double residuo = mean-peakEnergy;
            double fwhm = sigma*2.355;
            residuo_acumulado_Short = residuo_acumulado_Short + residuo;

            int bin1=histo->FindBin(xmin);
            int bin2=histo->FindBin(xmax);

            double integral = histo->Integral(bin1,bin2);

            GAGGS_Residuo_Short[i]->Fill(peakEnergy,residuo);
            GAGGS_FWHM_Short[i]->Fill(peakEnergy,fwhm);

            // Escribir fila en el TXT
            outFileShort << i+1        << "\t"
            << peakEnergy << "\t"
            << mean       << "\t" << err_mean  << "\t"
            << sigma      << "\t" << err_sigma << "\t"
            << slope      << "\t" << err_slope << "\t"
            << y_inter    << "\t" << err_y_inter<< "\t"
            << fwhm       << "\t" << err_sigma * 2.355<<"\t"
            << area       << "\t" << err_area << "\n";
        }

        // Residuo acumulado por GAGG como fila especial
        outFileShort << "# GAGG " << i+1 << " Residuo_Acumulado\t" << residuo_acumulado_Short << "\n";

    }
    outFileShort.close();

    Save_GAGGs(outputFile, GAGGS_Raw_hist_Total_Long, GAGGS_Raw_hist_Total_Short, GAGGS_Calib_hist_Total_Long, GAGGS_Calib_hist_Total_Short,
    GAGGS_Calib_Total_Long, GAGGS_Calib_Total_Short, GAGGS_Residuo_Long, GAGGS_Residuo_Short, GAGGS_FWHM_Long, GAGGS_FWHM_Short, 
    GAGGS_Coinci_Two_hist, GAGGS_Coinci_Two_Full_hist, GAGGS_Coinci_with_511_hist, GAGGS_Coinci_with_511_Full_hist, 
    GAGGS_Multi_2_hist, GAGGS_Coinci_2_hist, GAGGS_Short_vs_Long_hist, GAGGs_Check_Calib_hist, GAGGS_Add_back_hist_Total_Long, GAGGS_n_Coinci_hist,GAGGS_3hit_PairSum_vs_Third_hist);


    if(Substract_BKG){
        // Histogramas background calibrado Long
        std::vector<double> Calibration_Values_BKG = Select_Calibration_Values("BKG");
        const int nPeaks_BKG = Calibration_Values_BKG.size();

        char GAGGS_BKG_name_Long[8][100];
        char GAGGS_BKG_name_histo_Long[8][200];
        static TH1F* GAGGS_BKG_hist_Long[8] = {nullptr};

        char GAGGS_BKG_name_Short[8][100];
        char GAGGS_BKG_name_histo_Short[8][200];
        static TH1F* GAGGS_BKG_hist_Short[8] = {nullptr};

        for(int i = 0; i < 8; i++){
            sprintf(GAGGS_BKG_name_Long[i],"GAGG_%i BKG Long", i+1);
            sprintf(GAGGS_BKG_name_histo_Long[i], "GAGG_%i BKG Long; Energia (KeV); Counts", i+1);
            GAGGS_BKG_hist_Long[i] = new TH1F(GAGGS_BKG_name_Long[i], GAGGS_BKG_name_histo_Long[i], 5000/2.0, 0, 5000);

            sprintf(GAGGS_BKG_name_Short[i], "GAGG_%i BKG Short", i+1);
            sprintf(GAGGS_BKG_name_histo_Short[i], "GAGG_%i BKG Short; Energia (KeV); Counts", i+1);
            GAGGS_BKG_hist_Short[i] = new TH1F(GAGGS_BKG_name_Short[i],GAGGS_BKG_name_histo_Short[i], 5000/2.0, 0, 5000);

            for(int j = 0; j < nPeaks_BKG; j++){
                Draw_Calib_lines(GAGGS_BKG_hist_Long[i], Calibration_Values_BKG[j]);
                Draw_Calib_lines(GAGGS_BKG_hist_Short[i], Calibration_Values_BKG[j]);
            }

        }

        TH1F* GAGGS_BKG_Total_Long = new TH1F("GAGG_Total_BKG_Long","GAGG Total BKG Long; Energia (KeV); Counts", 5000, 0, 5000);
        TH1F* GAGGS_BKG_Total_Short = new TH1F("GAGG_Total_BKG_Short","GAGG Total BKG Short; Energia (KeV); Counts", 5000, 0, 5000);
        for(int j = 0; j < nPeaks_BKG; j++){
            Draw_Calib_lines(GAGGS_BKG_Total_Long, Calibration_Values_BKG[j]);
            Draw_Calib_lines(GAGGS_BKG_Total_Short, Calibration_Values_BKG[j]);
        }

        double a_MDPP16_Long_BKG[8], b_MDPP16_Long_BKG[8];
        double a_MDPP16_Short_BKG[8], b_MDPP16_Short_BKG[8];
        LoadMDPP16CalibParameters(BKG_Calib_filename, a_MDPP16_Long_BKG, b_MDPP16_Long_BKG, a_MDPP16_Short_BKG, b_MDPP16_Short_BKG);
        

        Long64_t nEntriesBKG = chainBKG->GetEntries();

        std::cout << "Number of BKG Entries: " << nEntriesBKG << std::endl;
        int printFreqBKG = nEntriesBKG / 100;

        int Loops_Time_Stamp_BKG = 0; //Por si no tengo Extended Timestamp
        int Last_Timestamp_BKG = 0;

        for(int i = 0; i < int(nEntriesBKG); i++){
            chainBKG->GetEvent(i);
            if(Last_Timestamp_BKG > EventTimeStamp_MDPP16[0]){
                Loops_Time_Stamp_BKG++;
            }
            Last_Timestamp_BKG = EventTimeStamp_MDPP16[0];

            for(int j = 0; j < mdpp16_counter; j++){
                for(Int_t k = 0; k < mul_MDPP16_Long[j]; k++){
                    int GAGG_Number = detCh_MDPP16_Long[j][k];
                    double CalibratedEnergy_MDPP16_Long = (a_MDPP16_Long_BKG[GAGG_Number]*(values_MDPP16_Long[j][k]+RandomR3.Uniform()-0.5)+b_MDPP16_Long_BKG[GAGG_Number]);
                    GAGGS_BKG_hist_Long[GAGG_Number]->Fill(CalibratedEnergy_MDPP16_Long);
                    GAGGS_BKG_Total_Long->Fill(CalibratedEnergy_MDPP16_Long);
                }
                for(Int_t k = 0; k < mul_MDPP16_Short[j]; k++){
                    int GAGG_Number = detCh_MDPP16_Short[j][k];
                    double CalibratedEnergy_MDPP16_Short = (a_MDPP16_Short_BKG[GAGG_Number]*(values_MDPP16_Short[j][k]+RandomR3.Uniform()-0.5)+b_MDPP16_Short_BKG[GAGG_Number]);
                    GAGGS_BKG_hist_Short[GAGG_Number]->Fill(CalibratedEnergy_MDPP16_Short);
                    GAGGS_BKG_Total_Short->Fill(CalibratedEnergy_MDPP16_Short);
                }
            }

            if(i % printFreqBKG == 0 || i == nEntriesBKG - 1){
                std::cout << " Building BKG histos: "
                        << (double(i)/nEntriesBKG)*100.0 << "% \r";
                std::cout.flush();
            }
        }

        chainBKG->GetEvent(nEntriesBKG);//Vamos a hacer que aparezca por pantalla el tiempo que ha estado midiendo.
        std::cout << "Time Stamp: "<< EventTimeStamp_MDPP16[0] << "\n";
        std::cout << "Extended Time Stamp: "<< ExtTimeStamp_MDPP16[0] << "\n";
        std::cout << "Loops: "<< Loops_Time_Stamp_BKG << "\n";
        double Time_BKG;
        Time_BKG = (EventTimeStamp_MDPP16[0] + Loops_Time_Stamp_BKG*pow(2,30))/16.0;
        double Time_BKG_hours = Time_BKG*pow(10,-6)/3600.0;
        std::cout<<"El Background ha medido: "<< Time_BKG_hours <<" horas"<<"\n";

        double ratio = Time / Time_BKG;
        std::cout << "Factor de normalización BKG: " << ratio << std::endl;

        for(int i = 0; i < 8; i++){
            GAGGS_BKG_hist_Long[i]->Scale(ratio);
            GAGGS_BKG_hist_Short[i]->Scale(ratio);
        }
        GAGGS_BKG_Total_Long->Scale(ratio);
        GAGGS_BKG_Total_Short->Scale(ratio);

        Save_BKG(outputFile, GAGGS_BKG_hist_Long, GAGGS_BKG_hist_Short, GAGGS_BKG_Total_Long, GAGGS_BKG_Total_Short);
    }
}