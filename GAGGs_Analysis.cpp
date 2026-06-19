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
#include <TStyle.h>
#include <TTimer.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include "TMatrixD.h"


// Incluir códigos extra
//#include "src/Calibracion_General.cpp"
#include "src/Sort_GAGGs.cpp"

bool Calibrar = 0;
bool Sort = 1;

int GAGGs_Analysis(){ 
    gStyle->SetOptStat(1111111); 
    std::vector<std::string> Runs_Analizar = { 
        "ROOT/IS690B_Calb_run003_251021_173736_part001.root"
    };

    int Amount_of_Runs = Runs_Analizar.size();

    for(int i = 0; i < Amount_of_Runs; i++){
        TChain* eventChain = new TChain("EventTree");

        //Cojo la ruta de los archivos a analizar
        std::string direction = Runs_Analizar[i]; 
        size_t rootPos = direction.find(".root"); 
        size_t backslashPos = direction.rfind('/', rootPos);
        std::string Name = direction.substr(backslashPos + 1, rootPos - backslashPos - 1);
        eventChain->Add(direction.c_str());
        std::cout << "Procesing "<< Name << ".root" << "\n";
        std::cout << "-------------------------- "<< "\n";

        //std::string isotope = "60Co";
        std::string isotope = "152Eu";
        //std::string isotope = "22Na";
        //std::string isotope = "137Cs";
        //std::string isotope = "207Bi";
        //std::string isotope = "226Ra";
        //std::string isotope = "BKG";

        if(Calibrar){
            std::vector<std::string> Calib_Nombre ={"_Calib.root"};
            std::string Appendix_Calib = Calib_Nombre[0];
            std::string outputFileName_Calib = "./Analized/"+ Name + Appendix_Calib; 
            TFile* outputFile_Calib = new TFile(outputFileName_Calib.c_str(),"recreate");
            std::cout << "Output Name: "<< outputFileName_Calib  << "\n";
            std::cout << "-------------------------- " << "\n"; 

            TStopwatch timer;
            timer.Start();
                    
            //Calibracion_General(eventChain, outputFile_Calib, Name, isotope);

            timer.Stop();
            timer.Print();
        }
    
        if(Sort){
            std::vector<std::string> Nombre ={"_Analized.root"}; 
            TChain* eventChainBKG = new TChain("EventTree");
            std::string Appendix = Nombre[0];
            std::string outputFileName_Analized = "./Analized/"+ Name + Appendix;
            TFile* outputFile_Analized = new TFile(outputFileName_Analized.c_str(),"recreate");
            std::cout << "Output Name: "<< outputFileName_Analized  << "\n";
            std::cout << "-------------------------- " << "\n";
            
            eventChainBKG->Add("ROOT/Practicas_Pablo_BKGrun001_260306_161317_part001.root");

            //bool Substract_BKG = true;
            bool Substract_BKG = false;

            std::string Calib_filename = "Calibracion/Calibration_IS690B_Calb_run003_251021_173736_part001.txt";
            std::string BKG_Calib_filename = "Calibracion/Calibration_Practicas_Pablo_22Na_Coincidencias_SetUp1_run001_260317_112856_part001.txt";
            //std::string BKG_Calib_filename = "Calibracion/Calibration_Prueba.dat";
            TStopwatch timer;
            timer.Start();

            Sort_GAGGs(eventChain, eventChainBKG, outputFile_Analized, Name, isotope, Calib_filename, Substract_BKG, BKG_Calib_filename);

            timer.Stop();
            timer.Print();        
        }                                                                                                                                                                                                                                                                           
    }

    //Para poder visualizar lo reactiva porque arriba lo hemos desactivado
    gROOT->SetBatch(kFALSE);
    
    // Tiene que devolver un entero y normalmente se pone 0 porque eso significa que todo fue bien
    
    return 0;
}

