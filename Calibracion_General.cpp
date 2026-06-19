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

#include "General.cpp"

// Si por ejemplo este código llamara a otros códigos/scripts (es común y lo harás en el futuro)
// Se escribe: 
// #include "ruta al codigo" (por ejemplo "src/Codigo2.cpp")

/////////////////////////////////////////
/////// Ir al final del Código /////////
////////////////////////////////////////

//Veamos como leemos los datos. Esta función es un void porque no devuelve nada solo analiza/cambia los paránmetros que le metemos
void Calibracion_General(TChain* chain, TFile* outputFile, std::string filename, std::string isotope 
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
    }

    //Defino 8 Histogramas uno por GAGG
    char GAGGS_Raw_name_Total_Long[8][100]; //Nombre en el archivo outputFile. Son caracteres. [8] indica el número de histogramas y [100] la longitud de la cadena de caracteres
    char GAGGS_Raw_name_histo_Total_Long[8][200]; //Nombre del histograma
    static TH1F* GAGGS_Raw_hist_Total_Long[8]= {nullptr}; //Histogramas 1 dimensional TH1 y la F es para floats. Lo inicializo a nullptr que significa puntero nulo para evitar fallos

    char GAGGS_Raw_name_Total_Short[8][100]; //Lo mismo pero para el short
    char GAGGS_Raw_name_histo_Total_Short[8][200]; 
    static TH1F* GAGGS_Raw_hist_Total_Short[8]= {nullptr};

    //Ahora hay que darle nombre a los histogramas
    for (int i = 0; i < 8 ; i++){
        sprintf(GAGGS_Raw_name_Total_Long[i],"GAGG_%i Raw Long", i+1); //Nombre del histograma en el archivo %i significa que ahi metes un entero
        sprintf(GAGGS_Raw_name_histo_Total_Long[i],"GAGG_%i Raw Long; Channel ; Counts", i+1); //Nombre del histogram ; eje X; eje Y
        GAGGS_Raw_hist_Total_Long[i] = new TH1F(GAGGS_Raw_name_Total_Long[i],GAGGS_Raw_name_histo_Total_Long[i],65536.0/4.0,0,65536); // Le metes los nombre, número de bines, bin inicial, bin final   
    
        sprintf(GAGGS_Raw_name_Total_Short[i],"GAGG_%i Raw Short", i+1); //Lo mismo pero con short
        sprintf(GAGGS_Raw_name_histo_Total_Short[i],"GAGG_%i Raw Short; Channel ; Counts", i+1); 
        GAGGS_Raw_hist_Total_Short[i] = new TH1F(GAGGS_Raw_name_Total_Short[i],GAGGS_Raw_name_histo_Total_Short[i],65536.0/4.0,0,65536); 
  
    }

    //Cogemos el número de entradas para hacernos una idea
    Long64_t nEntries = chain->GetEntries();
    cout << "Number of Entries: " << nEntries << endl;
    int printFrequency = nEntries / 100;

    //Bucle sobre las entradas
    for (int i = 0; i < int(nEntries); i++){
        //Cogemos el evento, super importante
        chain -> GetEvent(i);

        for (int j = 0; j < mdpp16_counter; j++){// Bucle por si hubiera mas MDPP16
            for (Int_t k = 0; k < mul_MDPP16_Long[j]; k++){ 
                int idx= detCh_MDPP16_Long[j][k];
                if (idx<0|| idx>=8) continue; //Nos quedamos sólo con los canales con GAGGs conectados
                GAGGS_Raw_hist_Total_Long[detCh_MDPP16_Long[j][k]]->Fill(values_MDPP16_Long[j][k]);//Values son los canales electrónicos, detCh los físicos (detectores)
            }
            for (Int_t k = 0; k < mul_MDPP16_Short[j]; k++){//Bucle de llenado short
                int idx= detCh_MDPP16_Long[j][k];
                if (idx<0|| idx>=8) continue; //Nos quedamos sólo con los canales con GAGGs conectados
                GAGGS_Raw_hist_Total_Short[detCh_MDPP16_Short[j][k]]->Fill(values_MDPP16_Short[j][k]);//Values son los canales electrónicos, detCh los físicos (detectores)
            }
        }
        if (i % printFrequency == 0 || i == nEntries - 1) {  // Imprimir cada 1% o en la última iteración
            double progress = (double(i) / nEntries) * 100.0;
            std::cout << " Building histos Raw: " << progress << "% c.\r";
            std::cout.flush();  // Asegurar que la salida se actualice en la consola
        }
    }
    
    //Pasamos a la parte de Calibración
    std::vector<double> Calibration_Values = Select_Calibration_Values(isotope);
    std::vector<double> Ranges_values = Select_Calibration_Ranges(isotope);
    const int nPeaks = Calibration_Values.size();

    std::vector<double> centroides_raw_Long(nPeaks,0.0);
    std::vector<double> err_centroides_Long(nPeaks,0.0);
    std::vector<double> centroides_raw_Short(nPeaks,0.0);
    std::vector<double> err_centroides_Short(nPeaks,0.0);
    
    double a_Long[8], err_a_Long[8];
    double b_Long[8], err_b_Long[8];
    double a_Short[8], err_a_Short[8];
    double b_Short[8], err_b_Short[8];

    for (int i = 0; i < 8;i++){//TSpectrum
        double noise_cut;

        TH1F* histo = GAGGS_Raw_hist_Total_Long[i];
        if (isotope=="226Ra"){
            if (i == 0){
                noise_cut=3000;
            }
            else if (i==1||i==2||i==3||i==5||i==7){
                noise_cut=2800;
            }
            else if (i==4){
                noise_cut=2050;
            }
            else if (i==6){
                noise_cut=2400;
            }

        }

        else {
            if (i == 4 || i == 5){
                noise_cut=950;//2000 para el Bi//canal a partir del que empezamos a buscar
            }
            else {
                noise_cut=1050;//2000 para el Bi
            }
            //Para la medida de europio de isolde
            //if (i==0|| i==1){
            //  noise_cut=2700;
            //}
            //if (i==2|| i==3 || i==5 || i==7){
                //noise_cut=2400;
            //}
            //if (i==4 ||i==6){
            //  noise_cut=2000;
            //}
        }
        int GAGG5;
        int GAGG7;
        std::vector<double> peaks;
        if (i==4){
            GAGG5=1;
            GAGG7=0;
        }
        else if (i==6){
            GAGG7=1;
            GAGG5=0;
        }
        else{
            GAGG5=0;
            GAGG7=0;
        }
        Buscar_Picos(histo, noise_cut, peaks, nPeaks, isotope,GAGG5,GAGG7);
        //std::cout << "Numero de picos encontrados: " << peaks.size() << "\n";

        //for (auto& p : peaks) std::cout << p << " ";
        //std::cout << std::endl;
        //std::cout<<"GAGG "<<i+1<<" Long"<<std::endl;

        centroides_raw_Long.clear();
        err_centroides_Long.clear();

        for (int j = 0; j < peaks.size(); j++){
            double peakValue = peaks[j];

            double amplitud,mean,sigma,slope,y_inter,area, cov_01;
            double err_amp, err_mean, err_sigma, err_slope, err_y_inter, err_area;

            Ajuste_Gauss_Lineal(histo, peakValue, Ranges_values[j],Ranges_values[j], amplitud, mean, sigma, slope, y_inter, area,err_amp, err_mean, err_sigma, err_slope, err_y_inter, err_area, cov_01);//Ajuste gaussiano+lineal
            centroides_raw_Long.push_back(mean);
            err_centroides_Long.push_back(err_mean);            
        }
        Ajuste_Lineal(centroides_raw_Long, Calibration_Values, err_centroides_Long, a_Long[i], b_Long[i], err_a_Long[i], err_b_Long[i]);
    }

    for (int i = 0; i < 8; i++){//TSpectrum para short
        double noise_cut;
        TH1F* histo = GAGGS_Raw_hist_Total_Short[i];
        if (i==4 || i==5){
            noise_cut=650;//canal a partir del que empezamos a buscar
        }
        else {
            noise_cut=800;
        }
        std::vector<double> peaks;
        int GAGG5;
        int GAGG7;
        if (i==4){
            GAGG5=1;
            GAGG7=0;
        }
        else if (i==6){
            GAGG7=1;
            GAGG5=0;
        }
        else{
            GAGG5=0;
            GAGG7=0;
        }
        Buscar_Picos(histo, noise_cut, peaks, nPeaks, isotope, GAGG5,GAGG7);
        //std::cout << "Numero de picos encontrados: " << peaks.size() << "\n";
        //std::cout<<"GAGG "<<i+1<<" Short"<<std::endl;
        centroides_raw_Short.clear();
        err_centroides_Short.clear();

        for (int j = 0; j < peaks.size(); j++){
            double peakValue = peaks[j];

            double amplitud,mean,sigma,slope,y_inter,area, cov_01;
            double err_amp, err_mean, err_sigma, err_slope, err_y_inter, err_area;

            Ajuste_Gauss_Lineal(histo, peakValue, Ranges_values[j],Ranges_values[j], amplitud, mean, sigma, slope, y_inter, area,err_amp, err_mean, err_sigma, err_slope, err_y_inter, err_area,cov_01);//Ajuste gaussiano+lineal
            centroides_raw_Short.push_back(mean);
            err_centroides_Short.push_back(err_mean);
        }
        //cout << " Croto " <<  centroides_raw_Short.si
        Ajuste_Lineal(centroides_raw_Short, Calibration_Values, err_centroides_Short, a_Short[i], b_Short[i], err_a_Short[i], err_b_Short[i]);
    }

    //Extraer y guardar calibración
    std::string ArchCalib = "Calibracion/Calibration_" + filename + ".txt";
    std::string ErrArchCalib = "Calibracion/Calibration_errors_" + filename + ".txt";
    
    // Guardar calibración manual en txt
    std::ofstream calibFile(ArchCalib);
    if (!calibFile.is_open()) {
        std::cerr << "Error: no se puede abrir el fichero de calibración" << ArchCalib <<  std::endl;
    } else {
        for (int i = 0; i < 8; i++) {
            calibFile << std::fixed << std::setprecision(8)
                    << i << "  "
                    << a_Long[i]  << "  " << b_Long[i]  << "  "
                    << a_Short[i] << "  " << b_Short[i]
                    << std::endl;
        }
        calibFile.close();
        std::cout << "Calibración guardada en: " << ArchCalib << "\n";
    }

    std::ofstream err_calibFile(ErrArchCalib);
    if (!err_calibFile.is_open()) {
        std::cerr << "Error: no se puede abrir el fichero de incertidumbres de calibración" << ErrArchCalib<< std::endl;
    } else {
        for (int i = 0; i < 8; i++) {
            err_calibFile << std::fixed << std::setprecision(8)
                    << i << "  "
                    << err_a_Long[i]  << "  " << err_b_Long[i]  << "  "
                    << err_a_Short[i] << "  " << err_b_Short[i]
                    << std::endl;
        }
        err_calibFile.close();
    }

    //Desactivo la interfaz gráfica para que vaya mas rapido y si algo esta abierto se cierre
    gROOT->SetBatch(kTRUE);

    //Creo una carpeta llamada GAGGs
    TDirectory* GAGGsDir =  outputFile->mkdir("GAGGs");

    //Dentro creo una carpeta para los histogramas que hemos creado
    TDirectory* GAGG_Raw_Long_Dir = GAGGsDir->mkdir("Raw Integration Long");
    TDirectory* GAGG_Raw_Short_Dir = GAGGsDir->mkdir("Raw Integration Short");

    //Meto los histogramas ahí
    for ( int i = 0; i < 8; i++){
        Save_Histo(GAGGS_Raw_hist_Total_Long[i],GAGG_Raw_Long_Dir);
        Save_Histo(GAGGS_Raw_hist_Total_Short[i],GAGG_Raw_Short_Dir);
    }
}
