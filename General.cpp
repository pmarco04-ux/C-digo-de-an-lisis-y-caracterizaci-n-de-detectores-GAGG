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

#include "General.hh"
#include <stdexcept>
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMatrixD.h"


// Si por ejemplo este código llamara a otros códigos/scripts (es común y lo harás en el futuro)
// Se escribe: 
// #include "ruta al codigo" (por ejemplo "src/Codigo2.cpp")

/////////////////////////////////////////
/////// Ir al final del Código /////////
////////////////////////////////////////

void LoadMDPP16CalibParameters(//La función es un void porque no devuelve ningún valor 
const std::string& inputCalibMDPP16,
double (&a_MDPP16_Long)[8],
double (&b_MDPP16_Long)[8],
double (&a_MDPP16_Short)[8],
double (&b_MDPP16_Short)[8])//Entre paréntesis van los parámetros de la función
{
    // Open the input file stream
    std::ifstream inputFile(inputCalibMDPP16);//Comprueba si se puede abrir el archivo.
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file: " << inputCalibMDPP16 << std::endl;//En caso de no poder abrirse, cerr deja de correr el código
        return;
    }

    // Leer el archivo línea por línea
    std::string lineMDPP16;
    while (std::getline(inputFile, lineMDPP16)) {
        std::istringstream issMDPP16(lineMDPP16);//Leemos línea a línea

        double a_MDPP16_value_Long = 0.0, b_MDPP16_value_Long = 0.0, a_MDPP16_value_Short = 0.0, b_MDPP16_value_Short = 0.0;
        int channel=-1;//Hay que inicializar todas las variables, con un valor que no coincida con el real

        // Extraer los valores de la línea
        if (!(issMDPP16  >> channel >> a_MDPP16_value_Long >> b_MDPP16_value_Long >> a_MDPP16_value_Short >> b_MDPP16_value_Short)) {
            std::cerr << "Error parsing line: " << lineMDPP16 << std::endl;
            continue; // Saltar a la siguiente línea en caso de error
        }

        // Asignar valores a las matrices
        a_MDPP16_Long[channel]= a_MDPP16_value_Long;
        b_MDPP16_Long[channel]= b_MDPP16_value_Long;
        a_MDPP16_Short[channel]= a_MDPP16_value_Short;
        b_MDPP16_Short[channel] = b_MDPP16_value_Short;
    }
}

//Implementación TSpectrum
void Buscar_Picos(TH1F* histo, double noise_cut, std::vector<double> &peaks, const int nPeaks, std::string isotope,int GAGG5=0, int GAGG7=0) {
    double xmax = 15000;
    gErrorIgnoreLevel = kFatal;
    TSpectrum spectrum(nPeaks);

    // Búsqueda normal sobre el histograma original (se dibuja bien)
    histo->GetXaxis()->SetRangeUser(noise_cut, xmax);
    int nfound = spectrum.Search(histo, 20, "", 0.05);
    double *xpeaks = spectrum.GetPositionX();
    peaks.clear();
    for (int i = 0; i < nfound; i++) {
        peaks.push_back(xpeaks[i]);
    }

    // Segunda búsqueda sobre un clon para no tocar el original
    if (isotope=="207Bi"){
        if (nfound < nPeaks) {
            TH1F* hClone = (TH1F*)histo->Clone("hClone");
            hClone->SetDirectory(0); // desvincularlo de cualquier fichero/canvas
            hClone->GetXaxis()->SetRangeUser(7000, 11000);

            TSpectrum spectrum2(nPeaks-nfound);
            int nfound2 = spectrum2.Search(hClone, 20, "nodraw", 0.05); // clon sin dibujar
            double *xpeaks2 = spectrum2.GetPositionX();
            for (int i = 0; i < nfound2; i++) {
                bool isDuplicate = false;
                for (double existing : peaks) {
                    if (std::abs(xpeaks2[i] - existing) < 100) {
                        isDuplicate = true;
                        break;
                    }
                }
                if (!isDuplicate)
                    peaks.push_back(xpeaks2[i]);
            }
            delete hClone;
        }
    }
    if (isotope=="226Ra"){
        for (int i=2;i<nPeaks;i++){
            TH1F* hClone = (TH1F*)histo->Clone("hClone");
            hClone->SetDirectory(0);
            if (GAGG5==1){
                if (i==2){
                    hClone->GetXaxis()->SetRangeUser(4400, 5800);
                }
                if (i==3){
                    hClone->GetXaxis()->SetRangeUser(6500, 7800);
                }
                if (i==4){
                    hClone->GetXaxis()->SetRangeUser(7800, 8600);
                }
                if (i==5){
                    hClone->GetXaxis()->SetRangeUser(8600, 10000);
                }
                if (i==6){
                    hClone->GetXaxis()->SetRangeUser(10700, 13000);
                }
                if (i==7){
                    hClone->GetXaxis()->SetRangeUser(13000, 15000);
                }
                if (i==8){
                    hClone->GetXaxis()->SetRangeUser(15000, 17500);
                }
            }
            else if (GAGG7==1){
                if (i==2){
                    hClone->GetXaxis()->SetRangeUser(5400, 6700);
                }
                if (i==3){
                    hClone->GetXaxis()->SetRangeUser(7950, 9300);
                }
                if (i==4){
                    hClone->GetXaxis()->SetRangeUser(9300, 10000);
                }
                if (i==5){
                    hClone->GetXaxis()->SetRangeUser(10300, 12000);
                }
                if (i==6){
                    hClone->GetXaxis()->SetRangeUser(12800, 15000);
                }
                if (i==7){
                    hClone->GetXaxis()->SetRangeUser(16000, 17500);
                }
                if (i==8){
                    hClone->GetXaxis()->SetRangeUser(18000, 20000);
                }
            }
            else{
                if (i==2){
                    hClone->GetXaxis()->SetRangeUser(6000, 8000);
                }
                if (i==3){
                    hClone->GetXaxis()->SetRangeUser(9000, 11000);
                }
                if (i==4){
                    hClone->GetXaxis()->SetRangeUser(11000, 12300);
                }
                if (i==5){
                    hClone->GetXaxis()->SetRangeUser(12300, 15000);
                }
                if (i==6){
                    hClone->GetXaxis()->SetRangeUser(15600, 17600);
                }
                if (i==7){
                    hClone->GetXaxis()->SetRangeUser(19300, 21700);
                }
                if (i==8){
                    hClone->GetXaxis()->SetRangeUser(21800, 24000);
                }
            }
            TSpectrum spectrum2(1);
            int nfound2=spectrum2.Search(hClone,20,"nodraw",0.06);
            double *xpeaks2=spectrum2.GetPositionX();
            bool isDuplicate = false;
            if (nfound2>0){
                for (double existing : peaks) {
                    if (std::abs(xpeaks2[0] - existing) < 100) {
                        isDuplicate = true;
                        break;
                    }
                }
                if (!isDuplicate)
                    peaks.push_back(xpeaks2[0]);
                delete hClone;
            }
        }
    }
        
    std::sort(peaks.begin(), peaks.end());
    histo->GetXaxis()->SetRangeUser(0, xmax); // restaurar rango original para el dibujo
    gErrorIgnoreLevel = kInfo;
    }

//Ajustes
void Ajuste_Lineal(std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y, double &a, double &b, double &err_a, double &err_b) {
    int n = (int)x.size();
    if (n != (int)y.size() || n != (int)sigma_y.size()) {
        std::cerr << "Ajuste_Lineal: vectores de distinto tamaño x: " << n << " Energias: "<< y.size() << " sigma " << sigma_y.size()<< std::endl;
        return;
    }

    double S=0, Sx=0, Sy=0, Sxx=0, Sxy=0;
    for (int i = 0; i < n; i++) {
        if (sigma_y[i] <= 0) {
            std::cerr << "Ajuste_Lineal: sigma <= 0 en punto " << i << ", se omite" << std::endl;
            continue;
        }
        double w = 1.0 / (sigma_y[i] * sigma_y[i]); // peso = 1/σ²
        S   += w;
        Sx  += w * x[i];
        Sy  += w * y[i];
        Sxx += w * x[i] * x[i];
        Sxy += w * x[i] * y[i];
    }

    double Delta = S * Sxx - Sx * Sx;
    if (Delta == 0) {
        std::cerr << "Ajuste_Lineal: Delta=0, ajuste imposible" << std::endl;
        return;
    }

    a     = (S   * Sxy - Sx * Sy)  / Delta;
    b     = (Sxx * Sy  - Sx * Sxy) / Delta;
    err_a = TMath::Sqrt(S   / Delta);
    err_b = TMath::Sqrt(Sxx / Delta);
}

void Ajuste_Gauss(TH1F* Histo, double media_tentativa, double rango_izq, double rango_der, Double_t &area, Double_t &mean, Double_t &sigma){
    TF1 *gauss = new TF1("gauss", "gaus", media_tentativa-rango_izq, media_tentativa+rango_der); //Defino la funcion de ajuste y el rango a analizar y donde
    gauss->SetLineColor(kRed); //Color del ajuste


    //Inicializo parámetors
    gauss->SetParameter(0, Histo->GetMaximum());
    gauss->SetParameter(1, media_tentativa);
    gauss->SetParameter(2, 5.0);
    gauss->SetParLimits(2,0.01,1e9);


    Histo->Fit(gauss,"RQ+");// Ajusto mi histograma, el RQ+ hace que se vea laf
    
    //Cojo los parametros
    Double_t A = gauss->GetParameter(0);
    mean = gauss->GetParameter(1);
    sigma = gauss->GetParameter(2);   
    double bin_width =  Histo->GetBinWidth(1);
    area = A * TMath::Sqrt(2 * TMath::Pi()) * sigma/bin_width; //Formula Wikipedia
    
    delete gauss; //Diria que lo puedes borrar sin problema y la linea queda
}

void Ajuste_Gauss_Lineal(TH1F* Histo, double media_tentativa, double rango_izq, double rango_der, double &amplitud, double &mean, double &sigma, double &slope, double &y_inter, double &area, 
    double &err_amplitud, double &err_mean, double &err_sigma, double &err_slope, double &err_y_inter, double &err_area, double &cov_01 ){

        double xmin = media_tentativa-rango_izq;
        double xmax = media_tentativa+rango_der;

        // Defino la funcion a ajustar y los parametros con su orden, es una lineal mas una gaussiana
        TF1 *gauss_lineal = new TF1("gauss_lineal", "pol1(0) + gaus(2)", xmin, xmax);

        double yPico = Histo->GetBinContent(Histo->FindBin(media_tentativa));
        double yIzq = Histo->GetBinContent(Histo->FindBin(xmin));
        double yDer = Histo->GetBinContent(Histo->FindBin(xmax));
        
        // Parámetros estimados
        
        double slopeEst = (yDer-yIzq)/(xmax-xmin);
        double fondoEst = (yIzq+yDer)/2.0;
        //double fondoEst = (yIzq-slopeEst*xmin);
        double ampEst = yPico-fondoEst;//amplitud sin fondo
        double sigmaEst = (xmax-xmin)/6.0;

        //Doy unos parámetros iniciales
        gauss_lineal->SetLineColor(kBlack); //Color del ajuste
        gauss_lineal->SetParameter(0, fondoEst);   // ordenada en el origen
        gauss_lineal->SetParameter(1, slopeEst);                   // pendiente inicial
        gauss_lineal->SetParameter(2, ampEst);   // amplitud gaussiana
        gauss_lineal->SetParameter(3, media_tentativa);       // media
        gauss_lineal->SetParameter(4, sigmaEst);// sigma inicial (ajusta según tu resolución)                 
        TFitResultPtr result= Histo->Fit(gauss_lineal, "RQS+");


        y_inter = gauss_lineal->GetParameter(0);
        slope       = gauss_lineal->GetParameter(1);
        amplitud  = gauss_lineal->GetParameter(2);
        mean        = gauss_lineal->GetParameter(3);
        sigma       =abs(gauss_lineal->GetParameter(4));

        err_y_inter  = gauss_lineal->GetParError(0);
        err_slope    = gauss_lineal->GetParError(1);
        err_amplitud = gauss_lineal->GetParError(2);
        err_mean     = gauss_lineal->GetParError(3);
        err_sigma    = gauss_lineal->GetParError(4);

        double bin_width =  Histo->GetBinWidth(1);
        area = amplitud * TMath::Sqrt(2 * TMath::Pi()) * sigma/bin_width;//Habrá que dividirla entre el width de los canales
        err_area = (TMath::Sqrt(2 * TMath::Pi()) *TMath::Sqrt(TMath::Power(sigma * err_amplitud, 2) +TMath::Power(amplitud * err_sigma, 2)))/bin_width;

        TMatrixD covarianza=result->GetCovarianceMatrix();
        cov_01=covarianza(0,1);

        delete gauss_lineal;
}

void Ajuste_Dos_Gauss(TH1F* Histo, double media_tentativa1, double media_tentativa2, double rango_izq, double rango_der, 
    double &amplitud1, double &amplitud2, double &mean1, double &mean2, double &sigma1, double &sigma2, double &area1, double &area2, 
    double &err_amplitud1, double &err_amplitud2, double &err_mean1, double &err_mean2, double &err_sigma1, double &err_sigma2, double &err_area1, double &err_area2){
        
        double xmin = media_tentativa1-rango_izq;
        double xmax = media_tentativa1+rango_der;

        //Defino la funcion a ajustar y los parametros con su orden, es una lineal mas una gaussiana
        TF1 *gauss_gauss = new TF1("gauss_gauss", "gaus(0) + gaus(3)", xmin, xmax);

        double yPico1 = Histo->GetBinContent(Histo->FindBin(media_tentativa1));
        double yPico2 = Histo->GetBinContent(Histo->FindBin(media_tentativa2));

        //double fondoEst = (yIzq+yDer)/2.0;
        double ampEst1 = yPico1;//amplitud sin fondo
        double sigmaEst1 = (xmax-xmin)/6.0;
        double ampEst2 = yPico2;//amplitud sin fondo
        double sigmaEst2 = (xmax-xmin)/6.0;
        
        //Doy unos parámetros iniciales
        gauss_gauss->SetLineColor(kBlack);
        gauss_gauss->SetParameter(0, ampEst1);   // amplitud gaussiana
        gauss_gauss->SetParameter(1, media_tentativa1); // media
        gauss_gauss->SetParameter(2, sigmaEst1); //sigma inicial (ajusta según tu resolución)                 
        gauss_gauss->SetParameter(3, ampEst2); // amplitud gaussiana
        gauss_gauss->SetParameter(4, media_tentativa2); // media
        gauss_gauss->SetParameter(5, sigmaEst2); //sigma inicial (ajusta según tu resolución) 

        Histo->Fit(gauss_gauss, "RQ+");
        amplitud1  = gauss_gauss->GetParameter(0);
        mean1        = gauss_gauss->GetParameter(1);
        sigma1       =abs(gauss_gauss->GetParameter(2));
        amplitud2  = gauss_gauss->GetParameter(3);
        mean2        = gauss_gauss->GetParameter(4);
        sigma2       =abs(gauss_gauss->GetParameter(5));

        err_amplitud1 = gauss_gauss->GetParError(0);
        err_mean1     = gauss_gauss->GetParError(1);
        err_sigma1    = gauss_gauss->GetParError(2);
        err_amplitud2 = gauss_gauss->GetParError(3);
        err_mean2     = gauss_gauss->GetParError(4);
        err_sigma2    = gauss_gauss->GetParError(5);

        double bin_width =  Histo->GetBinWidth(1);
        area1 = amplitud1 * TMath::Sqrt(2 * TMath::Pi()) * sigma1/bin_width; //Habrá que dividirla entre el width de los canales
        area2 = amplitud2 * TMath::Sqrt(2 * TMath::Pi()) * sigma2/bin_width; //Habrá que dividirla entre el width de los canales
        err_area1 = (TMath::Sqrt(2 * TMath::Pi()) *TMath::Sqrt(TMath::Power(sigma1 * err_amplitud1, 2) +TMath::Power(amplitud1 * err_sigma1, 2)))/bin_width;
        err_area2 = (TMath::Sqrt(2 * TMath::Pi()) *TMath::Sqrt(TMath::Power(sigma2 * err_amplitud2, 2) +TMath::Power(amplitud2 * err_sigma2, 2)))/bin_width;

        delete gauss_gauss;
}

void Ajuste_Gauss_Gauss_Lineal(TH1F* Histo, double media_tentativa1, double media_tentativa2, double rango_izq, double rango_der, 
    double &amplitud1, double &amplitud2 ,double &mean1, double &mean2, double &sigma1, double &sigma2, double &slope, double &y_inter, double &area1, double &area2, 
    double &err_amplitud1, double &err_amplitud2, double &err_mean1, double &err_mean2, double &err_sigma1, double &err_sigma2, double &err_slope, double &err_y_inter, double &err_area1, double &err_area2){

        double xmin = media_tentativa1-rango_izq;
        double xmax = media_tentativa1+rango_der;

        //Defino la funcion a ajustar y los parametros con su orden, es una lineal mas una gaussiana
        TF1 *gauss_gauss_lineal = new TF1("gauss_gauss_lineal", "pol1(0) + gaus(2) + gaus(5)", xmin, xmax);

        double yPico1 = Histo->GetBinContent(Histo->FindBin(media_tentativa1));
        double yPico2 = Histo->GetBinContent(Histo->FindBin(media_tentativa2));

        double yIzq = Histo->GetBinContent(Histo->FindBin(xmin));
        double yDer = Histo->GetBinContent(Histo->FindBin(xmax));
        //double fondoEst = (yIzq+yDer)/2.0;
        double slopeEst=(yDer-yIzq)/(xmax-xmin);
        double fondoEst = (yIzq-slopeEst*xmin);

        double ampEst1= yPico1-fondoEst;//amplitud sin fondo
        double sigmaEst1=(xmax-xmin)/6.0;
        double ampEst2= yPico2-fondoEst;//amplitud sin fondo
        double sigmaEst2=(xmax-xmin)/6.0;

        //Doy unos parámetros iniciales
        gauss_gauss_lineal->SetLineColor(kBlack);
        gauss_gauss_lineal->SetParameter(0, fondoEst);   // ordenada en el origen
        gauss_gauss_lineal->SetParameter(1, slopeEst);                   // pendiente inicial
        gauss_gauss_lineal->SetParameter(2, ampEst1);   // amplitud gaussiana
        gauss_gauss_lineal->SetParameter(3, media_tentativa1);       // media
        gauss_gauss_lineal->SetParameter(4, sigmaEst1);// sigma inicial (ajusta según tu resolución)                 
        gauss_gauss_lineal->SetParameter(5, ampEst2);   // amplitud gaussiana
        gauss_gauss_lineal->SetParameter(6, media_tentativa2);       // media
        gauss_gauss_lineal->SetParameter(7, sigmaEst2);
        Histo->Fit(gauss_gauss_lineal, "RQ+");

        y_inter = gauss_gauss_lineal->GetParameter(0);
        slope       = gauss_gauss_lineal->GetParameter(1);
        amplitud1  = gauss_gauss_lineal->GetParameter(2);
        mean1        = gauss_gauss_lineal->GetParameter(3);
        sigma1       =abs(gauss_gauss_lineal->GetParameter(4));
        amplitud2  = gauss_gauss_lineal->GetParameter(5);
        mean2        = gauss_gauss_lineal->GetParameter(6);
        sigma2       =abs(gauss_gauss_lineal->GetParameter(7));

        err_y_inter  = gauss_gauss_lineal->GetParError(0);
        err_slope    = gauss_gauss_lineal->GetParError(1);
        err_amplitud1 = gauss_gauss_lineal->GetParError(2);
        err_mean1     = gauss_gauss_lineal->GetParError(3);
        err_sigma1    = gauss_gauss_lineal->GetParError(4);
        err_amplitud2 = gauss_gauss_lineal->GetParError(5);
        err_mean2     = gauss_gauss_lineal->GetParError(6);
        err_sigma2    = gauss_gauss_lineal->GetParError(7);

        double bin_width =  Histo->GetBinWidth(1);
        area1 = amplitud1 * TMath::Sqrt(2 * TMath::Pi()) * sigma1/bin_width; //Habrá que dividirla entre el width de los canales
        area2 = amplitud2 * TMath::Sqrt(2 * TMath::Pi()) * sigma2/bin_width; //Habrá que dividirla entre el width de los canales

        err_area1 = (TMath::Sqrt(2 * TMath::Pi()) *TMath::Sqrt(TMath::Power(sigma1 * err_amplitud1, 2) +TMath::Power(amplitud1 * err_sigma1, 2)))/bin_width;
        err_area2 = (TMath::Sqrt(2 * TMath::Pi()) *TMath::Sqrt(TMath::Power(sigma2 * err_amplitud2, 2) +TMath::Power(amplitud2 * err_sigma2, 2)))/bin_width;

        delete gauss_gauss_lineal;
    }

void Efecto_Compton(double E, double &E_Borde){
    E_Borde=E*2*E/(2*E+511);
}

void BackScatter(double E, double &E_BS){
    E_BS=E/(1+2*E/511);
}
    

void Save_Histo(TObject* Histo, TDirectory* Directory_name) {
    Directory_name->cd();
    Histo->Write();
    delete Histo;
}


std::vector <double> Select_Calibration_Values(const std::string &isotope){
    if(isotope == "60Co"){
        return {1173.2, 1332.5};
    }
    else if(isotope == "152Eu"){
        return {344.278, 778.9, 964.057, 1099.82, 1408.013};//ISOLDE
        //return {244.78, 344.278,778.9, 964.057, 1099.82, 1408.013};//IEM
    }
    else if(isotope == "22Na"){
        return {511,1274.537};        
    }
    else if(isotope == "137Cs"){
        return {661.657};
    }
    else if(isotope == "207Bi"){
        return {569.70, 1063.7, 1770.2};
    }
    else if(isotope == "226Ra"){
        return {351,609,774.6,1123.8,1390,1760,2204,2447};
    }
    else if(isotope == "BKG"){
        return {609.3, 1460.8, 1764.5, 2614.5};
    }
    else{
        throw std::runtime_error("Isotope not found: " + isotope);
    }
}

std::vector <double> Select_Calibration_Ranges(const std::string &isotope){
    if(isotope == "60Co"){
        return {350, 300};
    }
    else if(isotope == "152Eu"){
        //return {400, 500, 450, 550, 700};//ISOLDE
        return {276, 336, 500, 300, 400, 600};//IEM
    }
    else if(isotope == "22Na"){
        return {400, 700};       
    }
    else if(isotope == "137Cs"){
        return {70};  //Cambiar por buenos
    }
    else if(isotope == "207Bi"){
        return {400, 600, 580}; //Cambiar por buenos 
    }
    else if(isotope == "226Ra"){
        return {300,900,700,700,800,1500,2000,1500};
    }
    else if(isotope == "BKG"){
        return {200, 300, 500, 700}; //Cambiar por buenos 
    }
    else{
        throw std::runtime_error("Isotope not found: " + isotope);
    }
}

std::vector <double> Select_Fit_Ranges(const std::string &isotope){
    if(isotope == "60Co"){
        return {90, 90};
    }
    else if(isotope == "152Eu"){
        return {50, 80, 60, 80, 110};//ISOLDE
        //return {45, 50, 80, 60, 80, 120};//IEM
    }
    else if(isotope == "22Na"){
        return {100, 100};       
    }
    else if(isotope == "137Cs"){
        return {110};  //Cambiar por buenos
    }
    else if(isotope == "207Bi"){
        return {110, 120, 130}; //Cambiar por buenos 
    }
    else if(isotope == "226Ra"){
        return {35,120,80,80,80,140,110,107};
    }
    else if(isotope == "BKG"){
        return {200, 200, 200, 200}; //Cambiar por buenos 
    }
    else{
        throw std::runtime_error("Isotope not found: " + isotope);
    }
}




void add_zero_line(TH2F* Histo){
    TLine *zero_line = new TLine(0, 0, 2000, 0);
    zero_line->SetLineColor(kRed);
    zero_line->SetLineWidth(2);
    Histo->GetListOfFunctions()->Add(zero_line);
}

void Draw_Calib_lines(TH1F* Histo, const double Calib_value){
    double x_Points[2] = {Calib_value,Calib_value};
    double y_Points[2] = {0.0,1000000.0};
    TPolyLine* Calib_line = new TPolyLine(2,x_Points,y_Points);
    Calib_line ->SetLineColor(kRed);
    Calib_line ->SetLineWidth(2);
    Calib_line ->SetLineStyle(2);
    Histo->GetListOfFunctions()->Add(Calib_line); 
}

void Draw_Compton_lines(TH1F* Histo, const double Calib_value){
    double E_Compton;
    Efecto_Compton(Calib_value, E_Compton);
    double x_Points[2] = {E_Compton,E_Compton};
    double y_Points[2] = {0.0,1000000.0};
    TPolyLine* Compton_line = new TPolyLine(2,x_Points,y_Points);
    Compton_line ->SetLineColor(kBlue);
    Compton_line ->SetLineWidth(2);
    Compton_line ->SetLineStyle(2);
    Histo->GetListOfFunctions()->Add(Compton_line); 
}

void Draw_BackScatter_lines(TH1F* Histo, const double Calib_value){
    double E_BS;
    BackScatter(Calib_value, E_BS);
    double x_Points[2] = {E_BS,E_BS};
    double y_Points[2] = {0.0,1000000.0};
    TPolyLine* BackScatter_line = new TPolyLine(2,x_Points,y_Points);
    BackScatter_line ->SetLineColor(kGreen);
    BackScatter_line ->SetLineWidth(2);
    BackScatter_line ->SetLineStyle(2);
    Histo->GetListOfFunctions()->Add(BackScatter_line); 
}

void Draw_Check_Calib_lines(TH2F* Histo, const double Calib_value){
    double x_Points[2] = {0.0,1000000.0};
    double y_Points[2] = {Calib_value,Calib_value};
    TPolyLine* Calib_line = new TPolyLine(2,x_Points,y_Points);
    Calib_line ->SetLineColor(kRed);
    Calib_line ->SetLineWidth(2);
    Calib_line ->SetLineStyle(2);
    Histo->GetListOfFunctions()->Add(Calib_line); 
}

void Save_GAGGs(TFile* outputFile,
    TH1F* GAGGS_Raw_hist_Total_Long[8],
    TH1F* GAGGS_Raw_hist_Total_Short[8],
    TH1F* GAGGS_Calib_hist_Total_Long[8],
    TH1F* GAGGS_Calib_hist_Total_Short[8],
    TH1F* GAGGS_Calib_Total_Long,
    TH1F* GAGGS_Calib_Total_Short,
    TH2F* GAGGS_Residuo_Long[8],
    TH2F* GAGGS_Residuo_Short[8],
    TH2F* GAGGS_FWHM_Long[8],
    TH2F* GAGGS_FWHM_Short[8],
    TH1F* GAGGS_Coinci_Two_hist[8],
    TH1F* GAGGS_Coinci_Two_Full_hist,
    TH1F* GAGGS_Coinci_with_511_hist[8],
    TH1F* GAGGS_Coinci_with_511_Full_hist,
    TH2F* GAGGS_Multi_2_hist,
    TH2F* GAGGS_Coinci_2_hist,
    TH2F* GAGGS_Short_vs_Long_hist[8],
    TH2F* GAGGs_Check_Calib_hist,
    TH1F* GAGGS_Add_back_hist_Total_Long[8], 
    TH1F* GAGGS_n_Coinci_hist, 
    TH2F* GAGGS_3hit_PairSum_vs_Third_hist)
    {
    
    gROOT->SetBatch(kTRUE);
    
    TDirectory* GAGGsDir =  outputFile->mkdir("GAGGs");
    TDirectory* GAGG_Raw_Long_Dir = GAGGsDir->mkdir("Raw Integration Long");
    TDirectory* GAGG_Raw_Short_Dir = GAGGsDir->mkdir("Raw Integration Short");
    TDirectory* GAGG_Calib_Long_Dir = GAGGsDir->mkdir("Calibrated Integration Long");
    TDirectory* GAGG_Calib_Short_Dir = GAGGsDir->mkdir("Calibrated Integration Short");
    TDirectory* Res_Long_Dir = GAGGsDir->mkdir("Residuos Long");
    TDirectory* Res_Short_Dir = GAGGsDir->mkdir("Residuos Short");
    TDirectory* FWHM_Long_Dir = GAGGsDir->mkdir("FWHM Long");
    TDirectory* FWHM_Short_Dir = GAGGsDir->mkdir("FWHM Short");
    TDirectory* Short_vs_Long_Dir = GAGGsDir->mkdir("Short vs Long");
    TDirectory* Coinci_GAGGs_Dir = GAGGsDir->mkdir("Coinci Two GAGGS");
    TDirectory* Coinci_GAGGs_511_Dir = GAGGsDir->mkdir("Coinci Two GAGGS one is 511");
    TDirectory* Total_Dir = GAGGsDir->mkdir("Total");
    TDirectory* Add_back_Dir = GAGGsDir->mkdir("Add back");
    for(int i = 0; i < 8; i++){
        Save_Histo(GAGGS_Raw_hist_Total_Long[i],GAGG_Raw_Long_Dir);
        Save_Histo(GAGGS_Raw_hist_Total_Short[i],GAGG_Raw_Short_Dir);

        Save_Histo(GAGGS_Calib_hist_Total_Long[i],GAGG_Calib_Long_Dir);
        Save_Histo(GAGGS_Calib_hist_Total_Short[i],GAGG_Calib_Short_Dir);

        Save_Histo(GAGGS_Residuo_Long[i],Res_Long_Dir);
        Save_Histo(GAGGS_Residuo_Short[i],Res_Short_Dir);

        Save_Histo(GAGGS_FWHM_Long[i],FWHM_Long_Dir);
        Save_Histo(GAGGS_FWHM_Short[i],FWHM_Short_Dir);

        Save_Histo(GAGGS_Coinci_Two_hist[i],Coinci_GAGGs_Dir);
        Save_Histo(GAGGS_Coinci_with_511_hist[i],Coinci_GAGGs_511_Dir);

        Save_Histo(GAGGS_Short_vs_Long_hist[i],Short_vs_Long_Dir);

        Save_Histo(GAGGS_Add_back_hist_Total_Long[i], Add_back_Dir);
    }

    Save_Histo(GAGGS_Calib_Total_Long,Total_Dir);
    Save_Histo(GAGGS_Calib_Total_Short,Total_Dir);
    Save_Histo(GAGGS_Coinci_Two_Full_hist,Total_Dir);
    Save_Histo(GAGGS_Coinci_with_511_Full_hist,Total_Dir);
    Save_Histo(GAGGS_Multi_2_hist,Total_Dir);
    Save_Histo(GAGGS_Coinci_2_hist,Total_Dir);
    Save_Histo(GAGGs_Check_Calib_hist,Total_Dir);
    Save_Histo(GAGGS_n_Coinci_hist,Total_Dir);
    Save_Histo(GAGGS_3hit_PairSum_vs_Third_hist,Total_Dir);
    std::cout << " GAGGs guardados " << "\n";
}

void Save_BKG(TFile* outputFile, 
    TH1F* GAGGS_BKG_hist_Long[8],
    TH1F* GAGGS_BKG_hist_Short[8],
    TH1F* GAGGS_BKG_Total_Long,
    TH1F* GAGGS_BKG_Total_Short
    ){
    gROOT->SetBatch(kTRUE);
    TDirectory* BKG_Dir =  outputFile->mkdir("BKG");
    TDirectory* GAGG_Calib_Long_Dir = BKG_Dir->mkdir("Calibrated Integration Long");
    TDirectory* GAGG_Calib_Short_Dir = BKG_Dir->mkdir("Calibrated Integration Short");
    TDirectory* Total_Dir = BKG_Dir->mkdir("Total");

    for(int i = 0; i < 8; i ++){
        Save_Histo(GAGGS_BKG_hist_Long[i],GAGG_Calib_Long_Dir);
        Save_Histo(GAGGS_BKG_hist_Short[i],GAGG_Calib_Short_Dir);
    }
    Save_Histo(GAGGS_BKG_Total_Long,Total_Dir);
    Save_Histo(GAGGS_BKG_Total_Short,Total_Dir);
    std::cout << " BKG guardado " << "\n";
}

int General(){ 
    return 0;
}

