#ifndef _IS_Model_H_
#define _IS_Model_H_

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

//condições iniciais do pedaço de tecido
//each 10 space equals 1 cm (1 × 10^−2 m)
const int    Xspace   = 10.0; // 10 mm = 1 cm
const int    Yspace   = 10.0; // 10 mm = 1 cm
const int    Zspace   = 10.0; // 10 mm = 1 cm
const int    SPACE    = Xspace * Yspace * Zspace; // 1 cm^3
const int    IC_SPACE = 0.4*(Xspace * Yspace * Zspace); //initial condition
const int    source   = 100*pow(10,0);
const double SCALE    = pow(10,-3);
const double MOL      = 6.02*pow(10,23);
const int    buffer   = 2;

class IS_Model{

  private:

    double A [buffer][Xspace][Yspace][Zspace];  //S. aureus Bacteria
    double MR [buffer][Xspace][Yspace][Zspace]; //Resting Macrophages
    double MA [buffer][Xspace][Yspace][Zspace];  //Activated Macrophages
    double F [buffer][Xspace][Yspace][Zspace];  //Antigens

    int simCase;
    int days;
    int points;
    double deltaT;  //intervalo de tempo
    double iterPerDay; //numero de iterações por dia
    double deltaX, deltaY, deltaZ;  //intervalo de espacamento (volume)

    int lnv;  //volume do linfonodo
    int bv; //blood vessels
    double tol;
    double source_mr;
    double migration_ma; //migração macrófagos
    double migration_f;   //migração antigens

    //Initial values of the coupled model (parameters)
    double m0;
    double a0;
    double th0;
    double b0;
    double p0;
    double f0;
    double t_estrela;
    double b_estrela;
    double p_estrela;
    double f_estrela;
    double m_estrela;
    //
    double MA_T, MA_L, MR_T; //concentration of active macrophages in the tissue,active macrophages in the lympho node,resting macrophages in the tissue
    double Th, B, P; //T-helper lymphocytes, B-lymphocytes, Plasma
    double F_T, F_L, A_T; //Antigens in the tissue,Antigens in the lympho node,S. aureus in the tissue
    //Diffusion coefficients.
    double d_a; 
    double d_mr;
    double d_ma; 
    double d_f;
    //Replication, decay, activation, and phagocytosis rates.
    double beta_A;
    double k_A;
    double m_A;
    double m_Mr;
    double m_Ma;
    double gamma_ma;
    double lambda_mr;
    double lambda_ma;
    double lambda_afmr;
    double lambda_afma;
    //Other coefficients used in the coupled model.
    double b_th;
    double b_p;
    double b_pb;
    double b_pp;
    double ro_t;
    double ro_b;
    double ro_p;
    double ro_f;
    double alpha_Ma;
    double alpha_t;
    double alpha_b;
    double alpha_p;
    double alpha_f;
    double alpha_mr;
    

    int saveFiles;
    char *dir;
    FILE* datamatlabA;
    FILE* datamatlabMr;
    FILE* datamatlabMa;
    FILE* datamatlabT;
    FILE* datamatlabB;
    FILE* datamatlabP;
    FILE* datamatlabF;
    FILE* datamatlabL;

    std::string Header();
    std::string Footer(long int t);
    int checkFile(FILE* theFile);
    int calcIntegral(double vec[][Xspace][Yspace][Zspace], double *V);
    int calcIntegral_lv(double vec[][Xspace][Yspace][Zspace], double *V);
    int calcIntegral_bv(double vec[][Xspace][Yspace][Zspace], double *V);
    void initialize();
    void update(double vec[][Xspace][Yspace][Zspace]);
    double laplacian(double vec[][Xspace][Yspace][Zspace], int x, int y, int z);
    int is_bvase(int x, int y, int z);
    int is_lnvase(int x, int y, int z);

  public:
    IS_Model(double simdefs[]);
    //IS_Model(int sCase, int sFile);
    //~IS_Model();
    void setSaveFiles(int sf);
    void setSimulationCase(int sc);
    int solve();

};

#endif
