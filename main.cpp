#include "IS_Model.h"

using namespace std;

int main(){
  // simulation definitions
  double simdefs[6] = {0,1,30,720,2,2}; //caso 0 - simulates coupled model,1 - saves all files(1 ponto por hora), 30 dias, 720 pontos
  //lnv = 2 - contact with lymph vessels given by function
  //bv = 2 - contact with blood vessels given by function
  IS_Model* model = new IS_Model(simdefs); //cria modelo
  //model parameters of interest

  model->solve();//solve
  return 0;
}
