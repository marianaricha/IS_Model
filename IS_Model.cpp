#include "IS_Model.h"

/******************************************************************************
 * 
 * IS_Model - Imunne System Model 
 * 
 * This software simmulates the behavior of the main cells of the immune
 * system both innate and acquired during antigen presentation.
 * 
 * Based on the works from Pigozzo (2011) and Marchuk (1997).
 * 
 * Date modified :30/08/2019 
 *                06/12/2013.
 * 
 * Recquires: 'IS_Model.h'.
 * 
 * Outputs : 
 *  
 *          '*.csv' files for time steps chosen for each PDE in the model.
 *          'L.dat' containing the averages of cells in the tissue.
 *          'T.dat', 'B.dat', 'P.dat' containing these cells concentrations 
 *                                    over time.
 * 
 * Use-me : 
 *  
 *          1. Instantiate an object of IS_Model class:
 * 
 *             IS_Model* model = new IS_Model();
 * 
 *          2. Call solve() method:
 *  
 *             model->solve();
 * 
 * 
 * Obs: Considers logistic growth of bacteria.
 * 
 *  
 ******************************************************************************/

using namespace std;

/**
 * Greetings from the simulator
 */
std::string IS_Model::Header(){
    int size = 40;
    std::string returnstring;
    std::ostringstream sstream;
    sstream << a0;
    std::string a0str = sstream.str();

    for(int i = 0; i < size; i++) returnstring += "*";
    returnstring += "\n* Begin of simulation *\n";
    returnstring += "*\n* Initial bacteria = "+a0str+".\n*\n";
    for(int i = 0; i < size; i++) returnstring += "*";
    returnstring += "\n";
    return returnstring;
}

/**
 * Goodbye message from the simulator
 */
std::string IS_Model::Footer(long int t){
    std::string returnstring;
    int size = 40;
    std::ostringstream sstream1, sstream2;
    sstream1 << t;
    sstream2 << days;
    std::string tstr = sstream1.str();
    std::string daysstr = sstream2.str();

    returnstring += "The end. \n ...of simulation!";
    returnstring += "\n"+tstr+" time steps for "+daysstr+" days.\n";
    for(int i = 0; i < size; i++) returnstring += "*";
    sstream1 << A_T;
    std::string atstr = sstream1.str();

    returnstring += "\nBacteria in the end : "+atstr+" \n";
    for(int i = 0; i < size; i++) returnstring += "*";
    returnstring += "\n";
    return returnstring;
}

/**
 * Setters
 */
void IS_Model::setSaveFiles(int sf){
    this->saveFiles = sf;
}
void IS_Model::setSimulationCase(int sc){
    this->simCase = sc;
}

/**
* Constructor set parameters
*/
IS_Model::IS_Model(double simdefs[]){

  /** 
   * 0 - simulates coupled model, 
   * 1 - simulates only antigen diffusion,
   * 2 - simulates only innate response,
   * 3 - complete model without diffusion.
   */  
  //simCase   = 0;  
   this->simCase   = simdefs[0];
  /** 
   * 0 - saves only edos files, 
   * 1 - saves all files.
   */
  //saveFiles = 1;
  this->saveFiles = simdefs[1];
  /**
   * number of days simulated
   */
  //days      = 30;
  this->days      = simdefs[2];
  /**
   * number of files saved
   */
  //points    = 720;
  this->points    = simdefs[3];
  /**
   * 0 - contact with lymph vessels only on one border, 
   * 1 - homogeneous contact with lymph vessels.
   * 2 - contact with lymph vessels given by function.
   */
  //lnv       = 2;
  this->lnv       = simdefs[4];
  /**
   * 0 - contact with blood vessels only on one border,
   * 1 - homogeneous contact with blood vessels,  
   * 2 - contact with blood vessels given by function.
   */
  //bv        = 2;
  this->bv        = simdefs[5];
  /**
   * output directory
   */
  this->dir       = (char *) "output/";
}

/**
* Set conditions and parameter values
*/
void IS_Model::initialize(){

  //TODO: Dúvida a respeito dos valores.
  /**
   * each (5/(pow(10,6)) = 0.0000002 days or 0,01728 secs
   */
  deltaT     = 0.001;
  /**
   * each 1000000 iterations represents 1 day
   */
  iterPerDay = 10000;
  /**
   * each deltaX represents 100 micrometers ((1 × 10^-6 m)), a cell
   * has 1000 cubic micrometers, each discretized space has 1.000.000 cubic micrometers,
   * leading to 1000 cells for each space
   */
  deltaX     = 0.1; //mm
  deltaY     = 0.1;
  deltaZ     = 0.1;
  
  tol        = pow(10,-6); //tolerância para quantidade de bacterias (prox 0)

  //Table 1: Initial values of the coupled model.
  a0         = 2.0;//1.7*pow(10,2);
  m0         = 0.0; //MA0
  f0         = 0.;//1.0*pow(10,1);//*MOL;
  th0        = 0.0;
  b0         = 0.0;
  p0         = 0.0;
  t_estrela  = 8.4*pow(10,-3);//*MOL;
  b_estrela  = 8.4*pow(10,-4);//*MOL;
  p_estrela  = 8.4*pow(10,-6);//*MOL;
  f_estrela  = 0.;//9.5*pow(10,-6);//*MOL;
  m_estrela  = 4.0;//2.3*pow(10,2);//*MOL //MR0

  
  MA_T    = 0.0; //concentration of active macrophages in the tissue
  MA_L    = 0.0; //active macrophages in the lympho node
  MR_T    = m_estrela; //resting macrophages in the tissue
  Th      = th0; //T-helper lymphocytes
  B       = b0;  //B-lymphocytes
  P       = p0;  //Plasma cells
  F_T     = f0;  //Antigens in the tissue
  F_L     = f0;  //f_estrela; Antigens in the lympho node
  A_T     = a0;  //bacteria S. aureus in the tissue

  //Table 2:diffusion coefficients
  d_a        = 0.00037;       //antigen diffusion (Haessler)
  d_mr       = 0.0432;        //resting macrophage diffusion (estimated)
  d_ma       = 0.3;           //active macrophage diffusion (estimated)
  d_f        = 0.016;         //antibody diffusion (estimated)

  //Table 3: Replication, decay, activation, and phagocytosis rates.
  beta_A     = 2.0;            //bacteria replication
  k_A        = 50.0;           //bacteria maximum capacity
  m_A        = 0.1;            //bacteria natural decay
  m_Mr       = 0.033;          //resting macrophage natural decay
  m_Ma       = 0.07;           //activated macrophage natural decay
  gamma_ma   = 8.30*pow(10,-2); //macrophage activation
  lambda_mr  = 5.98*pow(10,-3); //resting macrophage fagocitosis rate
  lambda_ma  = 5.98*pow(10,-2); //activated macrophage fagocitosis rate
  lambda_afmr= 1.66*pow(10,-3); //resting macrophage fagocitosis rate for opsonized antigen
  lambda_afma= 7.14*pow(10,-2); //activated macrophage fagocitosis rate for opsonized antigen
  
  //Table 4: Other coefficients used in the coupled model.
  alpha_Ma   = 0.001;//*pow(10,1);  //migration (estimated)
  alpha_t    = 0.01;           //th2 natural decay
  alpha_b    = 1.0;            //b natural decay
  alpha_p    = 5.0;            //plasma decay
  alpha_f    = 0.43;           //antibody migration
  alpha_mr   = 4.0;            //resting macrophage source coefficient
  b_th       = 1.7*pow(10,-2); //th2 stimuli
  b_p        = 1.*pow(10,5);   //th2 expenditure to stimulate b
  b_pb       = 6.02*pow(10,3); //b stimuli
  b_pp       = 2.3*pow(10,6);  //b stimuli while describes plasma cell
  ro_t       = 2.0;            //th2 descendents
  ro_b       = 16.0;           //b descendents
  ro_p       = 3.0;            //p descendents
  ro_f       = 5.1*pow(10,4);  //antibody release
  //VLN simplificado

  /**
   * Initial Conditions
   */
  for(int x = 0; x < Xspace; x++) {
    for(int y = 0; y < Yspace; y++) {
      for(int z = 0; z < Zspace; z++) {
        if (simCase == 3){ //no diffusion
          A[0][x][y][z] = a0;
	      }else{
          //bacteria only in the center of the cubic domain
          if ((x > (0.2*Xspace)&&( x < (0.7*Xspace)))
            && (y > (0.2*Yspace)&&( y < (0.7*Yspace)))
            && (z > (0.2*Zspace)&&( z < (0.7*Zspace)))) {
		          A[0][x][y][z] = a0;///IC_SPACE;
		      } else {
		           A[0][x][y][z] = 0.0;
          }	  
	      }
	      MR[0][x][y][z]  = m_estrela;
        MA[0][x][y][z]  = 0.0;
        F[0][x][y][z]   = f0;//SPACE;
      }
    }
  }
}

/**
 * Updates the current results to position 1
 * and sets the previous results to zero
 */
void IS_Model::update(double vec[][Xspace][Yspace][Zspace]){
  for(int x = 0; x < Xspace; x++) {
    for(int y = 0; y < Yspace; y++) {
      for(int z = 0; z < Zspace; z++) {
        vec[0][x][y][z] = vec[1][x][y][z];
        /**
         * In case it is necessary to keep the previous
	       * just comment the line below
	       */
	      //vec[1][x][y][z] = 0.;
      }
    }
  }
}

/**
 * Calculates the laplacian for given value and position
 */
double IS_Model::laplacian(double vec[][Xspace][Yspace][Zspace], int x, int y, int z){
	double resX = 0, resY = 0, resZ = 0;
	// same boundary condition to every equation
	if(x == 0) {
		resX = (vec[0][x+1][y][z] - vec[0][x][y][z])/(deltaX*deltaX);
	} else if(x == Xspace-1 ) {
		resX = (vec[0][x-1][y][z] - vec[0][x][y][z])/(deltaX*deltaX);
	} else {//dentro do dominio mas fora da extremidade
		resX = (vec[0][x+1][y][z] -2 * vec[0][x][y][z] + vec[0][x-1][y][z])/(deltaX*deltaX);
	}
	if(y == 0) {
		resY = (vec[0][x][y+1][z] - vec[0][x][y][z])/(deltaY*deltaY);
	} else if( y == Yspace-1) {
		resY = (vec[0][x][y-1][z] - vec[0][x][y][z])/(deltaY*deltaY);
	} else {
		resY = (vec[0][x][y+1][z] -2 * vec[0][x][y][z] + vec[0][x][y-1][z])/(deltaY*deltaY);
	}
	if(z == 0) {
		resZ = (vec[0][x][y][z+1] - vec[0][x][y][z])/(deltaZ*deltaZ);
	} else if( z == Zspace-1) {
		resZ = (vec[0][x][y][z-1] - vec[0][x][y][z])/(deltaZ*deltaZ);
	} else {
		resZ = (vec[0][x][y][z+1] -2 * vec[0][x][y][z] + vec[0][x][y][z-1])/(deltaZ*deltaZ);
	}
	return resX+resY+resZ;
}

/**
 * Tests if the file could be created and exits the simulation if
 * there was an error
 */
int IS_Model::checkFile(FILE* theFile){
  if (theFile==NULL){
    cout << "Error opening file!!!\n Make sure the path is correct! \n";
    return 1;
  }
  return 0;
}

/**
 * Calculates integrals of cells in the tissue and return the value as
 * a pointer 
 */
int IS_Model::calcIntegral(double vec[][Xspace][Yspace][Zspace], double *V){

  for(int x = 0; x < Xspace; x++) {
	for(int y = 0; y < Yspace; y++) {
	  for(int z = 0; z < Zspace; z++) {
	    if (vec[0][x][y][z]>0.0) *V += vec[0][x][y][z];
	  }
	}
  }
  if (*V > 0.0) *V = (*V/(SPACE)); else *V = 0.0;
  return 0;
}

/**
 * For activated macrophages the integral is calculated considering only 
 * the cells in contact with lymph vessels
 */
int IS_Model::calcIntegral_lv(double vec[][Xspace][Yspace][Zspace], double *V){

  for(int x = 0; x < Xspace; x++) {
    for(int y = 0; y < Yspace; y++) {
      for(int z = 0; z < Zspace; z++) {
	      if (((lnv==0)&&(x==0))||(lnv==1)||((lnv==2)&&(is_lnvase(x,y,z)))){
	        if (vec[0][x][y][z]>0.0) *V += vec[0][x][y][z];
	      }
      }
    }
  }
  if (*V > 0.0) *V = (*V/(SPACE)); else *V = 0.0;
    return 0;
}

/**
 * for antibodies consider only cells in contact with blood vessels
 */
int IS_Model::calcIntegral_bv(double vec[][Xspace][Yspace][Zspace], double *V){

  for(int x = 0; x < Xspace; x++) {
    for(int y = 0; y < Yspace; y++) {
      for(int z = 0; z < Zspace; z++) {
	      if (((bv==0)&&(x==0))||(bv==1)||((bv==2)&&(is_bvase(x,y,z)))){
	       if (vec[0][x][y][z]>0.0) *V += vec[0][x][y][z];
      	}
      }
    }
  }
  if (*V > 0.0) *V = (*V/(SPACE)); else *V = 0.0;
    return 0;
}

/**
 * returns 1 if the point is a blood vase and 0 if it is not
 */
int IS_Model::is_bvase(int x, int y, int z){
  //if(((x >= 0)&&(x <= 1))||((x>=4)&&(x<=5))||((x>=8)&&(x<=9)))
  //  if(((z >= 0)&&(z <= 1))||((z>=4)&&(z<=5))||((z>=8)&&(z<=9)))
 if(((x >= 0)&&(x <= 1))||((x>=8)&&(x<=9)))
    if(((z >= 0)&&(z <= 1))||((z>=8)&&(z<=9)))  
      return 1;
  
  return 0;     
}

/**
 * returns 1 if the point is a lymph vase and 0 if it is not
 */
int IS_Model::is_lnvase(int x, int y, int z){
  if(((x >= 2)&&(x <= 3))||((x>=6)&&(x<=7)))
    if(((z >= 0)&&(z <= 1))||((z>=4)&&(z<=5)))
      return 1;
  
  return 0;     
}


/******************************************************************************
* Solve model equations
*******************************************************************************/
int IS_Model::solve(){

  int i       = 0;
  long int t  = 0;

  //set initial conditions
  initialize();

  //print program header
  cout << Header();

  // memory allocation for filename
  char *fileName = (char *)malloc(20*sizeof(char));
  
  sprintf(fileName, "%s%s", dir, "L.dat");
  datamatlabL = fopen(fileName, "w");
    
  //check valid dir
  if (checkFile(datamatlabL)) return 1;

  sprintf(fileName, "%s%s", dir, "T.dat");
  datamatlabT = fopen(fileName, "w");

  sprintf(fileName, "%s%s", dir, "B.dat");
  datamatlabB = fopen(fileName, "w");

  sprintf(fileName, "%s%s", dir, "P.dat");
  datamatlabP = fopen(fileName, "w");

  /**
   * begin time loop
   */
  do{

    if (t == 0) cout << "Calculating...\n"; //????

    int value = ((int)iterPerDay*days)/points; //fora do if?

    if(t%value == 0) {

  	  cout << "Saving files : iteration ..."<< t << "\n";

      if (saveFiles){//check before saving all files
        sprintf(fileName,"%sA_%ld.csv", dir,t);
        datamatlabA = fopen(fileName, "w");
   	    sprintf(fileName,"%sMr_%ld.csv",dir,t);
   	    datamatlabMr = fopen(fileName, "w");
   	    sprintf(fileName,"%sMa_%ld.csv",dir,t);
   	    datamatlabMa = fopen(fileName, "w");
   	    sprintf(fileName,"%sF_%ld.csv",dir,t);
   	    datamatlabF = fopen(fileName, "w");
   	    fprintf(datamatlabT, "%ld %.2E \n", t, Th);
   	    fprintf(datamatlabB, "%ld %.2E \n", t, B);
   	    fprintf(datamatlabP, "%ld %.2E \n", t, P);

        fprintf(datamatlabL, "%ld %.2E %.2E %.2E %.2E %.2E %.2E\n", t, MA_T, F_T, MA_L, F_L, A_T, MR_T);
   	    for(int x = 0; x < Xspace; x++) {
  	      for(int y = 0; y < Yspace; y++) {
            for(int z = 0; z < Zspace; z++) {
              if( (x+1 == Xspace && y+1 == Yspace && z+1==Zspace) ) {
                fprintf(datamatlabA, "%d %d %d %E", x, y, z, A[0][x][y][z]);
     		        fprintf(datamatlabMr, "%d %d %d %E", x, y, z, MR[0][x][y][z]);
		            fprintf(datamatlabMa, "%d %d %d %E", x, y, z, MA[0][x][y][z]);
       		      fprintf(datamatlabF, "%d %d %d %E", x, y, z, F[0][x][y][z]);
              } else {
		            fprintf(datamatlabA, "%d %d %d %E\n", x, y, z, A[0][x][y][z]);
        	      fprintf(datamatlabMr, "%d %d %d %E\n", x, y, z, MR[0][x][y][z]);
                fprintf(datamatlabMa, "%d %d %d %E\n", x, y, z, MA[0][x][y][z]);
        	      fprintf(datamatlabF, "%d %d %d %E\n", x, y, z, F[0][x][y][z]);
	            }
	          }
    	    }
    	  }
        fclose(datamatlabA);
        fclose(datamatlabMr);
        fclose(datamatlabMa);
        fclose(datamatlabF);
      }else{
        fprintf(datamatlabT, "%ld %.2E \n", t, Th);
   	    fprintf(datamatlabB, "%ld %.2E \n", t, B);
   	    fprintf(datamatlabP, "%ld %.2E \n", t, P);
	      fprintf(datamatlabL, "%ld %.2E %.2E %.2E %.2E %.2E %.2E\n", t, MA_T, F_T, MA_L, F_L, A_T, MR_T);	
	
        //fprintf(datamatlabL, "%ld %.2E %.2E %.2E %.2E %.2E %.2E \n", t,
                //MA[0][0][0][0], F[0][0][0][0], MA_L, F_L, A[0][0][0][0],
                //MR[0][0][0][0]);
      }
    }
  
    //integral
    //cout << "Solve integrals. ";
    if (t > 0 && simCase!=3){ //with diffusion (0,1 e 2)
      MA_T = MR_T = F_T = A_T = 0.0;
      if (calcIntegral_lv(MA, &MA_T)!=0){
        cout << "Something went wrong with the integral!!! \n";
        return 1;
      }      
      calcIntegral(MR, &MR_T);
      calcIntegral_bv(F, &F_T);
      calcIntegral(A, &A_T);
    }

//*****************************************************************************
    //Complete model without diffusion (nao possui termo D*delta)
    if (simCase == 3){
        A[0][0][0][0] = ( beta_A*A[0][0][0][0]*(1-(A[0][0][0][0]/k_A))
		      - ( lambda_mr*MR[0][0][0][0]*A[0][0][0][0])
		      - (lambda_ma* MA[0][0][0][0]* A[0][0][0][0])
		      - (lambda_afma*F[0][0][0][0]*A[0][0][0][0]*MA[0][0][0][0])
		      - (lambda_afmr*F[0][0][0][0]*A[0][0][0][0]*MR[0][0][0][0])
		      - m_A * A[0][0][0][0]) * deltaT + A[0][0][0][0];
        //A*F

        MR[0][0][0][0] = ((- m_Mr * MR[0][0][0][0])
		       - (gamma_ma * MR[0][0][0][0] * A[0][0][0][0])       
		       + alpha_mr * (m_estrela - MR[0][0][0][0])
            ) * deltaT + MR[0][0][0][0];
  

        MA[0][0][0][0] = ((-m_Ma * MA[0][0][0][0])
		       + (gamma_ma * MR[0][0][0][0] * A[0][0][0][0])
		       - alpha_Ma * (MA_T - MA_L)
			     ) * deltaT + MA[0][0][0][0];

        F[0][0][0][0] = (- (lambda_afma * F[0][0][0][0]* A[0][0][0][0]*MA[0][0][0][0])
          - (lambda_afmr*F[0][0][0][0]*A[0][0][0][0]*MR[0][0][0][0])
		      - (alpha_f * (F_T - F_L))
			    )* deltaT + F[0][0][0][0];
        //FL-F?

        MA_L = ( alpha_Ma * (MA_T - MA_L)) * deltaT + MA_L;
        if (MA_L < 0.0) MA_L = 0.0;
        //VLV e VLN?

        Th = (b_th*(ro_t*Th*MA_L -Th*MA_L) -b_p*MA_L*Th*B
          + alpha_t*(t_estrela - Th)) * deltaT + Th;
        if (Th < 0.0) Th = t_estrela;

        B = (b_pb*(ro_b*Th*MA_L-Th*MA_L*B)
          + alpha_b*(b_estrela - B)) * deltaT + B;
        if (B < 0.0) B = b_estrela;

        P = (b_pp*(ro_p*Th*MA_L*B) + alpha_p*(p_estrela - P)) * deltaT + P;
        if (P < 0.0) P = p_estrela;

        F_L = (ro_f*P + alpha_f*(F_T-F_L)) * deltaT + F_L;
        if (F_L < 0.0) F_L = f_estrela;
        //FT-FL?
    }
//*****************************************************************************
    else{
      //Solve ODEs    
      if(simCase==0){
        MA_L = ( alpha_Ma * (MA_T - MA_L)) * deltaT + MA_L;
        if (MA_L < 0.0) MA_L = 0.0;

        Th = (b_th*(ro_t*Th*MA_L -Th*MA_L) -b_p*MA_L*Th*B
          + alpha_t*(t_estrela - Th)) * deltaT + Th;
        if (Th < 0.0) Th = t_estrela;

        B = (b_pb*(ro_b*Th*MA_L-Th*MA_L*B)
           + alpha_b*(b_estrela - B)) * deltaT + B;
        if (B < 0.0) B = b_estrela;

        P = (b_pp*(ro_p*Th*MA_L*B) + alpha_p*(p_estrela - P)) * deltaT + P;
        if (P < 0.0) P = p_estrela;

        F_L = (ro_f*P + alpha_f*(F_T-F_L)) * deltaT + F_L;
        if (F_L < 0.0) F_L = f_estrela;
    }

    //Solve PDEs
    for(int x = 0; x < Xspace; x++) {
      for(int y = 0; y < Yspace; y++) {
	      for(int z = 0; z < Zspace; z++) {
//*****************************************************************************
          //Simulates only antigen diffusion (nao tem lambdas)
	        if(simCase==1){
		
	          //Antigenos
	          A[1][x][y][z] = ( beta_A*A[0][x][y][z]*(1-(A[0][x][y][z]/k_A)) 
	            + (d_a * laplacian(A,x,y,z))
		          - m_A * A[0][x][y][z]) * deltaT + A[0][x][y][z];
              

	          if(A[1][x][y][z] != A[1][x][y][z]) {
	             cout << "A\t(NaN)-> i:"<<i<<"\t-> ("<< x << y << z << ")" << A[1][x][y][z] << "\n";
	            //  return 1;
	          }
//*****************************************************************************
	        //Simulates only innate response(equação completa)
	        }else if (simCase==2){

            //Antigenos
	          A[1][x][y][z] = ( beta_A*A[0][x][y][z]*(1-(A[0][x][y][z]/k_A))
	           - ( lambda_mr*MR[0][x][y][z]*A[0][x][y][z])
			       - (lambda_ma*MA[0][x][y][z]*A[0][x][y][z])
			       - m_A * A[0][x][y][z]
			       + (d_a * laplacian(A,x,y,z))
			       ) * deltaT + A[0][x][y][z];

	          if(A[1][x][y][z] != A[1][x][y][z]) {
	            cout << "A\t(NaN)-> i:"<<i<<"\t-> ("<< x << y << z << ")" << A[1][x][y][z]<< "\n";
	          //  return 1;
	          }

            //Macrophages     
            
	        /*****************************************************************
	        * Assuming: contact with blood vessels only on one border bv = 0
	        *           homogeneous contact with blood vessels bv = 1   
	        *           contact with blood vessels given by function bv = 2
          ******************************************************************/
	        source_mr = 0;            
	        if (((bv==0)&&(x==0))||(bv==1)||((bv==2)&&(is_bvase(x,y,z))))    
	          source_mr = alpha_mr * (m_estrela - MR[0][x][y][z]);
	        /*****************************************************************/             
	    
	        MR[1][x][y][z] = ((- m_Mr * MR[0][x][y][z])
	          - (gamma_ma * MR[0][x][y][z] * A[0][x][y][z])
			      + (d_mr * laplacian(MR,x,y,z))
			      + source_mr ) * deltaT + MR[0][x][y][z];
	   
	        if(MR[1][x][y][z] != MR[1][x][y][z] ) {
	          cout << "MR\t(NaN)-> i: " << i << "-> "<< MR[1][x][y][z] <<"(" << x << y << z << ")" << "\n";
	          // return 1;
	        }

	        MA[1][x][y][z] = ((-m_Ma * MA[0][x][y][z])
	         + (gamma_ma * MR[0][x][y][z] * A[0][x][y][z])
			     + (d_ma * laplacian(MA,x,y,z))
			     ) * deltaT + MA[0][x][y][z];
	   
	        if(MA[1][x][y][z] != MA[1][x][y][z] ) {
	          cout << "MA\t(NaN)-> i: " << i << "-> "<< MA[1][x][y][z] <<"(" << x << y << z << ")"<< "\n" ;
	          //  return 1;
	        }
//*****************************************************************************

	        //Simulates complete model
          }else if(simCase==0){
	          //Antigenos
            A[1][x][y][z] = ( beta_A*A[0][x][y][z]*(1-(A[0][x][y][z]/k_A))
	            - ( lambda_mr*MR[0][x][y][z]*A[0][x][y][z])
		          - ( lambda_ma * MA[0][x][y][z] * A[0][x][y][z])
		          - ( lambda_afma*F[0][x][y][z]*A[0][x][y][z]*MA[0][x][y][z])
              - ( lambda_afmr*F[0][x][y][z]*A[0][x][y][z]*MR[0][x][y][z])
		          - m_A * A[0][x][y][z]
		          + (d_a * laplacian(A,x,y,z))
		          ) * deltaT + A[0][x][y][z];
           
	          if(A[1][x][y][z] < tol) {A[1][x][y][z] = 0.0;} 
            
            if(A[1][x][y][z] != A[1][x][y][z]) {
              cout << "A\t(NaN)-> i:"<<i<<" -> ("<< x << y << z << ") ->" << A[1][x][y][z] << "\n";
	            //   return 1;
            }

            //Macrophages
            
	          /*****************************************************************
	          * Assuming: contact with blood vessels only on one border bv = 0
	          *           homogeneous contact with blood vessels bv = 1   
	          *           contact with blood vessels given by function bv = 2
	          ******************************************************************/
	          source_mr = 0;            
	          if (((bv==0)&&(x==0))||(bv==1)||((bv==2)&&(is_bvase(x,y,z))))    
	            source_mr = alpha_mr * (m_estrela - MR[0][x][y][z]);
	          /*****************************************************************/  
	    
	          MR[1][x][y][z] = ((- m_Mr * MR[0][x][y][z])
	           - (gamma_ma * MR[0][x][y][z] * A[0][x][y][z])       
	           + (d_mr * laplacian(MR,x,y,z))
		         + source_mr ) * deltaT + MR[0][x][y][z];

	          if(MR[1][x][y][z] != MR[1][x][y][z] ) {
	            cout << "MR\t(NaN)-> i: " << i << "-> "<< MR[1][x][y][z] <<"(" << x << y << z << ")" << "\n";
	            // return 1;
	          }

	          /*****************************************************************
	          * Assuming: contact with lymph vessels only on one border bv = 0
	          *           homogeneous contact with lymph vessels bv = 1   
	          *           contact with lymph vessels given by function bv = 2
	          ******************************************************************/
	          migration_ma = 0;            
	          if (((lnv==0)&&(x==0))||(lnv==1)||((lnv==2)&&(is_lnvase(x,y,z))))    
	            migration_ma = alpha_Ma * (MA_T - MA_L);
	          /*****************************************************************/ 
	    
	          MA[1][x][y][z] = ((-m_Ma * MA[0][x][y][z])
	            + (gamma_ma * MR[0][x][y][z] * A[0][x][y][z])
	            + (d_ma * laplacian(MA,x,y,z))
	            - migration_ma ) * deltaT + MA[0][x][y][z];
          
            if(MA[1][x][y][z] != MA[1][x][y][z] ) {
	            cout << "MA\t(NaN)-> i: " << i << "-> "<< MA[1][x][y][z] <<" -> (" << x << y << z << ")" << "\n";
	            //   return 1;
            }
	          //Antibody
	    
	          /*****************************************************************
	          * Assuming: contact with lymph vessels only on one border bv = 0
	          *           homogeneous contact with lymph vessels bv = 1   
	          *           contact with lymph vessels given by function bv = 2
	          ******************************************************************/
	          migration_f = 0;            
	          //if (((lnv==0)&&(x==0))||(lnv==1)||((lnv==2)&&(is_lnvase(x,y,z))))    
	          if (((bv==0)&&(x==0))||(bv==1)||((bv==2)&&(is_bvase(x,y,z))))
	            migration_f = (alpha_f * (F_T - F_L));
	          /*****************************************************************/ 
	    
	          F[1][x][y][z] = (
	            - ( lambda_afma * F[0][x][y][z] * A[0][x][y][z]*MA[0][x][y][z])
	            - ( lambda_afmr*F[0][x][y][z]*A[0][x][y][z]*MR[0][x][y][z])
			        - migration_f + (d_f * laplacian(F,x,y,z))
			        )* deltaT + F[0][x][y][z];

			      if(F[1][x][y][z] != F[1][x][y][z] ) {
	            cout << "F\t(NaN)-> i: " << i << "-> "<< F[1][x][y][z] <<"(" << x << y << z << ")"<< "\n" ;
	          //  return 1;
            }
	        }
        }
      }
    }

    //atualiza variaveis necessarias em cada caso
    update(A);
    if(simCase==2||simCase==0){
      update(MR);
      update(MA);
    }
    if(simCase==0){
      update(F);
    }
  }
  t++;

}while((t < (iterPerDay*days)) && (A_T > tol));

  cout << "teste\n" << Footer(t);
  return 0;
}
