/*
  VGCKIMO - Voltage-Gated Ion Channels Kinetic Modeling For Whole-Cell Voltage-Clamp Recordings
  Copyright (C) 2014 Laura L. Perissinotti (perissinotti@gmail.com), Pablo M. De Biase (pablodebiase@gmail.com)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <unistd.h>
#include <ctime>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <float.h> // levmar
#include "levmar.h"  //levmar
using namespace std;

//Universal Constants
const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
const double R = 8314.472;                              // J/mol*K
const double F = 96485.3415;
const double Ko = 5.4;  ////Bereki: external 5.4mM
const double Ki = 145; ///Bereki
const double T = 296;  //change to room temp (296) or 37 (310).
//Beating parameters
const double maxyy = 1000.0;
// Other Global
double i_init, i_inf;
// Structures
typedef struct Cell_param{
        double V, I_Kr, sum;
        double O, I, C1, C2, C3;
        double O_drug, I_drug, C1_drug, C2_drug, C3_drug;
        double Drug;
} Cell_param;
typedef struct protparam{
  double dt,t1,t2,t3,t4,t5,t6,v1,v2,v3,v4,v5,v6,drug; 
  int ptype;
  int e1; // for protocol type 7
  double ppeak;  // for protocol type 5
  int ndata; // for protocol type 5
  int nnum,nnumacum; // nnum[nps] ; nnumacum[nps+1]
  double ffg; //ffg[nps+1]
  double weight;
} protparam;
typedef struct threepoints{
      double i1,i2,i3,t1,t2,t3;
} threepoints;
typedef struct prnt{
      bool active;
      bool extra;
      string str;
} prnt; 
 // levmar
struct adata{
  double x[10000];
};

// Functions declaration
int main ();
string trim(string str);
void initrand();
double randdouble();
double randdouble(double min, double max);
double randgauss();
double randgauss(double media, double std);
void init_params(int which, Cell_param *Cell_ptr);
void Calculate_I_Kr(int which, Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);    //Potassium current
void curr1(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error); //Markov IKr for wt
void curr2(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error); //Markov IKr with drug
void curr3(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error); //Markov IKr with drug
void curr4(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error); //Markov IKr with drug
void curr5(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);
void curr6(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);
void curr7(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);
void curr8(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);
void curr9(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);
void curr10(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);
void curr11(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);
void curr12(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);
void curr13(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);
void curr14(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);
void curr15(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);
void curr16(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error);
void acum(threepoints *points3,double newi, double newt);
void calctau(int tit, int n, double *dat, double ipeak, double *tau, double *tau2, double *a, double *ar); //levmar
void inquprotprm(protparam *protprm);
void inquf(string lineout, double *fvar, double fdef); 
void inqui(string lineout, int *ivar, int idef); 
void inqus(string lineout, string *svar, string sdef);
void inqub(string lineout, bool *bvar, bool bdef);
double funk(double dx[], int *nn);
void fnk(double *p, double *x, int m, int n, void *data);
void jacfnk(double *p, double *jac, int m, int n, void *data);
void funkp(double dx[], int *nn, string str);
void funcval(int protnum, double *ffval);
void namefile(int protnum, int val, string &name, string extra = "");
string num2str(double num);
string int2str(int num);
void protocol(int protnum, int tit, double dx[],double *V,double *YY, bool *error);
void protocol0(int protnum, int tit, double dx[],double *V,double *YY, bool *error);
void protocol1(int protnum, int tit, double dx[],double *V,double *YY, bool *error);
void protocol2(int protnum, int tit, double dx[],double *V,double *YY, bool *error);
void protocol3(int protnum, int tit, double dx[],double *V,double *YY, bool *error);
void protocol4(int protnum, int tit, double dx[],double *V,double *YY, bool *error);
void protocol5(int protnum, int tit, double dx[],double *V,double *YY, bool *error);
void protocol6(int protnum, int tit, double dx[],double *V,double *YY, bool *error);
void protocol7(int protnum, int tit, double dx[],double *V,double *YY, bool *error);
void protocol8(int protnum, int tit, double dx[],double *V,double *YY, bool *error);
// to add more protocols change here
void timestamp ();
void convanal (int *nn, double dx[]);
Cell_param CellInit;
prnt prn;
const int nprot=7; // to add more protocols change here
//protocol type
string protname[nprot];
// protocols all {
int nps;
protparam *protprm;
// }
int curr;
int *nnumaux,*nnumacumaux;
//const int maxdata=1000;
//double xx[maxdata*nprot],ye[maxdata*nprot];
double *xx,*yy,*yd1,*yd2,*yd3;
double *xxaux, *yyaux;
int *trckp,td,tdaux;
int *frdx;
double *fldx;
int cnterr; // errors counter
bool *upbb,*lowbb,debugerrors,warn;
double *upbf,*lowbf;
double dta = 0.01; // delta for derivative purposes
//boundary penalty function penalty = odp*exp(bpr*(limit-x))
//at x=limit => penalty=odp
//at x=cbp+limit => penalty=cbp
const double obp = 10000.0; // In Boundary Penalty
//const double cbp = 0.001; // Close to Boundary Penalty
//double cbd; // Close to Boundary distance (same as h0)
//const double bpr = log(obp/cbp)/cbd ; // Boundary parameters rate
extern "C" {
  double praxis_ (double *t0, double *h0, int *n, int *prin, double x[], double f(double xx[],int *nn));
  void bax_ (int *nn, int *mm, double *b, double *a, double *x);
  void matmult_ (int *nn, int *mm, double *b, double *a, double *x);
}

int main () {
  int n;
  int m;
  double *dx,*mindx,*xxt,*y1t,*y2t,*y3t,*xtaux,*ytaux;
  double t0 = 0.000001;  //tolerance
  double h0 = 0.1;    //step size for minimization original=10.0
  double pr,minpr; 
  int prin = 3;
  int i,j,ln,lnt;
  bool *onoff,optim,doconvanal,guess,gauss;
  double gaussstd;
  int guessnum;
  string line;
  string *file;
  string fileaux;
  prn.active=false;
  prn.extra=false;
  CellInit.V = 0.0;
  CellInit.I_Kr = 0.0;
  CellInit.sum = 0.0;
  CellInit.O = 0.0;
  CellInit.I = 0.0;
  CellInit.C1 = 0.0;
  CellInit.C2 = 0.0;
  CellInit.C3 = 0.0;
  CellInit.O_drug = 0.0;
  CellInit.I_drug = 0.0;
  CellInit.C1_drug = 0.0;
  CellInit.C2_drug = 0.0;
  CellInit.C3_drug = 0.0;
  CellInit.Drug = 0.0;
  protname[0]="SSA-I";       //1
  protname[1]="ACT";         //6
  protname[2]="BLOCK";       //3
  protname[3]="DRUG";        //7
  protname[4]="DEACdouble";  //5
  protname[5]="DEACTable";   //8
  protname[6]="INAC03";      //4
  cnterr=0;
  cout << "\n\n==================================================\nVGC-Kimo v2015.04.23\nAuthors: Pablo M. De Biase & Laura L. Perissinotti\n==================================================\n\n\n";
  inqui ("Current Model   (1-16) [1]: ",&curr,1);
  cout << "Initial Parameters for Current Model: " << endl;
  inquf("  Cell I  [ 0.0 ] = ",&CellInit.I,0.0);
  inquf("  Cell C2 [ 0.0 ] = ",&CellInit.C2,0.0);
  inquf("  Cell C3 [ 1.0 ] = ",&CellInit.C3,1.0);
  inquf("  Cell O  [ 0.0 ] = ",&CellInit.O,0.0);
  if ( curr > 2 ) {
    inquf("  Cell Drug    [ 30.0 ] = ",&CellInit.Drug,30.0);
    inquf("  Cell I_Drug  [  0.0 ] = ",&CellInit.I_drug,0.0);
    inquf("  Cell C1_Drug [  0.0 ] = ",&CellInit.C1_drug,0.0);
    inquf("  Cell C2_Drug [  0.0 ] = ",&CellInit.C2_drug,0.0);
    inquf("  Cell C3_Drug [  0.0 ] = ",&CellInit.C3_drug,0.0);
    inquf("  Cell O_Drug  [  0.0 ] = ",&CellInit.O_drug,0.0);
    CellInit.C1 = 1 - ( CellInit.I + CellInit.C3 + CellInit.C2 + CellInit.O + CellInit.I_drug + CellInit.C3_drug + CellInit.C2_drug + CellInit.O_drug + CellInit.C1_drug ); //effect of the drug
    CellInit.sum=  CellInit.O + CellInit.C1 + CellInit.C2 + CellInit.C3 + CellInit.I + CellInit.O_drug + CellInit.C1_drug + CellInit.C2_drug + CellInit.C3_drug + CellInit.I_drug;
    }
  else {
    CellInit.C1 = 1 - ( CellInit.I + CellInit.C3 + CellInit.C2 + CellInit.O);
    CellInit.sum = CellInit.O + CellInit.C1 + CellInit.C2 + CellInit.C3 + CellInit.I;
    } 
  if      (curr==1) m=14;
  else if (curr==2) m=16;
  else { cerr << "Current Model not implemented" << endl; return 1; } 

  cout << "Number of Variables to for Current Model " << curr << " :  " << m << endl;
  fldx=new double[m]; 
  frdx=new int[m];
  onoff=new bool[m];
  for ( i=0 ; i<m ; i++ ) {
    cout << "Enter Initial Value for " << i << " [1.0] = ";
    inquf("",&fldx[i],1.0);
    }
  inqus ("Keep all Variables free (y/n)? [y]: ",&line,"y");
  n=0;
  for (i=0 ; i<m ; i++ ) {
    if ( line == "n" || line == "N" ) {
      cout << "Variable #" << i << " fixed (on/off) [off]: ";
      inqub ("",&onoff[i],false);  }
    else
      onoff[i]=false;
    if (!onoff[i]) {
       frdx[n]=i;
       n+=1; } 
    }
  dx=new double[n];
  upbb=new bool[n];
  lowbb=new bool[n];
  upbf=new double[n];
  lowbf=new double[n];
  for (i=0; i<n ; i++ ) { dx[i]=fldx[frdx[i]]; upbb[i]=false; lowbb[i]=false; }
  inqus ("Enable restrains for free variables (y/n)? [n]: ",&line,"n");
  if ( line == "y" || line == "Y" ) {
    for (i=0 ; i<n ; i++) {
      cout << "Variable #" << frdx[i] << ":      ";
      cout << "Lower boundary [-inf] = ";getline(cin, line); line=trim(line.substr(0, line.find("#"))); if ( line != "" ) { lowbf[i] = atof(line.c_str()); lowbb[i]=true;}
      if ( lowbb[i] )
        cout << lowbf[i] << "      |      ";
      else
        cout << "-inf         |          ";
      cout << "Upper boundary [+inf] = ";
      getline(cin, line); line=trim(line.substr(0, line.find("#"))); if ( line != "" ) { upbf[i] = atof(line.c_str()); upbb[i]=true;}
      if ( upbb[i] )
        cout << upbf[i] << endl;
      else
        cout << "+inf" << endl;
      }
    cout << " Defined boundaries: " << endl;
    for (i=0 ; i<n ; i++ ) {
      if ( lowbb[i] )
        cout << "     " << lowbf[i] << " < ";
      else
        cout << "-inf < ";
      cout << " Variable #" << frdx[i] << " < "; 
      if ( upbb[i] )
        cout << upbf[i] << endl;
      else
        cout << "+inf" << endl;
      }
    }
  inqui("Number of Protocols to use: ",&nps,1);
  if (nps<1) { cerr << "Must be higher than 0" << endl; return 1; }
  protprm=new protparam[nps+1];
  file=new string[nps];
  protprm[0].dt=0.01;
  protprm[0].t1=0.0;
  protprm[0].t2=0.0;
  protprm[0].t3=0.0;
  protprm[0].t4=0.0;
  protprm[0].t5=0.0;
  protprm[0].t6=0.0;
  protprm[0].v1=0.0;
  protprm[0].v2=0.0;
  protprm[0].v3=0.0;
  protprm[0].v4=0.0;
  protprm[0].v5=0.0;
  protprm[0].v6=0.0;
  protprm[0].drug=0.0;
  protprm[0].weight=1.0;
  protprm[0].e1=1;
  protprm[0].ppeak=1.0;
  protprm[0].ndata=100;
  for ( i=1; i<nps ; i++ ) protprm[i]=protprm[0];
  cout << "Protocol Types: " << endl;
  for ( i=0; i<nprot ; i++ ) cout << i << ".- " << protname[i] << endl ;

  for ( i=0; i<nps; i++ ) { 
    cout << "Protocol #"+int2str(i) << endl;
    inqui("   Set type number:  ",&protprm[i].ptype,1); 
    inqus("   Experimental Data Filename: ",&file[i],"file"+int2str(i));
    inquf("    dt = ",&protprm[i].dt,   protprm[i].dt);
    inquf("    t1 = ",&protprm[i].t1,   protprm[i].t1);
    inquf("    t2 = ",&protprm[i].t2,   protprm[i].t2);
    inquf("    t3 = ",&protprm[i].t3,   protprm[i].t3);
    inquf("    t4 = ",&protprm[i].t4,   protprm[i].t4);
   if (protprm[i].ptype==3) {
      inquf("    t5 = ",&protprm[i].t5,   protprm[i].t5);
      inquf("    t6 = ",&protprm[i].t6,   protprm[i].t6);
      }
    inquf("    v1 = ",&protprm[i].v1,   protprm[i].v1);
    inquf("    v2 = ",&protprm[i].v2,   protprm[i].v2);
    inquf("    v3 = ",&protprm[i].v3,   protprm[i].v3);
    inquf("    v4 = ",&protprm[i].v4,   protprm[i].v4);
    if (protprm[i].ptype==3) {
      inquf("    v5 = ",&protprm[i].v5,   protprm[i].v5);
      inquf("    v6 = ",&protprm[i].v6,   protprm[i].v6);
      }
    inquf("  drug = ",&protprm[i].drug, protprm[i].drug);
    inquf("weight = ",&protprm[i].weight, protprm[i].weight);
    if (protprm[i].ptype==3) {
      inqui("    e1 = ",&protprm[i].e1, protprm[i].e1);
      } 
    else if (protprm[i].ptype==4) {
      inquf("  Relation to Peak [0.95]= ",&protprm[i].ppeak,0.95);
      inqui("  Number of data points [100]= ",&protprm[i].ndata,100);
      }
    }
// Open Experimental Files
  protprm[0].nnumacum=0;
  xxt=new double[100000];
  y1t=new double[100000];
  y2t=new double[100000];
  y3t=new double[100000];
  for ( i=0; i<100000; i++ ) y1t[i]=0.0,y2t[i]=0.0,y3t[i]=0.0;
  lnt=0;
  for ( i=0 ; i<nps ; i++ ) { 
    protprm[i].nnum=0; 
    ifstream protfile(file[i].c_str());
    if (protfile.fail()) {
      cerr << file[i].c_str() << ": " << strerror(errno) << endl;
      return 1;
      }
    if (protfile.is_open()) {
      ln=0;
      while (!protfile.eof()) {
        protfile >> xxt[lnt];
        if (!protfile.eof()) protfile >> y1t[lnt];
        if (!protfile.eof() && protprm[i].ptype == 4) protfile >> y2t[lnt];  //deacdouble third
        if (!protfile.eof() && protprm[i].ptype == 4) protfile >> y3t[lnt];  //deacdouble fourth column
        if (!protfile.eof()) ln++;
        if (!protfile.eof()) lnt++;
        }
      protprm[i].nnum=ln;
      protfile.close();
      } 
      protprm[i+1].nnumacum=lnt; 
    }
  td=lnt;
  xx=new double[td];
  yd1=new double[td];
  yd2=new double[td];
  yd3=new double[td];
  yy=new double[td];
  trckp=new int [td];
  nnumaux=new int [td];
  nnumacumaux=new int [td+1];
  for ( i=0 ; i < td ; i++ ) {
    yy[i]=0.0;
    xx[i]=xxt[i]; 
    yd1[i]=y1t[i]; 
    yd2[i]=y2t[i]; // auxiliary
    yd3[i]=y3t[i]; //auxiliary
    }
  delete [] y2t;
  delete [] y3t;
  delete [] xxt;
  delete [] y1t;
  for ( i=0 ; i<nps ; i++ ) { for ( j=protprm[i].nnumacum ; j<protprm[i+1].nnumacum ; j++ ) trckp[j]=i; }
// auxiliar data for protocol types 7 & 8
  xtaux=new double[100000];
  ytaux=new double[100000];
  lnt=0;
  for ( i=0; i<td; i++ ) nnumaux[i]=0,nnumacumaux[i]=-1;
  for ( i=0 ; i<td ; i++ ) {
    j=trckp[i];
    if ( protprm[j].ptype == 3 || protprm[j].ptype == 5 ) {
      nnumacumaux[i]=lnt;
      cout << "Protocol " << j << " (Type: " << protprm[j].ptype << ") Data # " << i << " filename: ";
      inqus ("",&fileaux,"file.dat");
      ifstream ffile(fileaux.c_str());
      if (ffile.fail()) {
        cerr << fileaux.c_str() << ": " << strerror(errno) << endl;
        return 1;
        }
      if (ffile.is_open()) {
        ln=0;
        while (!ffile.eof()) {
          ffile >> xtaux[lnt];
          if (!ffile.eof()) ffile >> ytaux[lnt];
          if (!ffile.eof()) ln++,lnt++;
          }
        nnumaux[i]=ln;
        ffile.close();
        }
      nnumacumaux[i+1]=lnt;
      }
    }
    tdaux=lnt;
    xxaux=new double[tdaux];
    yyaux=new double[tdaux];
    for ( i=0 ; i < tdaux ; i++ ) {
      xxaux[i]=xtaux[i];
      yyaux[i]=ytaux[i];
    }
  delete [] xtaux;
  delete [] ytaux;
// to add more protocols change here
  inqub ("Report Warnings (on/off) [off]: ",&warn,false);
  cout << endl;
  inqub ("Debug Errors (on/off) [off]: ",&debugerrors,false);
  cout << endl;
  inqub ("Convergence Analysis (on/off) [on]: ",&doconvanal,true);
  cout << endl;
  inqub ("Optimize (on/off) [on]: ",&optim,true);
  cout << endl;
  inqub ("Generate initial Guess (on/off) [off]: ",&guess,false);
  cout << endl;
  inqub ("Print extra data (on/off) [off]: ",&prn.extra,false);
  cout << endl;
  timestamp ();
  cout << endl;
  if (guess) {
    mindx=new double[n];
    inqub("Use gaussian instead of uniform random number generator [off]: ",&gauss,false);
    if (gauss) inquf("Standard Deviation for Gaussian [0.1]: ",&gaussstd,0.1);
    inqui ("Number of Tries [1000]: ",&guessnum,1000);
    if (!gauss) {
      for (i=0; i<n ; i++ ) {
        if ( !lowbb[i] ) inquf("  Define Lower Boundary for Variable #"+int2str(frdx[i])+": ",&lowbf[i],1.0);
        if ( ! upbb[i] ) inquf("  Define Upper Boundary for Variable #"+int2str(frdx[i])+": ",&upbf[i],1.0);
        }
      }
    for (i=0; i<n; i++ ) mindx[i]=dx[i];
    minpr=funk(dx,&n);
    cout << "Initial funk = " << minpr << endl;
    initrand();
    for (i=0; i<guessnum; i++ ){
      if (gauss) { for (j=0; j<n; j++ ) dx[j]=randgauss(mindx[j],gaussstd); }
      else { for (j=0; j<n; j++ ) dx[j]=randdouble(lowbf[j],upbf[j]); }
      pr=funk(dx,&n);      // test values
      if (pr < minpr) {   // check if better values
        minpr=pr;      // store best val
        for (j=0; j<n; j++ ) mindx[j]=dx[j]; // store best params
        }
      }
    cout << "Best starting values giving funk = " << minpr << endl;
    for (i=0; i<n; i++ ) cout << "Variable # " << i << "    val = " << mindx[i] << endl;
    for (i=0; i<n; i++ ) dx[i]=mindx[i];
    delete [] mindx;
    timestamp ();
    cout << endl;
    }
  if (optim) {
    inquf("Minimization Tolerance [0.000001]: ",&t0,0.000001);
    inquf("Maximum step size move for minimization [0.1]: ",&h0,0.1);
    cout << endl << "Computing initial function value ..." << endl << endl;
    funkp(dx, &n, "ini");
    cout << "Optimizing parameters ..." << endl << endl;
    pr = praxis_(&t0, &h0, &n, &prin, dx, funk);
    cout << " Minimized:" << "\n";
    for ( i=0 ; i<m ; i++ ) {
      if (onoff[i]) line="                   fixed";
      else line="";
      cout << i << "   " << fldx[i] << line << "\n";
      }
    cout << "  FX = " << pr << "\n";
    cout << endl;
    timestamp ();
    cout << endl;
    }
  cout << endl;
  cout << "Computing final function value ..." << endl << endl;
  funkp(dx, &n,"end");
  cout << endl << "  Number of Errors: " << cnterr << endl;
  timestamp ();
  cout << endl;
  if (doconvanal) {
    convanal(&n,dx);
    timestamp ();
    cout << endl;
    }
  cout << " The End\n" ;
  delete [] dx;
  return 0;
}

//use this first function to seed the random number generator,
//call this before any of the other functions
void initrand()
{
    srand((unsigned)(time(0)));
}

//generates a psuedo-random double between 0.0 and 0.999...
double randdouble()
{
    return rand()/(double(RAND_MAX)+1);
}

//generates a psuedo-random double between min and max
double randdouble(double min, double max)
{
    if (min>max)
    {
        return randdouble()*(min-max)+max;
    }
    else
    {
        return randdouble()*(max-min)+min;
    }
}

double randgauss() {
  return cos(randdouble()*2.0*pi)*sqrt(-2.0*log(1.0-randdouble()));
}

double randgauss(double media, double std) {
  return randgauss()*std+media;
}

void convanal (int *nn, double dx[]) {
  int h,i,j,k,l,bb,xx;
  int nm=*nn;
  int nx=nm+nm*(nm+1)/2;
  int nb=nm*(nm-1)*4;
  int tnprot=nps+1;
  double x[nx];
  double b[nb],bt[nb];
  double a[nb*nx];
  double comp;
  double dxx[nm];
  double ffb[nb*tnprot];
  double ff0[tnprot]; 
  double df[nm],ddf[nm*nm];
  for (i=0; i<nb*nx; i++) a[i]=0.0; // initialize a
  for (i=0; i<nb*tnprot; i++) ffb[i]=0.0; // initialize ffb
  cout << endl << "Doing Convergence Analysis ..." << endl;
  ff0[0]=funk(dx,&nm); // Total Function Value
  for ( h=1 ; h<tnprot ; h++ ) ff0[h]=protprm[h-1].ffg;   // Compute Partials Function Values to add more protocols change here
  bb=0;
  for ( i=0 ; i<nm ; i++ ) dxx[i]=dx[i]; 
  xx=2*nm;
  cout << "Building Matrix ..." << endl;
  for ( i=0 ; i<nm-1 ; i++ ) {
    cout << i << "    "; 
    for ( j=i+1 ; j<nm ; j++ ) {
      cout << j << " ... "; 
      for ( k=-1 ; k<=1 ; k++ ) {
        for ( l=-1 ; l<=1 ; l++ ) {
          if ( !(k==0 && l==0) ) { 
            dxx[i]=dx[i]+k*dta;
            dxx[j]=dx[j]+l*dta;
            ffb[bb]=funk(dxx,&nm)-ff0[0];
            for ( h=1 ; h<tnprot ; h++ ) ffb[bb+h*nb]=protprm[h-1].ffg-ff0[h]; 
            for ( h=0 ; h<nx ; h++ ) a[bb+h*nb]=0.0;
            a[bb+i*nb]=k*dta;
            a[bb+j*nb]=l*dta;
            a[bb+(i+nm)*nb]=pow(k*dta,2)/2;
            a[bb+(j+nm)*nb]=pow(l*dta,2)/2;
            a[bb+xx*nb]=k*l*dta*dta;
            bb+=1;
          }
        }
      }
      xx+=1;
      dxx[j]=dx[j];
    }
    cout << endl;
    dxx[i]=dx[i];
  }
  for ( k=0 ; k<tnprot ; k++ ) {
    cout << endl << "For Function Partial Value (0 => Total Value): " << k << endl;
    cout << endl << "Solving linear equation ..." << endl;
    for ( h=0 ; h<nb ; h++ ) b[h]=ffb[h+k*nb];
    bax_ (&nb, &nx, b, a, x);
    for ( i=0 ; i<nm ; i++ ) df[i]=x[i];
    for ( i=0 ; i<nm ; i++ ) ddf[i + i*nm]=x[i+nm];
    xx=2*nm;
    for ( i=0 ; i<nm-1 ; i++ ) {
      for ( j=i+1 ; j<nm ; j++ ) {
        ddf[i+j*nm]=x[xx]; 
        ddf[j+i*nm]=x[xx]; 
        xx+=1;
      }
    }
    cout << endl << "Comparing results ..." << endl;
    matmult_(&nb, &nx, bt, a, x); 
    comp=0.0;
    for ( i=0 ; i<nb; i++ ) {
      comp+=pow(bt[i]-b[i],2);
 //    cout << b[i] << "   " << bt[i] << "   " << bt[i]-b[i] << endl;
      }
    cout << " Data RMSd= " << comp/nb << endl; 
    cout << "Variable     First derivative" << endl;
    for ( i=0 ; i<nm ; i++ ) cout << frdx[i] << "    " << df[i] << endl;
    cout << endl << endl << "Second derivative matrix" << endl;
    cout << "---    ";
    for ( i=0 ; i<nm ; i++ ) cout << frdx[i] << "    ";
    cout << endl;
    for ( j=0 ; j<nm ; j++ ) {
      cout << frdx[j] << "    ";
      for ( i=0 ; i<=j ; i++ ) cout << ddf[i+j*nm] << "    ";
      cout << endl;
      }
    }
  cout << endl << "  Number of Errors: " << cnterr << endl;
  timestamp ();
  cout << endl;
}

string trim(string str) {
  stringstream trimmer;
  trimmer << str;
  str.clear();
  trimmer >> str;
  return str;
}

void inquf(string lineout, double *fvar, double fdef) {
  string line;
  cout << lineout;
  getline(cin, line);
  line=trim(line.substr(0, line.find("#")));
  if ( line == "" ) { *fvar = fdef ; } 
  else { *fvar = atof(line.c_str()); }
  cout << *fvar << endl;
}

void inqui(string lineout, int *ivar, int idef) {
  string line;
  cout << lineout;
  getline(cin, line);
  line=trim(line.substr(0, line.find("#")));
  if ( line == "" ) { *ivar = idef ; }
  else { *ivar = atoi(line.c_str()); }
  cout << *ivar << endl;
}

void inqus(string lineout, string *svar, string sdef) {
  string line;
  cout << lineout;
  getline(cin, line);
  line=trim(line.substr(0, line.find("#")));
  if ( line == "" ) { *svar = sdef ; }
  else { *svar = line; }
  cout << *svar << endl;
} 

void inqub(string lineout, bool *bvar, bool bdef) {
  string line;
  cout << lineout;
  getline(cin, line);
  line=trim(line.substr(0, line.find("#")));
  if ( line == "off" ) { *bvar = false; }
  else if ( line == "on" ) { *bvar = true; }
  else { *bvar = bdef ; }
  if ( *bvar ) { cout << "on" << endl; }
  else { cout << "off" << endl; }
} 

void funcval(int protnum, double *ffval)
  {
  int i,ii,j;
  double iimax=1.0,ffv=0.0;
  ii=protnum;
  j=protprm[ii].ptype;
  if (j < 3) {
      iimax=yd1[protprm[ii].nnumacum]/yy[protprm[ii].nnumacum]; }
  // Compute RMS
  for (i=protprm[ii].nnumacum; i<protprm[ii+1].nnumacum; i++ ) {
    if (j == 3 || j == 4 || j == 5 ) ffv+=yy[i];
    else ffv+=pow(yy[i]*iimax-yd1[i],2);
    }
  *ffval=ffv;
  }

double funk(double dx[], int *nn) {
  int i,n,j;
  double ffn[nps],ff;
  bool error[td];
  n=*nn;
  ff=0.0;
  for (i=0 ; i<=nps ; i++ ) protprm[i].ffg=0.0;
  for (i=0; i<n ; i++ ) {
    if ( lowbb[i] ) { if ( dx[i] < lowbf[i] ) protprm[nps].ffg+=obp+obp*(lowbf[i]-dx[i]); }
    if ( upbb[i] )  { if ( dx[i] >  upbf[i] ) protprm[nps].ffg+=obp+obp*(dx[i]-upbf[i]); }
    }
  ff+=protprm[nps].ffg;
  if (protprm[nps].ffg==0.0) {
    for (i=0 ; i<n ; i++ ) fldx[frdx[i]]=dx[i];
    for ( i=0 ; i<nps ; i++ ) ffn[i]=0;
    #pragma omp parallel for default(shared) private(i,j)
    for ( i=0 ; i<td ; i++ ) {
      j=trckp[i];
      protocol(j,i,fldx,&xx[i],&yy[i],&error[i]);
      }
    for ( i=0 ; i<td ; i++) { if (error[i]) cnterr+=1; }
    for ( i=0 ; i<nps ; i++ ) {
      funcval(i, &ffn[i]);
      protprm[i].ffg=ffn[i]*protprm[i].weight/protprm[i].nnum; 
      ff+=protprm[i].ffg;
      }
    }
  return ff;
}

void funkp(double dx[], int *nn, string str) {
  int i,j,t;
  double ff;
  prn.active=true;
  prn.str=str;
  ff=funk(dx,nn);
  prn.active=false;
  for ( i=0 ; i<nps ; i++ ) {
    t=protprm[i].ptype;
    ofstream result_file((str+protname[t]+"_#"+int2str(i)+"_param.txt").c_str()   );    
    if ( t == 3 || t == 4 || t == 5) 
      for (j=protprm[i].nnumacum; j<protprm[i+1].nnumacum; j++) result_file << xx[j] << "   " << yy[j] << "   " << 0.0 << endl;
    else if (t<3) 
      for (j=protprm[i].nnumacum; j<protprm[i+1].nnumacum; j++) result_file << xx[j] << "   " << yy[j]*yd1[protprm[i].nnumacum]/yy[protprm[i].nnumacum] << "   " << yd1[j] << endl;
    else 
      for (j=protprm[i].nnumacum; j<protprm[i+1].nnumacum; j++) result_file << xx[j] << "   " << yy[j] << "   " << yd1[j] << endl;
    result_file.close();
    cout << "  f" << i << "(" << t << "):   " << protprm[i].ffg << endl; 
    } 
  cout << "  f" << nps << "(b):   " << protprm[nps].ffg << endl; 
  cout << "  ff   :   " << ff << endl;
  return;
}

string num2str(double num) {
  stringstream ss (stringstream::in | stringstream::out);
  double var;
  var=num;
  ss << var;
  return ss.str();
}

string int2str(int num) {
  stringstream ss (stringstream::in | stringstream::out);
  int var;
  var=num;
  ss << var;
  return ss.str();
}

void namefile(int protnum, int val, string &name, string extra) {
    name="";
    if (prn.active) name=prn.str+"_"+extra+"_"+int2str(protnum)+"_"+int2str(val)+".txt";
}

// Change here to add more protocols
void protocol(int protnum, int tit, double dx[], double *V, double *YY, bool *error) {
  int type=protprm[protnum].ptype;
       if (type == 0) protocol0(protnum,tit,dx,V,YY,error);
  else if (type == 1) protocol1(protnum,tit,dx,V,YY,error);
  else if (type == 2) protocol2(protnum,tit,dx,V,YY,error);
  else if (type == 3) protocol3(protnum,tit,dx,V,YY,error);
  else if (type == 4) protocol4(protnum,tit,dx,V,YY,error);
  else if (type == 5) protocol5(protnum,tit,dx,V,YY,error);
  else if (type == 6) protocol6(protnum,tit,dx,V,YY,error);
  else cout << "Invalid protocol number" << endl;
}

void protocol0(int protnum, int tit, double dx[], double *V,double *YY, bool *error)
{
  int tcy,p=protnum;
  int t1i,t2i,t3i,t4i;
  double t, dt, i_peak=0.0;
  threepoints ct;
  bool peakfound;
  Cell_param Cell;
  init_params(curr, &Cell);
  Cell.Drug=protprm[p].drug;
  dt = protprm[p].dt;
  string name;
  namefile(p,tit,name);
  ofstream file(name.c_str());
  t1i=(int)(protprm[p].t1/dt);
  t2i=(int)(protprm[p].t2/dt);
  t3i=(int)(protprm[p].t3/dt); // P2 converted to number of steps Jinqing=50 variation to -100  P2 here is same as P1 in bereki protocol
  t4i=(int)(protprm[p].t4/dt); // P3  final -100 mV
  ct.t1=0.0;
  ct.t2=0.0;
  ct.t3=0.0;
  ct.i1=0.0;
  ct.i2=0.0;
  ct.i3=0.0;
  t=0.0;
/* SSA-I Protocol
            P1(1s)   
           -----     P3 (0.1s)(current measurement)
          |     |--|---| 
          |     |--|   |
  rest----|     |--|   |----rest 
                |--|
                 P2 (5ms)
*/               
  t = dt;
  Cell.V = protprm[p].v1;
  for (tcy = 0; tcy < t1i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref1  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t = t + dt;
    }
  Cell.V = protprm[p].v2;  
  for (tcy = 0 ; tcy < t2i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref2  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t = t + dt;
    }
  Cell.V = *V;
  for (tcy = 0 ; tcy < t3i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref3  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t = t + dt;
    }
  tcy=0;
  peakfound=false;
  Cell.V = protprm[p].v4;  
  while (tcy<t4i && (!peakfound || prn.active) )
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref4  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
//    if (tcy==1) i_peak=ct.i2; 
//    if (tcy>0 && !peakfound) {if ( (ct.i2 > ct.i1 && ct.i2 > ct.i3) || (ct.i2 < ct.i1 && ct.i2 < ct.i3) ) i_peak=ct.i2,peakfound=true; }  // find peak
    if (tcy>0 && !peakfound) {
      if ( (ct.i2 > 0.0 && ct.i2 > i_peak) || (ct.i2 < 0.0 && ct.i2 < i_peak) ) { i_peak=ct.i2; }
      else if ( (ct.i2 > 0.0 && ct.i2 < i_peak) || (ct.i2 < 0.0 && ct.i2 > i_peak) ) { peakfound=true; }  // find peak
      }
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t+=dt;
    tcy++;
    }
  if (warn && !peakfound) cout << "Warning: Peak not found in protocol0. V = " << *V << endl;
  *YY=i_peak;
  if (prn.active) file.close();
  return;
}

void protocol1(int protnum, int tit, double dx[], double *V,double *YY, bool *error)
{
  int tcy,p=protnum;
  int tvar;
  int t1i,t3i;
//  int t2i,t4i;
  double t, dt, i_peak=0.0;
  threepoints ct;
  bool peakfound;
  Cell_param Cell;
  init_params(curr, &Cell);
  Cell.Drug=protprm[p].drug;
  dt = protprm[p].dt;
  string name;
  namefile(p,tit,name);
  ofstream file(name.c_str());
  t1i=(int)(protprm[p].t1/dt);
//  t2i=(int)(protprm[p].t2/dt);
  t3i=(int)(protprm[p].t3/dt); // P2 converted to number of steps Jinqing=50 variation to -100  P2 here is same as P1 in bereki protocol
//  t4i=(int)(protprm[p].t4/dt); // P3  final -100 mV
  tvar=(int)(*V/dt);
  ct.t1=0.0;
  ct.t2=0.0;
  ct.t3=0.0;
  ct.i1=0.0;
  ct.i2=0.0;
  ct.i3=0.0;
  t=0.0;
  t = dt;
  // WaitTime 
  Cell.V = protprm[p].v1;
  for (tcy = 0; tcy < t1i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref5  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t = t + dt;
    }
  // Potential Cycle Interval
  Cell.V = protprm[p].v2;  //P1
  for (tcy = 0 ; tcy < tvar; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref6  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t = t + dt;
    }
  Cell.V = protprm[p].v3;
  tcy=0;
  peakfound=false;
  while (tcy<t3i && (!peakfound || prn.active) )
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref7  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
     *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
//    if (tcy==1) i_peak=ct.i2; 
//    if (tcy>0 && !peakfound) {if ( (ct.i2 > ct.i1 && ct.i2 > ct.i3) || (ct.i2 < ct.i1 && ct.i2 < ct.i3) ) i_peak=ct.i2,peakfound=true; }  // find peak
    if (tcy>0 && !peakfound) {
      if ( (ct.i2 > 0.0 && ct.i2 > i_peak) || (ct.i2 < 0.0 && ct.i2 < i_peak) ) { i_peak=ct.i2; }
      else if ( (ct.i2 > 0.0 && ct.i2 < i_peak) || (ct.i2 < 0.0 && ct.i2 > i_peak) ) { peakfound=true; }  // find peak
      }
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t+=dt;
    tcy+=1;
    }
  if (warn && !peakfound) cout << "Warning: Peak not found in protocol1. t = " << *V << endl;
  *YY=i_peak;
  if (prn.active) file.close();
  return;
}

void protocol2(int protnum, int tit, double dx[], double *V,double *YY, bool *error)
{
  int tcy,p=protnum,t1i,t2i,t3i,t4i;
  double t, dt, i_peak=0.0;
  threepoints ct;
  const double t_hold = protprm[p].t2;  //P1 jinqing=-100
  const double waitTime = protprm[p].t1;  //***************You MUST hold for 30,000ms to get block of late current!!!!
  bool peakfound;
  Cell_param Cell;
  init_params(curr, &Cell);
  Cell.Drug=protprm[p].drug;
  dt = protprm[p].dt;
  string name;
  namefile(p,tit,name);
  ofstream file(name.c_str());
  t1i=(int)(waitTime/dt);
  t2i=(int)(t_hold/dt);
  t3i=(int)(protprm[p].t3/dt); // P2 converted to number of steps Jinqing=50 variation to -100  P2 here is same as P1 in bereki protocol
  t4i=(int)(protprm[p].t4/dt); // P3  final -100 mV
  Cell.Drug=*V;
  ct.t1=0.0,ct.t2=0.0,ct.t3=0.0,ct.i1=0.0,ct.i2=0.0,ct.i3=0.0;
// BLOCK
  t=0.0;
  t = dt;
  // First Cycle
  Cell.V = protprm[p].v1;
  for (tcy = 0; tcy < t1i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref8  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t+=dt;
    }
  // Second Cycle
  Cell.V = protprm[p].v2;
  for (tcy = 0; tcy < t2i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref9  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t = t + dt;
    }
  // Third Cycle, find peak
  Cell.V = protprm[p].v3;  //P1
  tcy=0;
  peakfound=false;
  while (tcy<t3i && (!peakfound || prn.active) )
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref10  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
     *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (tcy>0 && !peakfound) {
      if ( (ct.i2 > 0.0 && ct.i2 > i_peak) || (ct.i2 < 0.0 && ct.i2 < i_peak) ) { i_peak=ct.i2; }
      else if ( (ct.i2 > 0.0 && ct.i2 < i_peak) || (ct.i2 < 0.0 && ct.i2 > i_peak) ) { peakfound=true; }  // find peak
      }
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t+=dt;
    tcy+=1;
    }
  if (warn && !peakfound) cout << "Warning: Peak not found in protocol2. Drug = " << Cell.Drug << endl;
  *YY=i_peak;
  if ( prn.active ) { 
    // Fourth Cycle
    Cell.V = protprm[p].v4;
    for (tcy = 0; tcy < t4i; tcy++)
      {
      Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
      if ( *error ) {
        if (debugerrors) cerr << "  Ref11  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
        *YY=maxyy;
        return;
        }
      acum(&ct,Cell.I_Kr,t);
      if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; 
      t = t + dt;
      }
      file.close();
    }
  return;
}

void protocol3(int protnum, int tit, double dx[], double *V,double *YY, bool *error)
{
  int tcy,p=protnum,t1i,t2i,t3i,t4i,t5i,t6i;
  int i,ii,fi,e1; 
  double t, dt, i_peak=0.0,x7,y7,max;
  threepoints ct,cti;
  bool peakfound;
  Cell_param Cell;
  init_params(curr, &Cell);
  dt = protprm[p].dt;
  string name;
  namefile(p,tit,name);
  ofstream file(name.c_str());
  string nameextra;
  namefile(p,tit,nameextra,"extra");
  ofstream fileextra(nameextra.c_str());
  t1i=(int)(protprm[p].t1/dt);
  t2i=(int)(protprm[p].t2/dt);
  t3i=(int)(protprm[p].t3/dt);
  t4i=(int)(protprm[p].t4/dt); 
  t5i=(int)(protprm[p].t5/dt); 
  t6i=(int)(protprm[p].t6/dt); 
  if (nnumaux[tit] == 0 ) cerr << "Error: no data in protocol3" << endl; 
  ii=nnumacumaux[tit];
  fi=nnumacumaux[tit+1];
  e1=protprm[p].e1;
  Cell.Drug=protprm[p].drug;
  y7=0.0;
// DRUG
  ct.t1=0.0,ct.t2=0.0,ct.t3=0.0,ct.i1=0.0,ct.i2=0.0,ct.i3=0.0,cti=ct;
  t=0.0;
  t = dt;
  for ( i=0 ; i<e1 ; i++ ) {
    // First Cycle
    Cell.V = protprm[p].v1;
    for (tcy = 0; tcy < t1i; tcy++)
      {
      Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
      if ( *error ) {
        if (debugerrors) cerr << "  Ref12  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; //debug
        *YY=maxyy;
        return;
        }
      acum(&ct,Cell.I_Kr,t);
      if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
      t+=dt;
      }
    // Second Cycle
    Cell.V = protprm[p].v2;
    for (tcy = 0; tcy < t2i; tcy++)
      {
      Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
      if ( *error ) {
        if (debugerrors) cerr << "  Ref13  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; //debug
        *YY=maxyy;
        return;
        }
      acum(&ct,Cell.I_Kr,t);
      if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
      t+=dt;
      }
    // Third Cycle
    Cell.V = protprm[p].v3;
    for (tcy = 0; tcy < t3i; tcy++)
      {
      Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
      if ( *error ) {
        if (debugerrors) cerr << "  Ref14  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; //debug
        *YY=maxyy;
        return;
        }
      acum(&ct,Cell.I_Kr,t);
      if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
      t+=dt;
      }
    // Forth Cycle
    Cell.V = protprm[p].v4;  //P1
    peakfound=false;
    if (prn.active) file << " ii     max     x7     i_peak/max    xxaux[ii]     yyaux[ii]     y7 " << endl;
    for (tcy = 0; tcy < t4i; tcy++) {
      Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
      if ( *error ) {
        if (debugerrors) cerr << "  Ref15  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; //debug
        *YY=maxyy;
        return;
        }
      acum(&ct,Cell.I_Kr,t);
      if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
      if (tcy>0 && !peakfound) {
        if ( (ct.i2 > 0.0 && ct.i2 > i_peak) || (ct.i2 < 0.0 && ct.i2 < i_peak) ) { i_peak=ct.i2; }
        else if ( (ct.i2 > 0.0 && ct.i2 < i_peak) || (ct.i2 < 0.0 && ct.i2 > i_peak) ) { peakfound=true; }  // find peak
        }
//      if (tcy == 2) cti=ct;
//      if (tcy > 1 && !peakfound) {if ( (ct.i2 > ct.i1 && ct.i2 > ct.i3) || (ct.i2 < ct.i1 && ct.i2 < ct.i3) ) i_peak=ct.i2,peakfound=true; }  // find peak
//      if (fmod(tcy,1000)==0) cout << t << " " << Cell.I_Kr << endl; //debug
      t+=dt;
      }
    if (warn && !peakfound) cout << "Warning: Peak not found in protocol3 # "<< p << " . " << tit << endl;
//    if (!peakfound) i_peak=cti.i2;
    if (i==0) max=i_peak;
    x7=(t1i+t2i+t3i+t4i+t5i+t6i)*dt*i/1000.0;
    if (ii<fi) {
      if (x7>=xxaux[ii]) {
        y7+=pow(yyaux[ii]-i_peak/max,2);
        if (prn.active) file << ii << " " << max << " " << x7 << " " << i_peak/max <<  " " << xxaux[ii] << " " << yyaux[ii] << " " << y7 << endl;
        ii++;
        }
      }
    // Fifth Cycle
    Cell.V = protprm[p].v5;
    for (tcy = 0; tcy < t5i; tcy++)
      {
      Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
      if ( *error ) {
        if (debugerrors) cerr << "  Ref16  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; //debug
        *YY=maxyy;
        return;
        }
      acum(&ct,Cell.I_Kr,t);
      if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
      t+=dt;
      }
    // Sixth Cycle
    Cell.V = protprm[p].v6;
    for (tcy = 0; tcy < t6i; tcy++)
      {
      Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
      if ( *error ) {
        if (debugerrors) cerr << "  Ref17  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; //debug
        *YY=maxyy;
        return;
        }
      acum(&ct,Cell.I_Kr,t);
      if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
      t+=dt;
      }
  }
  if (prn.active) file.close();
  if (prn.extra) fileextra.close();
  *YY=y7;
  return;
}

void protocol4(int protnum, int tit, double dx[], double *V, double *YY, bool *error){
  int tcy,p=protnum;
  int t1i,t2i,t3i,t4i;
  double t0, t, dt, i_s, i_e, i_f, a, ar, tau, tau2, rms, ipeak=0.0, iat;
  double ppeak=protprm[p].ppeak;
  int ndata=protprm[p].ndata;
  threepoints ct;
  const double t_hold = protprm[p].t2; //P1=+50
  const double waitTime = protprm[p].t1;  //***************waitTime between secuence of pulses
  bool peakfound,peak75;
  Cell_param Cell;
  init_params(curr, &Cell);
  Cell.Drug=protprm[p].drug;
  dt = protprm[p].dt;
  int i,j,stp;
  double data[ndata*2];
  string name;
  namefile(p,tit,name);
  ofstream file(name.c_str());
  string nameextra;
  namefile(p,tit,nameextra,"extra");
  ofstream fileextra(nameextra.c_str());
  t1i=(int)(waitTime/dt);
  t2i=(int)(t_hold/dt);
  t3i=(int)(protprm[p].t3/dt); // 5ms converted to number of steps for P2=-100 (short pulse to recover from inactivation)
  t4i=(int)(protprm[p].t4/dt); // period of time when current is measured converted to nsteps, P3=from -100 to +60 mV
  ct.t1=0.0;
  ct.t2=0.0;
  ct.t3=0.0;
  ct.i1=0.0;
  ct.i2=0.0;
  ct.i3=0.0;
  i_inf=0.0;
  i_init=0.0;
  t = dt;
   // WaitTime 
  Cell.V = protprm[p].v1;
  for (tcy = 0; tcy < t1i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref18  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
    t = t + dt;
    }
  // Potential Cycle Interval
  Cell.V = protprm[p].v2; //P1
  for (tcy = 0; tcy < t2i; tcy++ ) {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref19  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
    t+=dt;
    }
  Cell.V = protprm[p].v3;
  for ( tcy = 0; tcy < t3i ; tcy++ ); {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref20  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
    t+=dt;
    }
  peakfound=false;
  peak75=false;
  tcy=0;
  stp=t4i/ndata;
  j=0;
  Cell.V=*V;
  t0=t;
  while ( tcy < t4i && j < ndata ) {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref21  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
    if  (!peakfound && tcy > 1) {if ( (ct.i2 > ct.i1 && ct.i2 > ct.i3) || (ct.i2 < ct.i1 && ct.i2 < ct.i3) ) j=ndata-1,peakfound=true,ipeak=ct.i2; }  // find peak
    if (peakfound && !peak75) { if (fabs(ct.i2) <= fabs(ppeak*ipeak)) peak75=true,stp=(t4i-tcy)/ndata,j=0; } // save 0.75 of the peak
    if (tcy == t4i+(j-ndata)*stp) data[2*j]=t-t0, data[2*j+1]=Cell.I_Kr, j++;
    t+=dt;
    tcy+=1;
    }
  if (j < ndata) cout << "Error: Data not captured" << endl;
  if (!peak75) {
    if (debugerrors) cerr << "Error: Peak not found in protocol4" << endl; 
    *YY=maxyy; }
  else {   
    a=0.0;
    ar=0.0;
    for (i=0;i<ndata;i++) a+=data[2*i+1]*(yd3[tit]*exp(-data[2*i]/yd1[tit])+(1.0-yd3[tit])*exp(-data[2*i]/yd2[tit]));
    for (i=0;i<ndata;i++) ar+=data[2*i+1]*data[2*i+1];
    iat=a/ar;
    *YY=0.0;
    for (i=0;i<ndata;i++) *YY+=pow(data[2*i+1]*iat-yd3[tit]*exp(-data[2*i]/yd1[tit])-(1.0-yd3[tit])*exp(-data[2*i]/yd2[tit]),2);
    *YY=*YY/ndata;
    }
  if ( isnan(*YY) ) {
    if (debugerrors) cerr << "  Ref22  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
    *error=true;
    *YY=maxyy;
    }
  if (prn.active) {
    calctau(tit,ndata,data,ipeak,&tau,&tau2,&a,&ar);
    rms=0.0;
    for (i=0;i<ndata;i++) rms+=pow(data[2*i+1]/a-ar*exp(-data[2*i]/tau)-(1.0-ar)*exp(-data[2*i]/tau2),2);
    rms=rms/ndata;
    file << "Parameters for equation: a1*exp(-t/tau1) + a2*exp(-t/tau2)" << endl;
    file << "                     or  a*ar*exp(-t/tau1) + a*(1-ar)*exp(-t/tau2)" << endl;
    file << "                   a*ar=a1 ,  a*(1-ar) = a2 ,  ar = a1/(a1+a2) ,  a = a1 + a2" << endl;
    file << endl;
    file << "   iat = " << iat << "   t0 = " << t0 << endl; 
    file << endl;
    file << " Experimental: " << endl;
    file << "   RMS(i_sim*iat vs i_exp) = " << *YY << "   tau1 = " << yd1[tit] << "   tau2 = " << yd2[tit] << "   a = " << 1.0 << "   ar = " << yd3[tit] << "   a1 = " << yd3[tit] << "   a2 = " << 1.0-yd3[tit] << endl;
    file << endl;
    file << " LevMar Fitting: " << endl;
    file << "   RMS(i_sim/a vs i_fit/a) = " << rms << "   tau1 = " << tau << "   tau2 = " << tau2 << "   a = " << a << "   ar = " << ar << "   a1 = " << a*ar << "   a2 = " << a*(1.0-ar) << endl;
    file << endl;
    file << " t-t0  ,  i_exp/iat  ,  i_sim  ,  i_fit  ,  ,  i_exp  ,  i_sim*iat  ,  i_fit*iat  ,  ,  i_exp/iat/a  ,  i_sim/a  ,  i_fit/a" << endl;
    for (i=0;i<ndata;i++) {
       t   = data[2*i];
       i_s = data[2*i+1];
       i_e = yd3[tit]*exp(-t/yd1[tit])+(1.0-yd3[tit])*exp(-t/yd2[tit]);
       i_f = a*ar*exp(-t/tau)+a*(1.0-ar)*exp(-t/tau2);
       file << t << " , " << i_e/iat << " , " << i_s  << " , " << i_f  << " ,  ,  " << i_e << " , " << i_s*iat << " , " << i_f*iat << " ,  ,  " << i_e/a/iat << " , " << i_s/a << " , " << i_f/a << " , " << endl;
      }
    }
  if (prn.active) file.close();
  if (prn.extra) fileextra.close();
  return;
}

void protocol5(int protnum, int tit, double dx[], double *V,double *YY, bool *error)
{
  int tcy,p=protnum;
  int t1i,t2i,t3i,t4i;
  int ii,fi,nnn;
  double t,t0,dt,ax,y2,e2,ye;
  Cell_param Cell;
  init_params(curr, &Cell);
  Cell.Drug=protprm[p].drug;
  dt = protprm[p].dt;
  string name;
  namefile(p,tit,name);
  ofstream file(name.c_str());
  string nameextra;
  namefile(p,tit,nameextra,"extra");
  ofstream fileextra(nameextra.c_str());
  t1i=(int)(protprm[p].t1/dt);
  t2i=(int)(protprm[p].t2/dt);
  t3i=(int)(protprm[p].t3/dt); // P2 converted to number of steps Jinqing=50 variation to -100  P2 here is same as P1 in bereki protocol
  t4i=(int)(protprm[p].t4/dt); // P3  final -100 mV
  ii=nnumacumaux[tit];
  fi=nnumacumaux[tit+1];
  nnn=nnumaux[tit];
  if (nnn == 0 ) cerr << "Error: no data in protocol5 #" << protnum << " . " << tit << endl;
  t=0.0;
  t = dt;
  // First Cycle
  Cell.V = protprm[p].v1;
  for (tcy = 0; tcy < t1i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref23  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; //debug
      *YY=maxyy;
      return;
      }
    if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
    t+=dt;
    }
  // Second Cycle
  Cell.V = protprm[p].v2;
  for (tcy = 0; tcy < t2i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) { 
      if (debugerrors) cerr << "  Ref24  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; //debug
      *YY=maxyy;
      return;
      } 
    if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
    t+=dt;
    }
  // Third Cycle
  Cell.V = protprm[p].v3;
  for (tcy = 0; tcy < t3i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref25  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; //debug
      *YY=maxyy;
      return;
      }
    if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
    t+=dt;
    }
  // Fourth Cycle
  Cell.V = *V;  //P1
  t0=0.0;
  y2=0.0,ye=0.0,e2=0.0;
  for (tcy = 0; tcy < t4i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref26  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; //debug
      *YY=maxyy;
      return;
      }
    if (prn.extra) if (fmod(tcy,1000)==0) fileextra << t << "  " << Cell.I_Kr << "  " << Cell.V << endl; // extraextra
    if (ii<fi) {
      if (t0>=xxaux[ii]) {
        y2+=yyaux[ii]*yyaux[ii];
        ye+=yyaux[ii]*Cell.I_Kr;
        e2+=Cell.I_Kr*Cell.I_Kr;
        if (prn.active) file << ii << " " << t0 << " " << Cell.I_Kr <<  " " << xxaux[ii] << " " << yyaux[ii] << endl;
        ii++;
      }
    }
    t+=dt;
    t0+=dt;
    }
  ax=ye/e2;
  *YY=(ax*ax*e2-2*ax*ye+y2)/(nnn*y2);
  if (prn.active) file << "nnn    ax    ye    y2    e2     ax*ax*e2-2*ax*ye+y2     YY " << endl;
  if (prn.active) file << nnn << " " << ax << " " << ye << " " << y2 << " " << e2 << " " << ax*ax*e2-2*ax*ye+y2 << " " << *YY  << endl;
  if (prn.active) file.close();
  if (prn.extra) fileextra.close();
  return;
}

void protocol6(int protnum, int tit, double dx[], double *V,double *YY, bool *error){
  int tcy,p=protnum;
  int t1i,t2i,t3i,t4i;
  double t, dt, i_peak=0.0, i_50ms=0.0;
  threepoints ct;
  bool peakfound=false,i50found=false;
  Cell_param Cell;
  init_params(curr, &Cell);
  Cell.Drug=protprm[p].drug;
  dt = protprm[p].dt;  // units in ms 
  int tcypeak = 0;
  string name;
  namefile(p,tit,name);
  ofstream file(name.c_str());
  t1i=(int)(protprm[p].t1/dt);
  t2i=(int)(protprm[p].t2/dt);
  t3i=(int)(protprm[p].t3/dt); // 5ms converted to number of steps for P2=-100 (short pulse to recover from inactivation)
  t4i=(int)(protprm[p].t4/dt); // period of time when current is measured converted to nsteps, P3=from -100 to +60 mV
  ct.t1=0.0;
  ct.t2=0.0;
  ct.t3=0.0;
  ct.i1=0.0;
  ct.i2=0.0;
  ct.i3=0.0;
  i_inf=0.0;
  i_init=0.0;
// Deactivation Protocol
//            P1     P3 1000ms
//           -----  ----
//          |500ms||----| 
//          |     ||----|
//  rest----|     ||----|----rest 
//                P2----
//                5ms
  t = dt;
  Cell.V = protprm[p].v1;
  for (tcy = 0; tcy < t1i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref27  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.active) { if(fmod(tcy,1000)==0) file << t << " , " << Cell.I_Kr << " , " << Cell.V << endl; }
    t = t + dt;
    }
  Cell.V = protprm[p].v2; 
  for (tcy = 0; tcy < t2i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref28  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.active) { if(fmod(tcy,1000)==0) file << t << " , " << Cell.I_Kr << " , " << Cell.V << endl; }
    t = t + dt;
    }
  Cell.V = protprm[p].v3; 
  for (tcy = 0; tcy < t3i; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref29  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.active) { if(fmod(tcy,1000)==0) file << t << " , " << Cell.I_Kr << " , " << Cell.V << endl; }
    t = t + dt;
    }
  tcy=0;
  peakfound=false;
  i50found=false;
  Cell.V = *V;
  while (tcy<t4i && (!i50found || prn.active) ) 
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cerr << "  Ref30  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (prn.active) { if(fmod(tcy,1000)==0) file << t << " , " << Cell.I_Kr << " , " << Cell.V << endl; }
    if (tcy>0 && !peakfound) {
      if ( (ct.i2 > 0.0 && ct.i2 > i_peak) || (ct.i2 < 0.0 && ct.i2 < i_peak) ) { i_peak=ct.i2,tcypeak=tcy-1; }
      else if ( (ct.i2 > 0.0 && ct.i2 < i_peak) || (ct.i2 < 0.0 && ct.i2 > i_peak) ) { peakfound=true; }  // find peak
      }
    if (tcy > 0 && peakfound && !i50found) { if (tcy == tcypeak+5000) i50found=true,i_50ms=Cell.I_Kr; } // save 0.75 of the peak
    t+=dt;
    tcy+=1;
    }
  if (warn && !peakfound) cout << "Warning: Peak not found in protocol6" << endl;
  if (warn && !i50found) cout << "Warning: i50 not found in protocol6" << endl;
  if (i50found) *YY=i_50ms/i_peak;
  else *YY=maxyy;
  if (prn.active) file.close();
  return;
}

// levmar
void fnk(double *p, double *x, int m, int n, void *data)
{
  register int i;
  double t;
  struct adata *dptr;
  dptr=(struct adata *)data;

  for ( i=0; i<n; i++) {
    t=dptr->x[i];
    x[i]=p[0]*(p[1]*exp(-t/p[2])+(1-p[1])*exp(-t/p[3])); }
}

void jacfnk(double *p, double *jac, int m, int n, void *data)
{
  register int i, j;
  double t,tmp1,tmp2;
  struct adata *dptr;
  dptr=(struct adata *)data;

  for ( i=j=0; i<n; i++) { 
    t=dptr->x[i];
    tmp1=exp(-t/p[2]);
    tmp2=exp(-t/p[3]);
    jac[j++]=p[1]*tmp1+(1.0-p[1])*tmp2;
    jac[j++]=p[0]*(tmp1-tmp2);
    jac[j++]=p[0]*p[1]*t*tmp1/(p[2]*p[2]);
    jac[j++]=p[0]*(1.0-p[1])*t*tmp2/(p[3]*p[3]);
    } 
}

void calctau(int tit, int n, double *dat, double ipeak, double *tau, double *tau2, double *a, double *ar){
  int m=4;
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  double p[m],lb[m],ub[m];
  double xd[n];
  struct adata data;
  int i,ret;
  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20; opts[4]=LM_DIFF_DELTA;
  for (i=0; i<m ; i++ ) lb[i]=-DBL_MAX,ub[i]=DBL_MAX;
  lb[1]=0,ub[1]=1.0;
 
  p[0]=ipeak;
  p[1]=yd3[tit];
  p[2]=yd1[tit];
  p[3]=yd2[tit];
  for (i=0; i<n; i++) {
    data.x[i]=dat[2*i];
    xd[i]=dat[2*i+1];
  }
  ret=dlevmar_bc_der(fnk, jacfnk, p, xd, m, n, lb, ub, NULL ,1000, opts, info, NULL, NULL, (void *)&data); // with analytic Jacobian
  *a=p[0];
  if ( p[2] < p[3] ) {
    *tau=p[2],*tau2=p[3]; 
    *ar=p[1]; }
  else {
    *tau=p[3],*tau2=p[2];
    *ar=1.0-p[1]; }
  return;
} 

void acum(threepoints *points3,double newi, double newt){
  points3->i1 = points3->i2;
  points3->i2 = points3->i3;
  points3->i3 = newi;
  points3->t1 = points3->t2;
  points3->t2 = points3->t3;
  points3->t3 = newt;
  return;
}

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}

void init_params(int which, Cell_param *Cell_ptr){
  if(which <= 2) {
    Cell_ptr->I = CellInit.I;
    Cell_ptr->C3 = CellInit.C3;
    Cell_ptr->C2 = CellInit.C2;
    Cell_ptr->O = CellInit.O;
    Cell_ptr->C1 = 1 - (Cell_ptr->I + Cell_ptr->C3 + Cell_ptr->C2 + Cell_ptr->O);
    Cell_ptr->sum = Cell_ptr->O+Cell_ptr->C1+Cell_ptr->C2+Cell_ptr->C3+Cell_ptr->I; }
  else {
    Cell_ptr->Drug = CellInit.Drug;
    Cell_ptr->I = CellInit.I;
    Cell_ptr->C3 = CellInit.C3;
    Cell_ptr->C2 = CellInit.C2;
    Cell_ptr->O = CellInit.O;
    Cell_ptr->I_drug = CellInit.I_drug;
    Cell_ptr->C1_drug = CellInit.C1_drug;
    Cell_ptr->C2_drug = CellInit.C2_drug;
    Cell_ptr->C3_drug = CellInit.C3_drug;
    Cell_ptr->O_drug = CellInit.O_drug;
    Cell_ptr->C1 = 1 - ( Cell_ptr->I + Cell_ptr->C3 + Cell_ptr->C2 + Cell_ptr->O + Cell_ptr->I_drug + Cell_ptr->C3_drug + Cell_ptr->C2_drug + Cell_ptr->O_drug + Cell_ptr->C1_drug ); //effect of the drug
    Cell_ptr->sum=  Cell_ptr->O + Cell_ptr->C1 + Cell_ptr->C2 + Cell_ptr->C3 + Cell_ptr->I+ Cell_ptr->O_drug + Cell_ptr->C1_drug + Cell_ptr->C2_drug + Cell_ptr->C3_drug + Cell_ptr->I_drug; }
}

void Calculate_I_Kr(int which, Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){	
  if(which == 1) 
    curr1(Cell_ptr, t, dt, dx, error);
  else if(which == 2)
    curr2(Cell_ptr, t, dt, dx, error);
  else if(which == 3)
    curr3(Cell_ptr, t, dt, dx, error);
  else if(which == 4)
    curr4(Cell_ptr, t, dt, dx, error);
  else if(which == 5)
    curr5(Cell_ptr, t, dt, dx, error);
  else if(which == 6)
    curr6(Cell_ptr, t, dt, dx, error);
  else if(which == 7)
    curr7(Cell_ptr, t, dt, dx, error);
  else if(which == 8)
    curr8(Cell_ptr, t, dt, dx, error);
  else if(which == 9)
    curr9(Cell_ptr, t, dt, dx, error);
  else if(which == 10)
    curr10(Cell_ptr, t, dt, dx, error);
  else if(which == 11)
    curr11(Cell_ptr, t, dt, dx, error);
  else if(which == 12)
    curr12(Cell_ptr, t, dt, dx, error);
  else if(which == 13)
    curr13(Cell_ptr, t, dt, dx, error);
  else if(which == 14)
    curr14(Cell_ptr, t, dt, dx, error);
  else if(which == 15)
    curr15(Cell_ptr, t, dt, dx, error);
  else if(which == 16)
    curr16(Cell_ptr, t, dt, dx, error);
  else {
    cout << "Invalid which" << endl;
//    exit(1);
    }
  return;
}

void curr1(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
  double E_Kr;
  double hh;
  hh = dt;
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0)); 	//var g_Kr_0: microS_per_nanoF {init: 0.024}; Ko dependence is included in the formulation of the current
	
//        const double Q10=3;    To account for correction in temperature

//        const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));

  E_Kr = (R*T/F)*log(Ko/Ki);			

//Markov Fink mechanism///////
//  C3 = C2 = C1 = O = I
/////////////////////////////
//rate cte definitions:
//C3_C2=ae
//C2_C3=be
//C2_C1=ain
//C1_C2=bin
//C1_O=aa
//O_C1=bb
//O_I=bi
//I_O=ai
////////////////////////////

/*
//Rate constant expresions Fink at 37 C

      double ae=exp(dx[0]*24.335+dx[1]*0.0112*Cell_ptr->V-dx[2]*25.914);
      double be=exp(dx[3]*13.688-dx[4]*0.0603*Cell_ptr->V-dx[5]*15.707);
      double ain=exp(dx[6]*22.746-dx[7]*25.914);
      double bin=exp(dx[8]*13.193-dx[9]*15.707);
      double ai=exp(dx[10]*30.061-dx[11]*0.0312*Cell_ptr->V-dx[12]*33.243);  // ai, bi steps different for the mutant.
      double bi=exp(dx[13]*30.016+dx[14]*0.0223*Cell_ptr->V-dx[15]*30.888)*pow((5.4/Ko),0.4);   //inactivation
      double aa=exp(dx[16]*22.098+dx[17]*0.0365*Cell_ptr->V-dx[18]*25.914);  //activation 
      double bb=exp(dx[19]*7.313-dx[20]*0.0399*Cell_ptr->V-dx[21]*15.707);  //deactivation

//Rate constants parameters at 23 C corrected for Jiqing Guo Data (NMS), and praxis opt but p1=1S
*/
  double ae=dx[0]*exp(24.288586641+dx[1]*0.0117*Cell_ptr->V-27.138936313); 
  double be=dx[2]*exp(13.637047252-dx[3]*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*exp(22.444099764-27.130772208);
  double bin=dx[5]*exp(13.110155065-16.420979385);
  double ai=dx[6]*exp(28.782084605-dx[7]*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*exp(29.690385428+dx[9]*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*exp(22.010640617+dx[11]*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*exp(7.266288229-dx[13]*0.0418*Cell_ptr->V-16.168847047);  //deactivation
/*
//Rate constants parameters at 23 C corrected for Jiqing Guo Data (NMS), and praxis opt done but p1=1S
  double ae=dx[0]*1.50063*exp(24.288586641+dx[1]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*0.598723*exp(13.637047252-dx[3]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*1.10345*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.819984*exp(13.110155065-16.420979385);
  double ai=dx[6]*1.02159*exp(28.782084605-dx[7]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*0.998552*exp(29.690385428+dx[9]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*1.86823*exp(22.010640617+dx[11]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*1.00124*exp(7.266288229-dx[13]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation

//Fink Rate constants at 23 aplying Q10 factor /To account for correction in temperature from Vandenberg et al. Am J Cell Physiol 2006
  const double Q10_a=2.1;   //same activation value for the three closed states
  const double Q10_d=1.7;   //deactivation
  const double Q10_i=2.5;   //inactivation
  const double Q10_r=2.6;   //recovery from inactivation

  const double Tfactor_a=1.0/(pow(Q10_a, (37.0-(T-273))/10.0));
  const double Tfactor_d=1.0/(pow(Q10_d, (37.0-(T-273))/10.0));
  const double Tfactor_i=1.0/(pow(Q10_i, (37.0-(T-273))/10.0));
  const double Tfactor_r=1.0/(pow(Q10_r, (37.0-(T-273))/10.0));

  double ae=Tfactor_a*exp(dx[0]*24.335+dx[1]*0.0112*Cell_ptr->V-dx[2]*25.914);
  double be=Tfactor_d*exp(dx[3]*13.688-dx[4]*0.0603*Cell_ptr->V-dx[5]*15.707);
  double ain=Tfactor_a*exp(dx[6]*22.746-dx[7]*25.914);
  double bin=Tfactor_d*exp(dx[8]*13.193-dx[9]*15.707);
  double ai=Tfactor_r*exp(dx[10]*30.061-dx[11]*0.0312*Cell_ptr->V-dx[12]*33.243);  // ai, bi steps different for the mutant.
  double bi=Tfactor_i*exp(dx[13]*30.016+dx[14]*0.0223*Cell_ptr->V-dx[15]*30.888)*pow((5.4/Ko),0.4);   //inactivation
  double aa=Tfactor_a*exp(dx[16]*22.098+dx[17]*0.0365*Cell_ptr->V-dx[18]*25.914);  //activation 
      double bb=Tfactor_d*exp(dx[19]*7.313-dx[20]*0.0399*Cell_ptr->V-dx[21]*15.707);  //deactivation


*/
//Initial contitions**********************************************************************************************************************

  // Calculation of k parameters for the Runge Kutta interations
  //Calculation of k1 values using the initial differential equation (set as f(t,V)
  double k1_O = hh*(ai*Cell_ptr->I + aa*Cell_ptr->C1 - Cell_ptr->O*(bi + bb));
  double k1_C1 = hh*(ain*Cell_ptr->C2 + bb*Cell_ptr->O - Cell_ptr->C1*(aa + bin));
  double k1_C2 = hh*(ae*Cell_ptr->C3 + bin*Cell_ptr->C1 - Cell_ptr->C2*(be + ain));
  double k1_C3 = hh*(be*Cell_ptr->C2 - ae*Cell_ptr->C3);
  double k1_I = hh*(bi*Cell_ptr->O - Cell_ptr->I*ai);
  double O_k1 = Cell_ptr->O + k1_O/2;
  //Value of O after the k1 iteration; used for calculaton of k2
  double C1_k1 = Cell_ptr->C1 + k1_C1/2;
  double C2_k1 = Cell_ptr->C2 + k1_C2/2;
  double C3_k1 = Cell_ptr->C3 + k1_C3/2;
  double I_k1 = Cell_ptr->I + k1_I/2;
  //calculation of K2 values
  double k2_O = hh*(ai*I_k1 + aa*C1_k1 - O_k1*(bi + bb));
  double k2_C1 = hh*(ain*C2_k1 + bb*O_k1 - C1_k1*(aa + bin));
  double k2_C2 = hh*(ae*C3_k1 + bin*C1_k1 - C2_k1*(be + ain));
  double k2_C3 = hh*(be*C2_k1 - ae*C3_k1);
  double k2_I = hh*(bi*O_k1 - I_k1*ai);
  double O_k2 = Cell_ptr->O + k2_O/2;
  double C1_k2 = Cell_ptr->C1 + k2_C1/2;
  double C2_k2 = Cell_ptr->C2 + k2_C2/2;
  double C3_k2 = Cell_ptr->C3 + k2_C3/2;
  double I_k2 = Cell_ptr->I + k2_I/2;
  //calculation of K3 values
  double k3_O = hh*(ai*I_k2 + aa*C1_k2 - O_k2*(bi + bb));
  double k3_C1 = hh*(ain*C2_k2 + bb*O_k2 - C1_k2*(aa + bin));
  double k3_C2 = hh*(ae*C3_k2 + bin*C1_k2 - C2_k2*(be + ain));
  double k3_C3 = hh*(be*C2_k2 - ae*C3_k2);
  double k3_I = hh*(bi*O_k2 - I_k2*ai);
  double O_k3 = Cell_ptr->O + k3_O;
  double C1_k3 = Cell_ptr->C1 + k3_C1;
  double C2_k3 = Cell_ptr->C2 + k3_C2;
  double C3_k3 = Cell_ptr->C3 + k3_C3;
  double I_k3 = Cell_ptr->I + k3_I;
  //calculation of K4 values
  double k4_O = hh*(ai*I_k3 + aa*C1_k3 - O_k3*(bi + bb));
  double k4_C1 = hh*(ain*C2_k3 + bb*O_k3 - C1_k3*(aa + bin));
  double k4_C2 = hh*(ae*C3_k3 + bin*C1_k3 - C2_k3*(be + ain));
  double k4_C3 = hh*(be*C2_k3 - ae*C3_k3);
  double k4_I = hh*(bi*O_k3 - I_k3*ai);
  //New computed values of the channel probabilities
  Cell_ptr->O = Cell_ptr->O + (k1_O + 2*k2_O + 2*k3_O + k4_O)/6;
  Cell_ptr->C1 = Cell_ptr->C1 + (k1_C1 + 2*k2_C1 + 2*k3_C1 + k4_C1)/6;
  Cell_ptr->C2 = Cell_ptr->C2 + (k1_C2 + 2*k2_C2 + 2*k3_C2 + k4_C2)/6;
  Cell_ptr->C3 = Cell_ptr->C3 + (k1_C3 + 2*k2_C3 + 2*k3_C3 + k4_C3)/6;
  Cell_ptr->I = Cell_ptr->I + (k1_I + 2*k2_I + 2*k3_I + k4_I)/6;
  
  Cell_ptr->sum=Cell_ptr->O+Cell_ptr->C1+Cell_ptr->C2+Cell_ptr->C3+Cell_ptr->I;

  if (fabs (Cell_ptr->sum -1.0) > 0.0001) {
//    cout << "Error in Wild Type Sum at t = " << t << ",\t sum = " << Cell_ptr->sum << endl;
    *error=true;
    Cell_ptr->I_Kr = 0.0;
    }
  else {
    *error=false;
    Cell_ptr->I_Kr = G_Kr*sqrt(Ko/5.4)*(Cell_ptr->O)*(Cell_ptr->V - E_Kr);
    }
  return;
}


void curr2(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
  double E_Kr;
  double hh;
  hh = dt;
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0));        //var g_Kr_0: microS_per_nanoF {init: 0.024}; Ko dependence is included in the formulation of the current

//        const double Q10=3;    To account for correction in temperature

//        const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));

  E_Kr = (R*T/F)*log(Ko/Ki);

// mechanism///////
//  C3 = C2 = C1 = O
//             \\ //
//               I
/////////////////////////////
//rate cte definitions:
//C3_C2=ae
//C2_C3=be
//C2_C1=ain
//C1_C2=bin
//C1_O=aa
//O_C1=bb
//O_I=bi
//I_O=ai
//I_C1=bo
//C1_I=ao
////////////////////////////
/*
//Rate constants parameters Silva & Rudy (2005) ??///////reformular
  double ae=dx[0]*3.02e-2*exp(dx[1]*1.48e-0*((Cell_ptr->V*96485.3415)/(8314.472*296)));
  double be=dx[2]*2.90e-3*exp(dx[3]*-9.78e-1*((Cell_ptr->V*96485.3415)/(8314.472*296)));
  double ain=dx[4]*2.17;
  double bin=dx[5]*1.08;
  double ai=dx[6]*5.45e-1*exp(dx[7]*-8.17e-1*((Cell_ptr->V*96485.3415)/(8314.472*296)))*(4.5/5.4);  // ai, bi steps different for the mutant.
  double bi=dx[8]*8.20e-1*exp(dx[9]*5.04e-1*((Cell_ptr->V*96485.3415)/(8314.472*296)))*pow((4.5/5.4),3);   //inactivation
  double aa=dx[10]*1.31e-2*exp(dx[11]*1.48e0*((Cell_ptr->V*96485.3415)/(8314.472*296)));  //activation
  double bb=dx[12]*3.3e-3*exp(dx[13]*-5.77e-1*((Cell_ptr->V*96485.3415)/(8314.472*296))); //deactivation  
  double ao=dx[14]*1.31e-2*exp(dx[15]*1.48e0*((Cell_ptr->V*96485.3415)/(8314.472*296)));
  double bo=ai*bb/bi;
*/
//Rates from MGWMN BJ review 2011  ///lau///
  double ae=dx[0]*0.0069*exp(dx[1]*0.0272*Cell_ptr->V);    //early act
  double be=dx[2]*0.0227*exp(dx[3]*-0.0431*Cell_ptr->V);   //late deact
  double ain=dx[4]*0.0266;                                     //rate limiting act
  double bin=dx[5]*0.1348;                                     //rate limiting deact
  double ai=dx[6]*0.0059*exp(dx[7]*-0.0443*Cell_ptr->V);  // recovery from inact
  double bi=dx[8]*0.0622*exp(dx[9]*0.0120*Cell_ptr->V);   //inactivation
  double aa=dx[10]*0.0218*exp(dx[11]*0.0262*Cell_ptr->V);  //activation
  double bb=dx[12]*0.0009*exp(dx[13]*-0.0269*Cell_ptr->V); //deactivation  
  double ao=dx[14]*1.29e-5*exp(dx[15]*2.71e-6*Cell_ptr->V);   //added
  double bo=(ao*ai*bb)/(aa*bi);                                        //added



//Initial contitions**********************************************************************************************************************

  // Calculation of k parameters for the Runge Kutta interations
  //Calculation of k1 values using the initial differential equation (set as f(t,V)
  double k1_O = hh*(ai*Cell_ptr->I + aa*Cell_ptr->C1 - Cell_ptr->O*(bi + bb));
  double k1_C1 = hh*(ain*Cell_ptr->C2 + bb*Cell_ptr->O + bo*Cell_ptr->I - Cell_ptr->C1*(aa + bin +ao));
  double k1_C2 = hh*(ae*Cell_ptr->C3 + bin*Cell_ptr->C1 - Cell_ptr->C2*(be + ain));
  double k1_C3 = hh*(be*Cell_ptr->C2 - ae*Cell_ptr->C3);
  double k1_I = hh*(bi*Cell_ptr->O + ao*Cell_ptr->C1 - Cell_ptr->I*(ai+bo));
  double O_k1 = Cell_ptr->O + k1_O/2;
  //Value of O after the k1 iteration; used for calculaton of k2
  double C1_k1 = Cell_ptr->C1 + k1_C1/2;
  double C2_k1 = Cell_ptr->C2 + k1_C2/2;
  double C3_k1 = Cell_ptr->C3 + k1_C3/2;
  double I_k1 = Cell_ptr->I + k1_I/2;
  //calculation of K2 values
  double k2_O = hh*(ai*I_k1 + aa*C1_k1 - O_k1*(bi + bb));
  double k2_C1 = hh*(ain*C2_k1 + bb*O_k1 +bo*I_k1 - C1_k1*(aa + bin +ao));
  double k2_C2 = hh*(ae*C3_k1 + bin*C1_k1 - C2_k1*(be + ain));
  double k2_C3 = hh*(be*C2_k1 - ae*C3_k1);
  double k2_I = hh*(bi*O_k1 + ao*C1_k1 - I_k1*(ai+bo));
  double O_k2 = Cell_ptr->O + k2_O/2;
  double C1_k2 = Cell_ptr->C1 + k2_C1/2;
  double C2_k2 = Cell_ptr->C2 + k2_C2/2;
  double C3_k2 = Cell_ptr->C3 + k2_C3/2;
  double I_k2 = Cell_ptr->I + k2_I/2;
  //calculation of K3 values
  double k3_O = hh*(ai*I_k2 + aa*C1_k2 - O_k2*(bi + bb));
  double k3_C1 = hh*(ain*C2_k2 + bb*O_k2 + bo*I_k2 - C1_k2*(aa + bin + ao));
  double k3_C2 = hh*(ae*C3_k2 + bin*C1_k2 - C2_k2*(be + ain));
  double k3_C3 = hh*(be*C2_k2 - ae*C3_k2);
  double k3_I = hh*(bi*O_k2 + ao*C1_k2 - I_k2*(ai + bo));
  double O_k3 = Cell_ptr->O + k3_O;
  double C1_k3 = Cell_ptr->C1 + k3_C1;
  double C2_k3 = Cell_ptr->C2 + k3_C2;
  double C3_k3 = Cell_ptr->C3 + k3_C3;
  double I_k3 = Cell_ptr->I + k3_I;
  //calculation of K4 values
  double k4_O = hh*(ai*I_k3 + aa*C1_k3 - O_k3*(bi + bb));
  double k4_C1 = hh*(ain*C2_k3 + bb*O_k3 + bo*I_k3 - C1_k3*(aa + bin + ao));
  double k4_C2 = hh*(ae*C3_k3 + bin*C1_k3 - C2_k3*(be + ain));
  double k4_C3 = hh*(be*C2_k3 - ae*C3_k3);
  double k4_I = hh*(bi*O_k3 + ao*C1_k3 - I_k3*(ai + bo));
  //New computed values of the channel probabilities
  Cell_ptr->O = Cell_ptr->O + (k1_O + 2*k2_O + 2*k3_O + k4_O)/6;
  Cell_ptr->C1 = Cell_ptr->C1 + (k1_C1 + 2*k2_C1 + 2*k3_C1 + k4_C1)/6;
  Cell_ptr->C2 = Cell_ptr->C2 + (k1_C2 + 2*k2_C2 + 2*k3_C2 + k4_C2)/6;
  Cell_ptr->C3 = Cell_ptr->C3 + (k1_C3 + 2*k2_C3 + 2*k3_C3 + k4_C3)/6;
  Cell_ptr->I = Cell_ptr->I + (k1_I + 2*k2_I + 2*k3_I + k4_I)/6;

  Cell_ptr->sum=Cell_ptr->O+Cell_ptr->C1+Cell_ptr->C2+Cell_ptr->C3+Cell_ptr->I;

  if (fabs (Cell_ptr->sum -1.0) > 0.0001) {
//    cout << "Error in Wild Type Sum at t = " << t << ",\t sum = " << Cell_ptr->sum << endl;
    *error=true;
    Cell_ptr->I_Kr = 0.0;
    }
  else {
    *error=false;
    Cell_ptr->I_Kr = G_Kr*sqrt(Ko/5.4)*(Cell_ptr->O)*(Cell_ptr->V - E_Kr);
    }
  return;
}
void curr3(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
return;
}
void curr4(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
return;
}
void curr5(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
return;
}
void curr6(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
  return;
}
void curr7(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
return;
}
void curr8(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
return;
}
void curr9(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
return;
}
void curr10(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
return;
}
void curr11(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
return;
}
void curr12(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
return;
}
void curr13(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
  return;
}
void curr14(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
  return;
}
void curr15(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
  return;
}
void curr16(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){
  return;
}



