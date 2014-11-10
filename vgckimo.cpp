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
#include <float.h> // levmar
#include "levmar.h"  //levmar
using namespace std;

//Universal Constants
const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
const double R = 8314.472;                              // J/mol*K
const double F = 96485.3415;
const double Ko = 5.4;  ////Bereki: external 5.4mM
const double Ki = 145; ///Bereki
const double V_duration = 3000;  //modify for drug, Lau
const double BCL_1 = 1000;
const double T = 296;  //change to room temp (296) or 37 (310).
//Beating parameters
const double V_o=-87.5;   ///V resting?
const double V_step= -110;
const double total_cycles = 300;
const double BCL_3 = 40;
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
  double thold,waittime,dt,tvc,tdec,v1,v2,v3;
} protparam;
typedef struct threepoints{
      double i1,i2,i3,t1,t2,t3;
} threepoints;
typedef struct prnt{
      bool active;
      string str;
} prnt; 
 // levmar
struct adata{
  double x[10000];
};

// Functions declaration
int main ();
string trim(string str);
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
void calctau(threepoints *pst, threepoints *pnd, double *tau); 
void calctau2(int tit, int n, double *dat, double ipeak, double *tau, double *tau2, double *a, double *ar); //levmar
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
void namefile(int protnum, double *valor, string name);
string num2str(double num);
void protocol(int protnum, int tit, double dx[],double *V,double *YY, bool *error);
void protocol0(double dx[],double *V,double *YY, bool *error);
void protocol1(double dx[],double *V,double *YY, bool *error);
void protocol2(double dx[],double *V,double *YY, bool *error);
void protocol3(double dx[],double *V,double *YY, bool *error);
void protocol4(double dx[],double *V,double *YY, bool *error);
void protocol5(int tit, double dx[],double *V,double *YY, bool *error);
void protocol6(double dx[],double *V,double *YY, bool *error);
// to add more protocols change here
void timestamp ();
void convanal (int *nn, double dx[]);
Cell_param CellInit;
prnt prn;
const int nprot=7; // to add more protocols change here
bool prot[nprot];   
int curr,nnum[nprot],nnumacum[nprot+1];
//const int maxdata=1000;
//double xx[maxdata*nprot],ye[maxdata*nprot];
double *xx,*yy,*ye,*yi;
int *trckp,td;
double ffg[nprot+1]; 
double ppeak;  // for protocol 5
int ndata; // for protocol 5
string protname[nprot];
protparam protprm[nprot];
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
  double *dx,*xxt,*yet,*yit;
  double t0 = 0.000001;  //tolerance
  double h0 = 0.1;    //step size for minimization original=10.0
  double pr; 
  int prin = 3;
  int i,j,ln,lnt;
  bool *onoff,optim,doconvanal;
  string file[nprot],line;
  prn.active=false;
  protprm[0].thold=500.0;
  protprm[0].waittime=10000.0;
  protprm[0].dt=0.01;
  protprm[0].tvc=5.0;
  protprm[0].tdec=1000.0;
  protprm[0].v1=-87.0;
  protprm[0].v2=50.0;
  protprm[0].v3=-100.0;
  for ( i=0; i<nprot ; i++ ) protprm[i]=protprm[0];
  protprm[0].thold=1000.0;
  protprm[0].tvc=1000.0;
  protprm[0].v2=-100.0;
  protprm[1].thold=1000.0;
  protprm[1].tdec=100.0;
  protprm[1].v3=40.0;
  protprm[6].tvc=500.0;
  protprm[6].tdec=500.0;
  protprm[6].thold=1000.0;
  protprm[6].v1=-80.0;
  protprm[6].v2=0.0;
  protprm[6].v3=-100.0;
  protprm[3].v1=-80.0;
  protprm[3].v2=20.0;
  protprm[3].v3=-100.0;
  protprm[3].thold=300000.0;
  protprm[3].tvc=100.0;
  protprm[3].tdec=100.0;
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
  protname[0]="SSA";
  protname[1]="SSI";
  protname[2]="DEAC";
  protname[3]="BLOCK";
  protname[4]="DEACExt";
  protname[5]="DEACdouble";
  protname[6]="ACT";
  cnterr=0;
  cout << "\n\n==================================================\nVGC-Kimo v2014.10.10\nAuthors: Pablo M. De Biase & Laura L. Perissinotti\n==================================================\n\n\n";
  inquf("Minimization Tolerance [0.000001]: ",&t0,0.000001);
  inquf("Maximum step size move for minimization [0.1]: ",&h0,0.1);
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
  inqui ("Number of Variables to Fit [14]: ",&m,14);
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
  for ( i=0; i<nprot ; i++ ) inqub (protname[i]+"   (on/off) [on]: ",&prot[i],true);
  for ( i=0; i<nprot ; i++ ) { if ( prot[i] == true ) inqus ("Experimental Data Filename for "+protname[i]+"  : ",&file[i],"file"+protname[i]); }
  cout << "Input parameters for protocols: " << endl;
  for ( i=0; i<nprot ; i++ ) { 
    if ( prot[i] == true ) {
      cout << "  " << protname[i] << ": " << endl;
      inquprotprm(&protprm[i]);
      }
    }
  cout << endl;
// Special Case 
  if ( prot[5] == true ) {
    cout << "For Protocol 5: " << endl;
    inquf("  Relation to Peak [0.95]= ",&ppeak,0.95);
    inqui("  Number of data points [100]= ",&ndata,100);
    }
// Open Experimental Files
  nnumacum[0]=0;
  xxt=new double[10000];
  yet=new double[10000];
  yit=new double[20000];
  for ( i=0; i<10000; i++ ) yit[i]=0.0;
  lnt=0;
  for ( i=0 ; i<nprot ; i++ ) { 
    nnum[i]=0; 
    if ( prot[i] == true ) {
      ifstream protfile(file[i].c_str());
      if (protfile.fail()) {
        cerr << file[i].c_str() << ": " << strerror(errno) << endl;
        return 1;
        }
      if (protfile.is_open()) {
        ln=0;
        while (!protfile.eof()) {
          protfile >> xxt[lnt];
          if (!protfile.eof()) protfile >> yet[lnt];
          if (!protfile.eof() && i == 5) {protfile >> yit[lnt]; protfile >> yit[lnt+10000]; } //deacdouble third and fourth column
          if (!protfile.eof()) ln++;
          if (!protfile.eof()) lnt++;
          }
        nnum[i]=ln;
        protfile.close();
        } 
      }
      nnumacum[i+1]=lnt; 
    }
  td=lnt;
  xx=new double[td];
  ye=new double[td];
  yi=new double[td*2];
  yy=new double[td];
  trckp=new int [td];
  for ( i=0 ; i < td ; i++ ) {
    yy[i]=0.0;
    xx[i]=xxt[i]; 
    ye[i]=yet[i]; 
    yi[2*i]=yit[i]; // auxiliary
    yi[2*i+1]=yit[i+10000]; //auxiliary
    }
  delete [] xxt;
  delete [] yet;
  delete [] yit;
  for ( i=0 ; i<nprot ; i++ ) { for ( j=nnumacum[i] ; j<nnumacum[i+1] ; j++ ) trckp[j]=i; }
// to add more protocols change here
  inqub ("Report Warnings (on/off) [off]: ",&warn,false);
  cout << endl;
  inqub ("Debug Errors (on/off) [off]: ",&debugerrors,false);
  cout << endl;
  inqub ("Convergence Analysis (on/off) [on]: ",&doconvanal,true);
  cout << endl;
  inqub ("Optimize (on/off) [on]: ",&optim,true);
  cout << endl;
  timestamp ();
  cout << endl;
  if (optim) {
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
  }
  cout << endl;
  cout << "Computing final function value ..." << endl << endl;
  funkp(dx, &n,"end");
  timestamp ();
  cout << endl << "  Number of Errors: " << cnterr << endl;
  if (doconvanal) convanal(&n,dx);
  cout << " The End\n" ;
  delete [] dx;
  return 0;
}

void convanal (int *nn, double dx[]) {
  int h,i,j,k,l,bb,xx;
  int nm=*nn;
  int nx=nm+nm*(nm+1)/2;
  int nb=nm*(nm-1)*4;
  int tnprot=nprot+1;
  bool tprot[tnprot];
  double x[nx];
  double b[nb],bt[nb];
  double a[nb*nx];
  double comp;
  double dxx[nm];
  double ffb[nb*tnprot];
  double ff0[tnprot]; 
  double df[nm],ddf[nm*nm];
  tprot[0]=true;
  for ( i=1 ; i<tnprot ; i++ ) tprot[i]=prot[i-1];

  cout << endl << "Doing Convergence Analysis ..." << endl;

  ff0[0]=funk(dx,&nm); // Total Function Value
  for ( h=1 ; h<tnprot ; h++ ) { if (tprot[h]) ff0[h]=ffg[h-1]; }  // Compute Partials Function Values to add more protocols change here
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
            for ( h=1 ; h<tnprot ; h++ ) { if (tprot[h]) ffb[bb+h*nb]=ffg[h-1]-ff0[h]; }
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
    if ( tprot[k] ) {
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
  //      cout << b[i] << "   " << bt[i] << "   " << bt[i]-b[i] << endl;
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
void inquprotprm(protparam *protprm) {
  inquf("    thold ["+num2str(protprm->thold)+"] = ",&protprm->thold, protprm->thold);
  inquf("    waittime ["+num2str(protprm->waittime)+"] = ",&protprm->waittime, protprm->waittime);
  inquf("    dt ["+num2str(protprm->dt)+"] = ",&protprm->dt, protprm->dt);
  inquf("    tvc ["+num2str(protprm->tvc)+"] = ",&protprm->tvc, protprm->tvc);
  inquf("    tdec ["+num2str(protprm->tdec)+"] = ",&protprm->tdec, protprm->tdec);
  inquf("    v1 ["+num2str(protprm->v1)+"] = ",&protprm->v1, protprm->v1);
  inquf("    v2 ["+num2str(protprm->v2)+"] = ",&protprm->v2, protprm->v2);
  inquf("    v3 ["+num2str(protprm->v3)+"] = ",&protprm->v3, protprm->v3);
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
  int i,ii;
  double imax,iimax=1.0,ffv=0.0;
  ii=protnum;
/*
  if (ii == 0) {
    // Find iMax
    imax=yy[nnumacum[ii]];
    for (i=nnumacum[ii]+1; i<nnumacum[ii+1]; i++ ) {
      if (fabs(yy[i]) > fabs(imax)) imax=yy[i];
      }
    iimax=1.0/imax;     }
*/
  if (ii < 2 || ii==3 || ii==6) {
      imax=yy[0];
      iimax=1.0/imax; }
  // Compute RMS
  for (i=nnumacum[ii]; i<nnumacum[ii+1]; i++ ) {
    if (ii == 5) ffv+=yy[i];
    else ffv+=pow(yy[i]*iimax-ye[i],2);
    }
  *ffval=ffv;
  }

double funk(double dx[], int *nn) {
  int i,n,j;
  double ffn[nprot],ff;
  bool error[td];
  n=*nn;
  ff=0.0;
  for (i=0 ; i<=nprot ; i++ ) ffg[i]=0.0;
  for (i=0; i<n ; i++ ) {
    if ( lowbb[i] ) { if ( dx[i] < lowbf[i] ) ffg[nprot]+=obp+obp*(lowbf[i]-dx[i]); }
    if ( upbb[i] )  { if ( dx[i] >  upbf[i] ) ffg[nprot]+=obp+obp*(dx[i]-upbf[i]); }
    }
  ff+=ffg[nprot];
  if (ffg[nprot]==0.0) {
    for (i=0 ; i<n ; i++ ) fldx[frdx[i]]=dx[i];
    for ( i=0 ; i<nprot ; i++ ) ffn[i]=0;
    #pragma omp parallel for default(shared) private(i,j)
    for ( i=0 ; i<td ; i++ ) {
      j=trckp[i];
      protocol(j,i,fldx,&xx[i],&yy[i],&error[i]);
      }
    for ( i=0 ; i<td ; i++) { if (error[i]) cnterr+=1; }
    for ( i=0 ; i<nprot ; i++ ) {
      if (prot[i]==true) {
        funcval(i, &ffn[i]);
	ffg[i]=ffn[i]/nnum[i]; 
        ff+=ffg[i];
        }
      }
    }
  return ff;
}

void funkp(double dx[], int *nn, string str) {
  int i,j;
  double ff;
  prn.active=true;
  prn.str=str;
  ff=funk(dx,nn);
  prn.active=false;
  for ( i=0 ; i<nprot ; i++ ) {
    if (prot[i]==true) {
      ofstream result_file((str+protname[i]+"_param.txt").c_str()   );    
      if ( i == 5 ) for (j=nnumacum[i]; j<nnumacum[i+1]; j++) result_file << xx[j] << "   " << yy[j] << "   " << 0.0 << endl;
      else if (i<2 || i==3 || i==6) for (j=nnumacum[i]; j<nnumacum[i+1]; j++) result_file << xx[j] << "   " << yy[j]/yy[nnumacum[i]] << "   " << ye[j] << endl;
      else          for (j=nnumacum[i]; j<nnumacum[i+1]; j++) result_file << xx[j] << "   " << yy[j] << "   " << ye[j] << endl;
      result_file.close();
      cout << "  f" << i << ":   " << ffg[i] << endl;    }
    } 
  cout << "  f" << nprot << ":   " << ffg[nprot] << endl; 
  cout << "  ff: " << ff << endl;
  return;
}

string num2str(double num) {
  stringstream ss (stringstream::in | stringstream::out);
  double var;
  var=num;
  ss << var;
  return ss.str();
}

void namefile(int protnum, double *valor, string *name) {
    double val;
    val=*valor;
    stringstream ss (stringstream::in | stringstream::out);
    ss << val;
    string test = ss.str();
    if (!prn.active) *name="";
    else *name=prn.str+test+"_curtime_"+protname[protnum]+".txt";
}

// Change here to add more protocols
void protocol(int protnum, int tit, double dx[], double *V, double *YY, bool *error) {
  if (protnum == 0) protocol0(dx,V,YY,error);
  else if (protnum == 1) protocol1(dx,V,YY,error);
  else if (protnum == 2) protocol2(dx,V,YY,error);
  else if (protnum == 3) protocol3(dx,V,YY,error);
  else if (protnum == 4) protocol4(dx,V,YY,error);
  else if (protnum == 5) protocol5(tit,dx,V,YY,error);
  else if (protnum == 6) protocol6(dx,V,YY,error);
  else cout << "Invalid protocol number" << endl;
}

void protocol0(double dx[], double *V,double *YY, bool *error)
{
  int tcy,p=0;
  int interval, waittimei;
  int tholdi,tvc,tdec;
  double t, dt, i_peak=0.0;
  threepoints ct;   
  const double t_hold = protprm[p].thold;  //P1 jinqing=-100
  const double waitTime = protprm[p].waittime;  //***************You MUST hold for 30,000ms to get block of late current!!!!
  bool peakfound; 
  Cell_param Cell;
  init_params(curr, &Cell);
  dt = protprm[p].dt;
  string name;
  namefile(p,V,&name);
  ofstream file(name.c_str());
  tholdi=(int)(t_hold/dt);
  tvc=(int)(protprm[p].tvc/dt); // P2 converted to number of steps Jinqing=50 variation to -100  P2 here is same as P1 in bereki protocol
  tdec=(int)(protprm[p].tdec/dt); // P3  final -100 mV
  interval=tvc + tdec + tholdi;
  waittimei=(int)(waitTime/dt);
  ct.t1=0.0;
  ct.t2=0.0;
  ct.t3=0.0;
  ct.i1=0.0;
  ct.i2=0.0;
  ct.i3=0.0;
  t=0.0;
// SSA Protocol---Jiqing
//              P2 1000ms (current measurment)   
//              ----
//             |----| 
//             |----|
//  rest----   |----|    ----rest 
//         |---|    |---|
//          P1        P3
  t = dt;
  // WaitTime 
  Cell.V = protprm[p].v1;  
  for (tcy = 0; tcy < waittimei; tcy++) 
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref1  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    t = t + dt;
    }
  // Potential Cycle Interval
  Cell.V = protprm[p].v2;  //P1
  tcy=0;
  peakfound=false;
  while (tcy<interval && (!peakfound || prn.active) )
    {
    if (tcy==tholdi) Cell.V = *V;   
    else if (tcy==tholdi+tvc) Cell.V = protprm[p].v3;  
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref2  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
     *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (tcy > tholdi+tvc+1 && !peakfound) {if ( (ct.i2 > ct.i1 && ct.i2 > ct.i3) || (ct.i2 < ct.i1 && ct.i2 < ct.i3) ) i_peak=ct.i2,peakfound=true; }  // find peak
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t+=dt;
    tcy+=1;
    }
//  if (warn && !peakfound) cout << "Warning: Peak not found in protocol0. V = " << Cell.V << endl;
  if (peakfound) *YY=i_peak;
  else *YY=ct.i3;
  if (prn.active) file.close();
  return;
} 

void protocol1(double dx[], double *V,double *YY, bool *error){
  int tcy,p=1;
  int interval, waittimei;
  int tholdi,tvc,tdec;
  double t, dt, i_peak=0.0;
  threepoints ct,cti;
  const double t_hold = protprm[p].thold;  
  const double waitTime = protprm[p].waittime;  //***************You MUST hold for 30,000ms to get block of late current!!!!
  bool peakfound; 
  Cell_param Cell;
  init_params(curr, &Cell);
  dt = protprm[p].dt;
  string name;
  namefile(p,V,&name);
  ofstream file(name.c_str());
  tholdi=(int)(t_hold/dt);
  tvc=(int)(protprm[p].tvc/dt); // 5ms converted to number of steps
  tdec=(int)(protprm[p].tdec/dt); // period of time when current is measured converted to nsteps
  interval=tvc + tdec + tholdi;
  waittimei=(int)(waitTime/dt);
  ct.t1=0.0;
  ct.t2=0.0;
  ct.t3=0.0;
  ct.i1=0.0;
  ct.i2=0.0;
  ct.i3=0.0;
  t=0.0;
/* SSI Protocol
            P1(1s)   
           -----     P3 (0.1s)(current measurement)
          |     |--|---| 
          |     |--|   |
  rest----|     |--|   |----rest 
                |--|
                 P2 (5ms)
*/               
  t = dt;
  // WaitTime 
  Cell.V = protprm[p].v1;
  for (tcy = 0; tcy < waittimei; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref3  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    t = t + dt;
    }
  // Potential Cycle Interval
  Cell.V = protprm[p].v2;  //P1
  peakfound=false;
  tcy=0;
  while (tcy<interval && !peakfound)
    {
    if (tcy==tholdi) {Cell.V = *V; } 
    else if (tcy==tholdi+tvc) {Cell.V = protprm[p].v3,cti=ct; }
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref4  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
     *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (tcy > tholdi+tvc+1 && !peakfound) {if ( (ct.i2 > ct.i1 && ct.i2 > ct.i3) || (ct.i2 < ct.i1 && ct.i2 < ct.i3) ) i_peak=ct.i2,peakfound=true; }  // find peak
    if (prn.active) { if (fmod(tcy,10000)==0) file << tcy << " , " << t << " , " << Cell.I_Kr << " ," << Cell.V << ", " << endl; }
    t+=dt;
    tcy+=1;
    }
  if (warn && !peakfound) cout << "Warning: Peak not found in protocol1" << endl;
  if (peakfound) *YY=i_peak;
  else *YY=ct.i3; 
  if (prn.active) file.close();
  return;
}

void protocol2(double dx[], double *V,double *YY, bool *error){
  int tcy,p=2;
  int interval, waittimei;
  int tholdi,tvc,tdec;
  double t, dt, tau, i_peak=0.0;
  threepoints ct,ct75,cti;
  const double t_hold = protprm[p].thold; //P1=+50
  const double waitTime = protprm[p].waittime;  //***************waitTime between secuence of pulses
  bool peakfound,peak75;
  Cell_param Cell;
  init_params(curr, &Cell);
  dt = protprm[p].dt;
  int abci[20][2],ii;
  float abcf[20][4];
  string name;
  namefile(p,V,&name);
  ofstream file(name.c_str());
  tholdi=(int)(t_hold/dt);
  tvc=(int)(protprm[p].tvc/dt); // 5ms converted to number of steps for P2=-100 (short pulse to recover from inactivation)
  tdec=(int)(protprm[p].tdec/dt); // period of time when current is measured converted to nsteps, P3=from -100 to +60 mV
  interval=tvc + tdec + tholdi;
  waittimei=(int)(waitTime/dt);
  ct.t1=0.0;
  ct.t2=0.0;
  ct.t3=0.0;
  ct.i1=0.0;
  ct.i2=0.0;
  ct.i3=0.0;
  ct75=ct;
  cti=ct;
  i_inf=0.0;
  i_init=0.0;
  tau=0.0;
/* Deactivation Protocol
            P1     P3 1000ms
           -----  ----
          |500ms||----| 
          |     ||----|
  rest----|     ||----|----rest 
                P2----
                5ms 
*/
  t = dt;
   // WaitTime 
  Cell.V = protprm[p].v1;
  for (tcy = 1; tcy < waittimei; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref5 " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    t = t + dt;
    }
  // Potential Cycle Interval
  Cell.V = protprm[p].v2; //P1
  peakfound=false;
  peak75=false;
  if (prn.active) ii=0;
  tcy=0;
  while ( tcy < interval )
    {
    if (tcy==tholdi) {Cell.V=protprm[p].v3;} //P2 
    else if (tcy==tholdi+tvc) {Cell.V =*V;} //P3 where current is measured
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref6 " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (tcy == tholdi+tvc+2) cti=ct;
    if (tcy > tholdi+tvc+1 && !peakfound) {if ( (ct.i2 > ct.i1 && ct.i2 > ct.i3) || (ct.i2 < ct.i1 && ct.i2 < ct.i3) ) i_peak=ct.i2,peakfound=true; }  // find peak
    if (tcy > tholdi+tvc+1 && peakfound && !peak75) { if (fabs(ct.i2) <= fabs(0.75*i_peak))  peak75=true,ct75=ct; } // save 0.75 of the peak
    if (prn.active) {
      if(fmod(tcy-50500,5000)==0 && Cell.V==-100 && tcy >= tholdi+tvc+2000){
        abci[ii][0]=tcy;
        abcf[ii][0]=t;
        abcf[ii][1]=Cell.I_Kr;
        abcf[ii][2]=Cell.V;
        ii++; }
      }
    t+=dt;
    tcy+=1;
    }
  if (warn && !peakfound) cout << "Warning: Peak not found in protocol3" << endl;
  if (!peak75) ct75=cti;
  calctau(&ct,&ct75,&tau);
  *YY=fabs(tau);
  if (prn.active) {
    for (tcy=0;tcy<ii;tcy++) { file << abci[tcy][0] << " , " << abcf[tcy][0] << " , " << abcf[tcy][1] << " ," << abcf[tcy][2] << endl; }
    }
  if ( isnan(*YY) ) {
    if (debugerrors) cout << "  Ref7 " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
    *error=true;
    *YY=maxyy;
    }
  if (prn.active) file.close();
  return;
}

void protocol3(double dx[], double *V,double *YY, bool *error)
{
  int tcy,p=3;
  int interval, waittimei;
  int tholdi,tvc,tdec;
  double t, dt, i_peak=0.0;
  threepoints ct,cti;
  const double t_hold = protprm[p].thold;  //P1 jinqing=-100
  const double waitTime = protprm[p].waittime;  //***************You MUST hold for 30,000ms to get block of late current!!!!
  bool peakfound;
  Cell_param Cell;
  init_params(curr, &Cell);
  dt = protprm[p].dt;
  string name;
  namefile(p,V,&name);
  ofstream file(name.c_str());
  tholdi=(int)(t_hold/dt);
  tvc=(int)(protprm[p].tvc/dt); // P2 converted to number of steps Jinqing=50 variation to -100  P2 here is same as P1 in bereki protocol
  tdec=(int)(protprm[p].tdec/dt); // P3  final -100 mV
  interval=tvc + tdec + tholdi;
  waittimei=(int)(waitTime/dt);
  Cell.Drug=*V;
  ct.t1=0.0;
  ct.t2=0.0;
  ct.t3=0.0;
  ct.i1=0.0;
  ct.i2=0.0;
  ct.i3=0.0;
  cti=ct;
// BLOCK
  t=0.0;
  t = dt;
  // WaitTime 
  Cell.V = protprm[p].v1;
  for (tcy = 0; tcy < waittimei; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref8  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    t = t + dt;
    }
  // Potential Cycle Interval
  Cell.V = protprm[p].v2;  //P1
  tcy=0;
  peakfound=false;
  while (tcy<interval && (!peakfound || prn.active) )
    {
    if (tcy==tholdi) Cell.V = protprm[p].v3;
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref9  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
     *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (tcy == tholdi+2) cti=ct;
    if (tcy > tholdi+1 && !peakfound) {if ( (ct.i2 > ct.i1 && ct.i2 > ct.i3) || (ct.i2 < ct.i1 && ct.i2 < ct.i3) ) i_peak=ct.i2,peakfound=true; }  // find peak
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t+=dt;
    tcy+=1;
    }
  if (warn && !peakfound) cout << "Warning: Peak not found in protocol0. V = " << Cell.V << endl;
  if (peakfound) *YY=i_peak;
  else *YY=cti.i3;
  if (prn.active) file.close();
  return;
}

void protocol4(double dx[], double *V,double *YY, bool *error){
  int tcy,p=4;
  int interval, waittimei;
  int tholdi,tvc,tdec;
  double t, dt, i_peak=0.0, i_50ms=0.0;
  threepoints ct;
  const double t_hold = protprm[p].thold; //P1=+50
  const double waitTime = protprm[p].waittime;  //***************waitTime between secuence of pulses
  bool peakfound=false,i50found=false;
  Cell_param Cell;
  init_params(curr, &Cell);
  dt = protprm[p].dt;  // units in ms 
  int tcypeak = 0;
  string name;
  namefile(p,V,&name);
  ofstream file(name.c_str());
  tholdi=(int)(t_hold/dt);
  tvc=(int)(protprm[p].tvc/dt); // 5ms converted to number of steps for P2=-100 (short pulse to recover from inactivation)
  tdec=(int)(protprm[p].tdec/dt); // period of time when current is measured converted to nsteps, P3=from -100 to +60 mV
  interval=tvc + tdec + tholdi;
  waittimei=(int)(waitTime/dt);
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
   // WaitTime 
  Cell.V = protprm[p].v1;
  for (tcy = 0; tcy < waittimei; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref10 " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    t = t + dt;
    }
  // Potential Cycle Interval
  Cell.V = protprm[p].v2; //P1
  tcy=0;
  peakfound=false;
  i50found=false;
  while (tcy<interval && !i50found) 
    {
    if (tcy==tholdi) Cell.V=protprm[p].v3; //P2 
    else if (tcy==tholdi+tvc) Cell.V=*V; //P3 where current is measured
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref11 " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if (tcy > tholdi+tvc+1 && !peakfound) {if ( (ct.i2 > ct.i1 && ct.i2 > ct.i3) || (ct.i2 < ct.i1 && ct.i2 < ct.i3) ) i_peak=ct.i2,peakfound=true,tcypeak=tcy-1; }  // find peak
    if (tcy > tholdi+tvc+1 && peakfound && !i50found) { if (tcy == tcypeak+5000) i50found=true,i_50ms=Cell.I_Kr; } // save 0.75 of the peak
    if (prn.active) { if(fmod(tcy,100)==0) file << t << " , " << Cell.I_Kr << " , " << Cell.V << " ," << *V << endl; }
    t+=dt;
    tcy+=1;
    }
  if (warn && !peakfound) cout << "Warning: Peak not found in protocol4" << endl;
  if (i50found) *YY=i_50ms/i_peak;
  else *YY=maxyy;
  if (prn.active) file.close();
  return;
}

void protocol5(int tit, double dx[], double *V, double *YY, bool *error){
  int tcy,p=5;
  int waittimei;
  int tholdi,tvc,tdec;
  double t0, t, dt, i_s, i_e, i_f, a, ar, tau, tau2, rms, ipeak=0.0, iat;
  threepoints ct;
  const double t_hold = protprm[p].thold; //P1=+50
  const double waitTime = protprm[p].waittime;  //***************waitTime between secuence of pulses
  bool peakfound,peak75;
  Cell_param Cell;
  init_params(curr, &Cell);
  dt = protprm[p].dt;
  int i,j,stp;
  double data[ndata*2];
  string name;
  namefile(p,V,&name);
  ofstream file(name.c_str());
  tholdi=(int)(t_hold/dt);
  tvc=(int)(protprm[p].tvc/dt); // 5ms converted to number of steps for P2=-100 (short pulse to recover from inactivation)
  tdec=(int)(protprm[p].tdec/dt); // period of time when current is measured converted to nsteps, P3=from -100 to +60 mV
  waittimei=(int)(waitTime/dt);
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
  for (tcy = 0; tcy < waittimei; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref12 " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    t = t + dt;
    }
  // Potential Cycle Interval
  Cell.V = protprm[p].v2; //P1
  for (tcy = 0; tcy < tholdi; tcy++ ) {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref13a " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    t+=dt;
    }
  Cell.V = protprm[p].v3;
  for ( tcy = 0; tcy < tvc ; tcy++ ); {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref13b " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    t+=dt;
    }
  peakfound=false;
  peak75=false;
  tcy=0;
  stp=tdec/ndata;
  j=0;
  Cell.V=*V;
  t0=t;
  while ( tcy < tdec && j < ndata ) {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref13c " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    if  (!peakfound && tcy > 1) {if ( (ct.i2 > ct.i1 && ct.i2 > ct.i3) || (ct.i2 < ct.i1 && ct.i2 < ct.i3) ) j=ndata-1,peakfound=true,ipeak=ct.i2; }  // find peak
    if (peakfound && !peak75) { if (fabs(ct.i2) <= fabs(ppeak*ipeak)) peak75=true,stp=(tdec-tcy)/ndata,j=0; } // save 0.75 of the peak
    if (tcy == tdec+(j-ndata)*stp) data[2*j]=t-t0, data[2*j+1]=Cell.I_Kr, j++;
    t+=dt;
    tcy+=1;
    }
  if (j < ndata) cout << "Error: Data not captured" << endl;
  if (!peak75) {
    if (debugerrors) cout << "Error: Peak not found in protocol5" << endl; 
    *YY=maxyy; }
  else {   
    a=0.0;
    ar=0.0;
    for (i=0;i<ndata;i++) a+=data[2*i+1]*(yi[2*tit+1]*exp(-data[2*i]/ye[tit])+(1.0-yi[2*tit+1])*exp(-data[2*i]/yi[2*tit]));
    for (i=0;i<ndata;i++) ar+=data[2*i+1]*data[2*i+1];
    iat=a/ar;
    *YY=0.0;
    for (i=0;i<ndata;i++) *YY+=pow(data[2*i+1]*iat-yi[2*tit+1]*exp(-data[2*i]/ye[tit])-(1.0-yi[2*tit+1])*exp(-data[2*i]/yi[2*tit]),2);
    *YY=*YY/ndata;
    }
  if ( isnan(*YY) ) {
    if (debugerrors) cout << "  Ref14 " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
    *error=true;
    *YY=maxyy;
    }
  if (prn.active) {
    calctau2(tit,ndata,data,ipeak,&tau,&tau2,&a,&ar);
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
    file << "   RMS(i_sim*iat vs i_exp) = " << *YY << "   tau1 = " << ye[tit] << "   tau2 = " << yi[2*tit] << "   a = " << 1.0 << "   ar = " << yi[2*tit+1] << "   a1 = " << yi[2*tit+1] << "   a2 = " << 1.0-yi[2*tit+1] << endl;
    file << endl;
    file << " LevMar Fitting: " << endl;
    file << "   RMS(i_sim/a vs i_fit/a) = " << rms << "   tau1 = " << tau << "   tau2 = " << tau2 << "   a = " << a << "   ar = " << ar << "   a1 = " << a*ar << "   a2 = " << a*(1.0-ar) << endl;
    file << endl;
    file << " t-t0  ,  i_exp/iat  ,  i_sim  ,  i_fit  ,  ,  i_exp  ,  i_sim*iat  ,  i_fit*iat  ,  ,  i_exp/iat/a  ,  i_sim/a  ,  i_fit/a" << endl;
    for (i=0;i<ndata;i++) {
       t   = data[2*i];
       i_s = data[2*i+1];
       i_e = yi[2*tit+1]*exp(-t/ye[tit])+(1.0-yi[2*tit+1])*exp(-t/yi[2*tit]);
       i_f = a*ar*exp(-t/tau)+a*(1.0-ar)*exp(-t/tau2);
       file << t << " , " << i_e/iat << " , " << i_s  << " , " << i_f  << " ,  ,  " << i_e << " , " << i_s*iat << " , " << i_f*iat << " ,  ,  " << i_e/a/iat << " , " << i_s/a << " , " << i_f/a << " , " << endl;
      }
    }
  if (prn.active) file.close();
  return;
}

void protocol6(double dx[], double *V,double *YY, bool *error)
{
  int tcy,p=6;
  int interval, waittimei;
  int tholdi,tvc,tdec,tvar;
  double t, dt, i_peak=0.0;
  threepoints ct;
  const double t_hold = protprm[p].thold;  //P1 jinqing=-100
  const double waitTime = protprm[p].waittime;  //***************You MUST hold for 30,000ms to get block of late current!!!!
  bool peakfound;
  Cell_param Cell;
  init_params(curr, &Cell);
  dt = protprm[p].dt;
  string name;
  namefile(p,V,&name);
  ofstream file(name.c_str());
  tholdi=(int)(t_hold/dt);
  tvc=(int)(protprm[p].tvc/dt); // P2 converted to number of steps Jinqing=50 variation to -100  P2 here is same as P1 in bereki protocol
  tdec=(int)(protprm[p].tdec/dt); // P3  final -100 mV
  interval=tvc + tdec + tholdi;
  waittimei=(int)(waitTime/dt);
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
  for (tcy = 0; tcy < waittimei; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref15  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    t = t + dt;
    }
  // Potential Cycle Interval
  Cell.V = protprm[p].v2;  //P1
  for (tcy = 0 ; tcy < tvar; tcy++)
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref16  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
      *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
    t = t + dt;
    }
  Cell.V = protprm[p].v3;
  tcy=0;
  peakfound=false;
  while (tcy<interval && (!peakfound || prn.active) )
    {
    Calculate_I_Kr(curr, &Cell, t, dt, dx, error);
    if ( *error ) {
      if (debugerrors) cout << "  Ref17  " << curr << "  " << t << "  " << Cell.I_Kr << "  " << Cell.V << "  " << *V << endl; //debug
     *YY=maxyy;
      return;
      }
    acum(&ct,Cell.I_Kr,t);
//    if (tcy==1) i_peak=ct.i2; 
//    if (tcy>0 && !peakfound) {if ( (ct.i2 > ct.i1 && ct.i2 > ct.i3) || (ct.i2 < ct.i1 && ct.i2 < ct.i3) ) i_peak=ct.i2,peakfound=true; }  // find peak
    if (tcy>0) {
      if ( (ct.i2 > 0.0 && ct.i2 > i_peak) || (ct.i2 < 0.0 && ct.i2 < i_peak) ) { i_peak=ct.i2; }
      else if ( (ct.i2 > 0.0 && ct.i2 < i_peak) || (ct.i2 < 0.0 && ct.i2 > i_peak) ) { peakfound=true; }  // find peak
      }
    if (prn.active) { if (fmod(tcy,1000)==0) file << tcy << " " << t << " " << Cell.I_Kr << " " << Cell.V << endl; }
    t+=dt;
    tcy+=1;
    }
  if (warn && !peakfound) cout << "Warning: Peak not found in protocol6. t = " << *V << endl;
  if (peakfound) *YY=i_peak;
  else *YY=ct.i3;
  if (prn.active) file.close();
  return;
}

void calctau(threepoints *pst, threepoints *pnd, double *tau){ 
  double ip1, ip2, t1, t2; 
  ip1=(pst->i3-pst->i1)/(pst->t3-pst->t1);
  ip2=(pnd->i3-pnd->i1)/(pnd->t3-pnd->t1);
  t1=pst->t2; 
  t2=pnd->t2; 
  *tau=(t2-t1)/log(ip1/ip2);
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

void calctau2(int tit, int n, double *dat, double ipeak, double *tau, double *tau2, double *a, double *ar){
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
  p[1]=yi[2*tit+1];
  p[2]=ye[tit];
  p[3]=yi[2*tit];
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
// there is an error
//void calciinit(threepoints *pst, threepoints *pnd, double tinit, double *iinit, double *iinf){
//  double k, ib, io, ip1, ip2, i1, i2, e1, e2, t1, t2,ie2e1;
//// i=(io-ib)*exp(-k*(t-t0)) + ib
//  ip1=(pst->i3-pst->i1)/(pst->t3-pst->t1);
//  ip2=(pnd->i3-pnd->i1)/(pnd->t3-pnd->t1);
//  t1=pst->t2;
//  t2=pnd->t2;
//  i1=pst->i2;
//  i2=pnd->i2;
//  k=-log(ip1/ip2)/(t1-t2);
//  e1=exp(-k*(t1-tinit));
//  e2=exp(-k*(t2-tinit));
//  ie2e1=1.0/(e2-e1);
//  ib=(i1*e2-i2*e1)*ie2e1;
//  io=(i2*(1.0-e1)-i1*(1.0-e2))*ie2e1;
//  cout << tinit << " " << ip1 << " " << ip2 << " " << t1 << " " << t2 << " " << i1 << " " << i2 << " " << k << " " << e1 << " " << e2 << " " << ib << " " << io << endl;
////  10505 -1.18644e-07 -0.0994513 11505 10505.3 0.01111 0.451418 0.0136425 1.18908e-06 0.996731 0.0111094 0.452862
//
//  *iinit=io;
//  *iinf=ib;
//  return;
//}


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
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0)); 	//var g_Kr_0: microS_per_nanoF {init: 0.024}; Lucia: Ko dependence is included in the formulation of the current
	
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

*/
/*
//Fink Rate constants parameters at 23 C Fit 16 params no Volt

      double ae=exp(dx[0]*24.2886+0.0117*Cell_ptr->V-dx[1]*27.1395);
      double be=exp(dx[2]*13.6416-0.0631*Cell_ptr->V-dx[3]*16.4494);
      double ain=exp(dx[4]*22.699-dx[5]*27.1395);
      double bin=exp(dx[6]*13.1473-dx[7]*16.4494);
      double ai=exp(dx[8]*30.0156-0.0327*Cell_ptr->V-dx[9]*34.8148);  // ai, bi steps different for the mutant.
      double bi=exp(dx[10]*29.9699+0.0233*Cell_ptr->V-dx[11]*32.3491)*pow((5.4/Ko),0.4);   //inactivation
      double aa=exp(dx[12]*22.0517+0.03822*Cell_ptr->V-dx[13]*27.1395);  //activation 
      double bb=exp(dx[14]*7.2663-0.0418*Cell_ptr->V-dx[15]*16.449);  //deactivation

//Fink Rate constants parameters at 23 C Fit 22 params no Volt
      double ae=exp(dx[0]*24.2886+dx[1]*0.0117*Cell_ptr->V-dx[2]*27.1395);
      double be=exp(dx[3]*13.6416-dx[4]*0.0631*Cell_ptr->V-dx[5]*16.4494);
      double ain=exp(dx[6]*22.699-dx[7]*27.1395);
      double bin=exp(dx[8]*13.1473-dx[9]*16.4494);
      double ai=exp(dx[10]*30.0156-dx[11]*0.0327*Cell_ptr->V-dx[12]*34.8148);  // ai, bi steps different for the mutant.
      double bi=exp(dx[13]*29.9699+dx[14]*0.0233*Cell_ptr->V-dx[15]*32.3491)*pow((5.4/Ko),0.4);   //inactivation
      double aa=exp(dx[16]*22.0517+dx[17]*0.03822*Cell_ptr->V-dx[18]*27.1395);  //activation 
      double bb=exp(dx[19]*7.2663-dx[20]*0.0418*Cell_ptr->V-dx[21]*16.449);  //deactivation


//Fink Rate constants parameters at 23 C Fit 14 params no Volt
      double ae=dx[0]*exp(24.2886+dx[1]*0.0117*Cell_ptr->V-27.1395);
      double be=dx[2]*exp(13.6416-dx[3]*0.0631*Cell_ptr->V-16.4494);
      double ain=dx[4]*exp(22.699-27.1395);
      double bin=dx[5]*exp(13.1473-16.4494);
      double ai=dx[6]*exp(30.0156-dx[7]*0.0327*Cell_ptr->V-34.8148);  // ai, bi steps different for the mutant.
      double bi=dx[8]*exp(29.9699+dx[9]*0.0233*Cell_ptr->V-32.3491)*pow((5.4/Ko),0.4);   //inactivation
      double aa=dx[10]*exp(22.0517+dx[11]*0.03822*Cell_ptr->V-27.1395);  //activation 
      double bb=dx[12]*exp(7.2663-dx[13]*0.0418*Cell_ptr->V-16.449);  //deactivation

//Rate constants parameters at 23 C corrected for Jiqing Data (NMS), and praxis opt but p1=1S
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
//Rate constants parameters at 23 C corrected for Jiqing Data (NMS), and praxis opt done but p1=1S
  double ae=dx[0]*1.50063*exp(24.288586641+dx[1]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*0.598723*exp(13.637047252-dx[3]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*1.10345*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.819984*exp(13.110155065-16.420979385);
  double ai=dx[6]*1.02159*exp(28.782084605-dx[7]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*0.998552*exp(29.690385428+dx[9]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*1.86823*exp(22.010640617+dx[11]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*1.00124*exp(7.266288229-dx[13]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation

//Rate constants parameters at 23 C corrected for L529I Data (NMS), and praxis opt done but p1=1S
  double ae=dx[0]*5.23371*exp(24.288586641+dx[1]*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*2.87555*exp(13.637047252-dx[3]*0.620343*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*0.6136*exp(22.444099764-27.130772208);
  double bin=dx[5]*1.22233*exp(13.110155065-16.420979385);
  double ai=dx[6]*1.88767*exp(28.782084605-dx[7]*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*1.54327*exp(29.690385428+dx[9]*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*7.98841*exp(22.010640617+dx[11]*1.10778*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*1.23079*exp(7.266288229-dx[13]*0.0418*Cell_ptr->V-16.168847047);  //deactivation


//Rate constants parameters at 23 C corrected for Jiqing Data (NMS), and praxis correction opt but p1=3S
  double ae=dx[0]*1.39429*exp(24.288586641+dx[1]*0.51998*0.0117*Cell_ptr->V-27.138936313); 
  double be=dx[2]*1.12436*exp(13.637047252-dx[3]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*0.80694*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.905579*exp(13.110155065-16.420979385);
  double ai=dx[6]*0.989603*exp(28.782084605-dx[7]*1.00013*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*1.00033*exp(29.690385428+dx[9]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*0.98968*exp(22.010640617+dx[11]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*0.99917*exp(7.266288229-dx[13]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation

//Fitting for NS + L529I (fit_7)
  double ae=1.43*dx[0]*exp(24.288586641+dx[1]*1.24*0.0117*Cell_ptr->V-27.138936313);
  double be=0.44*dx[2]*exp(13.637047252-dx[3]*0.0631*Cell_ptr->V-16.448205938);
  double ain=1.93*dx[4]*exp(22.444099764-27.130772208);
  double bin=1.99*dx[5]*exp(13.110155065-16.420979385);
  double ai=1.03*dx[6]*exp(28.782084605-0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=0.16*dx[7]*exp(29.690385428+1.50*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=1.49*dx[8]*exp(22.010640617+dx[9]*0.00095*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=1.25*dx[10]*exp(7.266288229-dx[11]*0.0418*Cell_ptr->V-16.168847047);  //deactivation
//ing for NS + WT
  double ae=dx[0]*exp(24.288586641+0.0117*Cell_ptr->V-27.138936313);
  double be=dx[1]*exp(13.637047252-0.0631*Cell_ptr->V-16.448205938);
  double ain=exp(22.444099764-27.130772208);
  double bin=exp(13.110155065-16.420979385);
  double ai=exp(28.782084605-0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=exp(29.690385428+0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[2]*exp(22.010640617+0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[3]*exp(7.266288229-0.0418*Cell_ptr->V-16.168847047);  //deactivation
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
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0));        //var g_Kr_0: microS_per_nanoF {init: 0.024}; Lucia: Ko dependence is included in the formulation of the current

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

//Rate constants parameters Silva & Rudy (2005)
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
  double E_Kr;
  double hh;
  hh = dt;
  
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0));       //var g_Kr_0: microS_per_nanoF {init: 0.024}; Lucia: Ko dependence is included in the formulation of the current

//        const double Q10=3;    To account for correction in temperature

//        const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));

  E_Kr = (R*T/F)*log(Ko/Ki);



//////Markov mechanism with NS drug: Option 1////////--lau--/////

//  C3 = C2 = C1 = O = I
//  |                 
//  C3D =C2D= C1D =OD= ID         

////////////////////////////////////////////////////
//rate cts. definitions:
//C3_C2=ae
//C2_C3=be
//C2_C1=ain
//C1_C2=bin
//C1_O=aa
//O_C1=bb
//O_I=bi
//I_O=ai
//ID_OD=aid
//OD_ID=bid
//C1D_OD=aad
//OD_CID=bbd
//C2D_C1D=aind
//C1D_C2D=bind
//C3D_C2D=aed
//C2D_C3D=bed
//*********Drug rates******
//C3_C3D=kC  ----> (microMolar*ms)-1)
//C3D_C3=rO  -----> (ms)-1)
//************************
//Drug equilibrium constants
//K_C=rC/kC
////////////////////////////

/*
//channel rates cts ----> 
//WT
      double ae=exp(24.288586641+0.0117*Cell_ptr->V-27.138936313); 
      double be=exp(13.637047252-0.0631*Cell_ptr->V-16.448205938);
      double ain=exp(22.444099764-27.130772208);
      double bin=exp(13.110155065-16.420979385);
      double ai=exp(28.782084605-0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
      double bi=exp(29.690385428+0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
      double aa=exp(22.010640617+0.03822*Cell_ptr->V-27.091932598);  //activation  
      double bb=exp(7.266288229-0.0418*Cell_ptr->V-16.168847047);  //deactivation
      double aed=dx[0]*exp(24.288586641+dx[1]*0.0117*Cell_ptr->V-27.138936313); 
      double bed=dx[2]*exp(13.637047252-0.0631*Cell_ptr->V-16.448205938);
      double aind=dx[3]*exp(22.444099764-27.130772208);
      double bind=dx[4]*exp(13.110155065-16.420979385);
      double aid=dx[5]*exp(28.782084605-dx[6]*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
      double bid=dx[7]*exp(29.690385428+dx[8]*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
      double aad=dx[9]*exp(22.010640617+dx[10]*0.03822*Cell_ptr->V-27.091932598);  //activation  
      double bbd=dx[11]*exp(7.266288229-0.0418*Cell_ptr->V-16.168847047);  //deactivation


//WT + drug (3S)

  double ae=dx[0]*1.39429*exp(24.288586641+dx[1]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*1.12436*exp(13.637047252-dx[3]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*0.80694*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.905579*exp(13.110155065-16.420979385);
  double ai=dx[6]*0.989603*exp(28.782084605-dx[7]*1.00013*0.0327*Cell_ptr->V-34.3968626); 
  double bi=dx[8]*1.00033*exp(29.690385428+dx[9]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*0.98968*exp(22.010640617+dx[11]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*0.99917*exp(7.266288229-dx[13]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.39429*exp(24.288586641+dx[15]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*1.12436*exp(13.637047252-dx[17]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*0.80694*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.905579*exp(13.110155065-16.420979385);
  double aid=dx[20]*0.989603*exp(28.782084605-dx[21]*1.00013*0.0327*Cell_ptr->V-34.3968626);  
  double bid=dx[22]*1.00033*exp(29.690385428+dx[23]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*0.98968*exp(22.010640617+dx[25]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*0.99917*exp(7.266288229-dx[27]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
*/
//WT + drug (1S)
  double ae=dx[0]*1.50063*exp(24.288586641+dx[1]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*0.598723*exp(13.637047252-dx[3]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*1.10345*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.819984*exp(13.110155065-16.420979385);
  double ai=dx[6]*1.02159*exp(28.782084605-dx[7]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*0.998552*exp(29.690385428+dx[9]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*1.86823*exp(22.010640617+dx[11]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*1.00124*exp(7.266288229-dx[13]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.50063*exp(24.288586641+dx[15]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*0.598723*exp(13.637047252-dx[17]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*1.10345*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.819984*exp(13.110155065-16.420979385);
  double aid=dx[20]*1.02159*exp(28.782084605-dx[21]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bid=dx[22]*0.998552*exp(29.690385428+dx[23]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*1.86823*exp(22.010640617+dx[25]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*1.00124*exp(7.266288229-dx[27]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation

  
  //drug rates cts ---->
  double kC=dx[28]; //microM-1s-1
  //double rC=dx[4];
  
  //relations based on K and microscopic reversibility
  double K_C=dx[29]; // (microM)//drug conscentration is 30 microM
  double rC=K_C*kC;  //s-1
  
  //double rC=kC*(rO*ae*ai*aa*bbd*bind*bed)/(kO*be*bin*bb*aad*aind*aed)
  
  //INitial conditions**********************************************************************

  // Calculation of k parameters for the Runge Kutta interations
  //Calculation of k1 values using the initial differential equation (set as f(t,V)
  double k1_O = hh*(ai*Cell_ptr->I + aa*Cell_ptr->C1 - Cell_ptr->O*(bi + bb));
  double k1_C1 = hh*(ain*Cell_ptr->C2 + bb*Cell_ptr->O - Cell_ptr->C1*(aa + bin));
  double k1_C2 = hh*(ae*Cell_ptr->C3 + bin*Cell_ptr->C1 - Cell_ptr->C2*(be + ain));
  double k1_C3 = hh*(be*Cell_ptr->C2 + rC*Cell_ptr->C3_drug - Cell_ptr->C3*(ae+kC*Cell_ptr->Drug));
  double k1_I = hh*(bi*Cell_ptr->O - ai*Cell_ptr->I);
  double k1_O_drug = hh*(aid*Cell_ptr->I_drug + aad*Cell_ptr->C1_drug - Cell_ptr->O_drug*(bid + bbd));
  double k1_C1_drug = hh*(aind*Cell_ptr->C2_drug + bbd*Cell_ptr->O_drug - Cell_ptr->C1_drug*(aad + bind));
  double k1_C2_drug = hh*(aed*Cell_ptr->C3_drug + bind*Cell_ptr->C1_drug - Cell_ptr->C2_drug*(bed + aind));
  double k1_C3_drug = hh*(bed*Cell_ptr->C2_drug +kC*Cell_ptr->Drug*Cell_ptr->C3 - Cell_ptr->C3_drug*(aed + rC ));
  double k1_I_drug= hh*(bid*Cell_ptr->O_drug - aid*Cell_ptr->I_drug);
  
  
  //Value of O after the k1 iteration; used for calculaton of k2
  double O_k1 = Cell_ptr->O + k1_O/2;
  double C1_k1 = Cell_ptr->C1 + k1_C1/2;
  double C2_k1 = Cell_ptr->C2 + k1_C2/2;
  double C3_k1 = Cell_ptr->C3 + k1_C3/2;
  double I_k1 = Cell_ptr->I + k1_I/2;
  double O_drug_k1 = Cell_ptr->O_drug + k1_O_drug/2;
  double C1_drug_k1 = Cell_ptr->C1_drug + k1_C1_drug/2;
  double C2_drug_k1 = Cell_ptr->C2_drug + k1_C2_drug/2;
  double C3_drug_k1 = Cell_ptr->C3_drug + k1_C3_drug/2;
  double I_drug_k1 = Cell_ptr->I_drug + k1_I_drug/2;
  
  
  //calculation of K2 values
  double k2_O = hh*(ai*I_k1 + aa*C1_k1 - O_k1*(bi + bb));
  double k2_C1 = hh*(ain*C2_k1 + bb*O_k1 - C1_k1*(aa + bin));
  double k2_C2 = hh*(ae*C3_k1 + bin*C1_k1 - C2_k1*(be + ain));
  double k2_C3 = hh*(be*C2_k1 +rC*C3_drug_k1 - C3_k1*(ae+ kC*Cell_ptr->Drug));
  double k2_I = hh*(bi*O_k1 - ai*I_k1);
  double k2_O_drug = hh*(aid*I_drug_k1 + aad*C1_drug_k1 - O_drug_k1*(bid + bbd));
  double k2_C1_drug = hh*(aind*C2_drug_k1 + bbd*O_drug_k1 - C1_drug_k1*(aad + bind));
  double k2_C2_drug = hh*(aed*C3_drug_k1 + bind*C1_drug_k1 - C2_drug_k1*(bed + aind));
  double k2_C3_drug = hh*(bed*C2_drug_k1 + kC*Cell_ptr->Drug*C3_k1 - C3_drug_k1*(aed + rC));
  double k2_I_drug= hh*(bid*O_drug_k1 - aid*I_drug_k1);
  
  
  double O_k2 = Cell_ptr->O + k2_O/2;
  double C1_k2 = Cell_ptr->C1 + k2_C1/2;
  double C2_k2 = Cell_ptr->C2 + k2_C2/2;
  double C3_k2 = Cell_ptr->C3 + k2_C3/2;
  double I_k2 = Cell_ptr->I + k2_I/2;
  double O_drug_k2 = Cell_ptr->O_drug + k2_O_drug/2;
  double C1_drug_k2 = Cell_ptr->C1_drug + k2_C1_drug/2;
  double C2_drug_k2 = Cell_ptr->C2_drug + k2_C2_drug/2;
  double C3_drug_k2 = Cell_ptr->C3_drug + k2_C3_drug/2;
  double I_drug_k2 = Cell_ptr->I_drug + k2_I_drug/2;
  
  //calculation of K3 values
  
  double k3_O = hh*(ai*I_k2 + aa*C1_k2 - O_k2*(bi + bb));
  double k3_C1 = hh*(ain*C2_k2 + bb*O_k2 - C1_k2*(aa + bin));
  double k3_C2 = hh*(ae*C3_k2 + bin*C1_k2 - C2_k2*(be + ain));
  double k3_C3 = hh*(be*C2_k2 +rC*C3_drug_k2 - C3_k2*(ae+ kC*Cell_ptr->Drug));
  double k3_I = hh*(bi*O_k2 - ai*I_k2);
  double k3_O_drug = hh*(aid*I_drug_k2 + aad*C1_drug_k2 - O_drug_k2*(bid + bbd));
  double k3_C1_drug = hh*(aind*C2_drug_k2 + bbd*O_drug_k2 - C1_drug_k2*(aad + bind));
  double k3_C2_drug = hh*(aed*C3_drug_k2 + bind*C1_drug_k2 - C2_drug_k2*(bed + aind));
  double k3_C3_drug = hh*(bed*C2_drug_k2 + kC*Cell_ptr->Drug*C3_k2 - C3_drug_k2*(aed + rC));
  double k3_I_drug= hh*(bid*O_drug_k2 - aid*I_drug_k2);
  
  double O_k3 = Cell_ptr->O + k3_O;
  double C1_k3 = Cell_ptr->C1 + k3_C1;
  double C2_k3 = Cell_ptr->C2 + k3_C2;
  double C3_k3 = Cell_ptr->C3 + k3_C3;
  double I_k3 = Cell_ptr->I + k3_I;
  double O_drug_k3 = Cell_ptr->O_drug + k3_O_drug;
  double C1_drug_k3 = Cell_ptr->C1_drug + k3_C1_drug;
  double C2_drug_k3 = Cell_ptr->C2_drug + k3_C2_drug;
  double C3_drug_k3 = Cell_ptr->C3_drug + k3_C3_drug;
  double I_drug_k3 = Cell_ptr->I_drug + k3_I_drug;
  
  //calculation of K4 values
  
  double k4_O = hh*(ai*I_k3 + aa*C1_k3 - O_k3*(bi + bb));
  double k4_C1 = hh*(ain*C2_k3 + bb*O_k3 - C1_k3*(aa + bin));
  double k4_C2 = hh*(ae*C3_k3 + bin*C1_k3 - C2_k3*(be + ain));
  double k4_C3 = hh*(be*C2_k3 +rC*C3_drug_k3 - C3_k3*(ae+ kC*Cell_ptr->Drug));
  double k4_I = hh*(bi*O_k3 - ai*I_k3);
  double k4_O_drug = hh*(aid*I_drug_k3 + aad*C1_drug_k3 - O_drug_k3*(bid + bbd));
  double k4_C1_drug = hh*(aind*C2_drug_k3 + bbd*O_drug_k3 - C1_drug_k3*(aad + bind));
  double k4_C2_drug = hh*(aed*C3_drug_k3 + bind*C1_drug_k3 - C2_drug_k3*(bed + aind));
  double k4_C3_drug = hh*(bed*C2_drug_k3 + kC*Cell_ptr->Drug*C3_k3 - C3_drug_k3*(aed + rC));
  double k4_I_drug= hh*(bid*O_drug_k3 - aid*I_drug_k3);
  
  //New computed values of the channel probabilities
  Cell_ptr->O = Cell_ptr->O + (k1_O + 2*k2_O + 2*k3_O + k4_O)/6;
  Cell_ptr->C1 = Cell_ptr->C1 + (k1_C1 + 2*k2_C1 + 2*k3_C1 + k4_C1)/6;
  Cell_ptr->C2 = Cell_ptr->C2 + (k1_C2 + 2*k2_C2 + 2*k3_C2 + k4_C2)/6;
  Cell_ptr->C3 = Cell_ptr->C3 + (k1_C3 + 2*k2_C3 + 2*k3_C3 + k4_C3)/6;
  Cell_ptr->I = Cell_ptr->I + (k1_I + 2*k2_I + 2*k3_I + k4_I)/6;
  Cell_ptr->O_drug = Cell_ptr->O_drug + (k1_O_drug + 2*k2_O_drug + 2*k3_O_drug + k4_O_drug)/6;
  Cell_ptr->C1_drug = Cell_ptr->C1_drug + (k1_C1_drug + 2*k2_C1_drug + 2*k3_C1_drug + k4_C1_drug)/6;
  Cell_ptr->C2_drug = Cell_ptr->C2_drug + (k1_C2_drug + 2*k2_C2_drug + 2*k3_C2_drug + k4_C2_drug)/6;
  Cell_ptr->C3_drug = Cell_ptr->C3_drug + (k1_C3_drug + 2*k2_C3_drug + 2*k3_C3_drug + k4_C3_drug)/6;
  Cell_ptr->I_drug = Cell_ptr->I_drug + (k1_I_drug + 2*k2_I_drug + 2*k3_I_drug + k4_I_drug)/6;
  
  Cell_ptr->sum = Cell_ptr->O + Cell_ptr->C1 + Cell_ptr->C2 + Cell_ptr->C3 + Cell_ptr->I+ Cell_ptr->O_drug + Cell_ptr->C1_drug + Cell_ptr->C2_drug + Cell_ptr->C3_drug + Cell_ptr->I_drug;

  if (fabs (Cell_ptr->sum -1.0) > 0.0001) {
//    cout << "Error in Wild Type Sum at t = " << t << ",\t sum = " << Cell_ptr->sum << endl;
    *error=true;
    Cell_ptr->I_Kr = 0.0;
    }
  else {
    *error=false;
    Cell_ptr->I_Kr = G_Kr*sqrt(Ko/5.4)*(Cell_ptr->O + Cell_ptr->O_drug)*(Cell_ptr->V - E_Kr);
    }
  return;
}

void curr4(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){

  double E_Kr;
  double hh;
  hh = dt;
  
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0));       //var g_Kr_0: microS_per_nanoF {init: 0.024}; Lucia: Ko dependence is included in the formulation of the current

//        const double Q10=3;    To account for correction in temperature
//        const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));

  E_Kr = (R*T/F)*log(Ko/Ki);

//////Markov mechanism with NS drug////////Option2 --lau--/////

//  C3 = C2 = C1 = O = I
//  |    |    |    |   
// C3D= C2D= C1D = OD= ID

////////////////////////////////////////////////////

//rate cts. definitions:
//C3_C2=ae
//C2_C3=be
//C2_C1=ain
//C1_C2=bin
//C1_O=aa
//O_C1=bb
//O_I=bi
//I_O=ai
//ID_OD=aid
//OD_ID=bid
//C1D_OD=aad
//OD_CID=bbd
//C2D_C1D=aind
//C1D_C2D=bind
//C3D_C2D=aed
//C2D_C3D=bed
//*********Drug rates******
//O_OD=kO  ----> (microMolar*ms)-1) 
//OD_O=rO  -----> (ms)-1)
//C3_C3DD=kC  ----> (microMolar*ms)-1)
//C3DD_C3=rO  -----> (ms)-1)
//************************
//Drug equilibrium constants
//K_O=rO/kO
//K_C=rC/kC
////////////////////////////

//channel rates cts ----> 


/*
//WT + drug (3S)
//
  double ae=dx[0]*1.39429*exp(24.288586641+dx[1]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*1.12436*exp(13.637047252-dx[3]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*0.80694*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.905579*exp(13.110155065-16.420979385);
  double ai=dx[6]*0.989603*exp(28.782084605-dx[7]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bi=dx[8]*1.00033*exp(29.690385428+dx[9]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*0.98968*exp(22.010640617+dx[11]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*0.99917*exp(7.266288229-dx[13]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.39429*exp(24.288586641+dx[15]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*1.12436*exp(13.637047252-dx[17]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*0.80694*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.905579*exp(13.110155065-16.420979385);
  double aid=dx[20]*0.989603*exp(28.782084605-dx[21]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bid=dx[22]*1.00033*exp(29.690385428+dx[23]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*0.98968*exp(22.010640617+dx[25]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*0.99917*exp(7.266288229-dx[27]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
*/
//wt + drug (1S)
  double ae=dx[0]*1.50063*exp(24.288586641+dx[1]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*0.598723*exp(13.637047252-dx[3]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*1.10345*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.819984*exp(13.110155065-16.420979385);
  double ai=dx[6]*1.02159*exp(28.782084605-dx[7]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*0.998552*exp(29.690385428+dx[9]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*1.86823*exp(22.010640617+dx[11]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*1.00124*exp(7.266288229-dx[13]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.50063*exp(24.288586641+dx[15]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*0.598723*exp(13.637047252-dx[17]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*1.10345*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.819984*exp(13.110155065-16.420979385);
  double aid=dx[20]*1.02159*exp(28.782084605-dx[21]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bid=dx[22]*0.998552*exp(29.690385428+dx[23]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*1.86823*exp(22.010640617+dx[25]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*1.00124*exp(7.266288229-dx[27]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation


//drug rates cts ---->
  double kO=dx[28];
//double rO=dx[1];
  double kC=dx[29];
//double rC=dx[4];

//relations based on K and microscopic reversibility
//  double K_C=10.0; // (microM)//drug conscentration is 30 microM
//  double K_O=3.0; //(microM)

  double K_C=dx[30]; // (microM)//drug conscentration is 30 microM
  double K_O=dx[31]; //(microM)

  double rC=K_C*kC;
  double rO=K_O*kO;
//double rC=kC*(rO*ae*ai*aa*bbd*bind*bed)/(kO*be*bin*bb*aad*aind*aed)

//INitial conditions**********************************************************************

  // Calculation of k parameters for the Runge Kutta interations
  //Calculation of k1 values using the initial differential equation (set as f(t,V)
  double k1_O = hh*(ai*Cell_ptr->I + aa*Cell_ptr->C1+ rO*Cell_ptr->O_drug - Cell_ptr->O*(bi + bb + kO*Cell_ptr->Drug));
  double k1_C1 = hh*(ain*Cell_ptr->C2 + bb*Cell_ptr->O + rC*Cell_ptr->C1_drug - Cell_ptr->C1*(aa + bin + kC*Cell_ptr->Drug));
  double k1_C2 = hh*(ae*Cell_ptr->C3 + bin*Cell_ptr->C1 + rC*Cell_ptr->C2_drug - Cell_ptr->C2*(be + ain + kC*Cell_ptr->Drug));
  double k1_C3 = hh*(be*Cell_ptr->C2 + rC*Cell_ptr->C3_drug - Cell_ptr->C3*(ae+kC*Cell_ptr->Drug));
  double k1_I = hh*(bi*Cell_ptr->O - ai*Cell_ptr->I);
  double k1_O_drug = hh*( aid*Cell_ptr->I_drug + aad*Cell_ptr->C1_drug + kO*Cell_ptr->Drug*Cell_ptr->O - Cell_ptr->O_drug*(bid + bbd+rO));  
  double k1_C1_drug = hh*(aind*Cell_ptr->C2_drug + bbd*Cell_ptr->O_drug + kC*Cell_ptr->Drug*Cell_ptr->C1 - Cell_ptr->C1_drug*(aad + bind +rC));
  double k1_C2_drug = hh*(aed*Cell_ptr->C3_drug + bind*Cell_ptr->C1_drug + kC*Cell_ptr->Drug*Cell_ptr->C2 - Cell_ptr->C2_drug*(bed + aind + rC));
  double k1_C3_drug = hh*(bed*Cell_ptr->C2_drug +kC*Cell_ptr->Drug*Cell_ptr->C3 - Cell_ptr->C3_drug*(aed + rC ));
  double k1_I_drug= hh*(bid*Cell_ptr->O_drug - aid*Cell_ptr->I_drug);

  //Value of O after the k1 iteration; used for calculaton of k2
  double O_k1 = Cell_ptr->O + k1_O/2;
  double C1_k1 = Cell_ptr->C1 + k1_C1/2;
  double C2_k1 = Cell_ptr->C2 + k1_C2/2;
  double C3_k1 = Cell_ptr->C3 + k1_C3/2;
  double I_k1 = Cell_ptr->I + k1_I/2;
  double O_drug_k1 = Cell_ptr->O_drug + k1_O_drug/2;
  double C1_drug_k1 = Cell_ptr->C1_drug + k1_C1_drug/2;
  double C2_drug_k1 = Cell_ptr->C2_drug + k1_C2_drug/2;
  double C3_drug_k1 = Cell_ptr->C3_drug + k1_C3_drug/2;
  double I_drug_k1 = Cell_ptr->I_drug + k1_I_drug/2;
  //calculation of K2 values
  double k2_O = hh*(ai*I_k1 + aa*C1_k1 + rO*O_drug_k1 - O_k1*(bi + bb+ kO*Cell_ptr->Drug));
  double k2_C1 = hh*(ain*C2_k1 + bb*O_k1 +rC*C1_drug_k1 - C1_k1*(aa + bin + kC*Cell_ptr->Drug));
  double k2_C2 = hh*(ae*C3_k1 + bin*C1_k1 +rC*C2_drug_k1 - C2_k1*(be + ain + kC*Cell_ptr->Drug));
  double k2_C3 = hh*(be*C2_k1 +rC*C3_drug_k1 - C3_k1*(ae+ kC*Cell_ptr->Drug));
  double k2_I = hh*(bi*O_k1 - ai*I_k1);
  double k2_O_drug = hh*( aid*I_drug_k1 + aad*C1_drug_k1 + kO*Cell_ptr->Drug*O_k1 - O_drug_k1*(bid + bbd + rO)); 
  double k2_C1_drug = hh*(aind*C2_drug_k1 + bbd*O_drug_k1 + kC*Cell_ptr->Drug*C1_k1 - C1_drug_k1*(aad + bind + rC));
  double k2_C2_drug = hh*(aed*C3_drug_k1 + bind*C1_drug_k1 + kC*Cell_ptr->Drug*C2_k1 - C2_drug_k1*(bed + aind + rC));
  double k2_C3_drug = hh*(bed*C2_drug_k1 + kC*Cell_ptr->Drug*C3_k1 - C3_drug_k1*(aed + rC));
  double k2_I_drug= hh*(bid*O_drug_k1 - aid*I_drug_k1); 
  double O_k2 = Cell_ptr->O + k2_O/2;
  double C1_k2 = Cell_ptr->C1 + k2_C1/2;
  double C2_k2 = Cell_ptr->C2 + k2_C2/2;
  double C3_k2 = Cell_ptr->C3 + k2_C3/2;
  double I_k2 = Cell_ptr->I + k2_I/2;
  double O_drug_k2 = Cell_ptr->O_drug + k2_O_drug/2;
  double C1_drug_k2 = Cell_ptr->C1_drug + k2_C1_drug/2;
  double C2_drug_k2 = Cell_ptr->C2_drug + k2_C2_drug/2;
  double C3_drug_k2 = Cell_ptr->C3_drug + k2_C3_drug/2;
  double I_drug_k2 = Cell_ptr->I_drug + k2_I_drug/2;
  //calculation of K3 values
  double k3_O = hh*(ai*I_k2 + aa*C1_k2 + rO*O_drug_k2 - O_k2*(bi + bb+ kO*Cell_ptr->Drug));
  double k3_C1 = hh*(ain*C2_k2 + bb*O_k2 +rC*C1_drug_k2 - C1_k2*(aa + bin + kC*Cell_ptr->Drug));
  double k3_C2 = hh*(ae*C3_k2 + bin*C1_k2 +rC*C2_drug_k2 - C2_k2*(be + ain + kC*Cell_ptr->Drug));
  double k3_C3 = hh*(be*C2_k2 +rC*C3_drug_k2 - C3_k2*(ae+ kC*Cell_ptr->Drug));
  double k3_I = hh*(bi*O_k2 - ai*I_k2);
  double k3_O_drug = hh*( aid*I_drug_k2 + aad*C1_drug_k2 + kO*Cell_ptr->Drug*O_k2 - O_drug_k2*(bid + bbd + rO));
  double k3_C1_drug = hh*(aind*C2_drug_k2 + bbd*O_drug_k2 + kC*Cell_ptr->Drug*C1_k2 - C1_drug_k2*(aad + bind + rC));
  double k3_C2_drug = hh*(aed*C3_drug_k2 + bind*C1_drug_k2 + kC*Cell_ptr->Drug*C2_k2 - C2_drug_k2*(bed + aind + rC));
  double k3_C3_drug = hh*(bed*C2_drug_k2 + kC*Cell_ptr->Drug*C3_k2 - C3_drug_k2*(aed + rC));
  double k3_I_drug= hh*(bid*O_drug_k2 - aid*I_drug_k2);
  double O_k3 = Cell_ptr->O + k3_O;
  double C1_k3 = Cell_ptr->C1 + k3_C1;
  double C2_k3 = Cell_ptr->C2 + k3_C2;
  double C3_k3 = Cell_ptr->C3 + k3_C3;
  double I_k3 = Cell_ptr->I + k3_I;
  double O_drug_k3 = Cell_ptr->O_drug + k3_O_drug;
  double C1_drug_k3 = Cell_ptr->C1_drug + k3_C1_drug;
  double C2_drug_k3 = Cell_ptr->C2_drug + k3_C2_drug;
  double C3_drug_k3 = Cell_ptr->C3_drug + k3_C3_drug;
  double I_drug_k3 = Cell_ptr->I_drug + k3_I_drug;
  //calculation of K4 values
  double k4_O = hh*(ai*I_k3 + aa*C1_k3 + rO*O_drug_k3 - O_k3*(bi + bb+ kO*Cell_ptr->Drug));
  double k4_C1 = hh*(ain*C2_k3 + bb*O_k3 +rC*C1_drug_k3 - C1_k3*(aa + bin +  kC*Cell_ptr->Drug));
  double k4_C2 = hh*(ae*C3_k3 + bin*C1_k3 +rC*C2_drug_k3 - C2_k3*(be + ain +  kC*Cell_ptr->Drug));
  double k4_C3 = hh*(be*C2_k3 +rC*C3_drug_k3 - C3_k3*(ae+ kC*Cell_ptr->Drug));
  double k4_I = hh*(bi*O_k3 - ai*I_k3);
  double k4_O_drug = hh*( aid*I_drug_k3 + aad*C1_drug_k3 + kO*Cell_ptr->Drug*O_k3 - O_drug_k3*(bid + bbd + rO));
  double k4_C1_drug = hh*(aind*C2_drug_k3 + bbd*O_drug_k3 + kC*Cell_ptr->Drug*C1_k3 - C1_drug_k3*(aad + bind + rC));
  double k4_C2_drug = hh*(aed*C3_drug_k3 + bind*C1_drug_k3 + kC*Cell_ptr->Drug*C2_k3 - C2_drug_k3*(bed + aind + rC));
  double k4_C3_drug = hh*(bed*C2_drug_k3 + kC*Cell_ptr->Drug*C3_k3 - C3_drug_k3*(aed + rC));
  double k4_I_drug= hh*(bid*O_drug_k3 - aid*I_drug_k3);
  
  
  //New computed values of the channel probabilities
  Cell_ptr->O = Cell_ptr->O + (k1_O + 2*k2_O + 2*k3_O + k4_O)/6;
  Cell_ptr->C1 = Cell_ptr->C1 + (k1_C1 + 2*k2_C1 + 2*k3_C1 + k4_C1)/6;
  Cell_ptr->C2 = Cell_ptr->C2 + (k1_C2 + 2*k2_C2 + 2*k3_C2 + k4_C2)/6;
  Cell_ptr->C3 = Cell_ptr->C3 + (k1_C3 + 2*k2_C3 + 2*k3_C3 + k4_C3)/6;
  Cell_ptr->I = Cell_ptr->I + (k1_I + 2*k2_I + 2*k3_I + k4_I)/6;
  Cell_ptr->O_drug = Cell_ptr->O_drug + (k1_O_drug + 2*k2_O_drug + 2*k3_O_drug + k4_O_drug)/6;
  Cell_ptr->C1_drug = Cell_ptr->C1_drug + (k1_C1_drug + 2*k2_C1_drug + 2*k3_C1_drug + k4_C1_drug)/6;
  Cell_ptr->C2_drug = Cell_ptr->C2_drug + (k1_C2_drug + 2*k2_C2_drug + 2*k3_C2_drug + k4_C2_drug)/6;
  Cell_ptr->C3_drug = Cell_ptr->C3_drug + (k1_C3_drug + 2*k2_C3_drug + 2*k3_C3_drug + k4_C3_drug)/6;
  Cell_ptr->I_drug = Cell_ptr->I_drug + (k1_I_drug + 2*k2_I_drug + 2*k3_I_drug + k4_I_drug)/6;
  
  Cell_ptr->sum = Cell_ptr->O + Cell_ptr->C1 + Cell_ptr->C2 + Cell_ptr->C3 + Cell_ptr->I+ Cell_ptr->O_drug + Cell_ptr->C1_drug + Cell_ptr->C2_drug + Cell_ptr->C3_drug + Cell_ptr->I_drug;


  if (fabs (Cell_ptr->sum -1.0) > 0.0001) {
//    cout << "Error in Wild Type Sum at t = " << t << ",\t sum = " << Cell_ptr->sum << endl;
    *error=true;
    Cell_ptr->I_Kr = 0.0;
    }
  else {
    *error=false;
    Cell_ptr->I_Kr = G_Kr*sqrt(Ko/5.4)*(Cell_ptr->O + Cell_ptr->O_drug)*(Cell_ptr->V - E_Kr);
    }
  return;
}

void curr5(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){

  double E_Kr;
  double hh;
  hh = dt;
  
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0));       //var g_Kr_0: microS_per_nanoF {init: 0.024}; Lucia: Ko dependence is included in the formulation of the current

//        const double Q10=3;    To account for correction in temperature

//        const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));

  E_Kr = (R*T/F)*log(Ko/Ki);



//////Markov mechanism with NS drug Option 3////////--lau--option1/////

//  C3 = C2 = C1 = O = I
//  |              |   
//  C3D =C2D= C1D =OD= ID         

////////////////////////////////////////////////////
//rate cts. definitions:
//C3_C2=ae
//C2_C3=be
//C2_C1=ain
//C1_C2=bin
//C1_O=aa
//O_C1=bb
//O_I=bi
//I_O=ai
//ID_OD=aid
//OD_ID=bid
//C1D_OD=aad
//OD_CID=bbd
//C2D_C1D=aind
//C1D_C2D=bind
//C3D_C2D=aed
//C2D_C3D=bed
//*********Drug rates******
//O_OD=kO  ----> (microMolar*ms)-1) 
//OD_O=rO  -----> (ms)-1)
//C3_C3DD=kC  ----> (microMolar*ms)-1)
//C3DD_C3=rO  -----> (ms)-1)
//************************
//Drug equilibrium constants
//K_O=rO/kO
//K_C=rC/kC
////////////////////////////

//channel rates cts ---->
/* 
//WT + drug (3S)

  double ae=dx[0]*1.39429*exp(24.288586641+dx[1]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*1.12436*exp(13.637047252-dx[3]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*0.80694*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.905579*exp(13.110155065-16.420979385);
  double ai=dx[6]*0.989603*exp(28.782084605-dx[7]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bi=dx[8]*1.00033*exp(29.690385428+dx[9]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*0.98968*exp(22.010640617+dx[11]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*0.99917*exp(7.266288229-dx[13]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.39429*exp(24.288586641+dx[15]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*1.12436*exp(13.637047252-dx[17]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*0.80694*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.905579*exp(13.110155065-16.420979385);
  double aid=dx[20]*0.989603*exp(28.782084605-dx[21]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bid=dx[22]*1.00033*exp(29.690385428+dx[23]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*0.98968*exp(22.010640617+dx[25]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*0.99917*exp(7.266288229-dx[27]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
*/
//WT+drug (1S)
  double ae=dx[0]*1.50063*exp(24.288586641+dx[1]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*0.598723*exp(13.637047252-dx[3]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*1.10345*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.819984*exp(13.110155065-16.420979385);
  double ai=dx[6]*1.02159*exp(28.782084605-dx[7]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*0.998552*exp(29.690385428+dx[9]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*1.86823*exp(22.010640617+dx[11]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*1.00124*exp(7.266288229-dx[13]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.50063*exp(24.288586641+dx[15]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*0.598723*exp(13.637047252-dx[17]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*1.10345*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.819984*exp(13.110155065-16.420979385);
  double aid=dx[20]*1.02159*exp(28.782084605-dx[21]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bid=dx[22]*0.998552*exp(29.690385428+dx[23]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*1.86823*exp(22.010640617+dx[25]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*1.00124*exp(7.266288229-dx[27]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation


//drug rates cts ---->
  double kO=dx[28];
//double rO=dx[1];
  double kC=dx[29];
//double rC=dx[4];


//relations based on K and microscopic reversibility
//double K_C=10.0; // (microM)//drug conscentration is 30 microM
//double K_O=3.0; //(microM)
  double K_C=dx[30]; // (microM)//drug conscentration is 30 microM
  double K_O=dx[31]; //(microM)
  double rC=K_C*kC;
  double rO=K_O*kO;


//double rC=kC*(rO*ae*ai*aa*bbd*bind*bed)/(kO*be*bin*bb*aad*aind*aed)



//INitial conditions**********************************************************************

  // Calculation of k parameters for the Runge Kutta interations
  //Calculation of k1 values using the initial differential equation (set as f(t,V)
  double k1_O = hh*(ai*Cell_ptr->I + aa*Cell_ptr->C1+ rO*Cell_ptr->O_drug - Cell_ptr->O*(bi + bb + kO*Cell_ptr->Drug));
  double k1_C1 = hh*(ain*Cell_ptr->C2 + bb*Cell_ptr->O - Cell_ptr->C1*(aa + bin));
  double k1_C2 = hh*(ae*Cell_ptr->C3 + bin*Cell_ptr->C1 - Cell_ptr->C2*(be + ain));
  double k1_C3 = hh*(be*Cell_ptr->C2 + rC*Cell_ptr->C3_drug - Cell_ptr->C3*(ae+kC*Cell_ptr->Drug));
  double k1_I = hh*(bi*Cell_ptr->O - ai*Cell_ptr->I);
  double k1_O_drug = hh*( aid*Cell_ptr->I_drug + aad*Cell_ptr->C1_drug + kO*Cell_ptr->Drug*Cell_ptr->O - Cell_ptr->O_drug*(bid + bbd+rO));  
  double k1_C1_drug = hh*(aind*Cell_ptr->C2_drug + bbd*Cell_ptr->O_drug - Cell_ptr->C1_drug*(aad + bind));
  double k1_C2_drug = hh*(aed*Cell_ptr->C3_drug + bind*Cell_ptr->C1_drug - Cell_ptr->C2_drug*(bed + aind));
  double k1_C3_drug = hh*(bed*Cell_ptr->C2_drug +kC*Cell_ptr->Drug*Cell_ptr->C3 - Cell_ptr->C3_drug*(aed + rC ));
  double k1_I_drug= hh*(bid*Cell_ptr->O_drug - aid*Cell_ptr->I_drug); 
  
  //Value of O after the k1 iteration; used for calculaton of k2
  double O_k1 = Cell_ptr->O + k1_O/2;
  double C1_k1 = Cell_ptr->C1 + k1_C1/2;
  double C2_k1 = Cell_ptr->C2 + k1_C2/2;
  double C3_k1 = Cell_ptr->C3 + k1_C3/2;
  double I_k1 = Cell_ptr->I + k1_I/2;
  double O_drug_k1 = Cell_ptr->O_drug + k1_O_drug/2;
  double C1_drug_k1 = Cell_ptr->C1_drug + k1_C1_drug/2;
  double C2_drug_k1 = Cell_ptr->C2_drug + k1_C2_drug/2;
  double C3_drug_k1 = Cell_ptr->C3_drug + k1_C3_drug/2;
  double I_drug_k1 = Cell_ptr->I_drug + k1_I_drug/2;
  
  
  //calculation of K2 values
  double k2_O = hh*(ai*I_k1 + aa*C1_k1 + rO*O_drug_k1 - O_k1*(bi + bb+ kO*Cell_ptr->Drug));
  double k2_C1 = hh*(ain*C2_k1 + bb*O_k1 - C1_k1*(aa + bin));
  double k2_C2 = hh*(ae*C3_k1 + bin*C1_k1 - C2_k1*(be + ain));
  double k2_C3 = hh*(be*C2_k1 +rC*C3_drug_k1 - C3_k1*(ae+ kC*Cell_ptr->Drug));
  double k2_I = hh*(bi*O_k1 - ai*I_k1);
  double k2_O_drug = hh*( aid*I_drug_k1 + aad*C1_drug_k1 + kO*Cell_ptr->Drug*O_k1 - O_drug_k1*(bid + bbd + rO)); 
  double k2_C1_drug = hh*(aind*C2_drug_k1 + bbd*O_drug_k1 - C1_drug_k1*(aad + bind));
  double k2_C2_drug = hh*(aed*C3_drug_k1 + bind*C1_drug_k1 - C2_drug_k1*(bed + aind));
  double k2_C3_drug = hh*(bed*C2_drug_k1 + kC*Cell_ptr->Drug*C3_k1 - C3_drug_k1*(aed + rC));
  double k2_I_drug= hh*(bid*O_drug_k1 - aid*I_drug_k1); 
  
  double O_k2 = Cell_ptr->O + k2_O/2;
  double C1_k2 = Cell_ptr->C1 + k2_C1/2;
  double C2_k2 = Cell_ptr->C2 + k2_C2/2;
  double C3_k2 = Cell_ptr->C3 + k2_C3/2;
  double I_k2 = Cell_ptr->I + k2_I/2;
  double O_drug_k2 = Cell_ptr->O_drug + k2_O_drug/2;
  double C1_drug_k2 = Cell_ptr->C1_drug + k2_C1_drug/2;
  double C2_drug_k2 = Cell_ptr->C2_drug + k2_C2_drug/2;
  double C3_drug_k2 = Cell_ptr->C3_drug + k2_C3_drug/2;
  double I_drug_k2 = Cell_ptr->I_drug + k2_I_drug/2;
  
  //calculation of K3 values
  
  double k3_O = hh*(ai*I_k2 + aa*C1_k2 + rO*O_drug_k2 - O_k2*(bi + bb+ kO*Cell_ptr->Drug));
  double k3_C1 = hh*(ain*C2_k2 + bb*O_k2 - C1_k2*(aa + bin));
  double k3_C2 = hh*(ae*C3_k2 + bin*C1_k2 - C2_k2*(be + ain));
  double k3_C3 = hh*(be*C2_k2 +rC*C3_drug_k2 - C3_k2*(ae+ kC*Cell_ptr->Drug));
  double k3_I = hh*(bi*O_k2 - ai*I_k2);
  double k3_O_drug = hh*( aid*I_drug_k2 + aad*C1_drug_k2 + kO*Cell_ptr->Drug*O_k2 - O_drug_k2*(bid + bbd + rO));
  double k3_C1_drug = hh*(aind*C2_drug_k2 + bbd*O_drug_k2 - C1_drug_k2*(aad + bind));
  double k3_C2_drug = hh*(aed*C3_drug_k2 + bind*C1_drug_k2 - C2_drug_k2*(bed + aind));
  double k3_C3_drug = hh*(bed*C2_drug_k2 + kC*Cell_ptr->Drug*C3_k2 - C3_drug_k2*(aed + rC));
  double k3_I_drug= hh*(bid*O_drug_k2 - aid*I_drug_k2);
  
  double O_k3 = Cell_ptr->O + k3_O;
  double C1_k3 = Cell_ptr->C1 + k3_C1;
  double C2_k3 = Cell_ptr->C2 + k3_C2;
  double C3_k3 = Cell_ptr->C3 + k3_C3;
  double I_k3 = Cell_ptr->I + k3_I;
  double O_drug_k3 = Cell_ptr->O_drug + k3_O_drug;
  double C1_drug_k3 = Cell_ptr->C1_drug + k3_C1_drug;
  double C2_drug_k3 = Cell_ptr->C2_drug + k3_C2_drug;
  double C3_drug_k3 = Cell_ptr->C3_drug + k3_C3_drug;
  double I_drug_k3 = Cell_ptr->I_drug + k3_I_drug;
  
  //calculation of K4 values
  
  double k4_O = hh*(ai*I_k3 + aa*C1_k3 + rO*O_drug_k3 - O_k3*(bi + bb+ kO*Cell_ptr->Drug));
  double k4_C1 = hh*(ain*C2_k3 + bb*O_k3 - C1_k3*(aa + bin));
  double k4_C2 = hh*(ae*C3_k3 + bin*C1_k3 - C2_k3*(be + ain));
  double k4_C3 = hh*(be*C2_k3 +rC*C3_drug_k3 - C3_k3*(ae+ kC*Cell_ptr->Drug));
  double k4_I = hh*(bi*O_k3 - ai*I_k3);
  double k4_O_drug = hh*( aid*I_drug_k3 + aad*C1_drug_k3 + kO*Cell_ptr->Drug*O_k3 - O_drug_k3*(bid + bbd + rO));
  double k4_C1_drug = hh*(aind*C2_drug_k3 + bbd*O_drug_k3 - C1_drug_k3*(aad + bind));
  double k4_C2_drug = hh*(aed*C3_drug_k3 + bind*C1_drug_k3 - C2_drug_k3*(bed + aind));
  double k4_C3_drug = hh*(bed*C2_drug_k3 + kC*Cell_ptr->Drug*C3_k3 - C3_drug_k3*(aed + rC));
  double k4_I_drug= hh*(bid*O_drug_k3 - aid*I_drug_k3);
  
  //New computed values of the channel probabilities
  Cell_ptr->O = Cell_ptr->O + (k1_O + 2*k2_O + 2*k3_O + k4_O)/6;
  Cell_ptr->C1 = Cell_ptr->C1 + (k1_C1 + 2*k2_C1 + 2*k3_C1 + k4_C1)/6;
  Cell_ptr->C2 = Cell_ptr->C2 + (k1_C2 + 2*k2_C2 + 2*k3_C2 + k4_C2)/6;
  Cell_ptr->C3 = Cell_ptr->C3 + (k1_C3 + 2*k2_C3 + 2*k3_C3 + k4_C3)/6;
  Cell_ptr->I = Cell_ptr->I + (k1_I + 2*k2_I + 2*k3_I + k4_I)/6;
  Cell_ptr->O_drug = Cell_ptr->O_drug + (k1_O_drug + 2*k2_O_drug + 2*k3_O_drug + k4_O_drug)/6;
  Cell_ptr->C1_drug = Cell_ptr->C1_drug + (k1_C1_drug + 2*k2_C1_drug + 2*k3_C1_drug + k4_C1_drug)/6;
  Cell_ptr->C2_drug = Cell_ptr->C2_drug + (k1_C2_drug + 2*k2_C2_drug + 2*k3_C2_drug + k4_C2_drug)/6;
  Cell_ptr->C3_drug = Cell_ptr->C3_drug + (k1_C3_drug + 2*k2_C3_drug + 2*k3_C3_drug + k4_C3_drug)/6;
  Cell_ptr->I_drug = Cell_ptr->I_drug + (k1_I_drug + 2*k2_I_drug + 2*k3_I_drug + k4_I_drug)/6;
  
  Cell_ptr->sum = Cell_ptr->O + Cell_ptr->C1 + Cell_ptr->C2 + Cell_ptr->C3 + Cell_ptr->I+ Cell_ptr->O_drug + Cell_ptr->C1_drug + Cell_ptr->C2_drug + Cell_ptr->C3_drug + Cell_ptr->I_drug;

  if (fabs (Cell_ptr->sum -1.0) > 0.0001) {
//    cout << "Error in Wild Type Sum at t = " << t << ",\t sum = " << Cell_ptr->sum << endl;
    *error=true;
    Cell_ptr->I_Kr = 0.0;
    }
  else {
    *error=false;
    Cell_ptr->I_Kr = G_Kr*sqrt(Ko/5.4)*(Cell_ptr->O + Cell_ptr->O_drug)*(Cell_ptr->V - E_Kr);
    }
  return;
}

void curr6(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){

  double E_Kr;
  double hh;
  hh = dt;
  
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0));       //var g_Kr_0: microS_per_nanoF {init: 0.024}; Lucia: Ko dependence is included in the formulation of the current

//        const double Q10=3;    To account for correction in temperature
//        const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));

  E_Kr = (R*T/F)*log(Ko/Ki);

//////Markov mechanism with NS drug////////Option4 --lau--/////

//  C3 = C2 = C1 = O = I
//  |    |         |   
// C3D= C2D= C1D = OD= ID

////////////////////////////////////////////////////

//rate cts. definitions:
//C3_C2=ae
//C2_C3=be
//C2_C1=ain
//C1_C2=bin
//C1_O=aa
//O_C1=bb
//O_I=bi
//I_O=ai
//ID_OD=aid
//OD_ID=bid
//C1D_OD=aad
//OD_CID=bbd
//C2D_C1D=aind
//C1D_C2D=bind
//C3D_C2D=aed
//C2D_C3D=bed
//*********Drug rates******
//O_OD=kO  ----> (microMolar*ms)-1) 
//OD_O=rO  -----> (ms)-1)
//C3_C3DD=kC  ----> (microMolar*ms)-1)
//C3DD_C3=rO  -----> (ms)-1)
//************************
//Drug equilibrium constants
//K_O=rO/kO
//K_C=rC/kC
////////////////////////////

//channel rates cts ----> 


/*
//WT + drug (3S)
//
  double ae=dx[0]*1.39429*exp(24.288586641+dx[1]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*1.12436*exp(13.637047252-dx[3]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*0.80694*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.905579*exp(13.110155065-16.420979385);
  double ai=dx[6]*0.989603*exp(28.782084605-dx[7]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bi=dx[8]*1.00033*exp(29.690385428+dx[9]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*0.98968*exp(22.010640617+dx[11]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*0.99917*exp(7.266288229-dx[13]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.39429*exp(24.288586641+dx[15]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*1.12436*exp(13.637047252-dx[17]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*0.80694*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.905579*exp(13.110155065-16.420979385);
  double aid=dx[20]*0.989603*exp(28.782084605-dx[21]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bid=dx[22]*1.00033*exp(29.690385428+dx[23]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*0.98968*exp(22.010640617+dx[25]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*0.99917*exp(7.266288229-dx[27]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
*/
//WT +drug (1S)
  double ae=dx[0]*1.50063*exp(24.288586641+dx[1]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*0.598723*exp(13.637047252-dx[3]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*1.10345*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.819984*exp(13.110155065-16.420979385);
  double ai=dx[6]*1.02159*exp(28.782084605-dx[7]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*0.998552*exp(29.690385428+dx[9]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*1.86823*exp(22.010640617+dx[11]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*1.00124*exp(7.266288229-dx[13]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.50063*exp(24.288586641+dx[15]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*0.598723*exp(13.637047252-dx[17]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*1.10345*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.819984*exp(13.110155065-16.420979385);
  double aid=dx[20]*1.02159*exp(28.782084605-dx[21]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bid=dx[22]*0.998552*exp(29.690385428+dx[23]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*1.86823*exp(22.010640617+dx[25]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*1.00124*exp(7.266288229-dx[27]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation


//drug rates cts ---->
  double kO=dx[28];
//double rO=dx[1];
  double kC=dx[29];
//double rC=dx[4];

//relations based on K and microscopic reversibility
//  double K_C=10.0; // (microM)//drug conscentration is 30 microM
//  double K_O=3.0; //(microM)

  double K_C=dx[30]; // (microM)//drug conscentration is 30 microM
  double K_O=dx[31]; //(microM)

  double rC=K_C*kC;
  double rO=K_O*kO;
//double rC=kC*(rO*ae*ai*aa*bbd*bind*bed)/(kO*be*bin*bb*aad*aind*aed)

//INitial conditions**********************************************************************

  // Calculation of k parameters for the Runge Kutta interations
  //Calculation of k1 values using the initial differential equation (set as f(t,V)
  double k1_O = hh*(ai*Cell_ptr->I + aa*Cell_ptr->C1+ rO*Cell_ptr->O_drug - Cell_ptr->O*(bi + bb + kO*Cell_ptr->Drug));
  double k1_C1 = hh*(ain*Cell_ptr->C2 + bb*Cell_ptr->O - Cell_ptr->C1*(aa + bin));
  double k1_C2 = hh*(ae*Cell_ptr->C3 + bin*Cell_ptr->C1 + rC*Cell_ptr->C2_drug - Cell_ptr->C2*(be + ain + kC*Cell_ptr->Drug));
  double k1_C3 = hh*(be*Cell_ptr->C2 + rC*Cell_ptr->C3_drug - Cell_ptr->C3*(ae+kC*Cell_ptr->Drug));
  double k1_I = hh*(bi*Cell_ptr->O - ai*Cell_ptr->I);
  double k1_O_drug = hh*( aid*Cell_ptr->I_drug + aad*Cell_ptr->C1_drug + kO*Cell_ptr->Drug*Cell_ptr->O - Cell_ptr->O_drug*(bid + bbd+rO));  
  double k1_C1_drug = hh*(aind*Cell_ptr->C2_drug + bbd*Cell_ptr->O_drug - Cell_ptr->C1_drug*(aad + bind));
  double k1_C2_drug = hh*(aed*Cell_ptr->C3_drug + bind*Cell_ptr->C1_drug + kC*Cell_ptr->Drug*Cell_ptr->C2 - Cell_ptr->C2_drug*(bed + aind + rC));
  double k1_C3_drug = hh*(bed*Cell_ptr->C2_drug +kC*Cell_ptr->Drug*Cell_ptr->C3 - Cell_ptr->C3_drug*(aed + rC ));
  double k1_I_drug= hh*(bid*Cell_ptr->O_drug - aid*Cell_ptr->I_drug);

  //Value of O after the k1 iteration; used for calculaton of k2
  double O_k1 = Cell_ptr->O + k1_O/2;
  double C1_k1 = Cell_ptr->C1 + k1_C1/2;
  double C2_k1 = Cell_ptr->C2 + k1_C2/2;
  double C3_k1 = Cell_ptr->C3 + k1_C3/2;
  double I_k1 = Cell_ptr->I + k1_I/2;
  double O_drug_k1 = Cell_ptr->O_drug + k1_O_drug/2;
  double C1_drug_k1 = Cell_ptr->C1_drug + k1_C1_drug/2;
  double C2_drug_k1 = Cell_ptr->C2_drug + k1_C2_drug/2;
  double C3_drug_k1 = Cell_ptr->C3_drug + k1_C3_drug/2;
  double I_drug_k1 = Cell_ptr->I_drug + k1_I_drug/2;
  //calculation of K2 values
  double k2_O = hh*(ai*I_k1 + aa*C1_k1 + rO*O_drug_k1 - O_k1*(bi + bb+ kO*Cell_ptr->Drug));
  double k2_C1 = hh*(ain*C2_k1 + bb*O_k1 - C1_k1*(aa + bin));
  double k2_C2 = hh*(ae*C3_k1 + bin*C1_k1 +rC*C2_drug_k1 - C2_k1*(be + ain + kC*Cell_ptr->Drug));
  double k2_C3 = hh*(be*C2_k1 +rC*C3_drug_k1 - C3_k1*(ae+ kC*Cell_ptr->Drug));
  double k2_I = hh*(bi*O_k1 - ai*I_k1);
  double k2_O_drug = hh*( aid*I_drug_k1 + aad*C1_drug_k1 + kO*Cell_ptr->Drug*O_k1 - O_drug_k1*(bid + bbd + rO)); 
  double k2_C1_drug = hh*(aind*C2_drug_k1 + bbd*O_drug_k1 - C1_drug_k1*(aad + bind));
  double k2_C2_drug = hh*(aed*C3_drug_k1 + bind*C1_drug_k1 + kC*Cell_ptr->Drug*C2_k1 - C2_drug_k1*(bed + aind + rC));
  double k2_C3_drug = hh*(bed*C2_drug_k1 + kC*Cell_ptr->Drug*C3_k1 - C3_drug_k1*(aed + rC));
  double k2_I_drug= hh*(bid*O_drug_k1 - aid*I_drug_k1); 
  double O_k2 = Cell_ptr->O + k2_O/2;
  double C1_k2 = Cell_ptr->C1 + k2_C1/2;
  double C2_k2 = Cell_ptr->C2 + k2_C2/2;
  double C3_k2 = Cell_ptr->C3 + k2_C3/2;
  double I_k2 = Cell_ptr->I + k2_I/2;
  double O_drug_k2 = Cell_ptr->O_drug + k2_O_drug/2;
  double C1_drug_k2 = Cell_ptr->C1_drug + k2_C1_drug/2;
  double C2_drug_k2 = Cell_ptr->C2_drug + k2_C2_drug/2;
  double C3_drug_k2 = Cell_ptr->C3_drug + k2_C3_drug/2;
  double I_drug_k2 = Cell_ptr->I_drug + k2_I_drug/2;
  //calculation of K3 values
  double k3_O = hh*(ai*I_k2 + aa*C1_k2 + rO*O_drug_k2 - O_k2*(bi + bb+ kO*Cell_ptr->Drug));
  double k3_C1 = hh*(ain*C2_k2 + bb*O_k2 - C1_k2*(aa + bin));
  double k3_C2 = hh*(ae*C3_k2 + bin*C1_k2 +rC*C2_drug_k2 - C2_k2*(be + ain + kC*Cell_ptr->Drug));
  double k3_C3 = hh*(be*C2_k2 +rC*C3_drug_k2 - C3_k2*(ae+ kC*Cell_ptr->Drug));
  double k3_I = hh*(bi*O_k2 - ai*I_k2);
  double k3_O_drug = hh*( aid*I_drug_k2 + aad*C1_drug_k2 + kO*Cell_ptr->Drug*O_k2 - O_drug_k2*(bid + bbd + rO));
  double k3_C1_drug = hh*(aind*C2_drug_k2 + bbd*O_drug_k2 - C1_drug_k2*(aad + bind));
  double k3_C2_drug = hh*(aed*C3_drug_k2 + bind*C1_drug_k2 + kC*Cell_ptr->Drug*C2_k2 - C2_drug_k2*(bed + aind + rC));
  double k3_C3_drug = hh*(bed*C2_drug_k2 + kC*Cell_ptr->Drug*C3_k2 - C3_drug_k2*(aed + rC));
  double k3_I_drug= hh*(bid*O_drug_k2 - aid*I_drug_k2);
  double O_k3 = Cell_ptr->O + k3_O;
  double C1_k3 = Cell_ptr->C1 + k3_C1;
  double C2_k3 = Cell_ptr->C2 + k3_C2;
  double C3_k3 = Cell_ptr->C3 + k3_C3;
  double I_k3 = Cell_ptr->I + k3_I;
  double O_drug_k3 = Cell_ptr->O_drug + k3_O_drug;
  double C1_drug_k3 = Cell_ptr->C1_drug + k3_C1_drug;
  double C2_drug_k3 = Cell_ptr->C2_drug + k3_C2_drug;
  double C3_drug_k3 = Cell_ptr->C3_drug + k3_C3_drug;
  double I_drug_k3 = Cell_ptr->I_drug + k3_I_drug;
  //calculation of K4 values
  double k4_O = hh*(ai*I_k3 + aa*C1_k3 + rO*O_drug_k3 - O_k3*(bi + bb+ kO*Cell_ptr->Drug));
  double k4_C1 = hh*(ain*C2_k3 + bb*O_k3 - C1_k3*(aa + bin));
  double k4_C2 = hh*(ae*C3_k3 + bin*C1_k3 +rC*C2_drug_k3 - C2_k3*(be + ain +  kC*Cell_ptr->Drug));
  double k4_C3 = hh*(be*C2_k3 +rC*C3_drug_k3 - C3_k3*(ae+ kC*Cell_ptr->Drug));
  double k4_I = hh*(bi*O_k3 - ai*I_k3);
  double k4_O_drug = hh*( aid*I_drug_k3 + aad*C1_drug_k3 + kO*Cell_ptr->Drug*O_k3 - O_drug_k3*(bid + bbd + rO));
  double k4_C1_drug = hh*(aind*C2_drug_k3 + bbd*O_drug_k3 - C1_drug_k3*(aad + bind));
  double k4_C2_drug = hh*(aed*C3_drug_k3 + bind*C1_drug_k3 + kC*Cell_ptr->Drug*C2_k3 - C2_drug_k3*(bed + aind + rC));
  double k4_C3_drug = hh*(bed*C2_drug_k3 + kC*Cell_ptr->Drug*C3_k3 - C3_drug_k3*(aed + rC));
  double k4_I_drug= hh*(bid*O_drug_k3 - aid*I_drug_k3);
  
  
  //New computed values of the channel probabilities
  Cell_ptr->O = Cell_ptr->O + (k1_O + 2*k2_O + 2*k3_O + k4_O)/6;
  Cell_ptr->C1 = Cell_ptr->C1 + (k1_C1 + 2*k2_C1 + 2*k3_C1 + k4_C1)/6;
  Cell_ptr->C2 = Cell_ptr->C2 + (k1_C2 + 2*k2_C2 + 2*k3_C2 + k4_C2)/6;
  Cell_ptr->C3 = Cell_ptr->C3 + (k1_C3 + 2*k2_C3 + 2*k3_C3 + k4_C3)/6;
  Cell_ptr->I = Cell_ptr->I + (k1_I + 2*k2_I + 2*k3_I + k4_I)/6;
  Cell_ptr->O_drug = Cell_ptr->O_drug + (k1_O_drug + 2*k2_O_drug + 2*k3_O_drug + k4_O_drug)/6;
  Cell_ptr->C1_drug = Cell_ptr->C1_drug + (k1_C1_drug + 2*k2_C1_drug + 2*k3_C1_drug + k4_C1_drug)/6;
  Cell_ptr->C2_drug = Cell_ptr->C2_drug + (k1_C2_drug + 2*k2_C2_drug + 2*k3_C2_drug + k4_C2_drug)/6;
  Cell_ptr->C3_drug = Cell_ptr->C3_drug + (k1_C3_drug + 2*k2_C3_drug + 2*k3_C3_drug + k4_C3_drug)/6;
  Cell_ptr->I_drug = Cell_ptr->I_drug + (k1_I_drug + 2*k2_I_drug + 2*k3_I_drug + k4_I_drug)/6;
  
  Cell_ptr->sum = Cell_ptr->O + Cell_ptr->C1 + Cell_ptr->C2 + Cell_ptr->C3 + Cell_ptr->I+ Cell_ptr->O_drug + Cell_ptr->C1_drug + Cell_ptr->C2_drug + Cell_ptr->C3_drug + Cell_ptr->I_drug;


  if (fabs (Cell_ptr->sum -1.0) > 0.0001) {
//    cout << "Error in Wild Type Sum at t = " << t << ",\t sum = " << Cell_ptr->sum << endl;
    *error=true;
    Cell_ptr->I_Kr = 0.0;
    }
  else {
    *error=false;
    Cell_ptr->I_Kr = G_Kr*sqrt(Ko/5.4)*(Cell_ptr->O + Cell_ptr->O_drug)*(Cell_ptr->V - E_Kr);
    }
  return;
}

void curr7(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){

  double E_Kr;
  double hh;
  hh = dt;
  
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0));       //var g_Kr_0: microS_per_nanoF {init: 0.024}; Lucia: Ko dependence is included in the formulation of the current

//        const double Q10=3;    To account for correction in temperature
//        const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));

  E_Kr = (R*T/F)*log(Ko/Ki);

//////Markov mechanism with NS drug////////Option5 --lau--/////

//  C3 = C2 = C1 = O = I
//  |         |    |   
// C3D= C2D= C1D = OD= ID

////////////////////////////////////////////////////

//rate cts. definitions:
//C3_C2=ae
//C2_C3=be
//C2_C1=ain
//C1_C2=bin
//C1_O=aa
//O_C1=bb
//O_I=bi
//I_O=ai
//ID_OD=aid
//OD_ID=bid
//C1D_OD=aad
//OD_CID=bbd
//C2D_C1D=aind
//C1D_C2D=bind
//C3D_C2D=aed
//C2D_C3D=bed
//*********Drug rates******
//O_OD=kO  ----> (microMolar*ms)-1) 
//OD_O=rO  -----> (ms)-1)
//C3_C3DD=kC  ----> (microMolar*ms)-1)
//C3DD_C3=rO  -----> (ms)-1)
//************************
//Drug equilibrium constants
//K_O=rO/kO
//K_C=rC/kC
////////////////////////////

//channel rates cts ----> 


/*
//WT + drug (3S)
//
  double ae=dx[0]*1.39429*exp(24.288586641+dx[1]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*1.12436*exp(13.637047252-dx[3]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*0.80694*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.905579*exp(13.110155065-16.420979385);
  double ai=dx[6]*0.989603*exp(28.782084605-dx[7]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bi=dx[8]*1.00033*exp(29.690385428+dx[9]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*0.98968*exp(22.010640617+dx[11]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*0.99917*exp(7.266288229-dx[13]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.39429*exp(24.288586641+dx[15]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*1.12436*exp(13.637047252-dx[17]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*0.80694*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.905579*exp(13.110155065-16.420979385);
  double aid=dx[20]*0.989603*exp(28.782084605-dx[21]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bid=dx[22]*1.00033*exp(29.690385428+dx[23]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*0.98968*exp(22.010640617+dx[25]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*0.99917*exp(7.266288229-dx[27]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
*/

//WT+drug (1S)
  double ae=dx[0]*1.50063*exp(24.288586641+dx[1]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*0.598723*exp(13.637047252-dx[3]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*1.10345*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.819984*exp(13.110155065-16.420979385);
  double ai=dx[6]*1.02159*exp(28.782084605-dx[7]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*0.998552*exp(29.690385428+dx[9]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*1.86823*exp(22.010640617+dx[11]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*1.00124*exp(7.266288229-dx[13]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.50063*exp(24.288586641+dx[15]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*0.598723*exp(13.637047252-dx[17]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*1.10345*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.819984*exp(13.110155065-16.420979385);
  double aid=dx[20]*1.02159*exp(28.782084605-dx[21]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bid=dx[22]*0.998552*exp(29.690385428+dx[23]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*1.86823*exp(22.010640617+dx[25]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*1.00124*exp(7.266288229-dx[27]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation


//drug rates cts ---->
  double kO=dx[28];
//double rO=dx[1];
  double kC=dx[29];
//double rC=dx[4];

//relations based on K and microscopic reversibility
//  double K_C=10.0; // (microM)//drug conscentration is 30 microM
//  double K_O=3.0; //(microM)

  double K_C=dx[30]; // (microM)//drug conscentration is 30 microM
  double K_O=dx[31]; //(microM)

  double rC=K_C*kC;
  double rO=K_O*kO;
//double rC=kC*(rO*ae*ai*aa*bbd*bind*bed)/(kO*be*bin*bb*aad*aind*aed)

//INitial conditions**********************************************************************

  // Calculation of k parameters for the Runge Kutta interations
  //Calculation of k1 values using the initial differential equation (set as f(t,V)
  double k1_O = hh*(ai*Cell_ptr->I + aa*Cell_ptr->C1+ rO*Cell_ptr->O_drug - Cell_ptr->O*(bi + bb + kO*Cell_ptr->Drug));
  double k1_C1 = hh*(ain*Cell_ptr->C2 + bb*Cell_ptr->O + rC*Cell_ptr->C1_drug - Cell_ptr->C1*(aa + bin + kC*Cell_ptr->Drug));
  double k1_C2 = hh*(ae*Cell_ptr->C3 + bin*Cell_ptr->C1 - Cell_ptr->C2*(be + ain));
  double k1_C3 = hh*(be*Cell_ptr->C2 + rC*Cell_ptr->C3_drug - Cell_ptr->C3*(ae+kC*Cell_ptr->Drug));
  double k1_I = hh*(bi*Cell_ptr->O - ai*Cell_ptr->I);
  double k1_O_drug = hh*( aid*Cell_ptr->I_drug + aad*Cell_ptr->C1_drug + kO*Cell_ptr->Drug*Cell_ptr->O - Cell_ptr->O_drug*(bid + bbd+rO));  
  double k1_C1_drug = hh*(aind*Cell_ptr->C2_drug + bbd*Cell_ptr->O_drug + kC*Cell_ptr->Drug*Cell_ptr->C1 - Cell_ptr->C1_drug*(aad + bind +rC));
  double k1_C2_drug = hh*(aed*Cell_ptr->C3_drug + bind*Cell_ptr->C1_drug - Cell_ptr->C2_drug*(bed + aind));
  double k1_C3_drug = hh*(bed*Cell_ptr->C2_drug +kC*Cell_ptr->Drug*Cell_ptr->C3 - Cell_ptr->C3_drug*(aed + rC ));
  double k1_I_drug= hh*(bid*Cell_ptr->O_drug - aid*Cell_ptr->I_drug);

  //Value of O after the k1 iteration; used for calculaton of k2
  double O_k1 = Cell_ptr->O + k1_O/2;
  double C1_k1 = Cell_ptr->C1 + k1_C1/2;
  double C2_k1 = Cell_ptr->C2 + k1_C2/2;
  double C3_k1 = Cell_ptr->C3 + k1_C3/2;
  double I_k1 = Cell_ptr->I + k1_I/2;
  double O_drug_k1 = Cell_ptr->O_drug + k1_O_drug/2;
  double C1_drug_k1 = Cell_ptr->C1_drug + k1_C1_drug/2;
  double C2_drug_k1 = Cell_ptr->C2_drug + k1_C2_drug/2;
  double C3_drug_k1 = Cell_ptr->C3_drug + k1_C3_drug/2;
  double I_drug_k1 = Cell_ptr->I_drug + k1_I_drug/2;
  //calculation of K2 values
  double k2_O = hh*(ai*I_k1 + aa*C1_k1 + rO*O_drug_k1 - O_k1*(bi + bb+ kO*Cell_ptr->Drug));
  double k2_C1 = hh*(ain*C2_k1 + bb*O_k1 +rC*C1_drug_k1 - C1_k1*(aa + bin + kC*Cell_ptr->Drug));
  double k2_C2 = hh*(ae*C3_k1 + bin*C1_k1 - C2_k1*(be + ain));
  double k2_C3 = hh*(be*C2_k1 +rC*C3_drug_k1 - C3_k1*(ae+ kC*Cell_ptr->Drug));
  double k2_I = hh*(bi*O_k1 - ai*I_k1);
  double k2_O_drug = hh*( aid*I_drug_k1 + aad*C1_drug_k1 + kO*Cell_ptr->Drug*O_k1 - O_drug_k1*(bid + bbd + rO)); 
  double k2_C1_drug = hh*(aind*C2_drug_k1 + bbd*O_drug_k1 + kC*Cell_ptr->Drug*C1_k1 - C1_drug_k1*(aad + bind + rC));
  double k2_C2_drug = hh*(aed*C3_drug_k1 + bind*C1_drug_k1 - C2_drug_k1*(bed + aind));
  double k2_C3_drug = hh*(bed*C2_drug_k1 + kC*Cell_ptr->Drug*C3_k1 - C3_drug_k1*(aed + rC));
  double k2_I_drug= hh*(bid*O_drug_k1 - aid*I_drug_k1); 
  double O_k2 = Cell_ptr->O + k2_O/2;
  double C1_k2 = Cell_ptr->C1 + k2_C1/2;
  double C2_k2 = Cell_ptr->C2 + k2_C2/2;
  double C3_k2 = Cell_ptr->C3 + k2_C3/2;
  double I_k2 = Cell_ptr->I + k2_I/2;
  double O_drug_k2 = Cell_ptr->O_drug + k2_O_drug/2;
  double C1_drug_k2 = Cell_ptr->C1_drug + k2_C1_drug/2;
  double C2_drug_k2 = Cell_ptr->C2_drug + k2_C2_drug/2;
  double C3_drug_k2 = Cell_ptr->C3_drug + k2_C3_drug/2;
  double I_drug_k2 = Cell_ptr->I_drug + k2_I_drug/2;
  //calculation of K3 values
  double k3_O = hh*(ai*I_k2 + aa*C1_k2 + rO*O_drug_k2 - O_k2*(bi + bb+ kO*Cell_ptr->Drug));
  double k3_C1 = hh*(ain*C2_k2 + bb*O_k2 +rC*C1_drug_k2 - C1_k2*(aa + bin + kC*Cell_ptr->Drug));
  double k3_C2 = hh*(ae*C3_k2 + bin*C1_k2 - C2_k2*(be + ain));
  double k3_C3 = hh*(be*C2_k2 +rC*C3_drug_k2 - C3_k2*(ae+ kC*Cell_ptr->Drug));
  double k3_I = hh*(bi*O_k2 - ai*I_k2);
  double k3_O_drug = hh*( aid*I_drug_k2 + aad*C1_drug_k2 + kO*Cell_ptr->Drug*O_k2 - O_drug_k2*(bid + bbd + rO));
  double k3_C1_drug = hh*(aind*C2_drug_k2 + bbd*O_drug_k2 + kC*Cell_ptr->Drug*C1_k2 - C1_drug_k2*(aad + bind + rC));
  double k3_C2_drug = hh*(aed*C3_drug_k2 + bind*C1_drug_k2 - C2_drug_k2*(bed + aind));
  double k3_C3_drug = hh*(bed*C2_drug_k2 + kC*Cell_ptr->Drug*C3_k2 - C3_drug_k2*(aed + rC));
  double k3_I_drug= hh*(bid*O_drug_k2 - aid*I_drug_k2);
  double O_k3 = Cell_ptr->O + k3_O;
  double C1_k3 = Cell_ptr->C1 + k3_C1;
  double C2_k3 = Cell_ptr->C2 + k3_C2;
  double C3_k3 = Cell_ptr->C3 + k3_C3;
  double I_k3 = Cell_ptr->I + k3_I;
  double O_drug_k3 = Cell_ptr->O_drug + k3_O_drug;
  double C1_drug_k3 = Cell_ptr->C1_drug + k3_C1_drug;
  double C2_drug_k3 = Cell_ptr->C2_drug + k3_C2_drug;
  double C3_drug_k3 = Cell_ptr->C3_drug + k3_C3_drug;
  double I_drug_k3 = Cell_ptr->I_drug + k3_I_drug;
  //calculation of K4 values
  double k4_O = hh*(ai*I_k3 + aa*C1_k3 + rO*O_drug_k3 - O_k3*(bi + bb+ kO*Cell_ptr->Drug));
  double k4_C1 = hh*(ain*C2_k3 + bb*O_k3 +rC*C1_drug_k3 - C1_k3*(aa + bin +  kC*Cell_ptr->Drug));
  double k4_C2 = hh*(ae*C3_k3 + bin*C1_k3 - C2_k3*(be + ain));
  double k4_C3 = hh*(be*C2_k3 +rC*C3_drug_k3 - C3_k3*(ae+ kC*Cell_ptr->Drug));
  double k4_I = hh*(bi*O_k3 - ai*I_k3);
  double k4_O_drug = hh*( aid*I_drug_k3 + aad*C1_drug_k3 + kO*Cell_ptr->Drug*O_k3 - O_drug_k3*(bid + bbd + rO));
  double k4_C1_drug = hh*(aind*C2_drug_k3 + bbd*O_drug_k3 + kC*Cell_ptr->Drug*C1_k3 - C1_drug_k3*(aad + bind + rC));
  double k4_C2_drug = hh*(aed*C3_drug_k3 + bind*C1_drug_k3 - C2_drug_k3*(bed + aind));
  double k4_C3_drug = hh*(bed*C2_drug_k3 + kC*Cell_ptr->Drug*C3_k3 - C3_drug_k3*(aed + rC));
  double k4_I_drug= hh*(bid*O_drug_k3 - aid*I_drug_k3);
  
  
  //New computed values of the channel probabilities
  Cell_ptr->O = Cell_ptr->O + (k1_O + 2*k2_O + 2*k3_O + k4_O)/6;
  Cell_ptr->C1 = Cell_ptr->C1 + (k1_C1 + 2*k2_C1 + 2*k3_C1 + k4_C1)/6;
  Cell_ptr->C2 = Cell_ptr->C2 + (k1_C2 + 2*k2_C2 + 2*k3_C2 + k4_C2)/6;
  Cell_ptr->C3 = Cell_ptr->C3 + (k1_C3 + 2*k2_C3 + 2*k3_C3 + k4_C3)/6;
  Cell_ptr->I = Cell_ptr->I + (k1_I + 2*k2_I + 2*k3_I + k4_I)/6;
  Cell_ptr->O_drug = Cell_ptr->O_drug + (k1_O_drug + 2*k2_O_drug + 2*k3_O_drug + k4_O_drug)/6;
  Cell_ptr->C1_drug = Cell_ptr->C1_drug + (k1_C1_drug + 2*k2_C1_drug + 2*k3_C1_drug + k4_C1_drug)/6;
  Cell_ptr->C2_drug = Cell_ptr->C2_drug + (k1_C2_drug + 2*k2_C2_drug + 2*k3_C2_drug + k4_C2_drug)/6;
  Cell_ptr->C3_drug = Cell_ptr->C3_drug + (k1_C3_drug + 2*k2_C3_drug + 2*k3_C3_drug + k4_C3_drug)/6;
  Cell_ptr->I_drug = Cell_ptr->I_drug + (k1_I_drug + 2*k2_I_drug + 2*k3_I_drug + k4_I_drug)/6;
  
  Cell_ptr->sum = Cell_ptr->O + Cell_ptr->C1 + Cell_ptr->C2 + Cell_ptr->C3 + Cell_ptr->I+ Cell_ptr->O_drug + Cell_ptr->C1_drug + Cell_ptr->C2_drug + Cell_ptr->C3_drug + Cell_ptr->I_drug;


  if (fabs (Cell_ptr->sum -1.0) > 0.0001) {
//    cout << "Error in Wild Type Sum at t = " << t << ",\t sum = " << Cell_ptr->sum << endl;
    *error=true;
    Cell_ptr->I_Kr = 0.0;
    }
  else {
    *error=false;
    Cell_ptr->I_Kr = G_Kr*sqrt(Ko/5.4)*(Cell_ptr->O + Cell_ptr->O_drug)*(Cell_ptr->V - E_Kr);
    }
  return;
}
void curr8(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){

  double E_Kr;
  double hh;
  hh = dt;
  
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0));       //var g_Kr_0: microS_per_nanoF {init: 0.024}; Lucia: Ko dependence is included in the formulation of the current

//        const double Q10=3;    To account for correction in temperature
//        const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));

  E_Kr = (R*T/F)*log(Ko/Ki);


//////Markov mechanism with NS drug////////Option6 --lau--/////

//  C3 = C2 = C1 = O = I
//            |    |   
// C3D= C2D= C1D = OD = ID

////////////////////////////////////////////////////

//rate cts. definitions:
//C3_C2=ae
//C2_C3=be
//C2_C1=ain
//C1_C2=bin
//C1_O=aa
//O_C1=bb
//O_I=bi
//I_O=ai
//C1D_OD=aad
//OD_CID=bbd
//C2D_C1D=aind
//C1D_C2D=bind
//C3D_C2D=aed
//C2D_C3D=bed
//*********Drug rates******
//O_OD=kO  ----> (microMolar*ms)-1) 
//OD_O=rO  -----> (ms)-1)
//C3_C3DD=kC  ----> (microMolar*ms)-1)
//C3DD_C3=rO  -----> (ms)-1)
//************************
//Drug equilibrium constants
//K_O=rO/kO
//K_C=rC/kC
////////////////////////////

//channel rates cts ----> 


/*
//WT + drug (3S)
//
  double ae=dx[0]*1.39429*exp(24.288586641+dx[1]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*1.12436*exp(13.637047252-dx[3]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*0.80694*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.905579*exp(13.110155065-16.420979385);
  double ai=dx[6]*0.989603*exp(28.782084605-dx[7]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bi=dx[8]*1.00033*exp(29.690385428+dx[9]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*0.98968*exp(22.010640617+dx[11]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*0.99917*exp(7.266288229-dx[13]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.39429*exp(24.288586641+dx[15]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*1.12436*exp(13.637047252-dx[17]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*0.80694*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.905579*exp(13.110155065-16.420979385);
  double aad=dx[20]*0.98968*exp(22.010640617+dx[21]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[22]*0.99917*exp(7.266288229-dx[23]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
*/

//WT+ drug (1S)
  double ae=dx[0]*1.50063*exp(24.288586641+dx[1]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*0.598723*exp(13.637047252-dx[3]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*1.10345*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.819984*exp(13.110155065-16.420979385);
  double ai=dx[6]*1.02159*exp(28.782084605-dx[7]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*0.998552*exp(29.690385428+dx[9]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*1.86823*exp(22.010640617+dx[11]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*1.00124*exp(7.266288229-dx[13]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.50063*exp(24.288586641+dx[15]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*0.598723*exp(13.637047252-dx[17]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*1.10345*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.819984*exp(13.110155065-16.420979385);
  double aid=dx[20]*1.02159*exp(28.782084605-dx[21]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bid=dx[22]*0.998552*exp(29.690385428+dx[23]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*1.86823*exp(22.010640617+dx[25]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*1.00124*exp(7.266288229-dx[27]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation


//drug rates cts ---->
  double kO=dx[28];
//double rO=dx[1];
  double kC=dx[29];
//double rC=dx[4];

//relations based on K and microscopic reversibility
//  double K_C=10.0; // (microM)//drug conscentration is 30 microM
//  double K_O=3.0; //(microM)

  double K_C=dx[30]; // (microM)//drug conscentration is 30 microM
  double K_O=dx[31]; //(microM)

  double rC=K_C*kC;
  double rO=K_O*kO;
//double rC=kC*(rO*ae*ai*aa*bbd*bind*bed)/(kO*be*bin*bb*aad*aind*aed)

//INitial conditions**********************************************************************

  // Calculation of k parameters for the Runge Kutta interations
  //Calculation of k1 values using the initial differential equation (set as f(t,V)
  double k1_O = hh*(ai*Cell_ptr->I + aa*Cell_ptr->C1+ rO*Cell_ptr->O_drug - Cell_ptr->O*(bi + bb + kO*Cell_ptr->Drug));
  double k1_C1 = hh*(ain*Cell_ptr->C2 + bb*Cell_ptr->O + rC*Cell_ptr->C1_drug - Cell_ptr->C1*(aa + bin + kC*Cell_ptr->Drug));
  double k1_C2 = hh*(ae*Cell_ptr->C3 + bin*Cell_ptr->C1 - Cell_ptr->C2*(be + ain));
  double k1_C3 = hh*(be*Cell_ptr->C2 - Cell_ptr->C3*ae);
  double k1_I = hh*(bi*Cell_ptr->O - ai*Cell_ptr->I);
  double k1_O_drug = hh*(aid*Cell_ptr->I_drug + aad*Cell_ptr->C1_drug + kO*Cell_ptr->Drug*Cell_ptr->O - Cell_ptr->O_drug*(bid + bbd + rO));  
  double k1_C1_drug = hh*(aind*Cell_ptr->C2_drug + bbd*Cell_ptr->O_drug + kC*Cell_ptr->Drug*Cell_ptr->C1 - Cell_ptr->C1_drug*(aad + bind +rC));
  double k1_C2_drug = hh*(aed*Cell_ptr->C3_drug + bind*Cell_ptr->C1_drug - Cell_ptr->C2_drug*(bed + aind));
  double k1_C3_drug = hh*(bed*Cell_ptr->C2_drug - Cell_ptr->C3_drug*aed);
  double k1_I_drug= hh*(bid*Cell_ptr->O_drug - aid*Cell_ptr->I_drug);

  //Value of O after the k1 iteration; used for calculaton of k2
  double O_k1 = Cell_ptr->O + k1_O/2;
  double C1_k1 = Cell_ptr->C1 + k1_C1/2;
  double C2_k1 = Cell_ptr->C2 + k1_C2/2;
  double C3_k1 = Cell_ptr->C3 + k1_C3/2;
  double I_k1 = Cell_ptr->I + k1_I/2;
  double O_drug_k1 = Cell_ptr->O_drug + k1_O_drug/2;
  double C1_drug_k1 = Cell_ptr->C1_drug + k1_C1_drug/2;
  double C2_drug_k1 = Cell_ptr->C2_drug + k1_C2_drug/2;
  double C3_drug_k1 = Cell_ptr->C3_drug + k1_C3_drug/2;
  double I_drug_k1 = Cell_ptr->I_drug + k1_I_drug/2;

  //calculation of K2 values
  double k2_O = hh*(ai*I_k1 + aa*C1_k1 + rO*O_drug_k1 - O_k1*(bi + bb+ kO*Cell_ptr->Drug));
  double k2_C1 = hh*(ain*C2_k1 + bb*O_k1 +rC*C1_drug_k1 - C1_k1*(aa + bin + kC*Cell_ptr->Drug));
  double k2_C2 = hh*(ae*C3_k1 + bin*C1_k1 - C2_k1*(be + ain ));
  double k2_C3 = hh*(be*C2_k1 - C3_k1*ae);
  double k2_I = hh*(bi*O_k1 - ai*I_k1);
  double k2_O_drug = hh*( aid*I_drug_k1 + aad*C1_drug_k1 + kO*Cell_ptr->Drug*O_k1 - O_drug_k1*(bid + bbd + rO)); 
  double k2_C1_drug = hh*(aind*C2_drug_k1 + bbd*O_drug_k1 + kC*Cell_ptr->Drug*C1_k1 - C1_drug_k1*(aad + bind + rC));
  double k2_C2_drug = hh*(aed*C3_drug_k1 + bind*C1_drug_k1 - C2_drug_k1*(bed + aind ));
  double k2_C3_drug = hh*(bed*C2_drug_k1 - C3_drug_k1*aed);
  double k2_I_drug= hh*(bid*O_drug_k1 - aid*I_drug_k1);
  double O_k2 = Cell_ptr->O + k2_O/2;
  double C1_k2 = Cell_ptr->C1 + k2_C1/2;
  double C2_k2 = Cell_ptr->C2 + k2_C2/2;
  double C3_k2 = Cell_ptr->C3 + k2_C3/2;
  double I_k2 = Cell_ptr->I + k2_I/2;
  double O_drug_k2 = Cell_ptr->O_drug + k2_O_drug/2;
  double C1_drug_k2 = Cell_ptr->C1_drug + k2_C1_drug/2;
  double C2_drug_k2 = Cell_ptr->C2_drug + k2_C2_drug/2;
  double C3_drug_k2 = Cell_ptr->C3_drug + k2_C3_drug/2;
  double I_drug_k2 = Cell_ptr->I_drug + k2_I_drug/2;
  //calculation of K3 values
  double k3_O = hh*(ai*I_k2 + aa*C1_k2 + rO*O_drug_k2 - O_k2*(bi + bb+ kO*Cell_ptr->Drug));
  double k3_C1 = hh*(ain*C2_k2 + bb*O_k2 +rC*C1_drug_k2 - C1_k2*(aa + bin + kC*Cell_ptr->Drug));
  double k3_C2 = hh*(ae*C3_k2 + bin*C1_k2 - C2_k2*(be + ain));
  double k3_C3 = hh*(be*C2_k2 - C3_k2*ae);
  double k3_I = hh*(bi*O_k2 - ai*I_k2);
  double k3_O_drug = hh*( aid*I_drug_k2 + aad*C1_drug_k2 + kO*Cell_ptr->Drug*O_k2 - O_drug_k2*(bid + bbd + rO));
  double k3_C1_drug = hh*(aind*C2_drug_k2 + bbd*O_drug_k2 + kC*Cell_ptr->Drug*C1_k2 - C1_drug_k2*(aad + bind + rC));
  double k3_C2_drug = hh*(aed*C3_drug_k2 + bind*C1_drug_k2 - C2_drug_k2*(bed + aind));
  double k3_C3_drug = hh*(bed*C2_drug_k2 - C3_drug_k2*aed);
  double k3_I_drug= hh*(bid*O_drug_k2 - aid*I_drug_k2);
  double O_k3 = Cell_ptr->O + k3_O;
  double C1_k3 = Cell_ptr->C1 + k3_C1;
  double C2_k3 = Cell_ptr->C2 + k3_C2;
  double C3_k3 = Cell_ptr->C3 + k3_C3;
  double I_k3 = Cell_ptr->I + k3_I;
  double O_drug_k3 = Cell_ptr->O_drug + k3_O_drug;
  double C1_drug_k3 = Cell_ptr->C1_drug + k3_C1_drug;
  double C2_drug_k3 = Cell_ptr->C2_drug + k3_C2_drug;
  double C3_drug_k3 = Cell_ptr->C3_drug + k3_C3_drug;
  double I_drug_k3 = Cell_ptr->I_drug + k3_I_drug;
  //calculation of K4 values
  double k4_O = hh*(ai*I_k3 + aa*C1_k3 + rO*O_drug_k3 - O_k3*(bi + bb+ kO*Cell_ptr->Drug));
  double k4_C1 = hh*(ain*C2_k3 + bb*O_k3 +rC*C1_drug_k3 - C1_k3*(aa + bin +  kC*Cell_ptr->Drug));
  double k4_C2 = hh*(ae*C3_k3 + bin*C1_k3 - C2_k3*(be + ain ));
  double k4_C3 = hh*(be*C2_k3 - C3_k3*ae);
  double k4_I = hh*(bi*O_k3 - ai*I_k3);
  double k4_O_drug = hh*( aid*I_drug_k3 + aad*C1_drug_k3 + kO*Cell_ptr->Drug*O_k3 - O_drug_k3*(bid + bbd + rO));
  double k4_C1_drug = hh*(aind*C2_drug_k3 + bbd*O_drug_k3 + kC*Cell_ptr->Drug*C1_k3 - C1_drug_k3*(aad + bind + rC));
  double k4_C2_drug = hh*(aed*C3_drug_k3 + bind*C1_drug_k3 - C2_drug_k3*(bed + aind));
  double k4_C3_drug = hh*(bed*C2_drug_k3 - C3_drug_k3*aed);
  double k4_I_drug= hh*(bid*O_drug_k3 - aid*I_drug_k3);

  
  //New computed values of the channel probabilities
  Cell_ptr->O = Cell_ptr->O + (k1_O + 2*k2_O + 2*k3_O + k4_O)/6;
  Cell_ptr->C1 = Cell_ptr->C1 + (k1_C1 + 2*k2_C1 + 2*k3_C1 + k4_C1)/6;
  Cell_ptr->C2 = Cell_ptr->C2 + (k1_C2 + 2*k2_C2 + 2*k3_C2 + k4_C2)/6;
  Cell_ptr->C3 = Cell_ptr->C3 + (k1_C3 + 2*k2_C3 + 2*k3_C3 + k4_C3)/6;
  Cell_ptr->I = Cell_ptr->I + (k1_I + 2*k2_I + 2*k3_I + k4_I)/6;
  Cell_ptr->O_drug = Cell_ptr->O_drug + (k1_O_drug + 2*k2_O_drug + 2*k3_O_drug + k4_O_drug)/6;
  Cell_ptr->C1_drug = Cell_ptr->C1_drug + (k1_C1_drug + 2*k2_C1_drug + 2*k3_C1_drug + k4_C1_drug)/6;
  Cell_ptr->C2_drug = Cell_ptr->C2_drug + (k1_C2_drug + 2*k2_C2_drug + 2*k3_C2_drug + k4_C2_drug)/6;
  Cell_ptr->C3_drug = Cell_ptr->C3_drug + (k1_C3_drug + 2*k2_C3_drug + 2*k3_C3_drug + k4_C3_drug)/6;
  Cell_ptr->I_drug = Cell_ptr->I_drug + (k1_I_drug + 2*k2_I_drug + 2*k3_I_drug + k4_I_drug)/6;
  
  Cell_ptr->sum = Cell_ptr->O + Cell_ptr->C1 + Cell_ptr->C2 + Cell_ptr->C3 + Cell_ptr->I+ Cell_ptr->O_drug + Cell_ptr->C1_drug + Cell_ptr->C2_drug + Cell_ptr->C3_drug + Cell_ptr->I_drug;


  if (fabs (Cell_ptr->sum -1.0) > 0.0001) {
//    cout << "Error in Wild Type Sum at t = " << t << ",\t sum = " << Cell_ptr->sum << endl;
    *error=true;
    Cell_ptr->I_Kr = 0.0;
    }
  else {
    *error=false;
    Cell_ptr->I_Kr = G_Kr*sqrt(Ko/5.4)*(Cell_ptr->O + Cell_ptr->O_drug)*(Cell_ptr->V - E_Kr);
    }
  return;
}


void curr9(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){

  double E_Kr;
  double hh;
  hh = dt;
  
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0));       //var g_Kr_0: microS_per_nanoF {init: 0.024}; Lucia: Ko dependence is included in the formulation of the current

//        const double Q10=3;    To account for correction in temperature
//        const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));

  E_Kr = (R*T/F)*log(Ko/Ki);

//////Markov mechanism with NS drug////////Option7 --lau--/////

//  C3 = C2 = C1 = O = I
//                 |   
// C3D= C2D= C1D = OD=ID

////////////////////////////////////////////////////

//rate cts. definitions:
//C3_C2=ae
//C2_C3=be
//C2_C1=ain
//C1_C2=bin
//C1_O=aa
//O_C1=bb
//O_I=bi
//I_O=ai
//C1D_OD=aad
//OD_CID=bbd
//C2D_C1D=aind
//C1D_C2D=bind
//C3D_C2D=aed
//C2D_C3D=bed
//*********Drug rates******
//O_OD=kO  ----> (microMolar*ms)-1) 
//OD_O=rO  -----> (ms)-1)
//C3_C3DD=kC  ----> (microMolar*ms)-1)
//C3DD_C3=rO  -----> (ms)-1)
//************************
//Drug equilibrium constants
//K_O=rO/kO
//K_C=rC/kC
////////////////////////////

//channel rates cts ----> 



//WT+ drug (1S)
  double ae=dx[0]*1.50063*exp(24.288586641+dx[1]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*0.598723*exp(13.637047252-dx[3]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*1.10345*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.819984*exp(13.110155065-16.420979385);
  double ai=dx[6]*1.02159*exp(28.782084605-dx[7]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*0.998552*exp(29.690385428+dx[9]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*1.86823*exp(22.010640617+dx[11]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*1.00124*exp(7.266288229-dx[13]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.50063*exp(24.288586641+dx[15]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*0.598723*exp(13.637047252-dx[17]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*1.10345*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.819984*exp(13.110155065-16.420979385);
  double aid=dx[20]*1.02159*exp(28.782084605-dx[21]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bid=dx[22]*0.998552*exp(29.690385428+dx[23]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*1.86823*exp(22.010640617+dx[25]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*1.00124*exp(7.266288229-dx[27]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation


//drug rates cts ---->
  double kO=dx[28];
//double rO=dx[1];
//  double kC=dx[29];
//double rC=dx[4];

//relations based on K and microscopic reversibility
//  double K_C=10.0; // (microM)//drug conscentration is 30 microM
//  double K_O=3.0; //(microM)

//  double K_C=dx[30]; // (microM)//drug conscentration is 30 microM
  double K_O=dx[29]; //(microM)

//  double rC=K_C*kC;
  double rO=K_O*kO;
//double rC=kC*(rO*ae*ai*aa*bbd*bind*bed)/(kO*be*bin*bb*aad*aind*aed)

//INitial conditions**********************************************************************

  // Calculation of k parameters for the Runge Kutta interations
  //Calculation of k1 values using the initial differential equation (set as f(t,V)
  double k1_O = hh*(ai*Cell_ptr->I + aa*Cell_ptr->C1+ rO*Cell_ptr->O_drug - Cell_ptr->O*(bi + bb + kO*Cell_ptr->Drug));
  double k1_C1 = hh*(ain*Cell_ptr->C2 + bb*Cell_ptr->O - Cell_ptr->C1*(aa + bin));
  double k1_C2 = hh*(ae*Cell_ptr->C3 + bin*Cell_ptr->C1 - Cell_ptr->C2*(be + ain));
  double k1_C3 = hh*(be*Cell_ptr->C2 - Cell_ptr->C3*ae);
  double k1_I = hh*(bi*Cell_ptr->O - ai*Cell_ptr->I);
  double k1_O_drug = hh*(aid*Cell_ptr->I_drug + aad*Cell_ptr->C1_drug + kO*Cell_ptr->Drug*Cell_ptr->O - Cell_ptr->O_drug*(bid + bbd + rO));  
  double k1_C1_drug = hh*(aind*Cell_ptr->C2_drug + bbd*Cell_ptr->O_drug - Cell_ptr->C1_drug*(aad + bind ));
  double k1_C2_drug = hh*(aed*Cell_ptr->C3_drug + bind*Cell_ptr->C1_drug - Cell_ptr->C2_drug*(bed + aind));
  double k1_C3_drug = hh*(bed*Cell_ptr->C2_drug - Cell_ptr->C3_drug*aed);
  double k1_I_drug= hh*(bid*Cell_ptr->O_drug - aid*Cell_ptr->I_drug);

  //Value of O after the k1 iteration; used for calculaton of k2
  double O_k1 = Cell_ptr->O + k1_O/2;
  double C1_k1 = Cell_ptr->C1 + k1_C1/2;
  double C2_k1 = Cell_ptr->C2 + k1_C2/2;
  double C3_k1 = Cell_ptr->C3 + k1_C3/2;
  double I_k1 = Cell_ptr->I + k1_I/2;
  double O_drug_k1 = Cell_ptr->O_drug + k1_O_drug/2;
  double C1_drug_k1 = Cell_ptr->C1_drug + k1_C1_drug/2;
  double C2_drug_k1 = Cell_ptr->C2_drug + k1_C2_drug/2;
  double C3_drug_k1 = Cell_ptr->C3_drug + k1_C3_drug/2;
  double I_drug_k1 = Cell_ptr->I_drug + k1_I_drug/2;

  //calculation of K2 values
  double k2_O = hh*(ai*I_k1 + aa*C1_k1 + rO*O_drug_k1 - O_k1*(bi + bb+ kO*Cell_ptr->Drug));
  double k2_C1 = hh*(ain*C2_k1 + bb*O_k1 - C1_k1*(aa + bin ));
  double k2_C2 = hh*(ae*C3_k1 + bin*C1_k1 - C2_k1*(be + ain ));
  double k2_C3 = hh*(be*C2_k1 - C3_k1*ae);
  double k2_I = hh*(bi*O_k1 - ai*I_k1);
  double k2_O_drug = hh*( aid*I_drug_k1 + aad*C1_drug_k1 + kO*Cell_ptr->Drug*O_k1 - O_drug_k1*(bid + bbd + rO)); 
  double k2_C1_drug = hh*(aind*C2_drug_k1 + bbd*O_drug_k1 - C1_drug_k1*(aad + bind));
  double k2_C2_drug = hh*(aed*C3_drug_k1 + bind*C1_drug_k1 - C2_drug_k1*(bed + aind ));
  double k2_C3_drug = hh*(bed*C2_drug_k1 - C3_drug_k1*aed);
  double k2_I_drug= hh*(bid*O_drug_k1 - aid*I_drug_k1);
  double O_k2 = Cell_ptr->O + k2_O/2;
  double C1_k2 = Cell_ptr->C1 + k2_C1/2;
  double C2_k2 = Cell_ptr->C2 + k2_C2/2;
  double C3_k2 = Cell_ptr->C3 + k2_C3/2;
  double I_k2 = Cell_ptr->I + k2_I/2;
  double O_drug_k2 = Cell_ptr->O_drug + k2_O_drug/2;
  double C1_drug_k2 = Cell_ptr->C1_drug + k2_C1_drug/2;
  double C2_drug_k2 = Cell_ptr->C2_drug + k2_C2_drug/2;
  double C3_drug_k2 = Cell_ptr->C3_drug + k2_C3_drug/2;
  double I_drug_k2 = Cell_ptr->I_drug + k2_I_drug/2;
  //calculation of K3 values
  double k3_O = hh*(ai*I_k2 + aa*C1_k2 + rO*O_drug_k2 - O_k2*(bi + bb+ kO*Cell_ptr->Drug));
  double k3_C1 = hh*(ain*C2_k2 + bb*O_k2 - C1_k2*(aa + bin ));
  double k3_C2 = hh*(ae*C3_k2 + bin*C1_k2 - C2_k2*(be + ain));
  double k3_C3 = hh*(be*C2_k2 - C3_k2*ae);
  double k3_I = hh*(bi*O_k2 - ai*I_k2);
  double k3_O_drug = hh*( aid*I_drug_k2 + aad*C1_drug_k2 + kO*Cell_ptr->Drug*O_k2 - O_drug_k2*(bid + bbd + rO));
  double k3_C1_drug = hh*(aind*C2_drug_k2 + bbd*O_drug_k2 - C1_drug_k2*(aad + bind ));
  double k3_C2_drug = hh*(aed*C3_drug_k2 + bind*C1_drug_k2 - C2_drug_k2*(bed + aind));
  double k3_C3_drug = hh*(bed*C2_drug_k2 - C3_drug_k2*aed);
  double k3_I_drug= hh*(bid*O_drug_k2 - aid*I_drug_k2);
  double O_k3 = Cell_ptr->O + k3_O;
  double C1_k3 = Cell_ptr->C1 + k3_C1;
  double C2_k3 = Cell_ptr->C2 + k3_C2;
  double C3_k3 = Cell_ptr->C3 + k3_C3;
  double I_k3 = Cell_ptr->I + k3_I;
  double O_drug_k3 = Cell_ptr->O_drug + k3_O_drug;
  double C1_drug_k3 = Cell_ptr->C1_drug + k3_C1_drug;
  double C2_drug_k3 = Cell_ptr->C2_drug + k3_C2_drug;
  double C3_drug_k3 = Cell_ptr->C3_drug + k3_C3_drug;
  double I_drug_k3 = Cell_ptr->I_drug + k3_I_drug;
  //calculation of K4 values
  double k4_O = hh*(ai*I_k3 + aa*C1_k3 + rO*O_drug_k3 - O_k3*(bi + bb+ kO*Cell_ptr->Drug));
  double k4_C1 = hh*(ain*C2_k3 + bb*O_k3 - C1_k3*(aa + bin ));
  double k4_C2 = hh*(ae*C3_k3 + bin*C1_k3 - C2_k3*(be + ain ));
  double k4_C3 = hh*(be*C2_k3 - C3_k3*ae);
  double k4_I = hh*(bi*O_k3 - ai*I_k3);
  double k4_O_drug = hh*( aid*I_drug_k3 + aad*C1_drug_k3 + kO*Cell_ptr->Drug*O_k3 - O_drug_k3*(bid + bbd + rO));
  double k4_C1_drug = hh*(aind*C2_drug_k3 + bbd*O_drug_k3 - C1_drug_k3*(aad + bind ));
  double k4_C2_drug = hh*(aed*C3_drug_k3 + bind*C1_drug_k3 - C2_drug_k3*(bed + aind));
  double k4_C3_drug = hh*(bed*C2_drug_k3 - C3_drug_k3*aed);
  double k4_I_drug= hh*(bid*O_drug_k3 - aid*I_drug_k3);

  
  //New computed values of the channel probabilities
  Cell_ptr->O = Cell_ptr->O + (k1_O + 2*k2_O + 2*k3_O + k4_O)/6;
  Cell_ptr->C1 = Cell_ptr->C1 + (k1_C1 + 2*k2_C1 + 2*k3_C1 + k4_C1)/6;
  Cell_ptr->C2 = Cell_ptr->C2 + (k1_C2 + 2*k2_C2 + 2*k3_C2 + k4_C2)/6;
  Cell_ptr->C3 = Cell_ptr->C3 + (k1_C3 + 2*k2_C3 + 2*k3_C3 + k4_C3)/6;
  Cell_ptr->I = Cell_ptr->I + (k1_I + 2*k2_I + 2*k3_I + k4_I)/6;
  Cell_ptr->O_drug = Cell_ptr->O_drug + (k1_O_drug + 2*k2_O_drug + 2*k3_O_drug + k4_O_drug)/6;
  Cell_ptr->C1_drug = Cell_ptr->C1_drug + (k1_C1_drug + 2*k2_C1_drug + 2*k3_C1_drug + k4_C1_drug)/6;
  Cell_ptr->C2_drug = Cell_ptr->C2_drug + (k1_C2_drug + 2*k2_C2_drug + 2*k3_C2_drug + k4_C2_drug)/6;
  Cell_ptr->C3_drug = Cell_ptr->C3_drug + (k1_C3_drug + 2*k2_C3_drug + 2*k3_C3_drug + k4_C3_drug)/6;
  Cell_ptr->I_drug = Cell_ptr->I_drug + (k1_I_drug + 2*k2_I_drug + 2*k3_I_drug + k4_I_drug)/6;
  
  Cell_ptr->sum = Cell_ptr->O + Cell_ptr->C1 + Cell_ptr->C2 + Cell_ptr->C3 + Cell_ptr->I+ Cell_ptr->O_drug + Cell_ptr->C1_drug + Cell_ptr->C2_drug + Cell_ptr->C3_drug + Cell_ptr->I_drug;


  if (fabs (Cell_ptr->sum -1.0) > 0.0001) {
//    cout << "Error in Wild Type Sum at t = " << t << ",\t sum = " << Cell_ptr->sum << endl;
    *error=true;
    Cell_ptr->I_Kr = 0.0;
    }
  else {
    *error=false;
    Cell_ptr->I_Kr = G_Kr*sqrt(Ko/5.4)*(Cell_ptr->O + Cell_ptr->O_drug)*(Cell_ptr->V - E_Kr);
    }
  return;
}


void curr10(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){

  double E_Kr;
  double hh;
  hh = dt;
  
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0));       //var g_Kr_0: microS_per_nanoF {init: 0.024}; Lucia: Ko dependence is included in the formulation of the current

//        const double Q10=3;    To account for correction in temperature
//        const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));

  E_Kr = (R*T/F)*log(Ko/Ki);

//////Markov mechanism with NS drug////////Option4 --lau--/////

//  C3 = C2 = C1 = O = I
//       |         |   
// C3D= C2D= C1D = OD= ID

////////////////////////////////////////////////////

//rate cts. definitions:
//C3_C2=ae
//C2_C3=be
//C2_C1=ain
//C1_C2=bin
//C1_O=aa
//O_C1=bb
//O_I=bi
//I_O=ai
//ID_OD=aid
//OD_ID=bid
//C1D_OD=aad
//OD_CID=bbd
//C2D_C1D=aind
//C1D_C2D=bind
//C3D_C2D=aed
//C2D_C3D=bed
//*********Drug rates******
//O_OD=kO  ----> (microMolar*ms)-1) 
//OD_O=rO  -----> (ms)-1)
//C3_C3DD=kC  ----> (microMolar*ms)-1)
//C3DD_C3=rO  -----> (ms)-1)
//************************
//Drug equilibrium constants
//K_O=rO/kO
//K_C=rC/kC
////////////////////////////

//channel rates cts ----> 


/*
//WT + drug (3S)
//
  double ae=dx[0]*1.39429*exp(24.288586641+dx[1]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*1.12436*exp(13.637047252-dx[3]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*0.80694*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.905579*exp(13.110155065-16.420979385);
  double ai=dx[6]*0.989603*exp(28.782084605-dx[7]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bi=dx[8]*1.00033*exp(29.690385428+dx[9]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*0.98968*exp(22.010640617+dx[11]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*0.99917*exp(7.266288229-dx[13]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.39429*exp(24.288586641+dx[15]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*1.12436*exp(13.637047252-dx[17]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*0.80694*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.905579*exp(13.110155065-16.420979385);
  double aid=dx[20]*0.989603*exp(28.782084605-dx[21]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bid=dx[22]*1.00033*exp(29.690385428+dx[23]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*0.98968*exp(22.010640617+dx[25]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*0.99917*exp(7.266288229-dx[27]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
*/
//WT +drug (1S)
  double ae=dx[0]*1.50063*exp(24.288586641+dx[1]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*0.598723*exp(13.637047252-dx[3]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*1.10345*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.819984*exp(13.110155065-16.420979385);
  double ai=dx[6]*1.02159*exp(28.782084605-dx[7]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*0.998552*exp(29.690385428+dx[9]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*1.86823*exp(22.010640617+dx[11]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*1.00124*exp(7.266288229-dx[13]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.50063*exp(24.288586641+dx[15]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*0.598723*exp(13.637047252-dx[17]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*1.10345*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.819984*exp(13.110155065-16.420979385);
  double aid=dx[20]*1.02159*exp(28.782084605-dx[21]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bid=dx[22]*0.998552*exp(29.690385428+dx[23]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*1.86823*exp(22.010640617+dx[25]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*1.00124*exp(7.266288229-dx[27]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation


//drug rates cts ---->
  double kO=dx[28];
//double rO=dx[1];
  double kC=dx[29];
//double rC=dx[4];

//relations based on K and microscopic reversibility
//  double K_C=10.0; // (microM)//drug conscentration is 30 microM
//  double K_O=3.0; //(microM)

  double K_C=dx[30]; // (microM)//drug conscentration is 30 microM
  double K_O=dx[31]; //(microM)

  double rC=K_C*kC;
  double rO=K_O*kO;
//double rC=kC*(rO*ae*ai*aa*bbd*bind*bed)/(kO*be*bin*bb*aad*aind*aed)

//INitial conditions**********************************************************************

  // Calculation of k parameters for the Runge Kutta interations
  //Calculation of k1 values using the initial differential equation (set as f(t,V)
  double k1_O = hh*(ai*Cell_ptr->I + aa*Cell_ptr->C1+ rO*Cell_ptr->O_drug - Cell_ptr->O*(bi + bb + kO*Cell_ptr->Drug));
  double k1_C1 = hh*(ain*Cell_ptr->C2 + bb*Cell_ptr->O - Cell_ptr->C1*(aa + bin));
  double k1_C2 = hh*(ae*Cell_ptr->C3 + bin*Cell_ptr->C1 + rC*Cell_ptr->C2_drug - Cell_ptr->C2*(be + ain + kC*Cell_ptr->Drug));
  double k1_C3 = hh*(be*Cell_ptr->C2 - Cell_ptr->C3*ae);
  double k1_I = hh*(bi*Cell_ptr->O - ai*Cell_ptr->I);
  double k1_O_drug = hh*( aid*Cell_ptr->I_drug + aad*Cell_ptr->C1_drug + kO*Cell_ptr->Drug*Cell_ptr->O - Cell_ptr->O_drug*(bid + bbd+rO));  
  double k1_C1_drug = hh*(aind*Cell_ptr->C2_drug + bbd*Cell_ptr->O_drug - Cell_ptr->C1_drug*(aad + bind));
  double k1_C2_drug = hh*(aed*Cell_ptr->C3_drug + bind*Cell_ptr->C1_drug + kC*Cell_ptr->Drug*Cell_ptr->C2 - Cell_ptr->C2_drug*(bed + aind + rC));
  double k1_C3_drug = hh*(bed*Cell_ptr->C2_drug - Cell_ptr->C3_drug*aed);
  double k1_I_drug= hh*(bid*Cell_ptr->O_drug - aid*Cell_ptr->I_drug);

  //Value of O after the k1 iteration; used for calculaton of k2
  double O_k1 = Cell_ptr->O + k1_O/2;
  double C1_k1 = Cell_ptr->C1 + k1_C1/2;
  double C2_k1 = Cell_ptr->C2 + k1_C2/2;
  double C3_k1 = Cell_ptr->C3 + k1_C3/2;
  double I_k1 = Cell_ptr->I + k1_I/2;
  double O_drug_k1 = Cell_ptr->O_drug + k1_O_drug/2;
  double C1_drug_k1 = Cell_ptr->C1_drug + k1_C1_drug/2;
  double C2_drug_k1 = Cell_ptr->C2_drug + k1_C2_drug/2;
  double C3_drug_k1 = Cell_ptr->C3_drug + k1_C3_drug/2;
  double I_drug_k1 = Cell_ptr->I_drug + k1_I_drug/2;
  //calculation of K2 values
  double k2_O = hh*(ai*I_k1 + aa*C1_k1 + rO*O_drug_k1 - O_k1*(bi + bb+ kO*Cell_ptr->Drug));
  double k2_C1 = hh*(ain*C2_k1 + bb*O_k1 - C1_k1*(aa + bin));
  double k2_C2 = hh*(ae*C3_k1 + bin*C1_k1 +rC*C2_drug_k1 - C2_k1*(be + ain + kC*Cell_ptr->Drug));
  double k2_C3 = hh*(be*C2_k1 - C3_k1*ae);
  double k2_I = hh*(bi*O_k1 - ai*I_k1);
  double k2_O_drug = hh*( aid*I_drug_k1 + aad*C1_drug_k1 + kO*Cell_ptr->Drug*O_k1 - O_drug_k1*(bid + bbd + rO)); 
  double k2_C1_drug = hh*(aind*C2_drug_k1 + bbd*O_drug_k1 - C1_drug_k1*(aad + bind));
  double k2_C2_drug = hh*(aed*C3_drug_k1 + bind*C1_drug_k1 + kC*Cell_ptr->Drug*C2_k1 - C2_drug_k1*(bed + aind + rC));
  double k2_C3_drug = hh*(bed*C2_drug_k1 - C3_drug_k1*aed);
  double k2_I_drug= hh*(bid*O_drug_k1 - aid*I_drug_k1); 
  double O_k2 = Cell_ptr->O + k2_O/2;
  double C1_k2 = Cell_ptr->C1 + k2_C1/2;
  double C2_k2 = Cell_ptr->C2 + k2_C2/2;
  double C3_k2 = Cell_ptr->C3 + k2_C3/2;
  double I_k2 = Cell_ptr->I + k2_I/2;
  double O_drug_k2 = Cell_ptr->O_drug + k2_O_drug/2;
  double C1_drug_k2 = Cell_ptr->C1_drug + k2_C1_drug/2;
  double C2_drug_k2 = Cell_ptr->C2_drug + k2_C2_drug/2;
  double C3_drug_k2 = Cell_ptr->C3_drug + k2_C3_drug/2;
  double I_drug_k2 = Cell_ptr->I_drug + k2_I_drug/2;
  //calculation of K3 values
  double k3_O = hh*(ai*I_k2 + aa*C1_k2 + rO*O_drug_k2 - O_k2*(bi + bb+ kO*Cell_ptr->Drug));
  double k3_C1 = hh*(ain*C2_k2 + bb*O_k2 - C1_k2*(aa + bin));
  double k3_C2 = hh*(ae*C3_k2 + bin*C1_k2 +rC*C2_drug_k2 - C2_k2*(be + ain + kC*Cell_ptr->Drug));
  double k3_C3 = hh*(be*C2_k2 - C3_k2*ae);
  double k3_I = hh*(bi*O_k2 - ai*I_k2);
  double k3_O_drug = hh*( aid*I_drug_k2 + aad*C1_drug_k2 + kO*Cell_ptr->Drug*O_k2 - O_drug_k2*(bid + bbd + rO));
  double k3_C1_drug = hh*(aind*C2_drug_k2 + bbd*O_drug_k2 - C1_drug_k2*(aad + bind));
  double k3_C2_drug = hh*(aed*C3_drug_k2 + bind*C1_drug_k2 + kC*Cell_ptr->Drug*C2_k2 - C2_drug_k2*(bed + aind + rC));
  double k3_C3_drug = hh*(bed*C2_drug_k2 - C3_drug_k2*aed);
  double k3_I_drug= hh*(bid*O_drug_k2 - aid*I_drug_k2);
  double O_k3 = Cell_ptr->O + k3_O;
  double C1_k3 = Cell_ptr->C1 + k3_C1;
  double C2_k3 = Cell_ptr->C2 + k3_C2;
  double C3_k3 = Cell_ptr->C3 + k3_C3;
  double I_k3 = Cell_ptr->I + k3_I;
  double O_drug_k3 = Cell_ptr->O_drug + k3_O_drug;
  double C1_drug_k3 = Cell_ptr->C1_drug + k3_C1_drug;
  double C2_drug_k3 = Cell_ptr->C2_drug + k3_C2_drug;
  double C3_drug_k3 = Cell_ptr->C3_drug + k3_C3_drug;
  double I_drug_k3 = Cell_ptr->I_drug + k3_I_drug;
  //calculation of K4 values
  double k4_O = hh*(ai*I_k3 + aa*C1_k3 + rO*O_drug_k3 - O_k3*(bi + bb+ kO*Cell_ptr->Drug));
  double k4_C1 = hh*(ain*C2_k3 + bb*O_k3 - C1_k3*(aa + bin));
  double k4_C2 = hh*(ae*C3_k3 + bin*C1_k3 +rC*C2_drug_k3 - C2_k3*(be + ain +  kC*Cell_ptr->Drug));
  double k4_C3 = hh*(be*C2_k3 - C3_k3*ae);
  double k4_I = hh*(bi*O_k3 - ai*I_k3);
  double k4_O_drug = hh*( aid*I_drug_k3 + aad*C1_drug_k3 + kO*Cell_ptr->Drug*O_k3 - O_drug_k3*(bid + bbd + rO));
  double k4_C1_drug = hh*(aind*C2_drug_k3 + bbd*O_drug_k3 - C1_drug_k3*(aad + bind));
  double k4_C2_drug = hh*(aed*C3_drug_k3 + bind*C1_drug_k3 + kC*Cell_ptr->Drug*C2_k3 - C2_drug_k3*(bed + aind + rC));
  double k4_C3_drug = hh*(bed*C2_drug_k3 - C3_drug_k3*aed);
  double k4_I_drug= hh*(bid*O_drug_k3 - aid*I_drug_k3);
  
  
  //New computed values of the channel probabilities
  Cell_ptr->O = Cell_ptr->O + (k1_O + 2*k2_O + 2*k3_O + k4_O)/6;
  Cell_ptr->C1 = Cell_ptr->C1 + (k1_C1 + 2*k2_C1 + 2*k3_C1 + k4_C1)/6;
  Cell_ptr->C2 = Cell_ptr->C2 + (k1_C2 + 2*k2_C2 + 2*k3_C2 + k4_C2)/6;
  Cell_ptr->C3 = Cell_ptr->C3 + (k1_C3 + 2*k2_C3 + 2*k3_C3 + k4_C3)/6;
  Cell_ptr->I = Cell_ptr->I + (k1_I + 2*k2_I + 2*k3_I + k4_I)/6;
  Cell_ptr->O_drug = Cell_ptr->O_drug + (k1_O_drug + 2*k2_O_drug + 2*k3_O_drug + k4_O_drug)/6;
  Cell_ptr->C1_drug = Cell_ptr->C1_drug + (k1_C1_drug + 2*k2_C1_drug + 2*k3_C1_drug + k4_C1_drug)/6;
  Cell_ptr->C2_drug = Cell_ptr->C2_drug + (k1_C2_drug + 2*k2_C2_drug + 2*k3_C2_drug + k4_C2_drug)/6;
  Cell_ptr->C3_drug = Cell_ptr->C3_drug + (k1_C3_drug + 2*k2_C3_drug + 2*k3_C3_drug + k4_C3_drug)/6;
  Cell_ptr->I_drug = Cell_ptr->I_drug + (k1_I_drug + 2*k2_I_drug + 2*k3_I_drug + k4_I_drug)/6;
  
  Cell_ptr->sum = Cell_ptr->O + Cell_ptr->C1 + Cell_ptr->C2 + Cell_ptr->C3 + Cell_ptr->I+ Cell_ptr->O_drug + Cell_ptr->C1_drug + Cell_ptr->C2_drug + Cell_ptr->C3_drug + Cell_ptr->I_drug;


  if (fabs (Cell_ptr->sum -1.0) > 0.0001) {
//    cout << "Error in Wild Type Sum at t = " << t << ",\t sum = " << Cell_ptr->sum << endl;
    *error=true;
    Cell_ptr->I_Kr = 0.0;
    }
  else {
    *error=false;
    Cell_ptr->I_Kr = G_Kr*sqrt(Ko/5.4)*(Cell_ptr->O + Cell_ptr->O_drug)*(Cell_ptr->V - E_Kr);
    }
  return;
}

 
void curr11(Cell_param *Cell_ptr, double t, double dt, double dx[], bool *error){

  double E_Kr;
  double hh;
  hh = dt;
  
  const double G_Kr = 0.024*(T/35.0-(55.0/7.0));       //var g_Kr_0: microS_per_nanoF {init: 0.024}; Lucia: Ko dependence is included in the formulation of the current

//        const double Q10=3;    To account for correction in temperature
//        const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));

  E_Kr = (R*T/F)*log(Ko/Ki);

//////Markov mechanism with NS drug////////Option4 --lau--/////

//  C3 = C2 = C1 = O = I
//  |    |            
// C3D= C2D= C1D = OD= ID

////////////////////////////////////////////////////

//rate cts. definitions:
//C3_C2=ae
//C2_C3=be
//C2_C1=ain
//C1_C2=bin
//C1_O=aa
//O_C1=bb
//O_I=bi
//I_O=ai
//ID_OD=aid
//OD_ID=bid
//C1D_OD=aad
//OD_CID=bbd
//C2D_C1D=aind
//C1D_C2D=bind
//C3D_C2D=aed
//C2D_C3D=bed
//*********Drug rates******
//O_OD=kO  ----> (microMolar*ms)-1) 
//OD_O=rO  -----> (ms)-1)
//C3_C3DD=kC  ----> (microMolar*ms)-1)
//C3DD_C3=rO  -----> (ms)-1)
//************************
//Drug equilibrium constants
//K_O=rO/kO
//K_C=rC/kC
////////////////////////////

//channel rates cts ----> 


/*
//WT + drug (3S)
//
  double ae=dx[0]*1.39429*exp(24.288586641+dx[1]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*1.12436*exp(13.637047252-dx[3]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*0.80694*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.905579*exp(13.110155065-16.420979385);
  double ai=dx[6]*0.989603*exp(28.782084605-dx[7]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bi=dx[8]*1.00033*exp(29.690385428+dx[9]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*0.98968*exp(22.010640617+dx[11]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*0.99917*exp(7.266288229-dx[13]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.39429*exp(24.288586641+dx[15]*0.51998*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*1.12436*exp(13.637047252-dx[17]*1.09902*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*0.80694*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.905579*exp(13.110155065-16.420979385);
  double aid=dx[20]*0.989603*exp(28.782084605-dx[21]*1.00013*0.0327*Cell_ptr->V-34.3968626);
  double bid=dx[22]*1.00033*exp(29.690385428+dx[23]*0.995342*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*0.98968*exp(22.010640617+dx[25]*0.997978*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*0.99917*exp(7.266288229-dx[27]*1.00052*0.0418*Cell_ptr->V-16.168847047);  //deactivation
*/
//WT +drug (1S)
  double ae=dx[0]*1.50063*exp(24.288586641+dx[1]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double be=dx[2]*0.598723*exp(13.637047252-dx[3]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double ain=dx[4]*1.10345*exp(22.444099764-27.130772208);
  double bin=dx[5]*0.819984*exp(13.110155065-16.420979385);
  double ai=dx[6]*1.02159*exp(28.782084605-dx[7]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bi=dx[8]*0.998552*exp(29.690385428+dx[9]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aa=dx[10]*1.86823*exp(22.010640617+dx[11]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bb=dx[12]*1.00124*exp(7.266288229-dx[13]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation
  double aed=dx[14]*1.50063*exp(24.288586641+dx[15]*1.81284*0.0117*Cell_ptr->V-27.138936313);
  double bed=dx[16]*0.598723*exp(13.637047252-dx[17]*1.16347*0.0631*Cell_ptr->V-16.448205938);
  double aind=dx[18]*1.10345*exp(22.444099764-27.130772208);
  double bind=dx[19]*0.819984*exp(13.110155065-16.420979385);
  double aid=dx[20]*1.02159*exp(28.782084605-dx[21]*1.00008*0.0327*Cell_ptr->V-34.3968626);  // ai, bi steps different for the mutant.
  double bid=dx[22]*0.998552*exp(29.690385428+dx[23]*0.991727*0.0233*Cell_ptr->V-31.764996563)*pow((5.4/Ko),0.4);   //inactivation
  double aad=dx[24]*1.86823*exp(22.010640617+dx[25]*1.14629*0.03822*Cell_ptr->V-27.091932598);  //activation  
  double bbd=dx[26]*1.00124*exp(7.266288229-dx[27]*0.999428*0.0418*Cell_ptr->V-16.168847047);  //deactivation


//drug rates cts ---->
  double kC=dx[28];

//relations based on K and microscopic reversibility
//  double K_C=10.0; // (microM)//drug conscentration is 30 microM
//  double K_O=3.0; //(microM)

  double K_C=dx[29]; // (microM)//drug conscentration is 30 microM
//  double K_O=dx[31]; //(microM)

  double rC=K_C*kC;
//double rC=kC*(rO*ae*ai*aa*bbd*bind*bed)/(kO*be*bin*bb*aad*aind*aed)

//INitial conditions**********************************************************************

  // Calculation of k parameters for the Runge Kutta interations
  //Calculation of k1 values using the initial differential equation (set as f(t,V)
  double k1_O = hh*(ai*Cell_ptr->I + aa*Cell_ptr->C1- Cell_ptr->O*(bi + bb));
  double k1_C1 = hh*(ain*Cell_ptr->C2 + bb*Cell_ptr->O - Cell_ptr->C1*(aa + bin));
  double k1_C2 = hh*(ae*Cell_ptr->C3 + bin*Cell_ptr->C1 + rC*Cell_ptr->C2_drug - Cell_ptr->C2*(be + ain + kC*Cell_ptr->Drug));
  double k1_C3 = hh*(be*Cell_ptr->C2 + rC*Cell_ptr->C3_drug - Cell_ptr->C3*(ae+kC*Cell_ptr->Drug));
  double k1_I = hh*(bi*Cell_ptr->O - ai*Cell_ptr->I);
  double k1_O_drug = hh*( aid*Cell_ptr->I_drug + aad*Cell_ptr->C1_drug - Cell_ptr->O_drug*(bid + bbd));  
  double k1_C1_drug = hh*(aind*Cell_ptr->C2_drug + bbd*Cell_ptr->O_drug - Cell_ptr->C1_drug*(aad + bind));
  double k1_C2_drug = hh*(aed*Cell_ptr->C3_drug + bind*Cell_ptr->C1_drug + kC*Cell_ptr->Drug*Cell_ptr->C2 - Cell_ptr->C2_drug*(bed + aind + rC));
  double k1_C3_drug = hh*(bed*Cell_ptr->C2_drug +kC*Cell_ptr->Drug*Cell_ptr->C3 - Cell_ptr->C3_drug*(aed + rC ));
  double k1_I_drug= hh*(bid*Cell_ptr->O_drug - aid*Cell_ptr->I_drug);

  //Value of O after the k1 iteration; used for calculaton of k2
  double O_k1 = Cell_ptr->O + k1_O/2;
  double C1_k1 = Cell_ptr->C1 + k1_C1/2;
  double C2_k1 = Cell_ptr->C2 + k1_C2/2;
  double C3_k1 = Cell_ptr->C3 + k1_C3/2;
  double I_k1 = Cell_ptr->I + k1_I/2;
  double O_drug_k1 = Cell_ptr->O_drug + k1_O_drug/2;
  double C1_drug_k1 = Cell_ptr->C1_drug + k1_C1_drug/2;
  double C2_drug_k1 = Cell_ptr->C2_drug + k1_C2_drug/2;
  double C3_drug_k1 = Cell_ptr->C3_drug + k1_C3_drug/2;
  double I_drug_k1 = Cell_ptr->I_drug + k1_I_drug/2;
  //calculation of K2 values
  double k2_O = hh*(ai*I_k1 + aa*C1_k1 - O_k1*(bi + bb));
  double k2_C1 = hh*(ain*C2_k1 + bb*O_k1 - C1_k1*(aa + bin));
  double k2_C2 = hh*(ae*C3_k1 + bin*C1_k1 +rC*C2_drug_k1 - C2_k1*(be + ain + kC*Cell_ptr->Drug));
  double k2_C3 = hh*(be*C2_k1 +rC*C3_drug_k1 - C3_k1*(ae+ kC*Cell_ptr->Drug));
  double k2_I = hh*(bi*O_k1 - ai*I_k1);
  double k2_O_drug = hh*( aid*I_drug_k1 + aad*C1_drug_k1 - O_drug_k1*(bid + bbd)); 
  double k2_C1_drug = hh*(aind*C2_drug_k1 + bbd*O_drug_k1 - C1_drug_k1*(aad + bind));
  double k2_C2_drug = hh*(aed*C3_drug_k1 + bind*C1_drug_k1 + kC*Cell_ptr->Drug*C2_k1 - C2_drug_k1*(bed + aind + rC));
  double k2_C3_drug = hh*(bed*C2_drug_k1 + kC*Cell_ptr->Drug*C3_k1 - C3_drug_k1*(aed + rC));
  double k2_I_drug= hh*(bid*O_drug_k1 - aid*I_drug_k1); 
  double O_k2 = Cell_ptr->O + k2_O/2;
  double C1_k2 = Cell_ptr->C1 + k2_C1/2;
  double C2_k2 = Cell_ptr->C2 + k2_C2/2;
  double C3_k2 = Cell_ptr->C3 + k2_C3/2;
  double I_k2 = Cell_ptr->I + k2_I/2;
  double O_drug_k2 = Cell_ptr->O_drug + k2_O_drug/2;
  double C1_drug_k2 = Cell_ptr->C1_drug + k2_C1_drug/2;
  double C2_drug_k2 = Cell_ptr->C2_drug + k2_C2_drug/2;
  double C3_drug_k2 = Cell_ptr->C3_drug + k2_C3_drug/2;
  double I_drug_k2 = Cell_ptr->I_drug + k2_I_drug/2;
  //calculation of K3 values
  double k3_O = hh*(ai*I_k2 + aa*C1_k2 - O_k2*(bi + bb));
  double k3_C1 = hh*(ain*C2_k2 + bb*O_k2 - C1_k2*(aa + bin));
  double k3_C2 = hh*(ae*C3_k2 + bin*C1_k2 +rC*C2_drug_k2 - C2_k2*(be + ain + kC*Cell_ptr->Drug));
  double k3_C3 = hh*(be*C2_k2 +rC*C3_drug_k2 - C3_k2*(ae+ kC*Cell_ptr->Drug));
  double k3_I = hh*(bi*O_k2 - ai*I_k2);
  double k3_O_drug = hh*( aid*I_drug_k2 + aad*C1_drug_k2 - O_drug_k2*(bid + bbd));
  double k3_C1_drug = hh*(aind*C2_drug_k2 + bbd*O_drug_k2 - C1_drug_k2*(aad + bind));
  double k3_C2_drug = hh*(aed*C3_drug_k2 + bind*C1_drug_k2 + kC*Cell_ptr->Drug*C2_k2 - C2_drug_k2*(bed + aind + rC));
  double k3_C3_drug = hh*(bed*C2_drug_k2 + kC*Cell_ptr->Drug*C3_k2 - C3_drug_k2*(aed + rC));
  double k3_I_drug= hh*(bid*O_drug_k2 - aid*I_drug_k2);
  double O_k3 = Cell_ptr->O + k3_O;
  double C1_k3 = Cell_ptr->C1 + k3_C1;
  double C2_k3 = Cell_ptr->C2 + k3_C2;
  double C3_k3 = Cell_ptr->C3 + k3_C3;
  double I_k3 = Cell_ptr->I + k3_I;
  double O_drug_k3 = Cell_ptr->O_drug + k3_O_drug;
  double C1_drug_k3 = Cell_ptr->C1_drug + k3_C1_drug;
  double C2_drug_k3 = Cell_ptr->C2_drug + k3_C2_drug;
  double C3_drug_k3 = Cell_ptr->C3_drug + k3_C3_drug;
  double I_drug_k3 = Cell_ptr->I_drug + k3_I_drug;
  //calculation of K4 values
  double k4_O = hh*(ai*I_k3 + aa*C1_k3 - O_k3*(bi + bb));
  double k4_C1 = hh*(ain*C2_k3 + bb*O_k3 - C1_k3*(aa + bin));
  double k4_C2 = hh*(ae*C3_k3 + bin*C1_k3 +rC*C2_drug_k3 - C2_k3*(be + ain +  kC*Cell_ptr->Drug));
  double k4_C3 = hh*(be*C2_k3 +rC*C3_drug_k3 - C3_k3*(ae+ kC*Cell_ptr->Drug));
  double k4_I = hh*(bi*O_k3 - ai*I_k3);
  double k4_O_drug = hh*( aid*I_drug_k3 + aad*C1_drug_k3 - O_drug_k3*(bid + bbd));
  double k4_C1_drug = hh*(aind*C2_drug_k3 + bbd*O_drug_k3 - C1_drug_k3*(aad + bind));
  double k4_C2_drug = hh*(aed*C3_drug_k3 + bind*C1_drug_k3 + kC*Cell_ptr->Drug*C2_k3 - C2_drug_k3*(bed + aind + rC));
  double k4_C3_drug = hh*(bed*C2_drug_k3 + kC*Cell_ptr->Drug*C3_k3 - C3_drug_k3*(aed + rC));
  double k4_I_drug= hh*(bid*O_drug_k3 - aid*I_drug_k3);
  
  
  //New computed values of the channel probabilities
  Cell_ptr->O = Cell_ptr->O + (k1_O + 2*k2_O + 2*k3_O + k4_O)/6;
  Cell_ptr->C1 = Cell_ptr->C1 + (k1_C1 + 2*k2_C1 + 2*k3_C1 + k4_C1)/6;
  Cell_ptr->C2 = Cell_ptr->C2 + (k1_C2 + 2*k2_C2 + 2*k3_C2 + k4_C2)/6;
  Cell_ptr->C3 = Cell_ptr->C3 + (k1_C3 + 2*k2_C3 + 2*k3_C3 + k4_C3)/6;
  Cell_ptr->I = Cell_ptr->I + (k1_I + 2*k2_I + 2*k3_I + k4_I)/6;
  Cell_ptr->O_drug = Cell_ptr->O_drug + (k1_O_drug + 2*k2_O_drug + 2*k3_O_drug + k4_O_drug)/6;
  Cell_ptr->C1_drug = Cell_ptr->C1_drug + (k1_C1_drug + 2*k2_C1_drug + 2*k3_C1_drug + k4_C1_drug)/6;
  Cell_ptr->C2_drug = Cell_ptr->C2_drug + (k1_C2_drug + 2*k2_C2_drug + 2*k3_C2_drug + k4_C2_drug)/6;
  Cell_ptr->C3_drug = Cell_ptr->C3_drug + (k1_C3_drug + 2*k2_C3_drug + 2*k3_C3_drug + k4_C3_drug)/6;
  Cell_ptr->I_drug = Cell_ptr->I_drug + (k1_I_drug + 2*k2_I_drug + 2*k3_I_drug + k4_I_drug)/6;
  
  Cell_ptr->sum = Cell_ptr->O + Cell_ptr->C1 + Cell_ptr->C2 + Cell_ptr->C3 + Cell_ptr->I+ Cell_ptr->O_drug + Cell_ptr->C1_drug + Cell_ptr->C2_drug + Cell_ptr->C3_drug + Cell_ptr->I_drug;


  if (fabs (Cell_ptr->sum -1.0) > 0.0001) {
//    cout << "Error in Wild Type Sum at t = " << t << ",\t sum = " << Cell_ptr->sum << endl;
    *error=true;
    Cell_ptr->I_Kr = 0.0;
    }
  else {
    *error=false;
    Cell_ptr->I_Kr = G_Kr*sqrt(Ko/5.4)*(Cell_ptr->O + Cell_ptr->O_drug)*(Cell_ptr->V - E_Kr);
    }
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


