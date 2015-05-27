// Copyright (C) 2014 CCMA@PSU Maximilian Metti, Xiaozhe Hu


//*************************************************
//
//    Still need to improve the Convergence
//    criterion; the current implementation OF
//    relative residual will lead to oversolving
//
//**************************************************

#include <iostream>
#include <fstream>
#include <dolfin.h>
#include <sys/time.h>
#include <string.h>
#include "./include/two_species.h"
#include "./include/energynorm.h"
#include "./include/dissipation.h"
#include "./include/PoissonCellMarker.h"
#include "./include/GradientRecovery.h"
extern "C"
{
#include "fasp.h"
#include "fasp_functs.h"
    INT fasp_solver_bdcsr_krylov_block_3(block_dCSRmat *A,
                                       dvector *b,
                                       dvector *x,
                                       itsolver_param *itparam,
                                       AMG_param *amgparam,
                                       dCSRmat *A_diag);
#define FASP_BSR     ON  /** use BSR format in fasp */
}

double Dx;
double Dy;
double Dz;
double Lx;
double Ly;
double Lz;
double T;
int Nx;
int Ny;
int Nz;
int Nt;

double epsVal;
double DpVal;
double DnVal;

double tol;
uint   maxit;
double mu;
double adaptTol;



using namespace std;
using namespace dolfin;


//////////////////////////////////
//                              //
//      Function Definitions    //
//      of Initial Guesses      //
//                              //
//////////////////////////////////

//  Initial Cation Number Density Profile
class Cation : public Expression
{
    void eval(Array<double>& values, const Array<double>& x) const
    {
        double log10  = 2.30258509299;
        values[0] = log( 1.0*(x[0]+Lx)/(2.*Lx) + 0.1*(x[0]-Lx)/(-2.*Lx) );
        //exp( log10*(x[0]-Lx)/(2.*Lx));
    }
};

//  Initial Anion Number Density Profile
class Anion : public Expression
{
    void eval(Array<double>& values, const Array<double>& x) const
    {
        double log10  = 2.30258509299;
        values[0] = log( 0.1*(x[0]+Lx)/(2.*Lx) + 1.0*(x[0]-Lx)/(-2.*Lx) );
        //exp(-log10*(x[0]+Lx)/(2.*Lx));
    }
};

//  Initial Electro-Static Potential Energy
class ESPotential : public Expression
{
    void eval(Array<double>& values, const Array<double>& x) const
    {
        double log10  = 2.30258509299;
        double alpha  = (3.6*Lx*Lx-epsVal*log10*log10)/(epsVal*Lx*log10*log10);
        
        double p = exp( log10*(x[0]-Lx)/(2.*Lx));
        double n = exp(-log10*(x[0]+Lx)/(2.*Lx));
        
        values[0] = -x[0]/Lx;//alpha*x[0] - 4.*Lx*Lx*(p-n)/(epsVal*log10*log10);
    }
};






//////////////////////////
//                      //
//  Dirichlet BCs and   //
//  exact solutions     //
//                      //
//////////////////////////

// To solve PNP with Dirichlet BCs
class cationSource : public Expression
{
    void eval(Array<double>& values, const Array<double>& x) const
    {
        double log10  = 2.30258509299;
        double alpha  = (3.6*Lx*Lx-epsVal*log10*log10)/(epsVal*Lx*log10*log10);
        double p = exp( log10*(x[0]-Lx)/(2.*Lx));
        double n = exp(-log10*(x[0]+Lx)/(2.*Lx));
        
        
        values[0] = DpVal*( log10*log10*p/(4.*Lx*Lx)
                           + (log10*p/(2.*Lx)) * (alpha-(2.*Lx)*(p+n)/(epsVal*log10))
                           + p * (-(p-n)/epsVal) );
    }
};

class anionSource : public Expression
{
    void eval(Array<double>& values, const Array<double>& x) const
    {
        double log10  = 2.30258509299;
        double alpha  = (3.6*Lx*Lx-epsVal*log10*log10)/(epsVal*Lx*log10*log10);
        double p = exp( log10*(x[0]-Lx)/(2.*Lx));
        double n = exp(-log10*(x[0]+Lx)/(2.*Lx));
        
        values[0] = DnVal*( log10*log10*n/(4.*Lx*Lx)
                           - (-log10*n/(2.*Lx)) * (alpha-(2.*Lx)*(p+n)/(epsVal*log10))
                           - n * (-(p-n)/epsVal) );
    }
};

// p(0)=0.1;p(1)=1;
class dirichlet_exact_p : public Expression
{
    void eval(Array<double>& values, const Array<double>& x) const
    {
        double log10  = 2.30258509299;
        
        values[0] = exp(log10*(x[0]-Lx)/(2.*Lx));
    }
};

// n(0)=1;n(1)=0.1;
class dirichlet_exact_n : public Expression
{
    void eval(Array<double>& values, const Array<double>& x) const
    {
        double log10  = 2.30258509299;
        
        values[0] = exp(-log10*(x[0]+Lx)/(2.*Lx));
    }
};

// phi(0)=1;phi(1)=-1;
class dirichlet_exact_phi : public Expression
{
    void eval(Array<double>& values, const Array<double>& x) const
    {
        double log10  = 2.30258509299;
        double alpha  = (3.6*Lx*Lx-epsVal*log10*log10)/(epsVal*Lx*log10*log10);
        
        double p = exp( log10*(x[0]-Lx)/(2.*Lx));
        double n = exp(-log10*(x[0]+Lx)/(2.*Lx));
        
        //values[0] = alpha*x[0] - 4.*Lx*Lx*(p-n)/(epsVal*log10*log10);
        values[0] = (exp( x[0])-exp(-x[0]))/(exp(Lx)-exp(-Lx));
    }
};






//////////////////////
//                  //
//  Dirichlet BCs   //
//                  //
//////////////////////


// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return on_boundary && (x[0] < -Lx+2.*DOLFIN_EPS or x[0] > Lx-2.*DOLFIN_EPS);
    }
};



// Sub domain for homogeneous channel wall
class dielectricChannel : public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        bool toppatches = ((   (x[0] < -10./3.+DOLFIN_EPS) or (fabs(x[0]+5./6.) < 5./6.+DOLFIN_EPS)
                            or (fabs(x[0]-15./6.) < 5./6.+DOLFIN_EPS))
                           and x[2] > Lz - DOLFIN_EPS  );
        
        bool bottompatches = ((   (x[0] > 10./3.-DOLFIN_EPS) or (fabs(x[0]-5./6.) < 5./6.+DOLFIN_EPS)
                               or (fabs(x[0]+15./6.) < 5./6.+DOLFIN_EPS))
                              and x[2] < -Lz + DOLFIN_EPS  );
        
        return ( on_boundary && (toppatches or bottompatches) );
        //                &&  (x[1] < -Ly+DOLFIN_EPS or x[1] > Ly-DOLFIN_EPS
        //                   or   x[2] < -Lz+DOLFIN_EPS or x[2] > Lz-DOLFIN_EPS) );
    }
};
// Sub domain for homogeneous channel wall
class channelGate : public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return ( on_boundary && (x[2] < -Lz + DOLFIN_EPS or x[2] > Lz - DOLFIN_EPS) );
        //                && (x[1] < -Ly+DOLFIN_EPS or x[1] > Ly-DOLFIN_EPS
        //                  or  x[2] < -Lz+DOLFIN_EPS or x[2] > Lz-DOLFIN_EPS) );
    }
};




//////////////////////////
//                      //
//      Main Program    //
//                      //
//////////////////////////

int main()
{
    
    ////////////////////////////////////
    //                                //
    //    Setup the environment       //
    //    initialize the problem      //
    //                                //
    ////////////////////////////////////
    
    // Set linear algebra backend
    parameters["linear_algebra_backend"] = "uBLAS";
    
    printf(" \n-------------------------------------\n"); fflush(stdout);
    printf(" This code adaptively solves the steady \n"); fflush(stdout);
    printf(" Poisson-Nernst-Planck system of a \n"); fflush(stdout);
    printf(" single cation and single anion \n"); fflush(stdout);
    printf("-------------------------------------\n\n"); fflush(stdout);
    
    
    
    
    printf(" \n-------------------------------------\n"); fflush(stdout);
    printf(" Initializing the problem \n"); fflush(stdout);
    printf("-------------------------------------\n\n"); fflush(stdout);
    
    
    //***********************
    //  Read in parameters
    //***********************
    
    printf("Read in parameters for the solver and describing the PDE\n");
    fflush(stdout);
    char buffer[500];   // max number of char for each line
    int  val;
    ifstream expin;
    char paramRegime[128]; // output directory
    char meshIn[128];      // mesh input file
    char surfIn[128];      // mesh surfaces file
    char subdIn[128];      // mesh subdomains file
    
    char filenm[] = "./params/exp-names.dat";
    
    // if input file is not specified, use the default values
    if (filenm==NULL) {
        printf("### ERROR: No file specified for input params \n");
        exit(0);
    }
    
    FILE *fp = fopen(filenm,"r");
    if (fp==NULL) {
        printf("### ERROR: Could not open file %s...\n", filenm);
        fasp_chkerr(ERROR_OPEN_FILE, "fasp_param_input");
        
    }
    
    bool state = true;
    while ( state ) {
        int     ibuff;
        double  dbuff;
        char    sbuff[500];
        char   *fgetsPtr;
        
        val = fscanf(fp,"%s",buffer);
        if (val==EOF) break;
        if (val!=1){ state = false; break; }
        if (buffer[0]=='[' || buffer[0]=='%' || buffer[0]=='|') {
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
            continue;
        }
        
        // match keyword and scan for value
        if (strcmp(buffer,"outdir")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%s",sbuff);
            if (val!=1) { state = false; break; }
            strncpy(paramRegime,sbuff,128);
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"mesh_file")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%s",sbuff);
            if (val!=1) { state = false; break; }
            strncpy(meshIn,sbuff,128);
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }

        else if (strcmp(buffer,"surf_file")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%s",sbuff);
            if (val!=1) { state = false; break; }
            strncpy(surfIn,sbuff,128);
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"subd_file")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%s",sbuff);
            if (val!=1) { state = false; break; }
            strncpy(subdIn,sbuff,128);
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"x_length")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            Dx = dbuff; Lx = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"y_length")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            Dy = dbuff; Ly = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"z_length")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            Dz = dbuff; Lz = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"t_length")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            T = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"x_grid")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { state = false; break; }
            Nx = ibuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"y_grid")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { state = false; break; }
            Ny = ibuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"z_grid")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { state = false; break; }
            Nz = ibuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"t_grid")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { state = false; break; }
            Nt = ibuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        
        else if (strcmp(buffer,"nonlin_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            tol = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"nonlin_maxit")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { state = false; break; }
            maxit = ibuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"damp_factor")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            mu = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"adapt_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            adaptTol = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"rel_perm")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            epsVal = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"cat_diff")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            DpVal = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"ani_diff")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            DnVal = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else {
            state = false;
            printf(" Bad read-in \n\n"); fflush(stdout);
        }
        
        
    }
    
    fclose(fp);
    
    printf("\nSolver parameters \n"); fflush(stdout);
    printf("\t adaptivity tol           = %e \n",adaptTol); fflush(stdout);
    printf("\t nonlinear solver tol     = %e \n",tol); fflush(stdout);
    printf("\t nonlinear solver max it  = %d \n",maxit); fflush(stdout);
    printf("\t nonlinear solver damping = %e \n\n",mu); fflush(stdout);
    
    
    
    //*****************************
    //  Open files to write data
    //*****************************
    
    File meshfile("./output/adapt/mesh.pvd");
    File pfile("./output/adapt/pFinal.pvd");
    File nfile("./output/adapt/nFinal.pvd");
    File phifile("./output/adapt/phiFinal.pvd");
    File pErrorFile("./output/adapt/pError.pvd");
    File nErrorFile("./output/adapt/nError.pvd");
    File phiErrorFile("./output/adapt/phiError.pvd");
    
    
    //***********
    //  Domain
    //***********
    
    Mesh initMesh;
    printf("Initialize mesh: "); fflush(stdout);
    if ( strcmp(meshIn,"box")==0 ) {
      printf(" Domain set to [ %6.3f, %6.3f, %6.3f ] \n\n", Dx, Dy, Dz);
      printf(" Initial mesh is %d x %d x %d \n", Nx,Ny,Nz); fflush(stdout);
      BoxMesh bMesh(-Lx,-Ly,-Lz, Lx, Ly, Lz, Nx, Ny, Nz);
      initMesh = bMesh;
    }
    else {
      printf(" Reading in the mesh from %s \n", meshIn);
      Mesh rMesh(meshIn);
      initMesh = rMesh;
    }
    
    
    
    
    //************************
    //  Dimensional Analysis
    //************************
    
    printf("\nSet and rescale coefficients for the PDE \n"); fflush(stdout);
    printf(" Dp  = %e \n",  DpVal);
    printf(" Dn  = %e \n",  DnVal);
    printf(" eps = %e \n", epsVal); fflush(stdout);
    
    
    
    
    //************************
    //  Analytic Expressions
    //************************
    
    // Cation boundary interpolant
      //dirichlet_exact_p pp;
      Cation pp;
    // Anion boundary interpolant
      //dirichlet_exact_n nn;//
      Anion nn;
    // Electric Potential boundary interpolant
      dirichlet_exact_phi kk;//
      //ESPotential kk;
    
    
    
    
    
    
    //***************
    //  Time March
    //***************
    
    
    
    
    
    
    
    //**************
    //  Adaptivity
    //**************
    
    // Copy mesh and initialize CG space
    Mesh mesh0(initMesh);
    two_species::FunctionSpace W0(mesh0);
    uint totalRefines = 0;
    
    // Define initial guesses and residual
    Function initGuess(W0);
    Function initP(initGuess[0]);  initP.interpolate(pp);
    Function initN(initGuess[1]);  initN.interpolate(nn);
    Function initK(initGuess[2]);  initK.interpolate(kk);

    
    // Initialize adaptivity marker
    Function adaptFn0(W0);
    Function adaptFN(adaptFn0[2]);
    adaptFN.interpolate(initK);
    double b0_fasp = -1.0;
    bool converged = false;
    
    
    // Adapt mesh to current iterate
    for ( uint adaptInd=0; adaptInd<15; adaptInd++ ) {
    
    bool tooCoarse  = true;
    bool unrefined  = true;
    uint numRefines = 0;
    Mesh mesh(mesh0);
    
        
    // Coarsen until all elements are sufficiently refined
    printf("Refine mesh until recovered gradient is sufficiently accurate\n");
    fflush(stdout);
    while (tooCoarse) {
        
        // Gradient recovery
        printf(" Recovering gradient of electric potential... \n"); fflush(stdout);
          Mesh meshAdapt(mesh);
          GradientRecovery::FunctionSpace GR(meshAdapt);
          GradientRecovery::BilinearForm aGR(GR,GR);
          GradientRecovery::LinearForm LGR(GR);
          LGR.u = adaptFN;
          Function Dsoln(GR);
        //solve(aGR==LGR,Dsoln);
        
        
        //**************************************************
        //  Solve linear problem: interface with FASP
        //**************************************************
        
        Matrix adaptA; assemble(adaptA,aGR);
        Vector adaptb; assemble(adaptb,LGR);
        Function adaptSolu(GR);
        
        printf(" -------------------------------------\n"); fflush(stdout);
        printf(" Start interface to FASP  \n"); fflush(stdout);
        printf(" -------------------------------------\n\n"); fflush(stdout);
        
        printf(" Step 1: convert sparse matrix format and lump\n"); fflush(stdout);
        // Convert adaptA to CSR
        dCSRmat adaptA_fasp;
        
        unsigned int adaptnz = boost::tuples::get<3>(adaptA.data());
        int adaptrow = adaptA.size(0);
        int adaptcol = adaptA.size(1);
        int* adaptap = (int*)fasp_mem_calloc(adaptrow+1, sizeof(int));
        const size_t* ap_tmp = boost::tuples::get<0>(adaptA.data());
        for (int i=0; i<adaptrow+1; i++) {
            adaptap[i] = (int)ap_tmp[i];
        }
        int* adaptai = (int*)fasp_mem_calloc(adaptnz, sizeof(int));
        const size_t* ai_tmp = boost::tuples::get<1>(adaptA.data());
        for (int i=0; i<adaptnz; i++) {
            adaptai[i] = (int)ai_tmp[i];
        }
        double* adaptax = (double*)boost::tuples::get<2>(adaptA.data());
        
        // Lump matrix
        for ( uint rowInd=0; rowInd<adaptrow; rowInd++ ) {
            double  rowSum = 0.;
            int diagColInd = -1;
            for ( uint colInd=adaptap[rowInd]; colInd < adaptap[rowInd+1]; colInd++ ) {
                rowSum += adaptax[colInd];
                adaptax[colInd] = 0.0;
                if ( adaptai[colInd] == rowInd ) diagColInd = colInd;
            }
            adaptax[diagColInd] = rowSum;
            //printf(" adapt_A[%d,%d] = %e \n", rowInd, adaptai[diagColInd], adaptax[diagColInd]);
        }
        
        adaptA_fasp.row = adaptrow;
        adaptA_fasp.col = adaptcol;
        adaptA_fasp.nnz = adaptnz;
        adaptA_fasp.IA  = adaptap;
        adaptA_fasp.JA  = adaptai;
        adaptA_fasp.val = adaptax;
        
        // initialize RHS
        dvector adaptb_fasp;
        dvector b_fasp;
        adaptb_fasp.row = adaptb.size();
        adaptb_fasp.val = (double*)adaptb.data();
        
        // initialize solution
        dvector adaptsoluvec;
        fasp_dvec_alloc(adaptb_fasp.row, &adaptsoluvec);
        fasp_dvec_set(adaptb_fasp.row, &adaptsoluvec, 0.0);
        
        
         // Need solver for generic 3x3 Mass matrix
         //#if FASP_BSR
         // convert CSR to BSR
         dBSRmat adaptA_fasp_bsr = fasp_format_dcsr_dbsr(&adaptA_fasp, 3);
         
         // free CSR matrix
         //fasp_dcsr_free(&A_fasp);

         
         printf(" Step 2: initialize solver parameters\n"); fflush(stdout);
         // initialize solver parameters
         input_param     inpar;  // parameters from input files
         itsolver_param  itpar;  // parameters for itsolver
         AMG_param       amgpar; // parameters for AMG
         ILU_param       ilupar; // parameters for ILU
         //#endif
         
         // read in parameters from a input file
         //#if FASP_BSR
         char inputfile[] = "./params/bsr.dat";
         fasp_param_input(inputfile, &inpar);
         fasp_param_init(&inpar, &itpar, &amgpar, &ilupar, NULL);
        
        
        printf(" Step 3: solve the linear system\n"); fflush(stdout);
        // solve
        int status=FASP_SUCCESS;
        //fasp_param_amg_print(&amgpar);
        
        //#if FASP_BSR
        status = fasp_solver_dbsr_krylov_diag(&adaptA_fasp_bsr, &adaptb_fasp, &adaptsoluvec, &itpar);
        //status = fasp_solver_dbsr_krylov_amg(&adaptA_fasp_bsr, &adaptb_fasp, &adaptsoluvec, &itpar, &amgpar);
        //#endif
        
        
        
        if (status<0) {
            printf("\n### WARNING: Solver failed! Exit status = %d.\n\n", status); fflush(stdout);
        }
        else {
            printf("\nSolver finished successfully!\n\n"); fflush(stdout);
        }
        
        
        
        printf(" Step 4: convert solution back\n"); fflush(stdout);
        // Convert solution vector to FE solution
        double * adaptSolVal = adaptSolu.vector()->data();
        //#if FASP_BSR
        for(std::size_t i=0; i<adaptsoluvec.row; ++i) {
            adaptSolVal[i] = adaptsoluvec.val[i];
        }
        //#endif
        Dsoln = adaptSolu;
        
        
        
        
        printf(" Step 5: free memory\n"); fflush(stdout);
        // Free memory
        //#if FASP_BSR
        //fasp_dbsr_free(&adaptA_fasp_bsr);
        //#endif
        //fasp_dvec_free(&b_fasp);
        fasp_dvec_free(&adaptsoluvec);
        free(adaptap); free(adaptai);
        
        printf(" -------------------------------------\n"); fflush(stdout);
        printf(" End of interface to FASP \n"); fflush(stdout);
        printf(" -------------------------------------\n\n"); fflush(stdout);
        

        
        // Estimate error
        printf(" Estimating error... "); fflush(stdout);
          PoissonCellMarker::FunctionSpace DG(meshAdapt);
          PoissonCellMarker::LinearForm errForm(DG);
          errForm.Du = Dsoln;
          errForm.u  = adaptFN;
        
        // Mark elements for refinement
        printf(" Marking elements for refinement \n"); fflush(stdout);
          Vector errVec;
          assemble(errVec,errForm);
          uint refineCount=0;
          MeshFunction<bool> cellMark(meshAdapt,3,false);
          for ( uint errVecInd=0; errVecInd<errVec.size(); errVecInd++) {
            if ( errVec[errVecInd] > adaptTol ) {
                refineCount++;
                cellMark.values()[errVecInd] = true;
            }
          }
        
        // Refine marked elemetns
        if ( refineCount>0 ) {
            printf(" Refine mesh\n "); fflush(stdout);
            mesh = adapt(meshAdapt,cellMark);
            numRefines++;
            unrefined  = false;
        }
        
        // No elements marked for refinement
        else {
            printf(" No elements marked for refinement... adaptivity complete\n");
            fflush(stdout);
            tooCoarse = false;
        }
        
    }
        
    // Solve complete
    if ( unrefined && converged ) {
        printf(" Converged to a solution \n\n"); fflush(stdout);
        printf(" Review output files \n\n"); fflush(stdout);
        return 0;
    }
    
    // Count iterations
    totalRefines++;
    printf(" Adaptivity iteration %d\n\n\n",totalRefines); fflush(stdout);
    
    
    
    
        
    
    //***********************************
    //   Finite Element Space and Forms
    //***********************************
        
    printf("\nDiscretize the PNP system \n"); fflush(stdout);
    
    // Finite element space
    printf(" Define PNP finite elements \n"); fflush(stdout);
    two_species::FunctionSpace W(mesh);
    
    // Define variational forms
    printf(" Define PNP variational forms \n\n"); fflush(stdout);
    two_species::BilinearForm a(W, W);
    two_species::LinearForm L(W);
    
    
        
        
    //*********************************************
    //  Mark subdomains and impose Dirichlet B.C.
    //*********************************************
        
    // Define channel walls
    printf("Define subdomains \n"); fflush(stdout);
    dielectricChannel sChannel;
    channelGate gChannel;
    FacetFunction<std::size_t> channelBndry(mesh);
    channelBndry.set_all(1);
    //gChannel.mark(channelBndry,3);
    //sChannel.mark(channelBndry,2);
    a.ds = channelBndry; L.ds = channelBndry;
    
    
    // Define Dirichlet boundary
    printf(" Define Dirichlet boundary conditions \n\n"); fflush(stdout);
    Constant DirBC(0.0, 0.0, 0.0);
    DirichletBoundary boundary;
    DirichletBC bc(W, DirBC, boundary);
    
    // Write marked mesh to file
    meshfile   << channelBndry;
        
        
        
        
        
    //***********************************
    // Interpolate Analytic Expressions
    //***********************************
        
    
    
    //***********************
    // Assign coefficients
    //***********************
    
    printf("Assign coefficients for variational forms \n\n"); fflush(stdout);
    
    // p-Nernst-Planck eqn
    Constant Dp(DpVal);  a.Dp = Dp; L.Dp = Dp;  // Diffusivity
    Constant qp(1.0);    a.qp = qp; L.qp = qp;  // Valency
    Constant mp(DpVal);  a.mp = mp; L.mp = mp;  // Mobility
    //cationSource pf;   L.pf = pf;               // Source terms
    Constant newf(0.0);  L.pf = newf;
    
    // n-Nernst-Planck eqn
    Constant Dn(DnVal);  a.Dn = Dn; L.Dn = Dn;  // Diffusivity
    Constant qn(-1.0);   a.qn = qn; L.qn = qn;  // Valency
    Constant mn(DnVal);  a.mn = mn; L.mn = mn;  // Mobility
    //anionSource  nf;   L.nf = nf;               // Source terms
    L.nf = newf;
    
    // Poisson eqn
    Constant eps(epsVal);
    a.eps = eps; L.eps = eps;                   // Permittivity
    Constant schrg( 1.0);//(rs_chrg);
    L.s1 = schrg;                               // Channel surface charge
    Constant gchrg(-1.0);//(rg_chrg);
    L.s2 = gchrg;                               // Gate surface charge
    L.fc = newf;                                // Fixed charge
    
    
    
    
        
    //***********************************
    //  Interpolate solution iterate
    //***********************************
        
    printf("Define initial guess for Newton iteration \n\n"); fflush(stdout);
    Function iterate(W);
    Function concsoln(W);
    
    Function pIterate(iterate[0]);
    pIterate.interpolate(initP);
    
    Function nIterate(iterate[1]);
    nIterate.interpolate(initN);
    
    // Convert from log to density
    Function pconc(concsoln[0]);
    Function nconc(concsoln[1]);
    double * pconcval = pconc.vector()->data();
    double * pItval = pIterate.vector()->data();
    for(std::size_t i=0; i<pIterate.vector()->size(); ++i) {
        pconcval[i] = exp(pItval[i]);
    }
    double * nconcval = nconc.vector()->data();
    double * nItval = nIterate.vector()->data();
    for(std::size_t i=0; i<nIterate.vector()->size(); ++i) {
        nconcval[i] = exp(nItval[i]);
    }
    
        
        
    //************************************
    //  Solve for electrostatic potential
    //************************************
    printf("Solve for electric potential \n\n"); fflush(stdout);
    Function kIterate(iterate[2]);
    kIterate.interpolate(initK);
    
    
    
    //  Dissipation
    /*
    printf("Define dissipation-rate functional \n"); fflush(stdout);
    dissipation::Functional Diss(mesh);
    Diss.Dp  = Dp;
    Diss.Dn  = Dn;
    Diss.p   = pconc;
    Diss.n   = nconc;
    Diss.phi = kIterate;
    double diss;
    diss     = assemble(Diss);
    printf(" The initial dissipation-rate is %e \n\n", diss); fflush(stdout);
    */
    
    
    
    //***********************************
    // Define exact solutions and error
    //***********************************
    
    printf("Define exact solution and measure error \n"); fflush(stdout);
    dirichlet_exact_p   exactp;
    dirichlet_exact_n   exactn;
    dirichlet_exact_phi exactphi;
    Function exactSoln(W);
    Function exact_p(exactSoln[0]);   exact_p.interpolate(exactp);
    Function exact_n(exactSoln[1]);   exact_n.interpolate(exactn);
    Function exact_phi(exactSoln[2]); exact_phi.interpolate(exactphi);
    
    Function error_p(exact_p);     *(error_p.vector())   -= *(pconc.vector());
    Function error_n(exact_n);     *(error_n.vector())   -= *(nconc.vector());
    Function error_phi(exact_phi); *(error_phi.vector()) -= *(kIterate.vector());
    
    // H1-norm
    energynorm::Functional M(mesh); // H1-norm
    M.p   = error_p;
    M.n   = error_n;
    M.phi = error_phi;
    double h1norm_error   = pow(assemble(M),.5);
    double prev_h1;
    printf(" The initial H1-error is %e \n\n",h1norm_error);
    
    // Save initial error in VTK format
    /*File pInitErrorFile("pInitError.pvd");
    pInitErrorFile   << error_p;
    File nInitErrorFile("nInitError.pvd");
    nInitErrorFile   << error_n;
    File phiInitErrorFile("phiInitError.pvd");
    phiInitErrorFile << error_phi;
     */
    
    
    
    
    
    ////////////////////////
    //                    //
    //      Solver        //
    //                    //
    ////////////////////////
        
    printf("\n\n ----------------------------------------\n"); fflush(stdout);
    printf(" Initializing Newton solver \n"); fflush(stdout);
    printf(" ----------------------------------------\n\n"); fflush(stdout);
    
    
    //********************
    //   Solver objects
    //********************
    
    printf("\nInitialize solver objects \n"); fflush(stdout);
    
    // Counters
    int it=0, i;
    double prev_relR, normp, normn, normphi, normb_fasp;
    double relR_fasp = 1.0;
    bool   done = false;
    
    // Newton updates
    Function soln(W);
    Function p(soln[0]);
    Function n(soln[1]);
    Function phi(soln[2]);
    
    // Block GS updates
    Function solu(W);
    
    // Linear system
    Matrix A;
    Vector b;
    
    // Extract DoF indices
    printf(" Assemble solution indices \n"); fflush(stdout);
    std::vector<std::size_t> component(1);
    std::vector<dolfin::la_index> gidx_p;
    std::vector<dolfin::la_index> gidx_n;
    std::vector<dolfin::la_index> gidx_phi;
    const dolfin::la_index n0 = W.dofmap()->ownership_range().first;
    const dolfin::la_index n1 = W.dofmap()->ownership_range().second;
    const dolfin::la_index num_dofs = n1 - n0;
    component[0] = 0;
    boost::shared_ptr<GenericDofMap> dofmap_p   = W.dofmap()->extract_sub_dofmap(component,mesh);
    component[0] = 1;
    boost::shared_ptr<GenericDofMap> dofmap_n   = W.dofmap()->extract_sub_dofmap(component,mesh);
    component[0] = 2;
    boost::shared_ptr<GenericDofMap> dofmap_phi = W.dofmap()->extract_sub_dofmap(component,mesh);

    for ( CellIterator cell(mesh); !cell.end(); ++cell)
    {
        const std::vector<dolfin::la_index> cell_dofs_p   = dofmap_p->cell_dofs(cell->index());
        const std::vector<dolfin::la_index> cell_dofs_n   = dofmap_n->cell_dofs(cell->index());
        const std::vector<dolfin::la_index> cell_dofs_phi = dofmap_phi->cell_dofs(cell->index());
        for (std::size_t i = 0; i < cell_dofs_p.size(); ++i)
        {
            const std::size_t dof = cell_dofs_p[i];
            if (dof >= n0 && dof < n1)
                gidx_p.push_back(dof);
        }
        for (std::size_t i = 0; i < cell_dofs_n.size(); ++i)
        {
            const std::size_t dof = cell_dofs_n[i];
            if (dof >= n0 && dof < n1)
                gidx_n.push_back(dof);
        }
        for (std::size_t i = 0; i < cell_dofs_phi.size(); ++i)
        {
            const std::size_t dof = cell_dofs_phi[i];
            if (dof >= n0 && dof < n1)
                gidx_phi.push_back(dof);
        }
    }
    std::sort(gidx_p.begin(), gidx_p.end());
    std::sort(gidx_n.begin(), gidx_n.end());
    std::sort(gidx_phi.begin(), gidx_phi.end());
    // Remove duplicates
    gidx_p.erase(std::unique(gidx_p.begin(), gidx_p.end()), gidx_p.end());
    gidx_n.erase(std::unique(gidx_n.begin(), gidx_n.end()), gidx_n.end());
    gidx_phi.erase(std::unique(gidx_phi.begin(), gidx_phi.end()), gidx_phi.end());
    // write index to file
    /*
     
     printf(" Write index vectors to files \n"); fflush(stdout);
     FILE* pidxfile;
     pidxfile = fopen("pidx.dat", "w");
     fprintf(pidxfile, "%d \n",gidx_p.size());
     for(std::size_t i=0; i<gidx_p.size(); i++)
     fprintf(pidxfile, "%d \n",gidx_p[i]);
     fclose(pidxfile);
     
     FILE* nidxfile;
     nidxfile = fopen("nidx.dat", "w");
     fprintf(nidxfile, "%d \n",gidx_n.size());
     for(std::size_t i=0; i<gidx_n.size(); i++)
     fprintf(nidxfile, "%d \n",gidx_n.data()[i]);
     fclose(nidxfile);
     
     FILE* phiidxfile;
     phiidxfile = fopen("phiidx.dat", "w");
     fprintf(phiidxfile, "%d \n",gidx_phi.size());
     for(std::size_t i=0; i<gidx_phi.size(); i++)
     fprintf(phiidxfile, "%d \n",gidx_phi.data()[i]);
     fclose(phiidxfile);
     
     
     //  Write RHS to file
     printf(" Write RHS to file  \n"); fflush(stdout);
     dvector b_fasp;
     b_fasp.row = b.size();
     b_fasp.val = (double*)b.data();
     
     FILE* bfile;
     bfile = fopen("rhs.dat", "w");
     fprintf(bfile, "%d \n",b_fasp.row);
     for( i=0; i<b_fasp.row; i++)
     fprintf(bfile, "%f \n",b_fasp.val[i]);
     fclose(bfile);
     */

    
    // Estimate the initial residual
    printf(" Measure initial residual \n");
    L.pp = pIterate;
    L.nn = nIterate;
    L.kk = kIterate;
    assemble(b,L); bc.apply(b);
    
    dvector b_fasp;
    b_fasp.row = b.size();
    b_fasp.val = (double*)b.data();
    if (b0_fasp < 0.) {   // Initial residual
        b0_fasp = fasp_blas_array_norm2(b_fasp.row,b_fasp.val);
        printf(" The initial residual is %e \n\n", b0_fasp);
    }
    
    
    // Start Timer
	timeval tim ;
	gettimeofday(&tim, NULL) ;
	double runtime = tim.tv_sec+(tim.tv_usec/1000000.0) ;
    
    
    //********************
    //  Newton iteration
    //********************
    
    while (done==false) {
        
        // Save solution in VTK format
        printf(" Write solution to file \n\n"); fflush(stdout);
        pfile   << pconc;
        nfile   << nconc;
        phifile << kIterate;
        
        // Save final error in VTK format
        pErrorFile   << error_p;
        nErrorFile   << error_n;
        phiErrorFile << error_phi;
        
        // Update newton step
        it++;
        printf("Newton iteration %d\n", it);
        fflush(stdout);

        // Update PNP coefficients
        printf(" Construct Jacobian matrix \n"); fflush(stdout);
        a.pp = pIterate; L.pp = pIterate;
        a.nn = nIterate; L.nn = nIterate;
        a.kk = kIterate; L.kk = kIterate;
        assemble(A,a); bc.apply(A);
        
        
        
        //**************************************************
        //  Solve linear problem: interface with FASP
        //**************************************************
        
        printf(" -------------------------------------\n"); fflush(stdout);
        printf(" Start interface to FASP  \n"); fflush(stdout);
        printf(" -------------------------------------\n\n"); fflush(stdout);
        
        printf(" Step 1: convert sparse matrix format\n"); fflush(stdout);
        // Convert A to CSR
        dCSRmat A_fasp;
        
        unsigned int nz = boost::tuples::get<3>(A.data());
        int row = A.size(0);
        int col = A.size(1);
        int* ap = (int*)fasp_mem_calloc(row+1, sizeof(int));
        const size_t* ap_tmp = boost::tuples::get<0>(A.data());
        for (i=0; i<row+1; i++) {
            ap[i] = (int)ap_tmp[i];
        }
        int* ai = (int*)fasp_mem_calloc(nz, sizeof(int));
        const size_t* ai_tmp = boost::tuples::get<1>(A.data());
        for (i=0; i<nz; i++) {
            ai[i] = (int)ai_tmp[i];
        }
        double* ax = (double*)boost::tuples::get<2>(A.data());
        
        A_fasp.row = row;
        A_fasp.col = col;
        A_fasp.nnz = nz;
        A_fasp.IA  = ap;
        A_fasp.JA  = ai;
        A_fasp.val = ax;
        
#if FASP_BSR
        // convert CSR to BSR
        dBSRmat A_fasp_bsr = fasp_format_dcsr_dbsr(&A_fasp, 3);
#else
        // convert CSR to block CSR
        block_dCSRmat Abcsr;
        dCSRmat *A_diag;
        
        // get index
        INT nrow = A_fasp.row/3;
        ivector phi_idx;
        ivector n_idx;
        ivector p_idx;
        
        fasp_ivec_alloc(nrow, &phi_idx);
        fasp_ivec_alloc(nrow, &n_idx);
        fasp_ivec_alloc(nrow, &p_idx);
        
        for (i=0; i<nrow; i++){
            phi_idx.val[i] = 3*i;
            n_idx.val[i] = 3*i+1;
            p_idx.val[i] = 3*i+2;
        }
        
        // Assemble the matrix in block dCSR format
        Abcsr.brow = 3; Abcsr.bcol = 3;
        Abcsr.blocks = (dCSRmat **)calloc(9, sizeof(dCSRmat *));
        for (i=0; i<9 ;i++) {
            Abcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
        
        // A11
        fasp_dcsr_getblk(&A_fasp, phi_idx.val, phi_idx.val, nrow, nrow, Abcsr.blocks[0]);
        // A12
        fasp_dcsr_getblk(&A_fasp, phi_idx.val, n_idx.val, nrow, nrow, Abcsr.blocks[1]);
        // A13
        fasp_dcsr_getblk(&A_fasp, phi_idx.val, p_idx.val, nrow, nrow, Abcsr.blocks[2]);
        // A21
        fasp_dcsr_getblk(&A_fasp, n_idx.val, phi_idx.val, nrow, nrow, Abcsr.blocks[3]);
        // A22
        fasp_dcsr_getblk(&A_fasp, n_idx.val, n_idx.val, nrow, nrow, Abcsr.blocks[4]);
        // A23
        fasp_dcsr_getblk(&A_fasp, n_idx.val, p_idx.val, nrow, nrow, Abcsr.blocks[5]);
        // A31
        fasp_dcsr_getblk(&A_fasp, p_idx.val, phi_idx.val, nrow, nrow, Abcsr.blocks[6]);
        // A32
        fasp_dcsr_getblk(&A_fasp, p_idx.val, n_idx.val, nrow, nrow, Abcsr.blocks[7]);
        // A33
        fasp_dcsr_getblk(&A_fasp, p_idx.val, p_idx.val, nrow, nrow, Abcsr.blocks[8]);
        
 
        // setup diagonal blocks for the preconditioner
        A_diag = (dCSRmat *)fasp_mem_calloc(3, sizeof(dCSRmat));
        
        // first diagonal block
        A_diag[0].row = Abcsr.blocks[0]->row;
        A_diag[0].col = Abcsr.blocks[0]->col;
        A_diag[0].nnz = Abcsr.blocks[0]->nnz;
        A_diag[0].IA  = Abcsr.blocks[0]->IA;
        A_diag[0].JA  = Abcsr.blocks[0]->JA;
        A_diag[0].val = Abcsr.blocks[0]->val;
        
        // second diagonal block
        A_diag[1].row = Abcsr.blocks[4]->row;
        A_diag[1].col = Abcsr.blocks[4]->col;
        A_diag[1].nnz = Abcsr.blocks[4]->nnz;
        A_diag[1].IA  = Abcsr.blocks[4]->IA;
        A_diag[1].JA  = Abcsr.blocks[4]->JA;
        A_diag[1].val = Abcsr.blocks[4]->val;
        
        // third diagonal block
        A_diag[2].row = Abcsr.blocks[8]->row;
        A_diag[2].col = Abcsr.blocks[8]->col;
        A_diag[2].nnz = Abcsr.blocks[8]->nnz;
        A_diag[2].IA  = Abcsr.blocks[8]->IA;
        A_diag[2].JA  = Abcsr.blocks[8]->JA;
        A_diag[2].val = Abcsr.blocks[8]->val;
        
        
        //fasp_dcoo_write("aDiag0",&A_diag[0]);
        //fasp_dcoo_write("aDiag1",&A_diag[1]);
        //fasp_dcoo_write("aDiag2",&A_diag[2]);
        
        
        
        // convert right hand side
        dvector bbcsr;
        fasp_dvec_alloc(b_fasp.row, &bbcsr);
        for (i=0; i<nrow; i++){
            
            bbcsr.val[i]        = b_fasp.val[3*i];
            bbcsr.val[nrow+i]    = b_fasp.val[3*i+1];
            bbcsr.val[2*nrow+i]  = b_fasp.val[3*i+2];
            
        }
        
#endif
        
        
        // free CSR matrix
        //fasp_dcsr_free(&A_fasp);
        
        // initialize solution
        dvector soluvec;
        fasp_dvec_alloc(b_fasp.row, &soluvec);
        fasp_dvec_set(b_fasp.row, &soluvec, 0.0);
        
        printf(" Step 2: initialize solver parameters\n"); fflush(stdout);
        // initialize solver parameters
        input_param     inpar;  // parameters from input files
        itsolver_param  itpar;  // parameters for itsolver
        AMG_param       amgpar; // parameters for AMG
        ILU_param       ilupar; // parameters for ILU
        
        // read in parameters from a input file
#if FASP_BSR
        char inputfile[] = "./params/bsr.dat";
#else
        char inputfile[] = "./params/bcsr.dat";
#endif
        fasp_param_input(inputfile, &inpar);
        fasp_param_init(&inpar, &itpar, &amgpar, &ilupar, NULL);
        
        printf(" Step 3: solve the linear system\n"); fflush(stdout);
        // solve
        int status=FASP_SUCCESS;
        //fasp_param_amg_print(&amgpar);
        
#if FASP_BSR
        status = fasp_solver_dbsr_krylov_amg(&A_fasp_bsr, &b_fasp, &soluvec, &itpar, &amgpar);
#else
        status = fasp_solver_bdcsr_krylov_block_3(&Abcsr, &bbcsr, &soluvec, &itpar, &amgpar, A_diag);
#endif
        
        if (status<0) {
            printf("\n### WARNING: Solver failed! Exit status = %d.\n\n", status); fflush(stdout);
        }
        else {
            printf("\nSolver finished successfully!\n\n"); fflush(stdout);
        }
        
        printf(" Step 4: convert solution back\n"); fflush(stdout);
        
        // Convert solution vector to FE solution
        double * solval = solu.vector()->data();
#if FASP_BSR
        for(std::size_t i=0; i<soluvec.row; ++i) {
            solval[i] = soluvec.val[i];
        }
#else
        for(std::size_t i=0; i<nrow; ++i) {
            solval[3*i]    = soluvec.val[i];
            solval[3*i+1]  = soluvec.val[nrow+i];
            solval[3*i+2]  = soluvec.val[2*nrow+i];
        }
#endif
        
        printf(" Step 5: free memory\n"); fflush(stdout);
        // Free memory
#if FASP_BSR
        fasp_dbsr_free(&A_fasp_bsr);
#else
        fasp_bdcsr_free(&Abcsr);
        //for (std::size_t i=0; i<9; i++) fasp_dcsr_free(Abcsr.blocks[i]);
        //fasp_mem_free(&Abcsr);
        fasp_dvec_free(&bbcsr);
        fasp_ivec_free(&phi_idx);
        fasp_ivec_free(&n_idx);
        fasp_ivec_free(&p_idx);
#endif
        //fasp_dvec_free(&b_fasp);
        fasp_dvec_free(&soluvec);
        free(ap); free(ai);
        
        printf(" -------------------------------------\n"); fflush(stdout);
        printf(" End of interface to FASP \n"); fflush(stdout);
        printf(" -------------------------------------\n\n"); fflush(stdout);
        
        
        
        
        
        
        //***************************
        //      Update solution
        //***************************
        
        printf(" Update solution \n"); fflush(stdout);
        p = solu[0]; n = solu[1]; phi = solu[2];

        
        
        //  Backtrack line search
        Function update(W);
        Function updateP(update[0]);
        Function updateN(update[1]);
        Function updateK(update[2]);
        double testRelRes = 1.0+DOLFIN_EPS;
        double dampFactor = 2.0*mu;
        while ( testRelRes > 1.0-DOLFIN_EPS && dampFactor > 1.0e-2 ) {
            dampFactor  /= 2.0;
            updateP.interpolate(pIterate);
            updateN.interpolate(nIterate);
            updateK.interpolate(kIterate);
            *(p.vector())   *= dampFactor;
            *(n.vector())   *= dampFactor;
            *(phi.vector()) *= dampFactor;
            *(updateP.vector()) += *(p.vector());
            *(updateN.vector()) += *(n.vector());
            *(updateK.vector()) += *(phi.vector());
            
            L.pp = updateP;
            L.nn = updateN;
            L.kk = updateK;
            assemble(b,L); bc.apply(b);
            
            dvector newb_faspBT;
            newb_faspBT.row = b.size();
            newb_faspBT.val = (double*)b.data();
            double normb_faspBT = fasp_blas_array_norm2(newb_faspBT.row,newb_faspBT.val);
            testRelRes   = normb_faspBT/b0_fasp ;
            
            if ( testRelRes > 1.0-DOLFIN_EPS ) {
                printf(" Reducing the update by a factor of %e \n", dampFactor);
                fflush(stdout);
            }
        }
        if ( testRelRes > 1.0-DOLFIN_EPS ){
            printf(" BAD UPDATE!!! \n"); fflush(stdout);
        }
        
        // update solutions
        *(pIterate.vector()) += *(p.vector());
        *(nIterate.vector()) += *(n.vector());
        *(kIterate.vector()) += *(phi.vector());

        
        
        
        // Update residual
        printf(" Update variational forms \n"); fflush(stdout);
        L.pp = pIterate;
        L.nn = nIterate;
        L.kk = kIterate;
        assemble(b,L); bc.apply(b);
        
        // convert right hand side and measure relative residual
        printf(" Update the nonlinear residual \n");
        dvector newb_fasp;
        newb_fasp.row = b.size();
        newb_fasp.val = (double*)b.data();
        normb_fasp = fasp_blas_array_norm2(newb_fasp.row,newb_fasp.val);
        prev_relR = relR_fasp;
        relR_fasp = normb_fasp/b0_fasp ;
        printf(" After %d iteration(s) the relative residual is %e \n", it, relR_fasp);
        fflush(stdout);
        
        
        
        //************************
        //  Measure dissipation
        //************************
        
        // Convert from log to density
        double * pconcval = pconc.vector()->data();
        double * pItval = pIterate.vector()->data();
        for(std::size_t i=0; i<pIterate.vector()->size(); ++i) {
            pconcval[i] = exp(pItval[i]);
        }
        double * nconcval = nconc.vector()->data();
        double * nItval = nIterate.vector()->data();
        for(std::size_t i=0; i<nIterate.vector()->size(); ++i) {
            nconcval[i] = exp(nItval[i]);
        }
        
        /*
        Diss.p   = pconc;
        Diss.n   = nconc;
        Diss.phi = kIterate;
        diss     = assemble(Diss);
        printf(" The dissipation is %e \n\n", diss);
        fflush(stdout);
        */

        
        //****************
        // Update error
        //****************
        printf(" Update error \n");
        *(error_p.vector())   *= 0.0;
        *(error_n.vector())   *= 0.0;
        *(error_phi.vector()) *= 0.0;
        *(error_p.vector())   += *(exact_p.vector());
        *(error_n.vector())   += *(exact_n.vector());
        *(error_phi.vector()) += *(exact_phi.vector());
        *(error_p.vector())   -= *(pconc.vector());
        *(error_n.vector())   -= *(nconc.vector());
        *(error_phi.vector()) -= *(kIterate.vector());
        M.p   = error_p; M.n   = error_n; M.phi = error_phi;
        
        // compute and report the H1-error
        prev_h1 = h1norm_error;
        h1norm_error   = pow(assemble(M),.5);
        printf(" After %d iteration(s), the H1-error is %e ==> ", it, h1norm_error);
        if ( prev_h1>h1norm_error ) printf("Error Reduction\n\n");
        else                        printf("Error increase!!\n\n");
         
        
        
        
        
        
        //**************************
        //  Check for convergence
        //**************************
        
        if ( relR_fasp < tol ) {       // Converged
            printf("\nConvergence: The relative residual is below the desired tolerance. \n\n");
            fflush(stdout);
            done = true;
            converged = true;
            
            //	Stop Timer
            gettimeofday(&tim, NULL) ;
            runtime -= tim.tv_sec+(tim.tv_usec/1000000.0) ;
            runtime *= -1. ;
            printf(" Runtime:  %e\n\n", runtime);
            
            
            
            // Save solution in VTK format
            printf(" Write solution to file \n"); fflush(stdout);
            pfile   << pconc;
            nfile   << nconc;
            phifile << kIterate;
            
            // Save final error in VTK format
            printf(" Write error to file \n"); fflush(stdout);
            pErrorFile   << error_p;
            nErrorFile   << error_n;
            phiErrorFile << error_phi;
            
            
            //  Write simulation summary to file
            /*
             printf(" Output for main-log-steady.cpp \n");fflush(stdout);
             FILE *fout;
             strcat(paramRegime,"-log-main-steady");
             fout = fopen(paramRegime, "w");
             fprintf(fout,"Output for main-log-steady.cpp \n\n");
             fprintf(fout," Mesh: ( %d, %d, %d ) \n", Nx,Ny,Nz);
             fprintf(fout," eps = %e \n", epsVal);
             fprintf(fout," Dp  = %e \n", DpVal);
             fprintf(fout," Dn  = %e \n\n", DnVal);
             if (it>maxit)
             fprintf(fout," Status: No convergence \n");
             else
             fprintf(fout," Status: Convergence \n");
             fprintf(fout," H1-error:    %e \n", h1norm_error);
             fprintf(fout," Nonlin Res:  %e \n", normb_fasp);
             fprintf(fout," Rel Res:     %e \n", relR_fasp);
             fprintf(fout," Dissipation: %e \n", diss);
             fprintf(fout," Newton Iter: %d \n", it);
             fprintf(fout," Runtime:     %e\n\n", runtime);
             fclose(fout);
             */
            
            /*
            printf(" Review output files \n\n"); fflush(stdout);
            
            return 0;*/
            
            
        } else if ( it>maxit ) {                    // Too many iterations
            printf("\nDivergence: The relative residual has not converged after %d iterations. \n\n",maxit+1);
            fflush(stdout);
            done=true;
        }  else if ( relR_fasp >= prev_relR ) {     // Increasing residual
            printf("The relative residual has grown in this iteration... \n\n");fflush(stdout);
        }

    
    }
    
    
    printf(" -------------------------------------\n"); fflush(stdout);
    printf(" End of Newton Iteration \n"); fflush(stdout);
    printf(" -------------------------------------\n\n"); fflush(stdout);
    
        
        
        
        
    
    //********************
    //  Adaptivity pass
    //********************
    
    printf(" Update mesh \n\n"); fflush(stdout);
    mesh0   = mesh;
    boost::shared_ptr<const Mesh> adapted_mesh( new const Mesh(mesh0) );
    adaptFN = adapt(kIterate,adapted_mesh);
    initP   = adapt(pIterate,adapted_mesh);
    initN   = adapt(nIterate,adapted_mesh);
    initK   = adapt(kIterate,adapted_mesh);
    
    }
    
    
    printf(" Failed to converge below desired relative residual!!! \n\n"); fflush(stdout);
    
    printf(" Exiting... \n\n"); fflush(stdout);
    
    
   
    
    return 0;
}
