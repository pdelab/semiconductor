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
#include "./include/protein_four_species.h"
#include "./include/protein_poisson.h"
#include "./include/current_flux.h"
#include "./include/poisson_cell_marker.h"
#include "./include/gradient_recovery.h"
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

double Lx;
double Ly;
double Lz;

using namespace std;
using namespace dolfin;


//////////////////////////////////
//                              //
//      Function Definitions    //
//       of Initial Guesses     //
//                              //
//////////////////////////////////


//  Initial Sodium Number Density Profile
class LogIon : public Expression
{
public:
    LogIon(double ext_bulk, double int_bulk, double bc_dist, int bc_dir): Expression(),ext_ion_bulk(ext_bulk),int_ion_bulk(int_bulk),bc_distance(bc_dist),bc_direction(bc_dir) {}
    void eval(Array<double>& values, const Array<double>& x) const
    {
        values[0]  = log(ext_ion_bulk)*(x[bc_direction]+bc_distance/2.0)/(bc_distance);
        values[0] -= log(int_ion_bulk)*(x[bc_direction]-bc_distance/2.0)/(bc_distance);
    }
private:
    double ext_ion_bulk, int_ion_bulk, bc_distance;
    int bc_direction;
};


//  Voltage
class Voltage : public Expression
{
public:
    Voltage(double ext_volt, double int_volt, double bc_dist, int bc_dir): Expression(),ext_voltage(ext_volt),int_voltage(int_volt),bc_distance(bc_dist),bc_direction(bc_dir) {}
    void eval(Array<double>& values, const Array<double>& x) const
    {
        values[0]  = ext_voltage*(x[bc_direction]+bc_distance/2.0)/(bc_distance);
        values[0] -= int_voltage*(x[bc_direction]-bc_distance/2.0)/(bc_distance);
    }
private:
    double ext_voltage, int_voltage, bc_distance;
    int bc_direction;
};


//  Analytic Solution to Fixed Point-Charges
class CoulombPotentials : public Expression
{
public: CoulombPotentials(double perm,uint n_atoms,double* pt_charge_x,double* pt_charge_y,double* pt_charge_z,double* pt_charge_q): Expression(),eps(perm),N(n_atoms),x_coord(pt_charge_x),y_coord(pt_charge_y),z_coord(pt_charge_z),valence(pt_charge_q) {}
    void eval(Array<double>& values, const Array<double>& x) const
    {
        
        double PI= 3.1415961;
        double q;
        Array<double> y(3);
        values[0] = 0.0;
        
        for (uint atom_ind=0; atom_ind<N; atom_ind++) {
            y[0] = x_coord[atom_ind];
            y[1] = y_coord[atom_ind];
            y[2] = z_coord[atom_ind];
            q    = valence[atom_ind];
            
            values[0] += (q/(4.0*PI*eps))*( 1.0/sqrt(1.0e-8+(x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2])) );
        }
    }
private:
    double eps;
    uint N;
    double* x_coord;
    double* y_coord;
    double* z_coord;
    double* valence;
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

//  Dirichlet boundary condition
class DirBCval : public Expression
{
public:
    
    DirBCval() : Expression(5) {}
    
    void eval(Array<double>& values, const Array<double>& x) const
    {
        
        values[0] = 0.0;
        values[1] = 0.0;
        values[2] = 0.0;
        values[3] = 0.0;
        values[4] = 0.0;
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
    parameters["allow_extrapolation"] = true;
    
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
    char buffer[500];      // max number of char for each line
    int  val;
    ifstream expin;
    char outputfile[128];       // output directory
    char meshIn[128];            // mesh input file
    char surfIn[128];            // mesh surfaces file
    char subdIn[128];            // mesh subdomains file
    char point_charge_file[128]; // point charge position file

    // BoxMesh Parameters
    double Dx, Dy, Dz, T;
    int    Nx, Ny, Nz, Nt;

    // Newton Solver Parameters
    double tol;
    uint   maxit;
    double mu;
    double adaptTol;
    
    // Boundary conditions
    int    bc_direction;
    double bc_distance;
    double temperature;
    double lower_volt, upper_volt, delta_volt;
    double ext_na_bulk, int_na_bulk;
    double ext_k_bulk,  int_k_bulk;
    double ext_cl_bulk, int_cl_bulk;
    double ext_ca_bulk, int_ca_bulk;

    // PDE coefficients
    double solvent_perm, protein_perm;
    double na_diff, na_valency;
    double k_diff, k_valency;
    double cl_diff, cl_valency;
    double ca_diff, ca_valency;
    
    char filenm[] = "./params/iv-protein-exp-names.dat";
    
    // if input file is not specified, use the default values
    /*if (filenm=="") {
        printf("### ERROR: No file specified for input params \n");
        exit(0);
    }*/
    
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
            strncpy(outputfile,sbuff,128);
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
        
        else if (strcmp(buffer,"bc_direction")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { state = false; break; }
            bc_direction = ibuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"bc_distance")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            bc_distance = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"lower_volt")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            lower_volt = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"upper_volt")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            upper_volt = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }

        else if (strcmp(buffer,"delta_volt")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            delta_volt = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"ext_na_bulk")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            ext_na_bulk = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"int_na_bulk")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            int_na_bulk = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"ext_k_bulk")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            ext_k_bulk = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"int_k_bulk")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            int_k_bulk = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"ext_cl_bulk")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            ext_cl_bulk = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"int_cl_bulk")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            int_cl_bulk = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"ext_ca_bulk")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            ext_ca_bulk = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"int_ca_bulk")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            int_ca_bulk = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"temperature")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            temperature = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"point_charge_file")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%s",sbuff);
            if (val!=1) { state = false; break; }
            strncpy(point_charge_file,sbuff,128);
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"solvent_perm")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            solvent_perm = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"protein_perm")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            protein_perm = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"na_diff")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            na_diff = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"na_valency")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            na_valency = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"k_diff")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            k_diff = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"k_valency")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            k_valency = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"cl_diff")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            cl_diff = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"cl_valency")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            cl_valency = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"ca_diff")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            ca_diff = dbuff;
            fgetsPtr = fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"ca_valency")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                state = false; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { state = false; break; }
            ca_valency = dbuff;
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
    
    File meshfile("./output/charged_protein_adapt/mesh.pvd");
    File   Nafile("./output/charged_protein_adapt/NaFinal.pvd");
    File    Kfile("./output/charged_protein_adapt/KFinal.pvd");
    File   Cafile("./output/charged_protein_adapt/CaFinal.pvd");
    File   Clfile("./output/charged_protein_adapt/ClFinal.pvd");
    File  phifile("./output/charged_protein_adapt/phiFinal.pvd");
    File chargePotentialFile("./output/charged_protein_adapt/chargePhiFinal.pvd");

    
    
    //***********
    //  Domain
    //***********
    
    printf("Initialize mesh: "); fflush(stdout);
    Mesh initMesh;
    MeshFunction<std::size_t> surfaces;
    MeshFunction<std::size_t> subdoms;

    if ( strcmp(meshIn,"box")==0 ) {
      printf(" Domain set to [ %6.3f, %6.3f, %6.3f ] \n\n", Dx, Dy, Dz);
      printf(" Initial mesh is %d x %d x %d \n", Nx,Ny,Nz); fflush(stdout);
      BoxMesh bMesh(-Lx,-Ly,-Lz, Lx, Ly, Lz, Nx, Ny, Nz);
      initMesh = bMesh;
    } else {
      printf(" Reading in the mesh from %s \n", meshIn);
        Mesh rMesh(meshIn);
        initMesh = rMesh;
        printf(" Reading in the mesh surfaces from %s \n", surfIn);
        MeshFunction<std::size_t> surfReadIn(initMesh, surfIn);
        printf(" Reading in the mesh subdomains from %s \n", subdIn);
        MeshFunction<std::size_t>  subdReadIn(initMesh, subdIn);
          surfaces = surfReadIn;
          subdoms  = subdReadIn;
          //meshfile << surfaces;
          //meshfile << subdoms;
    }
    
    
    
    
    //************************
    //  Dimensional Analysis
    //************************
    
    // Reference values
    double k_B    = 1.38064880e-23;     // Boltzmann Constant (m^2 kg / s^2 K)
    double eps_0  = 8.85418782e-12 ;    // Vacuum Permittivity (s^4 A^2 / m^3 kg)
    double e_chrg = 1.60217657e-19 ;    // Elementary Positive Charge (A s)
    double n_avo  = 6.02214129e+23 ;    // Avogadro's number 1 / mol
    double V_ref  = 1.0e+3;             // voltage scale (V)
    double L_ref  = 1.0e-09;            // length scale (m)
    double p_ref  = 1.0e+0;             // ionic reference density (mol / m^3)
    double D_ref  = 1.334e-5;           // reference diffusivity (cm^2 / s)
    
    double debye, epsSolventVal, epsMolVal;
    debye = sqrt( (k_B*temperature*eps_0)/(e_chrg*e_chrg*n_avo*p_ref) ) / L_ref; // Debye Length
    
    printf(" The spatial scale is a nanometer \n"); fflush(stdout);
    printf(" The debye length is %e nanometers \n\n", debye); fflush(stdout);
    
    
    // PDE coefficients
    epsSolventVal = solvent_perm * debye;   // dim'less solvent permittivity
    epsMolVal     = protein_perm * debye;   // dim'less protein permittivity
    
    na_diff  = na_diff / D_ref;             // dim'less sodium diffusivity
    k_diff   = k_diff  / D_ref;             // dim'less potassium diffusivity
    cl_diff  = cl_diff / D_ref;             // dim'less chloride diffusivity
    ca_diff  = ca_diff / D_ref;             // dim'less calcium diffusivity
    
    printf(" Dimensionless parameters are given by: \n"); fflush(stdout);
    printf("    Solvent Permittivity:  \t %15.10f \n", epsSolventVal); fflush(stdout);
    printf("    Protein Permittivity:  \t %15.10f \n", epsMolVal); fflush(stdout);
    printf("    Sodium Diffusivity:    \t %15.10f \n", na_diff); fflush(stdout);
    printf("    Potassium Diffusivity: \t %15.10f \n", k_diff); fflush(stdout);
    printf("    Chloride Diffusivity:  \t %15.10f \n", cl_diff); fflush(stdout);
    printf("    Calcium Diffusivity:   \t %15.10f \n\n", ca_diff); fflush(stdout);
    fflush(stdout);
    
    
    // Boundary conditions
    double ext_voltage, int_voltage;
    int_na_bulk = int_na_bulk / p_ref;      // dim'less intracellular sodium bulk density
    ext_na_bulk = ext_na_bulk / p_ref;      // dim'less extracellular sodium bulk density
    int_k_bulk  = int_k_bulk  / p_ref;      // dim'less intracellular potassium bulk density
    ext_k_bulk  = ext_k_bulk  / p_ref;      // dim'less extracellular potassium bulk density
    int_cl_bulk = int_cl_bulk / p_ref;      // dim'less intracellular chloride bulk density
    ext_cl_bulk = ext_cl_bulk / p_ref;      // dim'less extracellular chloride bulk density
    int_ca_bulk = int_ca_bulk / p_ref;      // dim'less intracellular calcium bulk density
    ext_ca_bulk = ext_ca_bulk / p_ref;      // dim'less extracellular calcium bulk density
    
    
    double int_e_neutral;
    int_e_neutral = na_valency*int_na_bulk+k_valency*int_k_bulk+ca_valency*int_ca_bulk;
    if ( cl_valency*int_cl_bulk != -int_e_neutral ) {
        int_e_neutral = -int_e_neutral/cl_valency;
        printf(" Intracellular bulk density is not electronuetral\n ");
        printf("   Changing intracellular Chloride bulk from %e M to %e M\n", int_cl_bulk, int_e_neutral);
        int_cl_bulk = int_e_neutral;
    }
    
    double ext_e_neutral;
    ext_e_neutral = na_valency*ext_na_bulk+k_valency*ext_k_bulk+ca_valency*ext_ca_bulk;
    if ( cl_valency*ext_cl_bulk != -ext_e_neutral ) {
        ext_e_neutral = -ext_e_neutral/cl_valency;
        printf(" Intracellular bulk density is not electronuetral\n ");
        printf("   Changing intracellular Chloride bulk from %e M to %e M\n", ext_cl_bulk, ext_e_neutral);
        ext_cl_bulk = ext_e_neutral;
    }
    
    
    printf(" The ion bulk boundary conditions are: \n"); fflush(stdout);
    printf("    Extracellular Sodium Bulk Density:    %8.4f M\n", ext_na_bulk); fflush(stdout);
    printf("    Intracellular Sodium Bulk Density:    %8.4f M\n", int_na_bulk); fflush(stdout);
    printf("    Extracellular Potassium Bulk Density: %8.4f M\n", ext_k_bulk);  fflush(stdout);
    printf("    Intracellular Potassium Bulk Density: %8.4f M\n", int_k_bulk);  fflush(stdout);
    printf("    Extracellular Chloride Bulk Density:  %8.4f M\n", ext_cl_bulk); fflush(stdout);
    printf("    Intracellular Chloride Bulk Density:  %8.4f M\n", int_cl_bulk); fflush(stdout);
    printf("    Extracellular Calcium Bulk Density:   %8.4f M\n", ext_ca_bulk); fflush(stdout);
    printf("    Intracellular Calcium Bulk Density:   %8.4f M\n", int_ca_bulk); fflush(stdout);
    
    
    
    
    
    //************************
    //  Analytic Expressions
    //************************
    printf("\n Initializing analytic expressions\n"); fflush(stdout);
    
    // Log-ion boundary interpolant
    printf("   Interpolating ion bulk densities\n"); fflush(stdout);
      LogIon Na(ext_na_bulk,int_na_bulk,bc_distance,bc_direction);
      LogIon K(  ext_k_bulk, int_k_bulk,bc_distance,bc_direction);
      LogIon Cl(ext_cl_bulk,int_cl_bulk,bc_distance,bc_direction);
      LogIon Ca(ext_ca_bulk,int_ca_bulk,bc_distance,bc_direction);
    
    // Electric potential from fixed point charges
    uint   n_atoms;
    double *pt_charge_x;
    double *pt_charge_y;
    double *pt_charge_z;
    double *pt_charge_q;
    if ( strcmp(point_charge_file,"none")==0 ) {
        printf("   No point charges\n"); fflush(stdout);
        
    } else {
        printf("   Reading in point charge coordinates and valencies \n");
        char *pt_chrg_filenm = point_charge_file;
        
        if (pt_chrg_filenm==NULL) {
            printf("### ERROR: No file specified for input params \n");
            exit(0);
        }
        
        FILE *pt_chrg_fp = fopen(pt_chrg_filenm,"r");
        if (pt_chrg_fp==NULL) {
            printf("### ERROR: Could not open file %s...\n", pt_chrg_filenm);
            fasp_chkerr(ERROR_OPEN_FILE, "fasp_param_input");
            
        }
        
        bool state    = true;
        uint atom_ind = 0;
        while ( state ) {
            int     ibuff;
            double  dbuff;
            char    sbuff[500];
            char   *fgetsPtr;
            
            val = fscanf(pt_chrg_fp,"%s",buffer);
            if (val==EOF) break;
            if (val!=1){ state = false; break; }
            if (buffer[0]=='[' || buffer[0]=='%' || buffer[0]=='|') {
                fgetsPtr = fgets(buffer,500,pt_chrg_fp); // skip rest of line
                continue;
            }
            
            // match keyword and scan for value
            if (strcmp(buffer,"n_atoms")==0) {
                val = fscanf(pt_chrg_fp,"%s",buffer);
                if (val!=1 || strcmp(buffer,"=")!=0) {
                    state = false; break;
                }
                val = fscanf(pt_chrg_fp,"%d",&ibuff);
                if (val!=1) { state = false; break; }
                n_atoms = ibuff;
                fgetsPtr = fgets(buffer,500,pt_chrg_fp); // skip rest of line
                printf("     Expecting %d coordinates\n", n_atoms);fflush(stdout);
                pt_charge_x = (double*) malloc(n_atoms*sizeof(double));
                pt_charge_y = (double*) malloc(n_atoms*sizeof(double));
                pt_charge_z = (double*) malloc(n_atoms*sizeof(double));
                pt_charge_q = (double*) malloc(n_atoms*sizeof(double));
            }
            
            else if (strcmp(buffer,"charge")==0) {
                
                val = fscanf(pt_chrg_fp,"%lf",&dbuff);
                if (val!=1) { state = false; break; }
                //y[0] = dbuff;
                pt_charge_x[atom_ind] = dbuff;
                //printf("       read in x ");
                
                val = fscanf(pt_chrg_fp,"%lf",&dbuff);
                if (val!=1) { state = false; break; }
                //y[1] = dbuff;
                pt_charge_y[atom_ind] = dbuff;
                //printf(" y ");
                
                val = fscanf(pt_chrg_fp,"%lf",&dbuff);
                if (val!=1) { state = false; break; }
                //y[2] = dbuff;
                pt_charge_z[atom_ind] = dbuff;
                //printf(" z ");
                
                val = fscanf(pt_chrg_fp,"%lf",&dbuff);
                if (val!=1) { state = false; break; }
                //q    = dbuff;
                pt_charge_q[atom_ind] = dbuff;
                //printf(" q\n");
                
                printf("     Charge %3d puts a %5.2f charge at ( %10.3e, %10.3e, %10.3e )\n", atom_ind+1, pt_charge_q[atom_ind], pt_charge_x[atom_ind], pt_charge_y[atom_ind], pt_charge_z[atom_ind]);
                fflush(stdout);
                atom_ind++;                              // increment index
                
                fgetsPtr = fgets(buffer,500,pt_chrg_fp); // skip rest of line
            }
            
            else {
                state = false;
                printf(" Bad read-in \n\n"); fflush(stdout);
            }
            
        }
        fclose(pt_chrg_fp);
        
    }
    printf("\n"); fflush(stdout);
    CoulombPotentials coulomb(epsMolVal,n_atoms,pt_charge_x,pt_charge_y,pt_charge_z,pt_charge_q);
    
    
    
    //**********************
    //  Iterate of Voltages
    //**********************

    // Open outputfile.csv
    FILE *fout;
    fout = fopen(outputfile, "w");
    fprintf(fout,"current,voltage\n");


    for (double curveNumber = 0.0; curveNumber <= (upper_volt-lower_volt)+DOLFIN_EPS; curveNumber += delta_volt) { 
    
    // Electric potential boundary interpolant
    double volts = lower_volt + curveNumber;
    double current;
    printf(" Computing current at %8.4f mV \n", volts); fflush(stdout);
    int_voltage = -volts / 2.0;
    ext_voltage =  volts / 2.0;
    int_voltage = int_voltage / (V_ref*k_B*temperature/e_chrg);      // dim'less intracellular voltage
    ext_voltage = ext_voltage / (V_ref*k_B*temperature/e_chrg);      // dim'less intracellular voltage
    printf("   Interpolating voltage drop\n"); fflush(stdout);
      Voltage volt(ext_voltage,int_voltage,bc_distance,bc_direction);
    
    
    //***************
    //  Time March
    //***************
    
    
    
    
    
    
    
    //**************
    //  Adaptivity
    //**************
    
    // Copy mesh and initialize CG space
    Mesh mesh0(initMesh);
    protein_four_species::FunctionSpace W0(mesh0);
    uint totalRefines = 0;
    
    // Define initial guesses and residual
    printf("\n Define initial guesses\n"); fflush(stdout);
    Function initGuess(W0);
    Function initNa(initGuess[0]);  initNa.interpolate(Na);
    Function initK(initGuess[1]);   initK.interpolate(K);
    Function initCa(initGuess[2]);  initCa.interpolate(Ca);
    Function initCl(initGuess[3]);  initCl.interpolate(Cl);
    Function initPHI(initGuess[4]); initPHI.interpolate(volt);
    
    // Initialize adaptivity marker
    Function adaptFn0(W0);
    Function adaptFN(adaptFn0[2]);
    adaptFN.interpolate(initPHI);
    double b0_fasp = -1.0;
    bool converged = false;
    
    
    // Adapt mesh to current iterate
    //for ( uint adaptInd=0; adaptInd<15; adaptInd++ ) {
    
    bool tooCoarse  = true;
    bool unrefined  = true;
    uint numRefines = 0;
    Mesh mesh(mesh0);
    
        
    // Coarsen until all elements are sufficiently refined        
    if ( strcmp(meshIn,"box")==0 ) {
        printf("\n Refine mesh until recovered gradient is sufficiently accurate\n");
        fflush(stdout);
    } else {
        printf("\n Skipping refinement for given mesh\n"); fflush(stdout);
        tooCoarse=false;
    }
        
    while (tooCoarse) {
        
        // Gradient recovery
        printf(" Recovering gradient of electric potential... \n"); fflush(stdout);
          Mesh meshAdapt(mesh);
          gradient_recovery::FunctionSpace GR(meshAdapt);
          gradient_recovery::BilinearForm aGR(GR,GR);
          gradient_recovery::LinearForm LGR(GR);
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
          poisson_cell_marker::FunctionSpace DG(meshAdapt);
          poisson_cell_marker::LinearForm errForm(DG);
          errForm.Du = Dsoln;
          errForm.u  = adaptFN;
        
        // Mark elements for refinement
        printf(" Marking elements for refinement\n "); fflush(stdout);
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
            mesh = refine(meshAdapt,cellMark);
            numRefines++;
            unrefined  = false;
        }
        
        // No elements marked for refinement
        else {
            printf("No elements marked for refinement... adaptivity complete\n");
            fflush(stdout);
            tooCoarse = false;
        }
        
    }
        
    // Solve complete
    if ( unrefined && converged ) {
        printf(" Converged to a solution \n\n"); fflush(stdout);
        //printf(" Review output files \n\n"); fflush(stdout);
        //return 0;
    }
    
    // Adapt Subdomains and Surfaces
    std::shared_ptr<const Mesh> new_mesh( new const Mesh(mesh) );
    MeshFunction<std::size_t> adapted_surfaces;
    MeshFunction<std::size_t> adapted_subdoms;
    
    if ( strcmp(meshIn,"box")==0 ) {
        // Count iterations
        totalRefines++;
        printf(" Adaptivity iteration %d\n\n\n",totalRefines); fflush(stdout);
    } else {
        // Read in surfaces and subdomains
        printf(" Reading in the mesh surfaces from %s \n", surfIn);
        MeshFunction<std::size_t> surfReadIn(mesh, surfIn);
        printf(" Reading in the mesh subdomains from %s \n", subdIn);
        MeshFunction<std::size_t> subdReadIn(mesh, subdIn);
        adapted_surfaces = surfReadIn;
        adapted_subdoms  = subdReadIn;

        printf(" Writing mesh\n");
        meshfile << adapted_surfaces;
        meshfile << adapted_subdoms;
    }

    
    
        
    
    //***********************************
    //   Finite Element Space and Forms
    //***********************************
        
    printf("\nDiscretize the PNP system \n"); fflush(stdout);
    
    // Finite element space
    printf(" Define PNP finite elements \n"); fflush(stdout);
    protein_four_species::FunctionSpace W(mesh);
    
    // Define variational forms
    printf(" Define PNP variational forms \n\n"); fflush(stdout);
    protein_four_species::BilinearForm a(W, W);
    protein_four_species::LinearForm L(W);
    
    
        
        
    //*********************************************
    //  Mark subdomains and impose Dirichlet B.C.
    //*********************************************
        
    // Define subdomains and subdomains
    printf("Define subdomains \n"); fflush(stdout);
      a.dx = adapted_subdoms;
      L.dx = adapted_subdoms;
      L.dS = adapted_surfaces;
      /*dielectricChannel sChannel;
      channelGate gChannel;
      FacetFunction<std::size_t> channelBndry(mesh);
      channelBndry.set_all(1);
      gChannel.mark(channelBndry,3);
      sChannel.mark(channelBndry,2);
      a.ds = channelBndry; L.ds = channelBndry;*/
    
    
    // Define Dirichlet boundary conditions
    printf(" Define Dirichlet boundary condition \n\n"); fflush(stdout);
      DirBCval DirBC;
      DirichletBC bc_3(W,DirBC,adapted_surfaces,4);
      DirichletBC bc_4(W,DirBC,adapted_surfaces,5);

        
    
    
    //***********************
    // Assign coefficients
    //***********************
    
    printf("Assign coefficients for variational forms \n\n"); fflush(stdout);
    
    // Nernst-Planck eqns
        Constant Dna(na_diff);    a.Dna = Dna; L.Dna = Dna;  // Diffusivity
        Constant qna(na_valency); a.qna = qna; L.qna = qna;  // Valency
        Constant Dk(k_diff);      a.Dk  = Dk;  L.Dk  = Dk;   // Diffusivity
        Constant qk(k_valency);   a.qk  = qk;  L.qk  = qk;   // Valency
        Constant Dcl(cl_diff);    a.Dcl = Dcl; L.Dcl = Dcl;  // Diffusivity
        Constant qcl(cl_valency); a.qcl = qcl; L.qcl = qcl;  // Valency
        Constant Dca(ca_diff);    a.Dca = Dca; L.Dca = Dca;  // Diffusivity
        Constant qca(ca_valency); a.qca = qca; L.qca = qca;  // Valency
    
    // Poisson eqn
        Constant eps_s(epsSolventVal);
        Constant eps_m(epsMolVal);
        a.eps_s = eps_s; L.eps_s = eps_s;           // Solvent permittivity
        a.eps_m = eps_m; L.eps_m = eps_m;           // Molecule permittivity
    
    
    
    
        
    //***********************************
    //  Interpolate solution iterate
    //***********************************
        
    printf("Define initial guess for Newton iteration \n"); fflush(stdout);
    Function iterate(W);
    Function concsoln(W);
    
        Function NaIterate(iterate[0]);
        NaIterate.interpolate(initNa);
        
        Function KIterate(iterate[1]);
        KIterate.interpolate(initK);
        
        Function CaIterate(iterate[2]);
        CaIterate.interpolate(initCa);
        
        Function ClIterate(iterate[3]);
        ClIterate.interpolate(initCl);
        
        // Convert from log to density
        Function NaConc(concsoln[0]);
        Function KConc(concsoln[1]);
        Function CaConc(concsoln[2]);
        Function ClConc(concsoln[3]);
        double * NaConcval = NaConc.vector()->data();
        double * NaItval = NaIterate.vector()->data();
        for(std::size_t i=0; i<NaIterate.vector()->size(); ++i) {
            NaConcval[i] = exp(NaItval[i]);
        }
        double * KConcval = KConc.vector()->data();
        double * KItval = KIterate.vector()->data();
        for(std::size_t i=0; i<KIterate.vector()->size(); ++i) {
            KConcval[i] = exp(KItval[i]);
        }
        double * CaConcval = CaConc.vector()->data();
        double * CaItval = CaIterate.vector()->data();
        for(std::size_t i=0; i<CaIterate.vector()->size(); ++i) {
            CaConcval[i] = exp(CaItval[i]);
        }
        double * ClConcval = ClConc.vector()->data();
        double * ClItval = ClIterate.vector()->data();
        for(std::size_t i=0; i<ClIterate.vector()->size(); ++i) {
            ClConcval[i] = exp(ClItval[i]);
        }
    
        // Previous potential
        printf(" Interpolate voltage drop\n"); fflush(stdout);
        Function esIterate(iterate[4]);
        esIterate.interpolate(initPHI);
        
        





    //************************************
    //  Solve for electrostatic potential
    //************************************
    printf("Solve for electric potential \n"); fflush(stdout);

    // Fixed point charge potential
    printf("  Solve for potential from fixed point charges\n"); fflush(stdout);
        protein_poisson::FunctionSpace V_charge(mesh);
        protein_poisson::BilinearForm a_charge(V_charge,V_charge);
        protein_poisson::LinearForm L_charge(V_charge);
        
        // Assign subdomains
        a_charge.dx = adapted_subdoms;
        L_charge.dx = adapted_subdoms;
        
        // Assign Functions
        printf("  Interpolate Coulomb potential\n"); fflush(stdout);
        Function phi_h(V_charge);
        Function phi_charge(V_charge);
        phi_charge.interpolate(coulomb);
        chargePotentialFile << phi_charge;
        
        // Assign coefficients
        Constant FixedCharges(0.0);
        a_charge.eps = eps_m;
        L_charge.fc  = FixedCharges;
        L_charge.phi_coulomb = phi_charge;
        
        // Dirichlet BC
        DirichletBC bc_charge_1(   V_charge,phi_charge,adapted_surfaces,2);
        DirichletBC bc_charge_2(   V_charge,phi_charge,adapted_surfaces,3);
        DirichletBC bc_charge_3(   V_charge,phi_charge,adapted_surfaces,4);
        DirichletBC bc_charge_4(   V_charge,phi_charge,adapted_surfaces,5);
        DirichletBC bc_charge_5(   V_charge,phi_charge,adapted_surfaces,6);
        DirichletBC bc_charge_6(   V_charge,phi_charge,adapted_surfaces,7);
        //DirichletBC bc_charge_bulk(V_charge,phi_charge,adapted_subdoms, 3);
        std::vector<const DirichletBC*> bc_charge;
        bc_charge.push_back(&bc_charge_1);
        bc_charge.push_back(&bc_charge_2);
        bc_charge.push_back(&bc_charge_3);
        bc_charge.push_back(&bc_charge_4);
        bc_charge.push_back(&bc_charge_5);
        //bc_charge.push_back(&bc_charge_bulk);

        
        // FEniCS solve
        /*solve( a_charge==L_charge, phi_h, bc_charge );
        chargePotentialFile << phi_h;
        (phi_charge.vector()) -= *(phi_h.vector());
        chargePotentialFile << phi_charge;*/
        
        // **************************************************
        // Interface with FASP
        // **************************************************
        Matrix charge_A;
        assemble(charge_A,a_charge);
        bc_charge_1.apply(charge_A);
        bc_charge_2.apply(charge_A);
        bc_charge_3.apply(charge_A);
        bc_charge_4.apply(charge_A);
        bc_charge_5.apply(charge_A);
        bc_charge_6.apply(charge_A);
        //bc_charge_bulk.apply(charge_A);
        Vector charge_b;
        assemble(charge_b,L_charge);
        bc_charge_1.apply(charge_b);
        bc_charge_2.apply(charge_b);
        bc_charge_3.apply(charge_b);
        bc_charge_4.apply(charge_b);
        bc_charge_5.apply(charge_b);
        bc_charge_6.apply(charge_b);
        //bc_charge_bulk.apply(charge_b);
        Function charge_solu(V_charge);
        
        printf(" -------------------------------------\n"); fflush(stdout);
        printf(" Start interface to FASP\n"); fflush(stdout);
        printf(" -------------------------------------\n\n"); fflush(stdout);
        
        printf(" Step 1: convert sparse matrix format\n"); fflush(stdout);
        // Convert charge_A to CSR
        dCSRmat charge_A_fasp;
        
        unsigned int charge_nz = boost::tuples::get<3>(charge_A.data());
        int charge_row = charge_A.size(0);
        int charge_col = charge_A.size(1);
        int* charge_ap = (int*)fasp_mem_calloc(charge_row+1, sizeof(int));
        const size_t* charge_ap_tmp = boost::tuples::get<0>(charge_A.data());
        for (int i=0; i<charge_row+1; i++) {
            charge_ap[i] = (int)charge_ap_tmp[i];
        }
        int* charge_ai = (int*)fasp_mem_calloc(charge_nz, sizeof(int));
        const size_t* charge_ai_tmp = boost::tuples::get<1>(charge_A.data());
        for (int i=0; i<charge_nz; i++) {
            charge_ai[i] = (int)charge_ai_tmp[i];
        }
        double* charge_ax = (double*)boost::tuples::get<2>(charge_A.data());
        
        charge_A_fasp.row = charge_row;
        charge_A_fasp.col = charge_col;
        charge_A_fasp.nnz = charge_nz;
        charge_A_fasp.IA  = charge_ap;
        charge_A_fasp.JA  = charge_ai;
        charge_A_fasp.val = charge_ax;
        
        // convert CSR to BSR
        dBSRmat charge_A_fasp_bsr = fasp_format_dcsr_dbsr(&charge_A_fasp, 1);
        
        // initialize RHS
        dvector charge_b_fasp;
        charge_b_fasp.row = charge_b.size();
        charge_b_fasp.val = (double*)charge_b.data();
        
        // initialize solution
        dvector charge_soluvec;
        fasp_dvec_alloc(charge_b_fasp.row, &charge_soluvec);
        fasp_dvec_set(charge_b_fasp.row, &charge_soluvec, 0.0);
        
        printf(" Step 2: initialize solver parameters\n"); fflush(stdout);
        // initialize solver parameters
        input_param     inpar;  // parameters from input files
        itsolver_param  itpar;  // parameters for itsolver
        AMG_param       amgpar; // parameters for AMG
        ILU_param       ilupar; // parameters for ILU
        
        // read in parameters from a input file
        char inputfile[] = "./params/bsr.dat";
        fasp_param_input(inputfile, &inpar);
        fasp_param_init(&inpar, &itpar, &amgpar, &ilupar, NULL);
        
        printf(" Step 3: solve the linear system\n"); fflush(stdout);
        // solve
        int status=FASP_SUCCESS;
        //fasp_param_amg_print(&amgpar);
        status = fasp_solver_dbsr_krylov_amg(&charge_A_fasp_bsr, &charge_b_fasp, &charge_soluvec, &itpar, &amgpar);
        
        if (status<0) {
            printf("\n### WARNING: Solver failed! Exit status = %d.\n\n", status); fflush(stdout);
        }
        else {
            printf("\nSolver finished successfully!\n\n"); fflush(stdout);
        }
        
        printf(" Step 4: convert solution back\n"); fflush(stdout);
        
        // Convert solution vector to FE solution
        double * charge_solval = phi_h.vector()->data();
        for(std::size_t i=0; i<charge_soluvec.row; ++i) {
            charge_solval[i] = charge_soluvec.val[i];
        }
        
        printf(" Step 5: free memory\n"); fflush(stdout);
        // Free memory
        fasp_dbsr_free(&charge_A_fasp_bsr);
        //fasp_dvec_free(&b_fasp);
        fasp_dvec_free(&charge_soluvec);
        free(charge_ap); free(charge_ai);
        
        printf(" -------------------------------------\n"); fflush(stdout);
        printf(" End of interface to FASP \n"); fflush(stdout);
        printf(" -------------------------------------\n\n"); fflush(stdout);
        
        printf(" Solved for fixed charge potential \n"); fflush(stdout);
        printf(" Assign to surface potential \n"); fflush(stdout);
        *(phi_charge.vector()) -= *(phi_h.vector());
        L.phi_surf = phi_charge;
        
        printf(" Write solution to file\n"); fflush(stdout);
        
        // Write to file
        chargePotentialFile << phi_h;
        chargePotentialFile << phi_charge;

        
        // Current
        /*
        printf(" Compute current in channel \n"); fflush(stdout);
        current_flux::Functional Current(mesh);
        double current;
        Constant n_vector(1.0,0.0,0.0);
        Current.n_vec = n_vector;
        Current.dS  = adapted_surfaces;
        Current.Dna = Dna; Current.qna = qna;
        Current.Dk  = Dk;  Current.qk  = qk;
        Current.Dcl = Dcl; Current.qcl = qcl;
        Current.Dca = Dca; Current.qca = qca;
        Current.Na  = NaIterate;
        Current.K   = KIterate;
        Current.Cl  = ClIterate;
        Current.Ca  = CaIterate;
        Current.Es  = esIterate;
        current = assemble(Current);
        printf("   Current is computed as %e pA / nM^2 \n\n", current); fflush(stdout);
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
    /*Function soln(W);
    Function dNa(soln[0]);
    Function dK(soln[1]);
    Function dCa(soln[2]);
    Function dCl(soln[3]);
    Function dphi(soln[4]);*/
    
    // Block GS updates
    Function solu(W);
    
    // Linear system
    Matrix A;
    Vector b;
    
        // Extract DoF indices
        printf(" Assemble solution indices \n"); fflush(stdout);
        std::vector<std::size_t> component(1);
        std::vector<dolfin::la_index> gidx_Na;
        std::vector<dolfin::la_index> gidx_K;
        std::vector<dolfin::la_index> gidx_Ca;
        std::vector<dolfin::la_index> gidx_Cl;
        std::vector<dolfin::la_index> gidx_phi;
        const dolfin::la_index n0 = W.dofmap()->ownership_range().first;
        const dolfin::la_index n1 = W.dofmap()->ownership_range().second;
        const dolfin::la_index num_dofs = n1 - n0;
        component[0] = 0;
        std::shared_ptr<GenericDofMap> dofmap_Na  = W.dofmap()->extract_sub_dofmap(component,mesh);
        component[0] = 1;
        std::shared_ptr<GenericDofMap> dofmap_K   = W.dofmap()->extract_sub_dofmap(component,mesh);
        component[0] = 2;
        std::shared_ptr<GenericDofMap> dofmap_Ca  = W.dofmap()->extract_sub_dofmap(component,mesh);
        component[0] = 3;
        std::shared_ptr<GenericDofMap> dofmap_Cl  = W.dofmap()->extract_sub_dofmap(component,mesh);
        component[0] = 4;
        std::shared_ptr<GenericDofMap> dofmap_phi = W.dofmap()->extract_sub_dofmap(component,mesh);
        
        for ( CellIterator cell(mesh); !cell.end(); ++cell)
        {
            const std::vector<dolfin::la_index> cell_dofs_Na  = dofmap_Na->cell_dofs(cell->index());
            const std::vector<dolfin::la_index> cell_dofs_K   = dofmap_K->cell_dofs(cell->index());
            const std::vector<dolfin::la_index> cell_dofs_Ca  = dofmap_Ca->cell_dofs(cell->index());
            const std::vector<dolfin::la_index> cell_dofs_Cl  = dofmap_Cl->cell_dofs(cell->index());
            const std::vector<dolfin::la_index> cell_dofs_phi = dofmap_phi->cell_dofs(cell->index());
            for (std::size_t i = 0; i < cell_dofs_Na.size(); ++i)
            {
                const std::size_t dof = cell_dofs_Na[i];
                if (dof >= n0 && dof < n1)
                    gidx_Na.push_back(dof);
            }
            for (std::size_t i = 0; i < cell_dofs_K.size(); ++i)
            {
                const std::size_t dof = cell_dofs_K[i];
                if (dof >= n0 && dof < n1)
                    gidx_K.push_back(dof);
            }
            for (std::size_t i = 0; i < cell_dofs_Ca.size(); ++i)
            {
                const std::size_t dof = cell_dofs_Ca[i];
                if (dof >= n0 && dof < n1)
                    gidx_Ca.push_back(dof);
            }
            for (std::size_t i = 0; i < cell_dofs_Cl.size(); ++i)
            {
                const std::size_t dof = cell_dofs_Cl[i];
                if (dof >= n0 && dof < n1)
                    gidx_Cl.push_back(dof);
            }
            for (std::size_t i = 0; i < cell_dofs_phi.size(); ++i)
            {
                const std::size_t dof = cell_dofs_phi[i];
                if (dof >= n0 && dof < n1)
                    gidx_phi.push_back(dof);
            }
        }
        std::sort(gidx_Na.begin(), gidx_Na.end());
        std::sort(gidx_K.begin(), gidx_K.end());
        std::sort(gidx_Ca.begin(), gidx_Ca.end());
        std::sort(gidx_Cl.begin(), gidx_Cl.end());
        std::sort(gidx_phi.begin(), gidx_phi.end());
        // Remove duplicates
        gidx_Na.erase(std::unique(gidx_Na.begin(),   gidx_Na.end()),  gidx_Na.end());
        gidx_K.erase(std::unique(gidx_K.begin(),     gidx_K.end()),   gidx_K.end());
        gidx_Ca.erase(std::unique(gidx_Ca.begin(),   gidx_Ca.end()),  gidx_Ca.end());
        gidx_Cl.erase(std::unique(gidx_Cl.begin(),   gidx_Cl.end()),  gidx_Cl.end());
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
        L.NaNa = NaIterate;
        L.KK   = KIterate;
        L.CaCa = CaIterate;
        L.ClCl = ClIterate;
        L.EsEs = esIterate;
    assemble(b,L); bc_3.apply(b); bc_4.apply(b);
        //bc_membrane_0.apply(b); bc_membrane_1.apply(b); bc_membrane_2.apply(b); bc_membrane_3.apply(b); bc_protein_0.apply(b);
        //bc_protein_1.apply(b); bc_protein_2.apply(b); bc_protein_3.apply(b);
    
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
        Nafile  << NaConc;
        Kfile   << KConc;
        Cafile  << CaConc;
        Clfile  << ClConc;
        phifile << esIterate;
        
        
        
        
        // Update newton step
        it++;
        printf("Newton iteration %d\n", it);
        fflush(stdout);
        
        // Update PNP coefficients
        printf(" Construct Jacobian matrix \n"); fflush(stdout);
        a.NaNa = NaIterate;
        a.KK   = KIterate;
        a.CaCa = CaIterate;
        a.ClCl = ClIterate;
        a.EsEs = esIterate;
        assemble(A,a); bc_3.apply(A); bc_4.apply(A);
        //bc_membrane_0.apply(A); bc_membrane_1.apply(A); bc_membrane_2.apply(A); bc_membrane_3.apply(A);
        //bc_protein_0.apply(A); bc_protein_1.apply(A); bc_protein_2.apply(A); bc_protein_3.apply(A);
        
        
        
        // **************************************************
        // Interface with FASP
        // **************************************************
        printf(" -------------------------------------\n"); fflush(stdout);
        printf(" Start interface to FASP \n"); fflush(stdout);
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
        
        // Check for rows of zeros and add a unit diagonal entry
        bool nonzero_entry = false;
        for ( uint rowInd=0; rowInd<row; rowInd++ ) {
            
            // Check for nonzero entry
            nonzero_entry   = false;
            int diagColInd = -1;
            if ( ap[rowInd] < ap[rowInd+1] ) {
                for ( uint colInd=ap[rowInd]; colInd < ap[rowInd+1]; colInd++ ) {
                    if ( ax[colInd] != 0.0 )    nonzero_entry = true;
                    if ( ai[colInd] == rowInd ) diagColInd    = colInd;    // marks diagonal entry
                }
            }
            
            
            if ( diagColInd < 0 ) {
                printf(" ERROR: diagonal entry not allocated!!\n\n Exiting... \n \n"); fflush(stdout);
                printf("      for row %d\n",rowInd); fflush(stdout);
                for ( uint colInd=ap[rowInd]; colInd < ap[rowInd+1]; colInd++ ) {
                    printf("          %d\n",ai[colInd]);
                }
                return 0;
            }
            
            if ( nonzero_entry==false ) {
                //printf(" Row %d has only zeros! Setting diagonal entry to 1.0 \n", rowInd); fflush(stdout);
                ax[diagColInd] = 1.0;
                //b_fasp.val[rowInd] = -1.e+5;
            }
        }

        A_fasp.row = row;
        A_fasp.col = col;
        A_fasp.nnz = nz;
        A_fasp.IA  = ap;
        A_fasp.JA  = ai;
        A_fasp.val = ax;
        
#if FASP_BSR
        // convert CSR to BSR
        dBSRmat A_fasp_bsr = fasp_format_dcsr_dbsr(&A_fasp, 5);
        
        // output
        /*
         fasp_dbsr_write("A_fasp_bsr.dat", &A_fasp_bsr);
         fasp_dvec_write("b_fasp.dat", &b_fasp);
         */
#else
        // convert CSR to block CSR
        block_dCSRmat Abcsr;
        dCSRmat *A_diag;
        
        // get index
        INT nrow = A_fasp.row/5;
        ivector phi_idx;
        ivector Cl_idx;
        ivector Ca_idx;
        ivector K_idx;
        ivector Na_idx;
        
        fasp_ivec_alloc(nrow, &phi_idx);
        fasp_ivec_alloc(nrow, &Cl_idx);
        fasp_ivec_alloc(nrow, &Ca_idx);
        fasp_ivec_alloc(nrow, &K_idx);
        fasp_ivec_alloc(nrow, &Na_idx);
        
        for (i=0; i<nrow; i++){
            phi_idx.val[i] = 5*i;
            Cl_idx.val[i]  = 5*i+1;
            Ca_idx.val[i]  = 5*i+2;
            K_idx.val[i]   = 5*i+3;
            Na_idx.val[i]  = 5*i+4;
        }
        
        // Assemble the matrix in block dCSR format
        Abcsr.brow = 5; Abcsr.bcol = 5;
        Abcsr.blocks = (dCSRmat **)calloc(9, sizeof(dCSRmat *));
        for (i=0; i<9 ;i++) {
            Abcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
        
        // A11
        fasp_dcsr_getblk(&A_fasp, phi_idx.val, phi_idx.val, nrow, nrow, Abcsr.blocks[0]);
        // A12
        fasp_dcsr_getblk(&A_fasp, phi_idx.val, Cl_idx.val, nrow, nrow, Abcsr.blocks[1]);
        // A13
        fasp_dcsr_getblk(&A_fasp, phi_idx.val, Ca_idx.val, nrow, nrow, Abcsr.blocks[2]);
        // A14
        fasp_dcsr_getblk(&A_fasp, phi_idx.val, K_idx.val, nrow, nrow, Abcsr.blocks[3]);
        // A15
        fasp_dcsr_getblk(&A_fasp, phi_idx.val, Na_idx.val, nrow, nrow, Abcsr.blocks[4]);
        
        // A21
        fasp_dcsr_getblk(&A_fasp, Cl_idx.val, phi_idx.val, nrow, nrow, Abcsr.blocks[5]);
        // A22
        fasp_dcsr_getblk(&A_fasp, Cl_idx.val, Cl_idx.val, nrow, nrow, Abcsr.blocks[6]);
        // A23
        fasp_dcsr_getblk(&A_fasp, Cl_idx.val, Ca_idx.val, nrow, nrow, Abcsr.blocks[7]);
        // A24
        fasp_dcsr_getblk(&A_fasp, Cl_idx.val, K_idx.val,  nrow, nrow, Abcsr.blocks[8]);
        // A25
        fasp_dcsr_getblk(&A_fasp, Cl_idx.val, Na_idx.val, nrow, nrow, Abcsr.blocks[9]);
        
        // A31
        fasp_dcsr_getblk(&A_fasp, Ca_idx.val, phi_idx.val, nrow, nrow, Abcsr.blocks[10]);
        // A32
        fasp_dcsr_getblk(&A_fasp, Ca_idx.val, Cl_idx.val, nrow, nrow, Abcsr.blocks[11]);
        // A33
        fasp_dcsr_getblk(&A_fasp, Ca_idx.val, Ca_idx.val, nrow, nrow, Abcsr.blocks[12]);
        // A34
        fasp_dcsr_getblk(&A_fasp, Ca_idx.val, K_idx.val,  nrow, nrow, Abcsr.blocks[13]);
        // A35
        fasp_dcsr_getblk(&A_fasp, Ca_idx.val, Na_idx.val, nrow, nrow, Abcsr.blocks[14]);
        
        // A41
        fasp_dcsr_getblk(&A_fasp, K_idx.val, phi_idx.val, nrow, nrow, Abcsr.blocks[15]);
        // A42
        fasp_dcsr_getblk(&A_fasp, K_idx.val, Cl_idx.val, nrow, nrow, Abcsr.blocks[16]);
        // A43
        fasp_dcsr_getblk(&A_fasp, K_idx.val, Ca_idx.val, nrow, nrow, Abcsr.blocks[17]);
        // A44
        fasp_dcsr_getblk(&A_fasp, K_idx.val, K_idx.val,  nrow, nrow, Abcsr.blocks[18]);
        // A45
        fasp_dcsr_getblk(&A_fasp, K_idx.val, Na_idx.val, nrow, nrow, Abcsr.blocks[19]);
        
        // A51
        fasp_dcsr_getblk(&A_fasp, Na_idx.val, phi_idx.val, nrow, nrow, Abcsr.blocks[20]);
        // A52
        fasp_dcsr_getblk(&A_fasp, Na_idx.val, Cl_idx.val, nrow, nrow, Abcsr.blocks[21]);
        // A53
        fasp_dcsr_getblk(&A_fasp, Na_idx.val, Ca_idx.val, nrow, nrow, Abcsr.blocks[22]);
        // A54
        fasp_dcsr_getblk(&A_fasp, Na_idx.val, K_idx.val,  nrow, nrow, Abcsr.blocks[23]);
        // A55
        fasp_dcsr_getblk(&A_fasp, Na_idx.val, Na_idx.val, nrow, nrow, Abcsr.blocks[24]);
        
        
        // setup diagonal blocks for the preconditioner
        A_diag = (dCSRmat *)fasp_mem_calloc(5, sizeof(dCSRmat));
        
        // first diagonal block
        A_diag[0].row = Abcsr.blocks[0]->row;
        A_diag[0].col = Abcsr.blocks[0]->col;
        A_diag[0].nnz = Abcsr.blocks[0]->nnz;
        A_diag[0].IA  = Abcsr.blocks[0]->IA;
        A_diag[0].JA  = Abcsr.blocks[0]->JA;
        A_diag[0].val = Abcsr.blocks[0]->val;
        
        // second diagonal block
        A_diag[1].row = Abcsr.blocks[6]->row;
        A_diag[1].col = Abcsr.blocks[6]->col;
        A_diag[1].nnz = Abcsr.blocks[6]->nnz;
        A_diag[1].IA  = Abcsr.blocks[6]->IA;
        A_diag[1].JA  = Abcsr.blocks[6]->JA;
        A_diag[1].val = Abcsr.blocks[6]->val;
        
        // third diagonal block
        A_diag[2].row = Abcsr.blocks[12]->row;
        A_diag[2].col = Abcsr.blocks[12]->col;
        A_diag[2].nnz = Abcsr.blocks[12]->nnz;
        A_diag[2].IA  = Abcsr.blocks[12]->IA;
        A_diag[2].JA  = Abcsr.blocks[12]->JA;
        A_diag[2].val = Abcsr.blocks[12]->val;
        
        // fourth diagonal block
        A_diag[3].row = Abcsr.blocks[18]->row;
        A_diag[3].col = Abcsr.blocks[18]->col;
        A_diag[3].nnz = Abcsr.blocks[18]->nnz;
        A_diag[3].IA  = Abcsr.blocks[18]->IA;
        A_diag[3].JA  = Abcsr.blocks[18]->JA;
        A_diag[3].val = Abcsr.blocks[18]->val;
        
        // fifth diagonal block
        A_diag[4].row = Abcsr.blocks[24]->row;
        A_diag[4].col = Abcsr.blocks[24]->col;
        A_diag[4].nnz = Abcsr.blocks[24]->nnz;
        A_diag[4].IA  = Abcsr.blocks[24]->IA;
        A_diag[4].JA  = Abcsr.blocks[24]->JA;
        A_diag[4].val = Abcsr.blocks[24]->val;
        
        
        //fasp_dcoo_write("aDiag0",&A_diag[0]);
        //fasp_dcoo_write("aDiag1",&A_diag[1]);
        //fasp_dcoo_write("aDiag2",&A_diag[2]);
        
        
        
        // convert right hand side
        dvector bbcsr;
        fasp_dvec_alloc(b_fasp.row, &bbcsr);
        for (i=0; i<nrow; i++){
            bbcsr.val[i]        = b_fasp.val[5*i];
            bbcsr.val[nrow+i]   = b_fasp.val[5*i+1];
            bbcsr.val[2*nrow+i] = b_fasp.val[5*i+2];
            bbcsr.val[3*nrow+i] = b_fasp.val[5*i+3];
            bbcsr.val[4*nrow+i] = b_fasp.val[5*i+4];
            
        }
        
        // output the matrices
        /*
         fasp_dcoo_write("A.dat", &A_fasp);
         fasp_ivec_write("nidx.dat", &n_idx);
         fasp_ivec_write("phiidx.dat", &phi_idx);
         fasp_ivec_write("pidx.dat", &p_idx);
         fasp_dvec_write("rhs.dat", &bbcsr);
         
         getchar();
         */
        
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
            solval[5*i]    = soluvec.val[i];
            solval[5*i+1]  = soluvec.val[nrow+i];
            solval[5*i+2]  = soluvec.val[2*nrow+i];
            solval[5*i+3]  = soluvec.val[3*nrow+i];
            solval[5*i+4]  = soluvec.val[4*nrow+i];
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
        fasp_ivec_free(&Cl_idx);
        fasp_ivec_free(&Ca_idx);
        fasp_ivec_free(&K_idx);
        fasp_ivec_free(&Na_idx);
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
        Function NewtonUpdate(W);
        Function dNa(NewtonUpdate[0]);  dNa.interpolate(solu[0]);
        Function dK(NewtonUpdate[1]);   dK.interpolate(solu[1]);
        Function dCa(NewtonUpdate[2]);  dCa.interpolate(solu[2]);
        Function dCl(NewtonUpdate[3]);  dCl.interpolate(solu[3]);
        Function dphi(NewtonUpdate[4]); dphi.interpolate(solu[4]);
        
        // Update solutions
        /*
        *(NaIterate.vector()) += *(dNa.vector());
        *(KIterate.vector())  += *(dK.vector());
        *(CaIterate.vector()) += *(dCa.vector());
        *(ClIterate.vector()) += *(dCl.vector());
        *(esIterate.vector()) += *(dphi.vector());
        */
        
        /*
         // measure size of update
         *(p.vector())   *= mu;
         *(n.vector())   *= mu;
         *(phi.vector()) *= mu;
         normp   = p.vector()->norm("l2");
         normn   = n.vector()->norm("l2");
         normphi = phi.vector()->norm("l2");
         
         // update solutions
         *(pIterate.vector()) += *(p.vector());
         *(nIterate.vector()) += *(n.vector());
         *(esIterate.vector()) += *(phi.vector());
         */
        
        
        //  Backtrack line search
        printf("   Use backtracking to guarantee decreasing residual \n"); fflush(stdout);
        Function update(W);
        Function updateNA(update[0]);
        Function updateK(update[1]);
        Function updateCA(update[2]);
        Function updateCL(update[3]);
        Function updatePHI(update[4]);
        double testRelRes  = relR_fasp+DOLFIN_EPS;//1.0+DOLFIN_EPS;
        double dampFactor  = mu;
        bool   reducing  = true;   // loop boolean
        bool   reduced   = false;  // residual is reduced?
        bool   tooDamped = false;  // dampFactor too small?
        while ( reducing ) {//testRelRes > relR_fasp-DOLFIN_EPS && dampFactor > 1.e-4 ) {
            
            //  Compute update
            updateNA.interpolate(NaIterate);
            updateK.interpolate(KIterate);
            updateCA.interpolate(CaIterate);
            updateCL.interpolate(ClIterate);
            updatePHI.interpolate(esIterate);
            *(updateNA.vector())  += *(dNa.vector());
            *(updateK.vector())   += *(dK.vector());
            *(updateCA.vector())  += *(dCa.vector());
            *(updateCL.vector())  += *(dCl.vector());
            *(updatePHI.vector()) += *(dphi.vector());
            
            // Evaluate residual
            L.NaNa = updateNA;
            L.KK   = updateK;
            L.CaCa = updateCA;
            L.ClCl = updateCL;
            L.EsEs = updatePHI;
            assemble(b,L); bc_3.apply(b); bc_4.apply(b);
            //bc_membrane_0.apply(b); bc_membrane_1.apply(b); bc_membrane_2.apply(b); bc_membrane_3.apply(b);
            //bc_protein_0.apply(b); bc_protein_1.apply(b); bc_protein_2.apply(b); bc_protein_3.apply(b);
            
            // Compute relative residual for backtrack line search
            dvector newb_faspBT;
            newb_faspBT.row = b.size();
            newb_faspBT.val = (double*)b.data();
            double normb_faspBT = fasp_blas_array_norm2(newb_faspBT.row,newb_faspBT.val);
            testRelRes   = normb_faspBT/b0_fasp ;
            printf("   Residual after update is %e \n", testRelRes);
            fflush(stdout);
            
            // check for reduced residual or NaN
            if ( (testRelRes>relR_fasp-DOLFIN_EPS) || (testRelRes!=testRelRes) ) {
                printf("   Reducing the update by a factor of %e \n", dampFactor); fflush(stdout);
                dampFactor            *= 0.5;
                *(dNa.vector())       *= dampFactor;
                *(dK.vector())        *= dampFactor;
                *(dCa.vector())       *= dampFactor;
                *(dCl.vector())       *= dampFactor;
                *(dphi.vector())      *= dampFactor;
            } else {// Significant reduction
                reduced = true;
            }
            if (dampFactor<1.0e-7)   tooDamped = true;  // too much backtracking
            if (reduced || tooDamped) reducing = false; // end while
        }
        
        // If backtrack failed
        if ( testRelRes > relR_fasp-DOLFIN_EPS ){
            printf("   BAD UPDATE!!! \n"); fflush(stdout);
        }

        
        // Update solutions
        *(NaIterate.vector()) = *(updateNA.vector());
        *(KIterate.vector())  = *(updateK.vector());
        *(CaIterate.vector()) = *(updateCA.vector());
        *(ClIterate.vector()) = *(updateCL.vector());
        *(esIterate.vector()) = *(updatePHI.vector());
        


        // Update residual
        printf(" Update variational forms \n"); fflush(stdout);
        L.NaNa = NaIterate;
        L.KK   = KIterate;
        L.CaCa = CaIterate;
        L.ClCl = ClIterate;
        L.EsEs = esIterate;
        assemble(b,L); bc_3.apply(b); bc_4.apply(b);
        //bc_membrane_0.apply(b); bc_membrane_1.apply(b); bc_membrane_2.apply(b); bc_membrane_3.apply(b);
        //bc_protein_0.apply(b); bc_protein_1.apply(b); bc_protein_2.apply(b); bc_protein_3.apply(b);
        
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

        // Current
        printf(" Compute current in channel \n"); fflush(stdout);
        current_flux::Functional Current(mesh);
        Constant n_vector(1.0,0.0,0.0);
        Current.n_vec = n_vector;
        Current.dS  = adapted_surfaces;
        Current.Dna = Dna; Current.qna = qna;
        Current.Dk  = Dk;  Current.qk  = qk;
        Current.Dcl = Dcl; Current.qcl = qcl;
        Current.Dca = Dca; Current.qca = qca;
        Current.Na  = NaIterate;
        Current.K   = KIterate;
        Current.Cl  = ClIterate;
        Current.Ca  = CaIterate;
        Current.Es  = esIterate;
        current  = assemble(Current);
        current *= 1.0e-12*D_ref*p_ref*n_avo*e_chrg/L_ref;
        printf("   Current is computed as %e pA \n\n", current); fflush(stdout);
        
        
        
        //************************
        //  Measure dissipation
        //************************
        
        // Convert from log to density
        double * NaConcval = NaConc.vector()->data();
        double * NaItval = NaIterate.vector()->data();
        for(std::size_t i=0; i<NaIterate.vector()->size(); ++i) {
            NaConcval[i] = exp(NaItval[i]);
        }
        double * KConcval = KConc.vector()->data();
        double * KItval = KIterate.vector()->data();
        for(std::size_t i=0; i<KIterate.vector()->size(); ++i) {
            KConcval[i] = exp(KItval[i]);
        }
        double * CaConcval = CaConc.vector()->data();
        double * CaItval = CaIterate.vector()->data();
        for(std::size_t i=0; i<CaIterate.vector()->size(); ++i) {
            CaConcval[i] = exp(CaItval[i]);
        }
        double * ClConcval = ClConc.vector()->data();
        double * ClItval = ClIterate.vector()->data();
        for(std::size_t i=0; i<ClIterate.vector()->size(); ++i) {
            ClConcval[i] = exp(ClItval[i]);
        }
        
        
        // Save solution in VTK format
        /*printf(" Write solution to file \n"); fflush(stdout);
         Nafile  << NaConc;
         Kfile   << KConc;
         Cafile  << CaConc;
         Clfile  << ClConc;
         phifile << esIterate;
         
         
         Diss.Na = NaConc;
         Diss.K  = KConc;
         Diss.Ca = CaConc;
         Diss.Cl = ClConc;
         Diss.phi = esIterate;
         diss     = assemble(Diss);
         printf(" The dissipation is %e \n\n", diss);
         fflush(stdout);
         */
        
        
        
  
        
        
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
            Nafile  << NaConc;
            Kfile   << KConc;
            Cafile  << CaConc;
            Clfile  << ClConc;
            phifile << esIterate;
            
            //free(pt_charge_x); free(pt_charge_y); free(pt_charge_z); free(pt_charge_q);

            
            
            //  Write simulation summary to file
            /*
             printf(" Output for main-log-steady.cpp \n");fflush(stdout);
             FILE *fout;
             strcat(outputfile,"-log-main-steady");
             fout = fopen(outputfile, "w");
             fprintf(fout,"Output for main-log-steady.cpp \n\n");
             fprintf(fout," Mesh: ( %d, %d, %d ) \n", Nx,Ny,Nz);
             fprintf(fout," eps = %e \n", epsSolventVal);
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
            
            
        } else if ( it>=maxit ) {                    // Too many iterations
            printf("\nDivergence: The relative residual has not converged after %d iterations. \n\n",maxit+1);
            fflush(stdout);
            done=true;
        }  else if ( relR_fasp >= prev_relR ) {     // Increasing residual
            printf(" *** The relative residual has grown in this iteration *** \n\n");fflush(stdout);
        }

    
    }
    
    
    printf(" -------------------------------------\n"); fflush(stdout);
    printf(" End of Newton Iteration \n"); fflush(stdout);
    printf(" -------------------------------------\n\n"); fflush(stdout);
    
    printf(" The current at %8.4f mV is %e pA.\n\n", volts, current); fflush(stdout);    
    b0_fasp = -1.0;    
    
    fprintf(fout,"%e,%e\n", volts, current);
        
    
    //********************
    //  Adaptivity pass
    //********************
    
    /*printf(" Update mesh \n\n"); fflush(stdout);
    mesh0   = mesh;
    std::shared_ptr<const Mesh> adapted_mesh( new const Mesh(mesh0) );
    adaptFN = adapt(esIterate, adapted_mesh);
    initNa  = adapt(NaIterate, adapted_mesh);
    initK   = adapt(KIterate,  adapted_mesh);
    initCa  = adapt(CaIterate, adapted_mesh);
    initCl  = adapt(ClIterate, adapted_mesh);
    initPHI = adapt(esIterate, adapted_mesh);*/

    
    }
    
    
    //printf(" Failed to converge below desired relative residual!!! \n\n"); fflush(stdout);
    printf(" Done computing IV-curve... Exiting. \n\n"); fflush(stdout);
    fclose(fout);
    //}
    
    free(pt_charge_x); free(pt_charge_y); free(pt_charge_z); free(pt_charge_q);
   
    
    return 0;
}