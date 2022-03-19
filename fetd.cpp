// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.
//
//    ------------------------------------------------------------------
//    Maxwell Miniapp:  Simple Full-Wave Electromagnetic Simulation Code
//    ------------------------------------------------------------------
//
// This miniapp solves a simple 3D full-wave electromagnetic problem using the
// coupled, first-order equations:
//
//                 epsilon dE/dt = Curl 1/mu B - sigma E - J
//                         dB/dt = - Curl E
//
// The permittivity function is that of the vacuum with an optional dielectric
// sphere. The permeability function is that of the vacuum with an optional
// diamagnetic or paramagnetic spherical shell. The optional conductivity
// function is also a user-defined sphere.
//
// The optional current density is a pulse of current in the shape of a cylinder
// with a time dependence resembling the derivative of a Gaussian distribution.
//
// Boundary conditions can be 'natural' meaning zero tangential current,
// 'Dirichlet' which sets the time-derivative of the tangential components of E,
// or 'absorbing' (we use a simple Sommerfeld first order absorbing boundary
// condition).
//
// We discretize the electric field with H(Curl) finite elements (Nedelec edge
// elements) and the magnetic flux with H(Div) finite elements (Raviart-Thomas
// elements).
//
// The symplectic time integration algorithm used below is designed to conserve
// energy unless lossy materials or absorbing boundary conditions are used.
// When losses are expected, the algorithm uses an implicit method which
// includes the loss operators in the left hand side of the linear system.
//
// For increased accuracy the time integration order can be set to 2, 3, or 4
// (the default is 1st order).
//
// Compile with: make fetd

#include "cxxtimer.hpp"
#include "fetd_solver.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;
using namespace mfem::common;
using namespace mfem::electromagnetics;


double eps_freespace(const Vector &);
double epsilon(const Vector &x) { return eps_freespace(x); }

double mu_freespace(const Vector &);
double muInv(const Vector & x) { return 1.0/mu_freespace(x); }

double sigma_freespace(const Vector &);
double sigma(const Vector &x) { return sigma_freespace(x); }

void j_current_free(const Vector &x, double t, Vector &j);
void j_src(const Vector &x, double t, Vector &j) { j_current_free(x, t, j); }

// dE/dt Boundary Condition: The following function returns zero but any time
// dependent function could be used.
void dEdt_GausssianPulse_BCFunc(const Vector &x, double t, Vector &E);

// The following functions return zero but they could be modified to set initial
// conditions for the electric and magnetic fields
void EFieldFunc(const Vector &, Vector&);
void BFieldFunc(const Vector &, Vector&);

// Scale factor between input time units and seconds
static double tScale_ = 1e-9;  // Input time in nanosecond

int SnapTimeStep(double tmax, double dtmax, double & dt);

// Prints the program's logo to the given output stream
void display_banner(ostream & os);

int main(int argc, char *argv[])
{

   display_banner(cout);
   timer_start("the DGTD main function", '~'); // minutes
   // Parse command-line options.

   const char *mesh_file = "../ComputerCase_noLoad_empty_aperture_merged.g"; // cubit export setting: 3d, large file format

   int sOrder = 1;
   int tOrder = 1;
   bool visualization = false;
   bool visit = true;
   double dt;
   double dtsf = 0.95;
   double ti = 0.0;  // ns
   double ts = 0.1e-9;  // ns
   double tf = 10e-9; // ns

   Array<int> abcs;
   Array<int> dbcs_pec,dbcs_dE;

   int _abcs[] = {2}; // 2: abc
   int _dbcs_pec[] = {4}; // 4:pec
   int _dbcs_dE[] = {3}; // 3:inc

   abcs = Array<int>(_abcs,1);
   dbcs_pec = Array<int>(_dbcs_pec,1);
   dbcs_dE = Array<int>(_dbcs_dE,1);

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&sOrder, "-so", "--spatial-order", "Finite element order (polynomial degree).");
   args.AddOption(&tOrder, "-to", "--temporal-order", "Time integration order.");
   args.AddOption(&dtsf, "-sf", "--dt-safety-factor", "Used to reduce the time step below the upper bound.");
   args.AddOption(&ti, "-ti", "--initial-time", "Beginning of time interval to simulate (ns).");
   args.AddOption(&tf, "-tf", "--final-time", "End of time interval to simulate (ns).");
   args.AddOption(&ts, "-ts", "--snapshot-time", "Time between snapshots (ns).");
   
   args.AddOption(&abcs, "-abcs", "--absorbing-bc-surf", "Absorbing Boundary Condition Surfaces"); // ABC
   args.AddOption(&dbcs_pec, "-dbcs_pec", "--dirichlet-bc-surf-pec", "Dirichlet Boundary Condition Surfaces"); // PEC
   args.AddOption(&dbcs_dE, "-dbcs_dE", "--dirichlet-bc-surf-dE", "Dirichlet Boundary Condition Surfaces"); // inc

   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis", "--no-visualization", "Enable or disable GLVis visualization.");
   args.AddOption(&visit, "-visit", "--visit", "-no-visit", "--no-visualization", "Enable or disable VisIt visualization.");
   args.Parse();

   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);


   // Read the (serial) mesh from the given mesh file on all processors.  We can
   // handle triangular, quadrilateral, tetrahedral, hexahedral, surface and
   // volume meshes with the same code.
   MEM_USAGE();
   named_ifgzstream imesh(mesh_file);
   if (!imesh)
   {
      cerr << "\nCan not open mesh file: " << mesh_file << '\n' << endl;
      return 2;
   }
   Mesh mesh = Mesh(imesh, 1, 1);
   MEM_USAGE();
   // Create the Electromagnetic solver
   // relative dielectric constant， relative permeability， conductivity e constant: read from .m file
   // electric current J: zero
   // ABC: one of the side set
   // PEC: one of the side set
   // dEdt: one of the side set, should represent a gaussian pulse plane wave incident,
   // field plane: one of the side set, additional information to be addedl
   MaxwellSolver FETD(mesh,
                      sOrder,
                      epsilon,
                      muInv,
                      sigma,
                      j_src,
                      abcs,
                      dbcs_pec,
                      dbcs_dE,
                      dEdt_GausssianPulse_BCFunc);

   // Display the current number of DoFs in each finite element space
   FETD.PrintSizes();
   MEM_USAGE();
   // Set the initial conditions for both the electric and magnetic fields
   VectorFunctionCoefficient EFieldCoef(3,EFieldFunc);
   VectorFunctionCoefficient BFieldCoef(3,BFieldFunc);

   FETD.SetInitialEField(EFieldCoef);
   FETD.SetInitialBField(BFieldCoef);

   // Compute the energy of the initial fields
   double energy = FETD.GetEnergy();

   cout << "Energy(" << ti << "ns):  " << energy << "J" << endl;


   // Approximate the largest stable time step
   double dtmax = FETD.GetMaximumTimeStep();

   cout.setf(ios::scientific,ios_base::floatfield);
   cout.precision(10);
   cout << "Maximum Time Step:     " << dtmax / tScale_ << "ns" << endl;


   // Round down the time step so that tf-ti is an integer multiple of dt
   int nsteps = SnapTimeStep(tf-ti, dtsf * dtmax, dt);

   cout << "Number of Time Steps:  " << nsteps << endl;
   cout << "Time Step Size:        " << dt / tScale_ << "ns" << endl;


   // Create the ODE solver
   SIAVSolver siaSolver(tOrder);
   siaSolver.Init(FETD.GetNegCurl(), FETD);

   // Initialize GLVis visualization
   if (visualization)
   {
      FETD.InitializeGLVis();
   }

   // Initialize VisIt visualization
   VisItDataCollection visit_dc("FETD", &mesh);

   double t = ti;
   FETD.SetTime(t);

   if ( visit )
   {
      FETD.RegisterVisItFields(visit_dc);
   }

   // Write initial fields to disk for VisIt
   if ( visit )
   {
      FETD.WriteVisItFields(0);
   }

   // Send the initial condition by socket to a GLVis server.
   if (visualization)
   {
      FETD.DisplayToGLVis();
   }

   // The main time evolution loop. (Time Stepping)
   int it = 1;
   while (t < tf)
   {
      MEM_USAGE();
      timer_start("Step "+ to_string(it), 'm');
      // Run the simulation until a snapshot is needed
      siaSolver.Run(FETD.GetBField(), FETD.GetEField(), t, dt, max(t + dt, ti + ts * it));

      // Approximate the current energy if the fields
      energy = FETD.GetEnergy();

      cout << "Energy(" << t/tScale_ << "ns):  " << energy << "J" << endl;

      // Update local DoFs with current true DoFs
      FETD.SyncGridFuncs();

      // Write fields to disk for VisIt
      if ( visit )
      {
         FETD.WriteVisItFields(it);
      }

      // Send the solution by socket to a GLVis server.
      if (visualization)
      {
         FETD.DisplayToGLVis();
      }
      it++;
      timer_stop('m');
   }
   timer_stop('~');
   return 0;
}

// Print the Maxwell ascii logo to the given ostream
void display_banner(ostream & os)
{
   os << "███████╗███████╗████████╗██████╗ \n"
         "██╔════╝██╔════╝╚══██╔══╝██╔══██╗\n"
         "█████╗  █████╗     ██║   ██║  ██║\n"
         "██╔══╝  ██╔══╝     ██║   ██║  ██║\n"
         "██║     ███████╗   ██║   ██████╔╝\n"
         "╚═╝     ╚══════╝   ╚═╝   ╚═════╝ \n"
         "                                 "
      << endl;
}

// A sphere with constant permittivity.  The sphere has a radius, center, and
// permittivity specified on the command line and stored in ds_params_.
double eps_freespace(const Vector &x)
{
   return epsilon0_;
}

// A spherical shell with constant permeability.  The sphere has inner and outer
// radii, center, and relative permeability specified on the command line and
// stored in ms_params_.
double mu_freespace(const Vector &x)
{
   return mu0_;
}

// A sphere with constant conductivity.  The sphere has a radius, center, and
// conductivity specified on the command line and stored in ls_params_.
double sigma_freespace(const Vector &x)
{
   return 0.0;
}

// A cylindrical rod of current density.  The rod has two axis end points, a
// radius, a current amplitude in Amperes, a center time, and a width.  All of
// these parameters are stored in dp_params_.
void j_current_free(const Vector &x, double t, Vector &j) // x and t is the pde variable, j is the calculated
{

   j.SetSize(x.Size());
   j = 0.0;

//   Vector  v(x.Size());  // Normalized Axis vector
//   Vector xu(x.Size());  // x vector relative to the axis end-point
//
//   xu = x;
//
//   for (int i=0; i<x.Size(); i++)
//   {
//      xu[i] -= dp_params_[i];
//      v[i]   = dp_params_[x.Size()+i] - dp_params_[i]; // vector points from dp start to dp end
//   }
//
//   double h = v.Norml2();
//
//   if ( h == 0.0 )
//   {
//      return;
//   }
//   v /= h;
//
//   double r = dp_params_[2*x.Size()+0];
//   double a = dp_params_[2*x.Size()+1] * tScale_;
//   double b = dp_params_[2*x.Size()+2] * tScale_;
//   double c = dp_params_[2*x.Size()+3] * tScale_;
//
//   double xv = xu * v;
//
//   // Compute perpendicular vector from axis to x
//   xu.Add(-xv, v);
//
//   double xp = xu.Norml2();
//
//   if ( xv >= 0.0 && xv <= h && xp <= r )
//   {
//      j = v;
//   }
//
//   j *= a * (t - b) * exp(-0.5 * pow((t-b)/c, 2)) / (c * c); // gaussian pulse td shape
}

void
EFieldFunc(const Vector &x, Vector &E)
{
   E.SetSize(3);
   E = 0.0;
}

void
BFieldFunc(const Vector &x, Vector &B)
{
   B.SetSize(3);
   B = 0.0;
}
#define Pi                3.141592653589793239e+00
#define Vo                2.997924580000000000e+08
void
dEdt_GausssianPulse_BCFunc(const Vector &x, double t, Vector &dE)
{
   double Emagnitude = 1;
   double Frequency = 5e9;
   double to = 1.0e-09;
   double tau = 0.2e-09;

   dE.SetSize(3);
   double omega = 2.0  * Pi * Frequency;
   Vector  k(x.Size()); // k vector
   Vector  n(x.Size()); // normal vector
   double epol_[] = {0,0,1.0};
   Vector Epol = Vector(epol_, 3);

   double V_light = Vo;
   double SinModul = cos(omega*(t-to));
   double Exponent   = t - to  - (k * x / V_light);
   double IncidExcit = Emagnitude * exp(- (Exponent * Exponent) / (tau * tau));
   dE = (Epol*=(IncidExcit * SinModul));
}

int
SnapTimeStep(double tmax, double dtmax, double & dt)
{
   double dsteps = tmax/dtmax;

   int nsteps = pow(10,(int)ceil(log10(dsteps)));

   for (int i=1; i<=5; i++)
   {
      int a = (int)ceil(log10(dsteps/pow(5.0,i)));
      int nstepsi = (int)pow(5,i)*max(1,(int)pow(10,a));

      nsteps = min(nsteps,nstepsi);
   }

   dt = tmax / nsteps;

   return nsteps;
}
