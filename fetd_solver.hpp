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

#ifndef MFEM_MAXWELL_SOLVER
#define MFEM_MAXWELL_SOLVER

#include "fem_extras.hpp"
#include "mesh_extras.hpp"
#include "electromagnetics.hpp"
#include "debug.hpp"
#include <string>
#include <map>

using namespace std;
using namespace mfem;

namespace mfem
{

using common::ND_FESpace;
using common::RT_FESpace;
using common::DiscreteCurlOperator;

namespace electromagnetics
{

class MaxwellSolver : public TimeDependentOperator
{
public:
   MaxwellSolver(Mesh & mesh, int sOrder,
                 double (*eps     )(const Vector&),
                 double (*muInv   )(const Vector&),
                 double (*sigma   )(const Vector&),
                 void   (*j_src   )(const Vector&, double, Vector&),
                 Array<int> & abcs,
                 Array<int> & dbcs_pec,
                 Array<int> & dbcs_dE,
                 void   (*dEdt_bc )(const Vector&, double, Vector&));

   ~MaxwellSolver();

   int GetLogging() const { return logging_; }
   void SetLogging(int logging) { logging_ = logging; }

   int GetProblemSize();

   void PrintSizes();

   void SetInitialEField(VectorCoefficient & EFieldCoef);
   void SetInitialBField(VectorCoefficient & BFieldCoef);

   void Mult(const Vector &B, Vector &dEdt) const;

   void ImplicitSolve(const double dt, const Vector &x, Vector &k);

   double GetMaximumTimeStep() const;

   double GetEnergy() const;

   Operator & GetNegCurl() { return *NegCurl_; }

   Vector & GetEField() { return *E_; }
   Vector & GetBField() { return *B_; }

   void SyncGridFuncs();

   void RegisterVisItFields(VisItDataCollection & visit_dc);

   void WriteVisItFields(int it = 0);

   void InitializeGLVis();

   void DisplayToGLVis();

private:

   // This method alters mutable member data
   void setupSolver(const int idt, const double dt) const;

   void implicitSolve(const double dt, const Vector &x, Vector &k) const;

   int order_;
   int logging_;

   bool lossy_;

   double dtMax_;   // Maximum stable time step
   double dtScale_; // Used to scale dt before converting to an integer

   Mesh * mesh_;

   ND_FESpace * HCurlFESpace_;
   RT_FESpace * HDivFESpace_;

   BilinearForm * hDivMassMuInv_;
   BilinearForm * hCurlLosses_;
   MixedBilinearForm * weakCurlMuInv_;
   DiscreteCurlOperator * Curl_;

   GridFunction * e_;    // Electric Field (HCurl)
   GridFunction * b_;    // Magnetic Flux (HDiv)
   GridFunction * j_;    // Volumetric Current Density (HCurl)
   GridFunction * dedt_; // Time Derivative of Electric Field (HCurl)
   GridFunction * rhs_;  // Dual of displacement current, rhs vector (HCurl)
   LinearForm   * jd_;   // Dual of current density (HCurl)

   SparseMatrix * M1Losses_;
   SparseMatrix * M2MuInv_;
   SparseMatrix * NegCurl_;
   SparseMatrix * WeakCurlMuInv_;
   Vector * E_; // Current value of the electric field DoFs
   Vector * B_; // Current value of the magnetic flux DoFs
   mutable Vector * HD_; // Used in energy calculation
   mutable Vector * RHS_;

   Coefficient       * epsCoef_;    // Electric Permittivity Coefficient
   Coefficient       * muInvCoef_;  // Magnetic Permeability Coefficient
   Coefficient       * sigmaCoef_;  // Electric Conductivity Coefficient
   Coefficient       * etaInvCoef_; // Admittance Coefficient
   VectorCoefficient * eCoef_;      // Initial Electric Field
   VectorCoefficient * bCoef_;      // Initial Magnetic Flux
   VectorCoefficient * jCoef_;      // Time dependent current density
   VectorCoefficient * dEdtBCCoef_; // Time dependent boundary condition

   double (*eps_    )(const Vector&);
   double (*muInv_  )(const Vector&);
   double (*sigma_  )(const Vector&);
   void   (*j_src_  )(const Vector&, double, Vector&);

   // Array of 0's and 1's marking the location of absorbing surfaces
   Array<int> abc_marker_;

   // Array of 0's and 1's marking the location of Dirichlet boundaries
   Array<int> dbc_marker_pec; // pec
   Array<int> dbc_marker_dE; // dE
   void   (*dEdt_bc_)(const Vector&, double, Vector&);

   // Dirichlet degrees of freedom
   Array<int>   dbc_dofs_;

   // High order symplectic integration requires partial time steps of differing
   // lengths. If losses are present the system matrix includes a portion scaled
   // by the time step. Consequently, high order time integration requires
   // different system matrices. The following maps contain various objects that
   // depend on the time step.
   mutable std::map<int, BilinearForm    *> a1_;
   mutable std::map<int, SparseMatrix    *> A1_;
   mutable std::map<int, Coefficient     *> dtCoef_;
   mutable std::map<int, Coefficient     *> dtSigmaCoef_;
   mutable std::map<int, Coefficient     *> dtEtaInvCoef_;
   mutable std::map<int, SparseSmoother  *> diagScale_; // solver preconditioner GSS
   mutable std::map<int, IterativeSolver *> pcg_; // CGSolver

   // Data collection used to write VisIt files
   VisItDataCollection * visit_dc_;

   // Sockets used to communicate with GLVis
   std::map<std::string, socketstream*> socks_;
};

} // namespace FETD

} // namespace mfem




#endif // MFEM_MAXWELL_SOLVER
