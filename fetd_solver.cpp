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

#include "fetd_solver.hpp"

using namespace std;

namespace mfem
{

using namespace common;

namespace electromagnetics
{

// Used for combining scalar coefficients
double prodFunc(double a, double b) { return a * b; }

MaxwellSolver::MaxwellSolver(Mesh & mesh, int order,
                             double (*eps     )(const Vector&),
                             double (*muInv   )(const Vector&),
                             double (*sigma   )(const Vector&),
                             void   (*j_src   )(const Vector&, double, Vector&),
                             Array<int> & abcs,
                             Array<int> & dbcs_pec,
                             Array<int> & dbcs_dE,
                             void   (*dEdt_bc )(const Vector&, double, Vector&))
   : order_(order),
     logging_(1),
     dtMax_(-1.0),
     dtScale_(1.0e6),
     mesh_(&mesh),
     HCurlFESpace_(nullptr),
     HDivFESpace_(nullptr),
     hDivMassMuInv_(nullptr),
     hCurlLosses_(nullptr),
     weakCurlMuInv_(nullptr),
     E_(nullptr),
     B_(nullptr),
     HD_(nullptr),
     RHS_(nullptr),
     e_(nullptr),
     b_(nullptr),
     j_(nullptr),
     dedt_(nullptr),
     rhs_(nullptr),
     jd_(nullptr),
     Curl_(nullptr),
     epsCoef_(nullptr),
     muInvCoef_(nullptr),
     sigmaCoef_(nullptr),
     etaInvCoef_(nullptr),
     eCoef_(nullptr),
     bCoef_(nullptr),
     jCoef_(nullptr),
     dEdtBCCoef_(nullptr),
     eps_(eps),
     muInv_(muInv),
     sigma_(sigma),
     j_src_(j_src),
     dEdt_bc_(dEdt_bc),
     visit_dc_(nullptr)
{

   // Define compatible finite element spaces on the
   // mesh. Here we use arbitrary order H1, Nedelec, and Raviart-Thomas finite
   // elements.
   HCurlFESpace_ = new ND_FESpace(mesh_,order_,mesh_->Dimension());
   HDivFESpace_  = new RT_FESpace(mesh_,order_,mesh_->Dimension());

   this->height = HCurlFESpace_->GetTrueVSize();
   this->width  = HDivFESpace_->GetTrueVSize();

   // Check for absorbing materials or boundaries
   lossy_ = abcs.Size() > 0 || sigma_ != nullptr;

   // Require implicit handling of loss terms
   type = lossy_ ? IMPLICIT : EXPLICIT;

   // Electric permittivity
   if ( eps_ == nullptr )
   {
      epsCoef_ = new ConstantCoefficient(epsilon0_);
   }
   else
   {
      if ( logging_ > 0 )
      {
         cout << "Creating Permittivity Coefficient" << endl;
      }
      epsCoef_ = new FunctionCoefficient(eps_);
   }

   // Inverse of the magnetic permeability
   if ( muInv_ == nullptr )
   {
      muInvCoef_ = new ConstantCoefficient(1.0/mu0_);
   }
   else
   {
      if ( logging_ > 0 )
      {
         cout << "Creating Permeability Coefficient" << endl;
      }
      muInvCoef_ = new FunctionCoefficient(muInv_);
   }

   // Electric conductivity
   if ( sigma_ != nullptr )
   {
      if ( logging_ > 0 )
      {
         cout << "Creating Conductivity Coefficient" << endl;
      }
      sigmaCoef_ = new FunctionCoefficient(sigma_);
   }

   // Impedance of free space
   if ( abcs.Size() > 0 )
   {
      if ( logging_ > 0 )
      {
         cout << "Creating Admittance Coefficient" << endl;
      }

      AttrToMarker(mesh.bdr_attributes.Max(), abcs, abc_marker_);
      etaInvCoef_ = new ConstantCoefficient(sqrt(epsilon0_/mu0_));
   }

   // Electric Field Boundary Condition
   if ( logging_ > 0 )
   {
      cout << "Configuring Dirichlet BC" << endl;
   }

   DEBUG_INFO("mesh.bdr_attributes.Max(): "+ to_string(mesh.bdr_attributes.Max()));

   AttrToMarker(mesh.bdr_attributes.Max(), dbcs_pec, dbc_marker_pec);
   AttrToMarker(mesh.bdr_attributes.Max(), dbcs_dE, dbc_marker_dE);

   HCurlFESpace_->GetEssentialTrueDofs(dbc_marker_pec, dbc_dofs_);

   if ( dEdt_bc_ != nullptr )
   {
      dEdtBCCoef_ = new VectorFunctionCoefficient(3,dEdt_bc_);
   }
   else
   {
      Vector ebc(3); ebc = 0.0;
      dEdtBCCoef_ = new VectorConstantCoefficient(ebc);
   }


   // Bilinear Forms
   if ( logging_ > 0 )
   {
      cout << "Creating H(Div) Mass Operator" << endl;
   }
   hDivMassMuInv_ = new BilinearForm(HDivFESpace_);
   hDivMassMuInv_->AddDomainIntegrator(new VectorFEMassIntegrator(*muInvCoef_));

   if ( logging_ > 0 )
   {
      cout << "Creating Weak Curl Operator" << endl;
   }
   weakCurlMuInv_ = new MixedBilinearForm(HDivFESpace_,HCurlFESpace_); 
   weakCurlMuInv_->AddDomainIntegrator(new MixedVectorWeakCurlIntegrator(*muInvCoef_));

   //
   M1Losses_ = new SparseMatrix();
   M2MuInv_ = new SparseMatrix();
   NegCurl_ = new SparseMatrix();
   WeakCurlMuInv_ = new SparseMatrix();

   // Assemble Matrices
   Array<int> ess_tdof_list;

   hDivMassMuInv_->Assemble(); 
   weakCurlMuInv_->Assemble(); 

   hDivMassMuInv_->Finalize();
   weakCurlMuInv_->Finalize();

   if ( sigmaCoef_ || etaInvCoef_ )
   {
      if ( logging_ > 0 )
      {
         cout << "Creating H(Curl) Loss Operator" << endl;
      }
      hCurlLosses_  = new BilinearForm(HCurlFESpace_);
      if ( sigmaCoef_ )
      {
         if ( logging_ > 0 )
         {
            cout << "Adding domain integrator for conductive regions" << endl;
         }
         hCurlLosses_->AddDomainIntegrator(new VectorFEMassIntegrator(*sigmaCoef_));
      }
      if ( etaInvCoef_ )
      {
         if ( logging_ > 0 )
         {
            cout << "Adding boundary integrator for absorbing boundary" << endl;
         }
         hCurlLosses_->AddBoundaryIntegrator(new VectorFEMassIntegrator(*etaInvCoef_), abc_marker_);
      }
      hCurlLosses_->Assemble();
      hCurlLosses_->Finalize();

      ess_tdof_list.DeleteAll();
      if (mesh_->bdr_attributes.Size())
      {
         Array<int> ess_bdr(mesh_->bdr_attributes.Max());
         ess_bdr = 1;
         HCurlFESpace_->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
      }
      hCurlLosses_->FormSystemMatrix(ess_tdof_list, *M1Losses_); //FormSystemMatrix
   }

   // Create Linear Algebra Matrices
   ess_tdof_list.DeleteAll();
   if (mesh_->bdr_attributes.Size())
   {
      Array<int> ess_bdr(mesh_->bdr_attributes.Max());
      ess_bdr = 1;
      HDivFESpace_->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }
   hDivMassMuInv_->FormSystemMatrix(ess_tdof_list, *M2MuInv_); 

   //weakCurlMuInv_->FormSystemMatrix(ess_tdof_list, *WeakCurlMuInv_);

   if ( logging_ > 0 )
   {
      cout << "Creating discrete curl operator" << endl;
   }
   Curl_ = new DiscreteCurlOperator(HCurlFESpace_, HDivFESpace_);
   Curl_->Assemble();
   Curl_->Finalize();

   ess_tdof_list.DeleteAll();
   if (mesh_->bdr_attributes.Size())
   {
      Array<int> ess_bdr(mesh_->bdr_attributes.Max());
      ess_bdr = 1;
      HCurlFESpace_->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // Beware this modifies the matrix stored within the Curl_ object.
   *NegCurl_ = (Curl_->SpMat() *= (-1.0));

   // Build grid functions
   e_    = new GridFunction(HCurlFESpace_);  
   dedt_ = new GridFunction(HCurlFESpace_);  
   rhs_  = new GridFunction(HCurlFESpace_); 
   b_    = new GridFunction(HDivFESpace_); 

   E_ = new Vector();
   B_ = new Vector();
   HD_  = new Vector(HDivFESpace_->GetTrueVSize()); 
   RHS_ = new Vector(HCurlFESpace_->GetTrueVSize()); 

   // Initialize dedt to zero
   *dedt_ = 0.0;

   if ( j_src_)
   {
      if ( logging_ > 0 )
      {
         cout << "Creating Current Source" << endl;
      }
      jCoef_ = new VectorFunctionCoefficient(3,j_src_);

      j_  = new GridFunction(HCurlFESpace_);
      j_->ProjectCoefficient(*jCoef_);

      jd_ = new LinearForm(HCurlFESpace_);
      jd_->AddDomainIntegrator(new VectorFEDomainLFIntegrator(*jCoef_));
      jd_->Assemble();
   }
   dtMax_ = GetMaximumTimeStep();
}

MaxwellSolver::~MaxwellSolver()
{
   delete epsCoef_;
   delete muInvCoef_;
   delete etaInvCoef_;
   delete jCoef_;
   delete dEdtBCCoef_;

   delete E_;
   delete B_;
   delete HD_;
   delete RHS_;

   delete e_;
   delete b_;
   delete j_;
   delete dedt_;
   delete rhs_;
   delete jd_;

   delete Curl_;

   delete M1Losses_;
   delete M2MuInv_;
   delete NegCurl_;
   delete WeakCurlMuInv_;

   delete hDivMassMuInv_;
   delete hCurlLosses_;
   delete weakCurlMuInv_;

   delete HCurlFESpace_;
   delete HDivFESpace_;

   map<int, BilinearForm*>::iterator mit1;
   for (mit1=a1_.begin(); mit1!=a1_.end(); mit1++)
   {
      int i = mit1->first;
      delete pcg_[i];
      delete diagScale_[i];
      delete A1_[i];
      delete a1_[i];
   }

   map<int, Coefficient*>::iterator mit2;
   for (mit2=dtCoef_.begin(); mit2!=dtCoef_.end(); mit2++)
   {
      delete mit2->second;
   }
   for (mit2=dtSigmaCoef_.begin(); mit2!=dtSigmaCoef_.end(); mit2++)
   {
      delete mit2->second;
   }
   for (mit2=dtEtaInvCoef_.begin(); mit2!=dtEtaInvCoef_.end(); mit2++)
   {
      delete mit2->second;
   }

   map<string, socketstream*>::iterator mit3;
   for (mit3=socks_.begin(); mit3!=socks_.end(); mit3++)
   {
      delete mit3->second;
   }
}

int
MaxwellSolver::GetProblemSize()
{
   return HCurlFESpace_->GetTrueVSize();
}

void
MaxwellSolver::PrintSizes()
{
   int size_nd = HCurlFESpace_->GetTrueVSize();
   int size_rt = HDivFESpace_->GetTrueVSize();

  cout << "Number of H(Curl) unknowns: " << size_nd << endl;
  cout << "Number of H(Div)  unknowns: " << size_rt << endl << flush;

}

void
MaxwellSolver::SetInitialEField(VectorCoefficient & EFieldCoef)
{
   eCoef_ = &EFieldCoef;
   e_->ProjectCoefficient(EFieldCoef);
   e_->GetTrueDofs(*E_); 
   e_->SetFromTrueDofs(*E_); 
}

void
MaxwellSolver::SetInitialBField(VectorCoefficient & BFieldCoef)
{
   bCoef_ = &BFieldCoef;
   b_->ProjectCoefficient(BFieldCoef);
   b_->GetTrueDofs(*B_); 
   b_->SetFromTrueDofs(*B_); 
}

void
MaxwellSolver::Mult(const Vector &B, Vector &dEdt) const
{
   implicitSolve(0.0, B, dEdt);
}

void
MaxwellSolver::ImplicitSolve(double dt, const Vector &B, Vector &dEdt)
{
   implicitSolve(dt, B, dEdt);
}

void
MaxwellSolver::setupSolver(const int idt, const double dt) const
{
   if ( pcg_.find(idt) == pcg_.end() )
   {
      if ( logging_ > 0 )
      {
         cout << "Creating implicit operator for dt = " << dt << endl;
      }

      a1_[idt] = new BilinearForm(HCurlFESpace_);
      a1_[idt]->AddDomainIntegrator(new VectorFEMassIntegrator(epsCoef_));

      if ( idt != 0 )
      {
         dtCoef_[idt] = new ConstantCoefficient(0.5 * dt);
         if ( sigmaCoef_ )
         {
            dtSigmaCoef_[idt] = new TransformedCoefficient(dtCoef_[idt], sigmaCoef_, prodFunc);
            a1_[idt]->AddDomainIntegrator(new VectorFEMassIntegrator(dtSigmaCoef_[idt]));
         }
         if ( etaInvCoef_ )
         {
            dtEtaInvCoef_[idt] = new TransformedCoefficient(dtCoef_[idt], etaInvCoef_, prodFunc);
            a1_[idt]->AddBoundaryIntegrator(new VectorFEMassIntegrator(dtEtaInvCoef_[idt]), const_cast<Array<int>&>(abc_marker_));
         }
      }

      a1_[idt]->Assemble();
      a1_[idt]->Finalize();

      Array<int> ess_tdof_list;
      if (mesh_->bdr_attributes.Size())
      {
         Array<int> ess_bdr(mesh_->bdr_attributes.Max());
         ess_bdr = 1;
         HCurlFESpace_->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
      }
      A1_[idt] = new SparseMatrix();
      a1_[idt]->FormSystemMatrix(ess_tdof_list, *A1_[idt]);

      // Use a simple symmetric Gauss-Seidel preconditioner with PCG.
      diagScale_[idt] = new GSSmoother(*A1_[idt]);
      pcg_[idt] = new CGSolver();
      pcg_[idt]->SetRelTol(sqrt(1.0e-12));
      pcg_[idt]->SetAbsTol(sqrt(1.0e-12));
      pcg_[idt]->SetMaxIter(200);
      pcg_[idt]->SetPrintLevel(0);
      pcg_[idt]->SetOperator(*A1_[idt]);
      pcg_[idt]->SetPreconditioner(*diagScale_[idt]);
   }
}

void MaxwellSolver::implicitSolve(double dt, const Vector &B, Vector &dEdt) const
{
   int idt = hCurlLosses_ ? ((int)(dtScale_ * dt / dtMax_)) : 0;

   b_->SetFromTrueDofs(B);
   weakCurlMuInv_->Mult(*b_, *rhs_);

   if ( hCurlLosses_ )
   {
      e_->SetFromTrueDofs(*E_);
      hCurlLosses_->AddMult(*e_, *rhs_, -1.0);
   }

   if ( jd_ )
   {
      jCoef_->SetTime(t); // 't' is member data from mfem::TimeDependentOperator
      jd_->Assemble();
      *rhs_ -= *jd_;
   }

   if ( dEdtBCCoef_ )
   {
      dEdtBCCoef_->SetTime(t);
      dedt_->ProjectBdrCoefficientTangent(*dEdtBCCoef_,
                                          const_cast<Array<int>&>(dbc_marker_dE));
   }

   // Create objects and matrices for solving with the given time step
   setupSolver(idt, dt);

   // Apply essential BCs and determine true DoFs for the right hand side
   a1_[idt]->FormLinearSystem(dbc_dofs_, *dedt_, *rhs_, *A1_[idt], dEdt, *RHS_);

   // Solve for the time derivative of the electric field (true DoFs)
   pcg_[idt]->Mult(*RHS_, dEdt);

   // Distribute shared DoFs to relevant processors
   a1_[idt]->RecoverFEMSolution(dEdt, *rhs_, *dedt_);
}

void
MaxwellSolver::SyncGridFuncs()
{
   e_->SetFromTrueDofs(*E_);
   b_->SetFromTrueDofs(*B_);
}

double
MaxwellSolver::GetMaximumTimeStep() const
{
   if ( dtMax_ > 0.0 )
   {
      return dtMax_;
   }

   Vector * v0 = new Vector(HCurlFESpace_->GetTrueVSize());
   Vector * v1 = new Vector(HCurlFESpace_->GetTrueVSize());
   Vector * u0 = new Vector(HDivFESpace_->GetTrueVSize());

   v0->Randomize(1234);

   int iter = 0, nstep = 20;
   double dt0 = 1.0, dt1 = 1.0, change = 1.0, ptol = 0.001;

   // Create Solver assuming no loss operators
   setupSolver(0, 0.0);

   // Use power method to approximate the largest eigenvalue of the updatev operator.
   while ( iter < nstep && change > ptol )
   {
      double normV0 = InnerProduct(*v0,*v0);
      *v0 /= sqrt(normV0);

      NegCurl_->Mult(*v0,*u0);
      M2MuInv_->Mult(*u0,*HD_);
      NegCurl_->MultTranspose(*HD_,*RHS_);

      pcg_[0]->Mult(*RHS_,*v1);

      double lambda = InnerProduct(*v0,*v1);
      dt1 = 2.0/sqrt(lambda);
      change = fabs((dt1-dt0)/dt0);
      dt0 = dt1;

      if ( logging_ > 1 )
      {
         cout << iter << ":  " << dt0 << " " << change << endl;
      }

      std::swap(v0, v1);

      iter++;
   }

   return dt0;
}

double
MaxwellSolver::GetEnergy() const
{
   double energy = 0.0;

   A1_[0]->Mult(*E_,*RHS_);
   M2MuInv_->Mult(*B_,*HD_);

   energy = InnerProduct(*E_,*RHS_) + InnerProduct(*B_,*HD_);

   return 0.5 * energy;
}

void
MaxwellSolver::RegisterVisItFields(VisItDataCollection & visit_dc)
{
   visit_dc_ = &visit_dc;

   visit_dc.RegisterField("E", e_);
   //visit_dc.RegisterField("B", b_);
   if ( j_ )
   {
      //visit_dc.RegisterField("J", j_);
   }
}

void
MaxwellSolver::WriteVisItFields(int it)
{
   if ( visit_dc_ )
   {
      if ( logging_ > 1 )
      { cout << "Writing VisIt files ..." << flush; }

      if ( j_ )
      {
         jCoef_->SetTime(t);
         j_->ProjectCoefficient(*jCoef_);
      }

      visit_dc_->SetCycle(it);
      visit_dc_->SetTime(t);
      visit_dc_->Save();

      if ( logging_ > 1 ) { cout << " " << endl << flush; }
   }
}

void
MaxwellSolver::InitializeGLVis()
{
   if ( logging_ > 0 )
   { cout << "Opening GLVis sockets." << endl << flush; }

   socks_["E"] = new socketstream;
   socks_["E"]->precision(8);

   socks_["B"] = new socketstream;
   socks_["B"]->precision(8);

   if ( j_ )
   {
      socks_["J"] = new socketstream;
      socks_["J"]->precision(8);
   }

   if ( logging_ > 0 )
   { cout << "GLVis sockets open." << endl << flush; }
}

void
MaxwellSolver::DisplayToGLVis()
{
   if ( logging_ > 1 )
   { cout << "Sending data to GLVis ..." << flush; }

   char vishost[] = "localhost";
   int  visport   = 19916;

   int Wx = 0, Wy = 0; // window position
   int Ww = 350, Wh = 350; // window size
   int offx = Ww+10, offy = Wh+45; // window offsets

   VisualizeField(*socks_["E"], vishost, visport,
                  *e_, "Electric Field (E)", Wx, Wy, Ww, Wh);
   Wx += offx;

   VisualizeField(*socks_["B"], vishost, visport,
                  *b_, "Magnetic Flux Density (B)", Wx, Wy, Ww, Wh);

   if ( j_ )
   {
      Wx = 0;
      Wy += offy;

      jCoef_->SetTime(t); // Is member data from mfem::TimeDependentOperator
      j_->ProjectCoefficient(*jCoef_);

      VisualizeField(*socks_["J"], vishost, visport,
                     *j_, "Current Density (J)", Wx, Wy, Ww, Wh);
   }
   if ( logging_ > 1 ) { cout << " " << flush; }
}

} // namespace FETD

} // namespace mfem

