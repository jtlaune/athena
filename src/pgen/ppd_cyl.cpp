//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
// Created JTL 10/19/22
//========================================================================================
// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
#include <math.h>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"


namespace {
  Real gm1, Sig0, dslope, dfloor, R0, CS02;
  Real T_damp_in, T_damp_bdy, WDL1, WDL2, innerbdy, x1min;
}

Real DenProf(const Real rad);
Real VelProf(const Real rad);
// User source function
void DiskSourceFunction(MeshBlock *pmb, const Real time, const Real dt,
		  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
		  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
		  AthenaArray<Real> &cons_scalar);
// User-defined boundary conditions for disk simulations
void DiskCartInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
	       Real time, Real dt,
	       int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskCartOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
	       Real time, Real dt,
	       int il, int iu, int jl, int ju, int kl, int ku, int ngh);
Real splineKernel(Real x);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  EnrollUserExplicitSourceFunction(DiskSourceFunction);

  dfloor = pin->GetReal("hydro", "Sigma_floor");

  gm1 = pin->GetReal("problem", "GM_s");
  Sig0 = pin->GetReal("problem", "Sigma_0");
  dslope = pin->GetReal("problem", "delta");
  R0 = pin->GetReal("problem", "R0");
  CS02 = SQR(pin->GetReal("hydro", "iso_sound_speed"));

  WDL1 = pin->GetReal("problem", "WaveDampingLength_in");
  WDL2 = pin->GetReal("problem", "WaveDampingLength_out");
  innerbdy = pin->GetReal("problem", "innerbdy");
  x1min = pin->GetReal("mesh", "x1min");
  T_damp_bdy = pin->GetReal("problem", "T_damp_bdy");
  T_damp_in = pin->GetReal("problem", "T_damp_in");

  // enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiskCartInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiskCartOuterX1);
  }
  return;
}

Real DenProf(const Real rad) {
  return(std::max(Sig0*std::pow(rad/R0, dslope), dfloor));
}

Real VelProf(const Real rad) {
  return std::sqrt(dslope*CS02 + 1/rad);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real x1, x2, x3;
  Real Sig, vk, Cx, Cy;
  Real rstar, rstar_soft;
  for (int k=ks; k<=ke; ++k) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
	x1 = pcoord->x1v(i);

	rstar = x1;
	Sig = DenProf(rstar);
	vk = VelProf(rstar);

	phydro->u(IDN,k,j,i) = Sig;
	phydro->u(IM1,k,j,i) = 0.;
	phydro->u(IM2,k,j,i) = Sig*vk;
	phydro->u(IM3,k,j,i) = 0.;

      }
    }
  }
}

Real splineKernel(Real x) {
  Real W;
  if (x<0) {
    W = 1.;
  } else if (x <= 0.5) {
    W = 1. - 6.*x*x + 6.*x*x*x;
  } else if (x<=1) {
    W = 2.*(1-x)*(1-x)*(1-x);
  } else  {
    W = 0.;
  }
  return(W);
}

void DiskSourceFunction(MeshBlock *pmb, const Real time, const Real dt,
			const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
			const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
			AthenaArray<Real> &cons_scalar) {
  Real vk, Sig, Sig0, rstar, rstar_soft, x1, x2, x3;
  Real Fpr;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    x3 = pmb->pcoord->x3v(k);
    for (int j=pmb->js; j<=pmb->je; ++j) {
      x2 = pmb->pcoord->x2v(j);
      for (int i=pmb->is; i<=pmb->ie; ++i) {
	x1 = pmb->pcoord->x1v(i);

	rstar = x1;
	vk = VelProf(rstar);

	Sig0 = DenProf(rstar);
	Sig = prim(IDN,k,j,i);

	// force calculation
	Fpr = -gm1/rstar/rstar;
        cons(IM1, k, j, i) += dt * Sig * Fpr;
        cons(IM2, k, j, i) += 0;

	// wave damping regions
        if (rstar <= innerbdy) {
	  Real x = (rstar-x1min)/(innerbdy-x1min);
	  Real factor = splineKernel(x);
	  cons(IDN,k,j,i) += -factor*dt*(cons(IDN,k,j,i)-Sig0)/T_damp_in;
	  cons(IM1,k,j,i) += -factor*dt*(cons(IM1,k,j,i))/T_damp_in;
	  cons(IM2,k,j,i) += -factor*dt*(cons(IM2,k,j,i)-Sig0*vk)/T_damp_in;
	}
	if (rstar >= WDL1) {
	  Real x = 1-(rstar-WDL1)/(WDL2-WDL1);
	  Real factor = splineKernel(x);
	  cons(IDN,k,j,i) += -factor*dt*(cons(IDN,k,j,i)-Sig0)/T_damp_bdy;
	  cons(IM1,k,j,i) += -factor*dt*(cons(IM1,k,j,i))/T_damp_bdy;
	  cons(IM2,k,j,i) += -factor*dt*(cons(IM2,k,j,i)-Sig0*vk)/T_damp_bdy;
	}
      }
    }
  }
}

void DiskCartInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
		     Real time, Real dt,
		     int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real Sig, vk, rstar, rstar_soft, x1, x2, x3, Cx, Cy;
  for (int k=kl; k<=ku; ++k) {
    x3 = pco->x3v(k);
    for (int j=jl; j<=ju; ++j) {
      x2 = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
	x1 = pco->x1v(il-i);

	rstar = x1;
	Sig = DenProf(rstar);
	vk = VelProf(rstar);

	prim(IDN,k,j,il-i) = Sig;
	prim(IM1,k,j,il-i) = 0.;
	prim(IM2,k,j,il-i) = vk;
	prim(IM3,k,j,il-i) = 0.;
      }
    }
  }
}

void DiskCartOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
		     Real time, Real dt,
		     int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real Sig, vk, rstar, rstar_soft, x1, x2, x3, Cx, Cy;
  for (int k=kl; k<=ku; ++k) {
    x3 = pco->x3v(k);
    for (int j=jl; j<=ju; ++j) {
      x2 = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
	x1 = pco->x1v(iu+i);

	rstar = x1;
	Sig = DenProf(rstar);
	vk = VelProf(rstar);

	prim(IDN,k,j,iu+i) = Sig;
	prim(IM1,k,j,iu+i) = 0.;
	prim(IM2,k,j,iu+i) = vk;
	prim(IM3,k,j,iu+i) = 0.;
      }
    }
  }
}
