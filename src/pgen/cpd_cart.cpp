//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
// Edited by JTL 6/7/22
//========================================================================================
//! \file cpd.cpp
//! \brief Adapted from pgen/disk.cpp to do a global simulation
//! centered on the planet in a rotating frame.

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
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


void PrimaryGravity(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
		    AthenaArray<Real> &cons_scalar);
//void OrbiterGravity(MeshBlock *pmb, const Real time, const Real dt,
//              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
//              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
//		    AthenaArray<Real> &cons_scalar);
namespace {
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real PoverR(const Real rad, const Real phi, const Real z);
Real VelProfileCyl(const Real rad, const Real phi, const Real z);
// problem parameters which are useful to make global to this file
Real gm1, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas, gm_planet;
Real dfloor;
Real Omega0;
  
} // namespace

// User-defined boundary conditions for disk simulations
void DiskCartInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskCartOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskCartInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskCartOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

// Should be able to modify this so that the initial density &
// velocity profile are circum-primary rather than circum-orbiter. The
// orbiter is at the center of the coordinate system.
void Mesh::InitUserMeshData(ParameterInput *pin) {

  EnrollUserExplicitSourceFunction(PrimaryGravity);
  //EnrollUserExplicitSourceFunction(OrbiterGravity);
  // Get parameters for gravitatonal potential of central point mass
  //gm1 = pin->GetOrAddReal("problem","GM_s",0.0); // primary gravitational constant
  r0 = pin->GetOrAddReal("problem","r0",1.0);

  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);

  // Get parameters of initial pressure and cooling parameters
  if (NON_BAROTROPIC_EOS) {
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
  } else {
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }
  Real float_min = std::numeric_limits<float>::min();
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(float_min)));

  Omega0 = pin->GetOrAddReal("orbital_advection","Omega0",0.0);

  // enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiskCartInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiskCartOuterX1);
  }
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, DiskCartInnerX2);
  }
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, DiskCartOuterX2);
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad, phi, z;
  Real den, vel;
  Real x1, x2, x3;
  Real Cx, Cy;
  Real radp, phip;

  //OrbitalVelocityFunc &vK = porb->OrbitalVelocity;
  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        x1 = pcoord->x1v(i);
	// calculate relative to the primary
	radp = sqrt(x2*x2+(1+x1)*(1+x1));
	phip = std::atan2(x2, 1+x1);
        // compute initial conditions in cylindrical coordinates
	// relative to primary
        den = DenProfileCyl(radp,phip,z);
        vel = VelProfileCyl(radp,phip,z);
        //if (porb->orbital_advection_defined)
        //  vel -= vK(porb, x1, x2, x3);
        phydro->u(IDN,k,j,i) = den;
	Cx = x2/radp;
	Cy = -(1+x1)/radp;
	phydro->u(IM1,k,j,i) = den*vel*Cx;
	phydro->u(IM2,k,j,i) = den*vel*Cy;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          Real p_over_r = PoverR(radp,phip,z);
          phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))
				       +SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
	}
      }
    }
  }
  return;
}

namespace {
//----------------------------------------------------------------------------------------
//! transform to cylindrical coordinate FROM CARTESIAN

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  Real x, y;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x=pco->x1v(i);
    y=pco->x2v(i);
    rad=sqrt(x*x+y*y);
    phi=atan2(y, x);
    z=pco->x3v(k);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    rad=std::abs(pco->x1v(i)*std::sin(pco->x2v(j)));
    phi=pco->x3v(i);
    z=pco->x1v(i)*std::cos(pco->x2v(j));
  }
  return;
}

//----------------------------------------------------------------------------------------
//! computes density in cylindrical coordinates

Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real den;
  Real p_over_r = p0_over_r0;
  if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);
  Real denmid = rho0*std::pow(rad/r0,dslope);
  Real dentem = denmid*std::exp(1./p_over_r*(1./std::sqrt(SQR(rad)+SQR(z))-1./rad));
  den = dentem;
  return std::max(den,dfloor);
}

//----------------------------------------------------------------------------------------
//! computes pressure/density in cylindrical coordinates

Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  poverr = p0_over_r0*std::pow(rad/r0, pslope);
  return poverr;
}

//----------------------------------------------------------------------------------------
//! computes rotational velocity in cylindrical coordinates

Real VelProfileCyl(const Real rad, const Real phi, const Real z) {
  Real p_over_r = PoverR(rad, phi, z);
  Real vel = (dslope+pslope)*p_over_r/(1./rad) + (1.0+pslope)
             - pslope*rad/std::sqrt(rad*rad+z*z);
  // Omega=1
  vel = -std::sqrt(1./rad)*std::sqrt(vel)+rad;
  return vel;
}
} // namespace

void PrimaryGravity(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
		    AthenaArray<Real> &cons_scalar){
  Real l, x1, x2, x3, Fx, Fy, den;
  Real vx, vy;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    x3 = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      x2 = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        x1 = pmb->pcoord->x1v(i);
	den = prim(IDN,k,j,i);
	vx = prim(IVX,k,j,i);
	vy = prim(IVY,k,j,i);
	// implement primary source term here
	// GM* = 1
	// R = 1
	// Omega =1
	l = sqrt(x2*x2+(1+x1)*(1+x1));
	Fx = ((1+x1)/l)*(1/l/l-l)-2*vy;
	Fy = (x2/l)*(1/l/l-l)+2*vx;
	cons(IM1,k,j,i) -= dt*den*Fx;
	cons(IM2,k,j,i) -= dt*den*Fy;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskCartInnerX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad, phi, z;
  Real vel, x1, x2, x3, radp, phip, den, Cx, Cy;
  for (int k=kl; k<=ku; ++k) {
    x3 = pco->x3v(k);
    for (int j=1; j<=ngh; ++j) {
      x2 = pco->x2v(jl-j);
      for (int i=il; i<=iu; ++i) {
	x1 = pco->x1v(i);
	// calculate relative to the primary
	radp = sqrt(x2*x2+(1+x1)*(1+x1));
	phip = std::atan2(x2, 1+x1);
	Cx = x2/radp;
	Cy = -(1+x1)/radp;
        den = DenProfileCyl(radp,phip,z);
        prim(IDN,k,jl-j,i) = den;
        vel = VelProfileCyl(radp,phip,z);
	prim(IVX,k,jl-j,i) = vel*Cx;
	prim(IVY,k,jl-j,i) = vel*Cy;
        prim(IVZ,k,jl-j,i) = 0.0;
        if (NON_BAROTROPIC_EOS)
	  prim(IEN,k,jl-j,i) = PoverR(radp, phip, z)*prim(IDN,k,jl-j,i);
      }
    }
  }
}

void DiskCartOuterX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad, phi, z;
  Real vel, x1, x2, x3, radp, phip, den, Cx, Cy;
  for (int k=kl; k<=ku; ++k) {
    x3 = pco->x3v(k);
    for (int j=1; j<=ngh; ++j) {
      x2 = pco->x2v(ju+j);
      for (int i=il; i<=iu; ++i) {
	x1 = pco->x1v(i);
	// calculate relative to the primary
	radp = sqrt(x2*x2+(1+x1)*(1+x1));
	phip = std::atan2(x2, 1+x1);
	Cx = x2/radp;
	Cy = -(1+x1)/radp;
        den = DenProfileCyl(radp,phip,z);
        prim(IDN,k,ju+j,i) = den;
        vel = VelProfileCyl(radp,phip,z);
	prim(IVX,k,ju+j,i) = vel*Cx;
	prim(IVY,k,ju+j,i) = vel*Cy;
        prim(IVZ,k,ju+j,i) = 0.0;
        if (NON_BAROTROPIC_EOS)
	  prim(IEN,k,ju+j,i) = PoverR(radp, phip, z)*prim(IDN,k,ju+j,i);
      }
    }
  }
}

void DiskCartInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad, phi, z;
  Real vel, x1, x2, x3, radp, phip, den, Cx, Cy;
  for (int k=kl; k<=ku; ++k) {
    x3 = pco->x3v(k);
    for (int j=jl; j<=ju; ++j) {
      x2 = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
	x1 = pco->x1v(il-i);
	// calculate relative to the primary
	radp = sqrt(x2*x2+(1+x1)*(1+x1));
	phip = std::atan2(x2, 1+x1);
	Cx = x2/radp;
	Cy = -(1+x1)/radp;
        den = DenProfileCyl(radp,phip,z);
        prim(IDN,k,j,il-i) = den;
        vel = VelProfileCyl(radp,phip,z);
	prim(IVX,k,j,il-i) = vel*Cx;
	prim(IVY,k,j,il-i) = vel*Cy;
        prim(IVZ,k,j,il-i) = 0.0;
        if (NON_BAROTROPIC_EOS)
	  prim(IEN,k,j,il-i) = PoverR(radp, phip, z)*prim(IDN,k,j,il-i);
      }
    }
  }
}

void DiskCartOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad, phi, z;
  Real vel, x1, x2, x3, radp, phip, den, Cx, Cy;
  for (int k=kl; k<=ku; ++k) {
    x3 = pco->x3v(k);
    for (int j=jl; j<=ju; ++j) {
      x2 = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
	x1 = pco->x1v(iu);
	// calculate relative to the primary
	radp = sqrt(x2*x2+(1+x1)*(1+x1));
	phip = std::atan2(x2, 1+x1);
	Cx = x2/radp;
	Cy = -(1+x1)/radp;
        den = DenProfileCyl(radp,phip,z);
        prim(IDN,k,j,iu+i) = den;
        vel = VelProfileCyl(radp,phip,z);
	prim(IVX,k,j,iu+i) = vel*Cx;
	prim(IVY,k,j,iu+i) = vel*Cy;
        prim(IVZ,k,j,iu+i) = 0.0;
        if (NON_BAROTROPIC_EOS)
	  prim(IEN,k,j,iu+i) = PoverR(radp, phip, z)*prim(IDN,k,j,iu+i);
      }
    }
  }
}

//void OrbiterGravity(MeshBlock *pmb, const Real time, const Real dt,
//              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
//              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
//		    AthenaArray<Real> &cons_scalar){
//  Real l, x1, x2, x3, Fx, Fy, den;
//  for (int k = pmb->ks; k <= pmb->ke; ++k) {
//    x3 = pmb->pcoord->x3v(k);
//    for (int j = pmb->js; j <= pmb->je; ++j) {
//      x2 = pmb->pcoord->x2v(j);
//      for (int i = pmb->is; i <= pmb->ie; ++i) {
//        x1 = pmb->pcoord->x1v(i);
//	den = prim(IDN,k,j,i);
//	// implement primary source term here
//	// GM* = 1
//	// R = 1
//	l = sqrt(x2*x2+x1*x1);
//	Fx = (x1/l)/l/l; //std::pow(l,3);
//	Fy = (x2/l)/l/l; ///std::pow(l,3);
//	cons(IM1,k,j,i) -= gm1*dt*den*Fx;
//	cons(IM2,k,j,i) -= gm1*dt*den*Fy;
//      }
//    }
//  }
//  return;
//}
