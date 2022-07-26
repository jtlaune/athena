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
void CentrifugalForce(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
		      AthenaArray<Real> &cons_scalar);
namespace {
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real PoverR(const Real rad, const Real phi, const Real z);
Real VelProfileCyl(const Real rad, const Real phi, const Real z);
// problem parameters which are useful to make global to this file
Real gm0, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas, gm_planet;
Real dfloor;
Real Omega0;
  
} // namespace

// User-defined boundary conditions for disk simulations
void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);

// User-defined boundary conditions for disk simulations
void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
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
  EnrollUserExplicitSourceFunction(CentrifugalForce);
  // Get parameters for gravitatonal potential of central point mass
  gm0 = pin->GetOrAddReal("problem","GM_s",0.0); // primary gravitational constant
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
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiskInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiskOuterX1);
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real den, vel;
  Real x1, x2, x3;
  Real x, y, Cx, Cy;
  Real radp, phip;

  OrbitalVelocityFunc &vK = porb->OrbitalVelocity;
  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        x1 = pcoord->x1v(i);
        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
	x = rad*std::cos(phi);
	y = rad*std::sin(phi);

	// calculate relative to the primary
	radp = std::sqrt(std::pow(rad, 2) + 1 - 2*rad*std::cos(phi));
	phip = 3.14159265358979323846 - std::atan2(y, 1-x);

        // compute initial conditions in cylindrical coordinates
	// relative to primary
        den = DenProfileCyl(radp,phip,z);
        vel = VelProfileCyl(radp,phip,z);

        //if (porb->orbital_advection_defined)
        //  vel -= vK(porb, x1, x2, x3);
        phydro->u(IDN,k,j,i) = den;
        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {

	  Cx = y/radp;
	  Cy = -(1-x)/radp;

	  phydro->u(IM1,k,j,i) = den*vel*(std::sin(phi)*Cy + std::cos(phi)*Cx);
	  phydro->u(IM2,k,j,i) = den*vel*(std::cos(phi)*Cy - std::sin(phi)*Cx);
	  //phydro->u(IM1,k,j,i) = den*vel*Cx;
	  //phydro->u(IM2,k,j,i) = den*vel*Cy;
          //phydro->u(IM1,k,j,i) = 0.0;
	  //phydro->u(IM2,k,j,i) = den*vel;
          phydro->u(IM3,k,j,i) = 0.0;
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = den*vel;
        }

        if (NON_BAROTROPIC_EOS) {
          Real p_over_r = PoverR(rad,phi,z);
          phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  return;
}

namespace {
//----------------------------------------------------------------------------------------
//! transform to cylindrical coordinate

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
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
  vel = std::sqrt(1./rad)*std::sqrt(vel) - rad*Omega0;
  return vel;
}
} // namespace

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel, x, y, radp, phip, Cx, Cy, den;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,il-i,j,k);
	  x = rad*std::cos(phi);
	  y = rad*std::sin(phi);

	  // calculate relative to the primary
	  radp = std::sqrt(std::pow(rad, 2) + 1 - 2*rad*std::cos(phi));
	  phip = 3.14159265358979323846 - std::atan2(y, 1-x);

	  Cx = y/radp;
	  Cy = -(1-x)/radp;

          den = DenProfileCyl(radp,phip,z);
          prim(IDN,k,j,il-i) = den;
          vel = VelProfileCyl(radp,phip,z);

          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(il-i), pco->x2v(j), pco->x3v(k));
	  if (prim(IM1,k,j,il) > 0){prim(IM1,k,j,il-i) = 0.0;}
	  else {prim(IM1,k,j,il-i) = prim(IM1,k,j,il);}
          prim(IM2,k,j,il-i) = vel*(std::cos(phi)*Cy - std::sin(phi)*Cx);
          prim(IM3,k,j,il-i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,il-i,j,k);
          prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(il-i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,il-i) = 0.0;
          prim(IM2,k,j,il-i) = 0.0;
          prim(IM3,k,j,il-i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel, x, y, radp, phip, den, Cx, Cy;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
	  x = rad*std::cos(phi);
	  y = rad*std::sin(phi);

	  // calculate relative to the primary
	  radp = std::sqrt(std::pow(rad, 2) + 1 - 2*rad*std::cos(phi));
	  phip = 3.14159265358979323846 - std::atan2(y, 1-x);

	  Cx = y/radp;
	  Cy = -(1-x)/radp;

          den = DenProfileCyl(radp,phip,z);
          prim(IDN,k,j,iu+i) = den;
          vel = VelProfileCyl(radp,phip,z);
          //if (pmb->porb->orbital_advection_defined)
          //  vel -= vK(pmb->porb, pco->x1v(iu+i), pco->x2v(j), pco->x3v(k));
	  prim(IM1,k,j,iu+i) = vel*(std::sin(phi)*Cy + std::cos(phi)*Cx);
	  prim(IM2,k,j,iu+i) = vel*(std::cos(phi)*Cy - std::sin(phi)*Cx);
	  //prim(IM1,k,j,iu+i) = vel*Cx;
	  //prim(IM2,k,j,iu+i) = vel*Cy;
          prim(IM3,k,j,iu+i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
          prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(iu+i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,iu+i) = 0.0;
          prim(IM2,k,j,iu+i) = 0.0;
          prim(IM3,k,j,iu+i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
        }
      }
    }
  }
}

// Do this in cylindrical geometry, be sure to check for it. Check out
// pointmass.cpp:33 to see how they implement the source term for a
// point mass at the origin. Should be straightforward to modify,
// however in that example they only have to modify cons(IM1,k,j,i),
// i.e. the momentum along r.
void PrimaryGravity(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
		    AthenaArray<Real> &cons_scalar){
  Real l, x1, x2, x3, Fx, Fy, den;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    x3 = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      x2 = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
	den = prim(IDN,k,j,i);
        x1 = pmb->pcoord->x1v(i);
	// implement primary source term here
	// GM* = 1
	// R = 1
	l = std::sqrt(std::pow(x1, 2) + 1 - 2*x1*std::cos(x2));
	Fx = ((1-x1*std::cos(x2))/l)/l/l; //std::pow(l,3);
	Fy = (x1*std::sin(x2)/l)/l/l; ///std::pow(l,3);
	cons(IM1,k,j,i) -= dt*den*(std::sin(x2)*Fy + std::cos(x2)*Fx);
	cons(IM2,k,j,i) -= dt*den*(std::cos(x2)*Fy - std::sin(x2)*Fx);
      }
    }
  }
  return;
}

void CentrifugalForce(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
		    AthenaArray<Real> &cons_scalar){
  Real l, x1, x2, x3, Cx, Cy, den;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    x3 = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      x2 = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
	den = prim(IDN,k,j,i);
        x1 = pmb->pcoord->x1v(i);
	// implement primary source term here
	// GM* = 1
	// R = 1
	l = std::sqrt(std::pow(x1, 2) + 1 - 2*x1*std::cos(x2));
	Cx = 0.;//x1*(1-x1*std::cos(x2))/l;
	Cy = 0.;//x1*x1*std::sin(x2)/l;
	cons(IM1,k,j,i) += dt*den*(std::sin(x2)*Cy + std::cos(x2)*Cx);
	cons(IM2,k,j,i) += dt*den*(std::cos(x2)*Cy - std::sin(x2)*Cx);
      }
    }
  }
  return;
}
