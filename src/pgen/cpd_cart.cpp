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

Real GravForceCalc(MeshBlock *pmb, int iout);

void KeplerianWithRotationAndWaveDamping(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
		    AthenaArray<Real> &cons_scalar);

namespace {
  void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
  Real dDRdpdr(const Real rad);
  Real dDendr(const Real rad);
  Real DenProfileCyl(const Real rad, const Real phi, const Real z);
  Real PoverR(const Real rad, const Real phi, const Real z);
  Real VelProfileCyl(const Real rad, const Real phi, const Real z);
  // problem parameters which are useful to make global to this file
  Real gm1, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas, soft_length;
  Real dfloor, denl, wg, wt, OmegaP,x1min ,x1max ,x2min ,x2max ,WDL ,T_damp_prim, T_damp_bdy, innerbdy;
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

  EnrollUserExplicitSourceFunction(KeplerianWithRotationAndWaveDamping);
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, GravForceCalc, "Fgrav_x");
  EnrollUserHistoryOutput(1, GravForceCalc, "Fgrav_y");

  // Get parameters 
  gm1 = pin->GetReal("problem","GM_s"); 
  r0 = pin->GetReal("problem","r0");
  soft_length = pin->GetReal("problem","soft_length");
  denl = pin->GetReal("problem","den_gap");
  wg = pin->GetReal("problem","wg");
  wt = pin->GetReal("problem","wt");
  rho0 = pin->GetReal("problem","rho0");
  //dslope = pin->GetOrAddReal("problem","dslope",0.0);
  OmegaP = pin->GetReal("problem","Omega_P");
  innerbdy = pin->GetReal("problem", "innerbdy");
  dslope = pin->GetReal("problem", "dslope");

  // for wave damping in the source function
  x1min = pin->GetReal("mesh", "x1min");
  x1max = pin->GetReal("mesh", "x1max");
  x2min = pin->GetReal("mesh", "x2min");
  x2max = pin->GetReal("mesh", "x2max");
  WDL = pin->GetReal("mesh", "WaveDampingLength");
  T_damp_bdy = pin->GetReal("mesh", "T_damp_bdy");
  T_damp_prim = pin->GetReal("mesh", "T_damp_prim");

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
  Real radp, phip, r;

  //OrbitalVelocityFunc &vK = porb->OrbitalVelocity;
  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        x1 = pcoord->x1v(i);
	// calculate relative to the primary
	radp = sqrt(x2*x2+(r0+x1)*(r0+x1));
	phip = std::atan2(x2, r0+x1);
        // compute initial conditions in cylindrical coordinates
	// relative to primary
        den = DenProfileCyl(radp,phip,x3);
        vel = VelProfileCyl(radp,phip,x3);
        //if (porb->orbital_advection_defined)
        //  vel -= vK(porb, x1, x2, x3);
        phydro->u(IDN,k,j,i) = den;
	Cx = -x2/radp;
	Cy = (r0+x1)/radp;
	phydro->u(IM1,k,j,i) = den*vel*Cx;
	phydro->u(IM2,k,j,i) = den*vel*Cy;
        phydro->u(IM3,k,j,i) = 0.0;

	r = std::sqrt(x1*x1+x2*x2+x3*x3);
	
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
  //Real den;
  Real denmid;
  if (rad > innerbdy) {
    denmid = rho0*std::pow(rad/r0,dslope);
  } else {
    denmid = rho0*std::pow(innerbdy/r0,dslope);
  }
  //Real DRdp = std::abs(rad-r0);
  //Real p_over_r = p0_over_r0;
  //if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);

  //if (DRdp >= wg) {
  //  denmid = denl + (rho0-denl)/2*(2-exp((wg-DRdp)/wt));
  //}
  //else {
  //  denmid = denl + (rho0-denl)/2*exp((DRdp-wg)/wt);
  //}
  //Real dentem = denmid;
  //Real dentem = denmid;
  //den = dentem;
  return std::max(denmid,dfloor);
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
  Real dDRdpdr(const Real rad) {
    if (rad < r0) {
      return -1;
    } else {
      return 1;
    }
  }

  Real dDendr(const Real rad) {
    if (rad > innerbdy) {
      return(dslope*rho0*std::pow(rad/r0, dslope-1));
    } else {
      return(0.);
    }
    //Real DRdp = std::abs(rad-r0);
    //if (DRdp >= wg) {
    //  return 0.5*(rho0-denl)/wt*exp((-DRdp+wg)/wt)*dDRdpdr(rad);
    //} else {
    //  return 0.5*(rho0-denl)/wt*exp((DRdp-wg)/wt)*dDRdpdr(rad);
    //}
  }

Real VelProfileCyl(const Real rad, const Real phi, const Real z) {
  Real p_over_r = PoverR(rad, phi, z);
  Real den = DenProfileCyl(rad, phi, z);
  Real vel = 1/rad + rad/den*p_over_r*dDendr(rad);
  vel = std::sqrt(vel)-OmegaP*rad;
  return vel;
}
} // namespace

void KeplerianWithRotationAndWaveDamping(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
		    AthenaArray<Real> &cons_scalar){
  Real l, x1, x2, x3, Fx, Fy, den;
  Real vx, vy, rs, phip, vk;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    x3 = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      x2 = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        x1 = pmb->pcoord->x1v(i);
	den = prim(IDN,k,j,i);
	vx = prim(IVX,k,j,i);
	vy = prim(IVY,k,j,i);

	// primary gravity
	l = sqrt(x2*x2+(r0+x1)*(r0+x1));
	Fx = ((r0+x1)/l)*(1/l/l-OmegaP*OmegaP*l)-2*OmegaP*vy;
	Fy = (x2/l)*(1/l/l-OmegaP*OmegaP*l)+2*OmegaP*vx;
	cons(IM1,k,j,i) -= dt*den*Fx;
	cons(IM2,k,j,i) -= dt*den*Fy;

	// orbiter gravity
	rs = sqrt(x1*x1+x2*x2+soft_length*soft_length);
	Fx = (-gm1*x1/rs)/rs/rs;
        Fy = (-gm1*x2/rs)/rs/rs;
	cons(IM1,k,j,i) += dt*den*Fx;
	cons(IM2,k,j,i) += dt*den*Fy;
      }
    }
  }
  return;
}

Real GravForceCalc(MeshBlock *pmb, int iout) {
  Real tauX=0, tauY=0, den, x1, x2, x3, rs, vol, dF;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  for(int k=ks; k<=ke; k++) {
    x3 = pmb->pcoord->x3v(k);
    for(int j=js; j<=je; j++) {
      x2 = pmb->pcoord->x2v(j);
      for(int i=is; i<=ie; i++) {
	x1 = pmb->pcoord->x1v(i);
	den = pmb->phydro->u(IDN, k, j, i);
        vol = pmb->pcoord->GetCellVolume(k,j,i);
	rs = sqrt(x1*x1+x2*x2+soft_length*soft_length);
	dF = gm1*x1/rs/rs/rs*den*vol;
	tauX += (x1/rs)*dF;
	tauY += (x2/rs)*dF;
      }
    }
  }
  if (iout==0) {
    return tauX;
  } else if (iout == 1) {
    return tauY;
  } else {
    return 0.;
  }
}

void MeshBlock::UserWorkInLoop() {
  Real l, x1, x2, x3, den, radp, phip, Cx, Cy, vk, dt;
  for (int k=ks; k<=ke; k++) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; j++) {
      x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; i++) {
	x1 = pcoord->x1v(i);

	dt = pmy_mesh->dt;
	l = sqrt(x2*x2+(r0+x1)*(r0+x1));
	radp = sqrt(x2*x2+(r0+x1)*(r0+x1));
	phip = std::atan2(x2, r0+x1);
        den = DenProfileCyl(radp,phip,x3);
	Cx = -x2/l;
	Cy = (r0+x1)/l;
	den = DenProfileCyl(l,phip,x3);
	vk = VelProfileCyl(l,phip,x3);

	// wave damping at the inner primary
	if (l < innerbdy) {
	  phydro->u(IM1,k,j,i) += -dt*(phydro->u(IM1,k,j,i)-Cx*den*vk)/T_damp_prim;
	  phydro->u(IM2,k,j,i) += -dt*(phydro->u(IM2,k,j,i)-Cy*den*vk)/T_damp_prim;
	}
	// wave damping at the edges
	if ((x1 < x1min+WDL) or (x1>x1max-WDL) or (x2<x2min+WDL) or (x2>x2max-WDL)){
	  phydro->u(IM1,k,j,i) += -dt*(phydro->u(IM1,k,j,i)-Cx*den*vk)/T_damp_bdy;
	  phydro->u(IM2,k,j,i) += -dt*(phydro->u(IM2,k,j,i)-Cy*den*vk)/T_damp_bdy;
	}
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
	radp = sqrt(x2*x2+(r0+x1)*(r0+x1));
	phip = std::atan2(x2, r0+x1);
	Cx = -x2/radp;
	Cy = (r0+x1)/radp;
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
	radp = sqrt(x2*x2+(r0+x1)*(r0+x1));
	phip = std::atan2(x2, r0+x1);
	Cx = -x2/radp;
	Cy = (r0+x1)/radp;
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
	radp = sqrt(x2*x2+(r0+x1)*(r0+x1));
	phip = std::atan2(x2, r0+x1);
	Cx = -x2/radp;
	Cy = (r0+x1)/radp;
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
	radp = sqrt(x2*x2+(r0+x1)*(r0+x1));
	phip = std::atan2(x2, r0+x1);
	Cx = -x2/radp;
	Cy = (r0+x1)/radp;
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
