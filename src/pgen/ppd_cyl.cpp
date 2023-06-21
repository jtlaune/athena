//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code
// contributors Licensed under the 3-clause BSD License, see LICENSE file for
// details Created JTL 10/19/22
//========================================================================================
// C++ headers
#include <algorithm> // min
#include <cmath>     // sqrt
#include <cstdlib>   // srand
#include <cstring>   // strcmp()
#include <fstream>
#include <iostream> // endl
#include <limits>
#include <math.h>
#include <sstream>   // stringstream
#include <stdexcept> // runtime_error
#include <string>    // c_str()

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

namespace
{
  Real ovSig, lambda, dfloor, pfloor, R0, cs0, CS02, Omega0, nu_iso;
  Real T_damp_bdy, WDL1, WDL2, x1min, x1max;
  Real r_exclude, innerbdy, A, m;
  int nPtEval;
} // namespace

Real DenProf(const Real rad);
Real dDenProfdr(const Real rad);
Real AzimVelProf(const Real rad);
Real RadVelProf(const Real rad);
Real dPresProfdr(const Real rad);
// Real hProf(const Real rad);

// User source function
void DiskSourceFunction(MeshBlock *pmb, const Real time, const Real dt,
                        const AthenaArray<Real> &prim,
                        const AthenaArray<Real> &prim_scalar,
                        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                        AthenaArray<Real> &cons_scalar);
// User-defined boundary conditions for disk simulations
void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                 FaceField &b, Real time, Real dt, int il, int iu, int jl,
                 int ju, int kl, int ku, int ngh);
void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                 FaceField &b, Real time, Real dt, int il, int iu, int jl,
                 int ju, int kl, int ku, int ngh);
Real splineKernel(Real x);
// void MeshBlock::UserWorkInLoop();

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserExplicitSourceFunction(DiskSourceFunction);

  // Hydro parameters.
  dfloor = pin->GetReal("hydro", "Sigma_floor");
  pfloor = pin->GetReal("hydro", "P_floor");
  cs0 = pin->GetReal("hydro", "iso_sound_speed");
  CS02 = SQR(cs0);
  // Secondary parameters.
  // Disk parameters.
  ovSig = pin->GetReal("problem", "overlineSigma");
  lambda = pin->GetReal("problem", "lambda");
  R0 = pin->GetReal("problem", "R0");
  Omega0 = pin->GetReal("problem", "Omega0");
  nu_iso = pin->GetReal("problem", "nu_iso");
  // Boundary conditions.
  WDL1 = pin->GetReal("problem", "WaveDampingLength_in");
  WDL2 = pin->GetReal("problem", "WaveDampingLength_out");
  x1min = pin->GetReal("mesh", "x1min");
  x1max = pin->GetReal("mesh", "x1max");
  T_damp_bdy = pin->GetReal("problem", "T_damp_bdy");
  innerbdy = pin->GetReal("problem", "innerbdy");
  // Calculation parameters.
  r_exclude = pin->GetOrAddReal("problem", "gforce_r_exclude", 0);
  // Potential parameters
  innerbdy = pin->GetReal("problem", "innerbdy");
  m = pin->GetReal("problem", "m");
  A = pin->GetReal("problem", "A");

  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiskInnerX1);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiskOuterX1);
  return;
}

Real splineKernel(Real x)
{
  // Smooth transition f(x)=0 @ x=1 and f(x)=1 @ x=0.
  // Used for the wave damping zones at radial boundary (with appropriate
  // definitions for x).
  Real W;
  if (x < 0)
  {
    W = 1.;
  }
  else if (x <= 0.5)
  {
    W = 1. - 6. * x * x + 6. * x * x * x;
  }
  else if (x <= 1)
  {
    W = 2. * (1 - x) * (1 - x) * (1 - x);
  }
  else
  {
    W = 0.;
  }
  return (W);
}

/*
  These functions with "Prof" in their name are the unperturbed steady state (parameterized) background profiles.
  Currently, only supports cs=const (lambda=0, temperature=constant, global isothermal).
*/

Real DenProf(const Real rad)
{
  // Density profile Sigma(r) (unperturbed)
  // x1max and rad both in units R0
  return ovSig;
}

Real dDenProfdr(const Real rad)
{
  // Derivative of density profile d(Sigma(r))/dr (unperturbed)
  return 0;
}

Real dPresProfdr(const Real rad)
{
  Real dSigdr = dDenProfdr(rad);
  return std::pow(cs0, 2) * dSigdr;
}

Real RadVelProf(const Real rad)
{
  return -3 * nu_iso / 2 / (rad / R0);
}

Real AzimVelProf(const Real rad)
{
  // Velocity profile v(r)
  Real dPdr = dPresProfdr(rad);
  Real Sig = DenProf(rad);
  return std::sqrt(dPdr * rad / Sig + 1 / rad) - Omega0 * rad;
}

void DiskSourceFunction(MeshBlock *pmb, const Real time, const Real dt,
                        const AthenaArray<Real> &prim,
                        const AthenaArray<Real> &prim_scalar,
                        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                        AthenaArray<Real> &cons_scalar)
{
  Real vk, vr, vr0, vth, Sig, Sig0, rprim, rsecn;
  Real x1, x2, x3;
  Real Fpr, Cs;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
  {
    x3 = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j)
    {
      x2 = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i)
      {
        x1 = pmb->pcoord->x1v(i);

        rprim = x1;
        rsecn =
            std::sqrt(rprim * rprim + R0 * R0 - 2 * R0 * rprim * std::cos(x2));
        vk = AzimVelProf(rprim);
        Sig0 = DenProf(rprim);
        Sig = prim(IDN, k, j, i);

        // Primary gravity.
        // No softening because it is off-grid.
        Fpr = -1. / rprim / rprim;
        cons(IM1, k, j, i) += dt * Sig * Fpr;
        cons(IM2, k, j, i) += 0;

        // m component
        cons(IM2, k, j, i) += dt * Sig * A * m * std::sin(m * x2 - Omega0 * time) / rprim;

        // Centrifugal force.
        cons(IM1, k, j, i) += dt * Sig * Omega0 * Omega0 * rprim;

        // Coriolis force.
        vr = prim(IM1, k, j, i);
        vth = prim(IM2, k, j, i);
        cons(IM1, k, j, i) += 2 * dt * Sig * Omega0 * vth;
        cons(IM2, k, j, i) += -2 * dt * Sig * Omega0 * vr;

        // Wave damping regions.
        if ((rprim <= innerbdy) and (innerbdy != x1min))
        {
          Real x = (rprim - x1min) / (innerbdy - x1min);
          Real factor = splineKernel(x);
          cons(IDN, k, j, i) +=
              -factor * dt * (cons(IDN, k, j, i) - Sig0) / T_damp_bdy;
          cons(IM1, k, j, i) += -factor * dt * (cons(IM1, k, j, i) - Sig0 * vr0) / T_damp_bdy;
          cons(IM2, k, j, i) +=
              -factor * dt * (cons(IM2, k, j, i) - Sig0 * vk) / T_damp_bdy;
        }
        if ((rprim >= WDL1) and (WDL2 != WDL1))
        {
          Real x = 1 - (rprim - WDL1) / (WDL2 - WDL1);
          Real factor = splineKernel(x);
          vr0 = RadVelProf(rprim); //-1.5 * nu_iso / rprim;
          cons(IDN, k, j, i) +=
              -factor * dt * (cons(IDN, k, j, i) - Sig0) / T_damp_bdy;
          cons(IM1, k, j, i) +=
              -factor * dt * (cons(IM1, k, j, i) - Sig0 * vr0) / T_damp_bdy;
          cons(IM2, k, j, i) +=
              -factor * dt * (cons(IM2, k, j, i) - Sig0 * vk) / T_damp_bdy;
        }
      }
    }
  }
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real x1, x2, x3;
  Real Sig, vk, Cx, Cy;
  Real rprim, vr0;
  for (int k = ks; k <= ke; ++k)
  {
    x3 = pcoord->x3v(k);
    for (int j = js; j <= je; ++j)
    {
      x2 = pcoord->x2v(j);
      for (int i = is; i <= ie; ++i)
      {
        x1 = pcoord->x1v(i);

        rprim = x1;
        Sig = DenProf(rprim);
        vk = AzimVelProf(rprim);
        vr0 = RadVelProf(rprim); //-1.5 * nu_iso / rprim;

        phydro->u(IDN, k, j, i) = Sig;
        phydro->u(IM1, k, j, i) = Sig * vr0;
        phydro->u(IM2, k, j, i) = Sig * vk;
        phydro->u(IM3, k, j, i) = 0.;
      }
    }
  }
}

void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                 FaceField &b, Real time, Real dt, int il, int iu, int jl,
                 int ju, int kl, int ku, int ngh)
{
  Real x1, x2, x3;
  Real rprim, Sig, vk, vr;
  Real r_active = pco->x1v(il);
  for (int k = kl; k <= ku; ++k)
  {
    x3 = pco->x3v(k);
    for (int j = jl; j <= ju; ++j)
    {
      x2 = pco->x2v(j);
      for (int i = 1; i <= ngh; ++i)
      {
        x1 = pco->x1v(il - i);

        rprim = x1;
        Sig = prim(IDN, k, j, il); // DenProf(rprim);
        vr = prim(IVX, k, j, il);  // RadVelProf(rprim); //-1.5 * nu_iso / rprim;
        vk = prim(IVY, k, j, il);  // AzimVelProf(rprim);

        prim(IDN, k, j, il - i) = Sig; // prim(IDN, k, j, il);
        prim(IVY, k, j, il - i) = (vk + r_active * Omega0) * std::pow((x1 / r_active), -0.5) - rprim * Omega0;
        prim(IVZ, k, j, il - i) = 0.; // prim(IVZ, k, j, il);
        if (vr >= 0)
        {
          prim(IVX, k, j, il - i) = 0;
        }
        else
        {
          prim(IVX, k, j, il - i) = vr * std::pow((x1 / r_active), -1);
        }
        // dumb debugging for nlim=1, dcycle=1
        // if (i == 1)
        //{
        //  std::cout << rprim << "\n";
        //  std::cout << r_active << "\n";
        //  std::cout << prim(IDN, k, j, il - i) << "\n";
        //  std::cout << prim(IVX, k, j, il - i) << "\n";
        //  std::cout << prim(IVY, k, j, il) << "\n";
        //  std::cout << prim(IVY, k, j, il - i) << "\n";
        //  std::cout << "--------end step--------"
        //            << "\n";
        //}
      }
    }
  }
}

/*
  Set to initial conditions (wave damping is in disk source)
*/

void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                 FaceField &b, Real time, Real dt, int il, int iu, int jl,
                 int ju, int kl, int ku, int ngh)
{
  Real Sig, vk, rprim, x1, x2, x3, Cx, Cy, vr;
  Real r_active = pco->x1v(iu);
  for (int k = kl; k <= ku; ++k)
  {
    x3 = pco->x3v(k);
    for (int j = jl; j <= ju; ++j)
    {
      x2 = pco->x2v(j);
      for (int i = 1; i <= ngh; ++i)
      {
        x1 = pco->x1v(iu + i);

        rprim = x1;
        Sig = prim(IDN, k, j, iu); // DenProf(rprim);
        vr = prim(IVX, k, j, iu);  // RadVelProf(rprim); //-1.5 * nu_iso / rprim;
        vk = prim(IVY, k, j, iu);  // AzimVelProf(rprim);

        prim(IDN, k, j, iu + i) = Sig;
        prim(IVX, k, j, iu + i) = vr * std::pow((x1 / r_active), -1);
        prim(IVY, k, j, iu + i) = (vk + r_active * Omega0) * std::pow((x1 / r_active), -0.5) - rprim * Omega0;
        prim(IVZ, k, j, iu + i) = 0.;
      }
    }
  }
}