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

namespace {
Real gm1, Sig0, dslope, dfloor, R0, CS02, Omega0, soft_sat;
Real T_damp_in, T_damp_bdy, WDL1, WDL2, innerbdy, x1min, l_refine;
Real rSink, rEval;
int nPtEval;
} // namespace

Real DenProf(const Real rad);
Real VelProf(const Real rad);
int RefinementCondition(MeshBlock *pmb);
Real Measurements(MeshBlock *pmb, int iout);

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

void Mesh::InitUserMeshData(ParameterInput *pin) {
  EnrollUserExplicitSourceFunction(DiskSourceFunction);
  if (adaptive == true) {
    EnrollUserRefinementCondition(RefinementCondition);
  }
  AllocateUserHistoryOutput(5);
  EnrollUserHistoryOutput(0, Measurements, "Fs,grav_x=r");
  EnrollUserHistoryOutput(1, Measurements, "Fs,grav_y=th");
  EnrollUserHistoryOutput(2, Measurements, "AccretionEval");
  EnrollUserHistoryOutput(3, Measurements, "momxAccretion");
  EnrollUserHistoryOutput(4, Measurements, "momyAccretion");

  dfloor = pin->GetReal("hydro", "Sigma_floor");

  gm1 = pin->GetReal("problem", "GM_s");
  Sig0 = pin->GetReal("problem", "Sigma_0");
  dslope = pin->GetReal("problem", "delta");
  R0 = pin->GetReal("problem", "R0");
  Omega0 = pin->GetReal("problem", "Omega0");
  soft_sat = pin->GetReal("problem", "soft_sat");
  l_refine = pin->GetReal("problem", "l_refine");
  CS02 = SQR(pin->GetReal("hydro", "iso_sound_speed"));

  WDL1 = pin->GetReal("problem", "WaveDampingLength_in");
  WDL2 = pin->GetReal("problem", "WaveDampingLength_out");
  innerbdy = pin->GetReal("problem", "innerbdy");
  x1min = pin->GetReal("mesh", "x1min");
  T_damp_bdy = pin->GetReal("problem", "T_damp_bdy");
  T_damp_in = pin->GetReal("problem", "T_damp_in");

  rSink = pin->GetReal("problem", "sink_radius");
  rEval = pin->GetReal("problem", "eval_radius");
  nPtEval = pin->GetReal("problem", "N_eval_pts");

  // enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiskInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiskOuterX1);
  }
  return;
}

Real DenProf(const Real rad) {
  // Density profile Sigma(r)
  return (std::max(Sig0 * std::pow(rad / R0, dslope), dfloor));
}

Real VelProf(const Real rad) {
  // Velocity profile v(r)
  return std::sqrt(dslope * CS02 + 1 / rad) - Omega0 * rad;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real x1, x2, x3;
  Real Sig, vk, Cx, Cy;
  Real rprim;
  for (int k = ks; k <= ke; ++k) {
    x3 = pcoord->x3v(k);
    for (int j = js; j <= je; ++j) {
      x2 = pcoord->x2v(j);
      for (int i = is; i <= ie; ++i) {
        x1 = pcoord->x1v(i);

        rprim = x1;
        Sig = DenProf(rprim);
        vk = VelProf(rprim);

        phydro->u(IDN, k, j, i) = Sig;
        phydro->u(IM1, k, j, i) = 0.;
        phydro->u(IM2, k, j, i) = Sig * vk;
        phydro->u(IM3, k, j, i) = 0.;
      }
    }
  }
}

Real Measurements(MeshBlock *pmb, int iout) {
  // User-defined history function. Must calculate a sum of values
  // within each MeshBlock & athena turns it into a global sum automatically.
  Real F_x = 0, F_y = 0, Sig, x1, x2, x3;
  Real rs, vol, Fmag, rsecn, rsoft;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks,
      ke = pmb->ke;
  Real momxEvaldir, momyEvaldir, angEval, xEval, yEval;
  Real rf1, rf2, thf1, thf2;
  Real rpeval, theval;
  Real AccRate = 0;
  Real momxRate = 0;
  Real momyRate = 0;
  bool InsideCell;

  Real FxEvals[nPtEval];
  Real FyEvals[nPtEval];
  Real evalVals[nPtEval];
  for (int l = 0; l <= nPtEval; l++) {
    evalVals[l] = 0;
    FxEvals[l] = 0;
    FyEvals[l] = 0;
  }

  for (int k = ks; k <= ke; k++) {
    x3 = pmb->pcoord->x3v(k);
    for (int j = js; j <= je; j++) {
      x2 = pmb->pcoord->x2v(j);
      for (int i = is; i <= ie; i++) {
        x1 = pmb->pcoord->x1v(i);

        Sig = pmb->phydro->u(IDN, k, j, i);
        vol = pmb->pcoord->GetCellVolume(k, j, i);

        rsecn = std::sqrt(x1 * x1 + R0 * R0 - 2 * R0 * x1 * std::cos(x2));
        rsoft = std::sqrt(rsecn * rsecn + soft_sat * soft_sat);

        Fmag = gm1 / rsoft / rsoft / rsoft;
        // Fmag = gm1/rsecn/rsecn/rsecn;
        F_x += Sig * vol * Fmag * (std::cos(x2) * x1 - R0);
        F_y += Sig * vol * Fmag * x1 * std::sin(x2);

        // For now, doing nearest neighbor matching algorithm??
        for (int l = 0; l <= nPtEval; l++) {
          angEval = 2 * PI / nPtEval * l;

          // sample locus of circle at secondary with radius rEval
          xEval = rEval * std::cos(angEval) + R0;
          yEval = rEval * std::sin(angEval);

          rf1 = pmb->pcoord->x1f(i);
          rf2 = pmb->pcoord->x1f(i + 1);
          thf1 = pmb->pcoord->x2f(j);
          thf2 = pmb->pcoord->x2f(j + 1);

          InsideCell = false;
          rpeval = sqrt(xEval * xEval + yEval * yEval);
          theval = std::atan2(yEval, xEval);

          if ((rf1 <= rpeval) && (rpeval < rf2)) {
            if ((thf1 <= theval) && (theval < thf2)) {
              InsideCell = true;
            }
          }

          if (InsideCell) {
            // if (drEval < drvalEvals[l]) {
            //   dA = ((2*pi/nPtEval)*rEval)
            //   -Sig*((2*pi/nPtEval)*rEval)*(u1cos+u2sin)
            momxEvaldir = std::cos(x2) * (pmb->phydro->u(IM1, k, j, i)) -
                          std::sin(x2) * (pmb->phydro->u(IM2, k, j, i));
            momyEvaldir = std::sin(x2) * (pmb->phydro->u(IM1, k, j, i)) +
                          std::cos(x2) * (pmb->phydro->u(IM2, k, j, i));
            FxEvals[l] = -((2 * PI) / nPtEval / Sig) *
                         (momxEvaldir * momxEvaldir * std::cos(angEval) +
                          momxEvaldir * momyEvaldir * std::cos(angEval));
            FyEvals[l] = -((2 * PI) / nPtEval / Sig) *
                         (momyEvaldir * momxEvaldir * std::cos(angEval) +
                          momyEvaldir * momyEvaldir * std::cos(angEval));
            evalVals[l] = -((2 * PI) / nPtEval) * rEval *
                          ((momxEvaldir) * (std::cos(angEval)) +
                           (momyEvaldir) * (std::sin(angEval)));
          }
        }
      }
    }
  }

  for (int l = 0; l <= nPtEval; l++) {
    // some diagnostic-ing
    // remember this outputs for each meshblock so this will stop the code
    // every time if (drvalEvals[l] == 1) {
    //  std::cout << l;
    //  throw std::invalid_argument(
    //      " <- dr was not minimized at this rEval circle index");
    // }

    // evalVals[l] being different from its initialized value means it was
    // altered in this meshblock.
    if (evalVals[l] != 0) {
      momxRate += FxEvals[l];
      momyRate += FyEvals[l];
      AccRate += evalVals[l];
    }
  }

  if (iout == 0) {
    return F_x;
  } else if (iout == 1) {
    return F_y;
  } else if (iout == 2) {
    return AccRate;
  } else if (iout == 3) {
    return momxRate;
  } else if (iout == 4) {
    return momyRate;
  } else {
    return 0.;
  }
}

Real splineKernel(Real x) {
  // smooth transition f(x)=0 @ x=1 and f(x)=1 @ x=0
  Real W;
  if (x < 0) {
    W = 1.;
  } else if (x <= 0.5) {
    W = 1. - 6. * x * x + 6. * x * x * x;
  } else if (x <= 1) {
    W = 2. * (1 - x) * (1 - x) * (1 - x);
  } else {
    W = 0.;
  }
  return (W);
}

void Mesh::UserWorkInLoop() {
  Real rsecn, x1, x2, x3;
  for (int bn = 0; bn < nblocal; ++bn) {
    MeshBlock *pmb = my_blocks(bn);
    // LogicalLocation &loc = pmb->loc;
    // if (loc.level == nLevel) { // lowest level
    for (int k = pmb->ks; k <= pmb->ke; k++) {
      x3 = pmb->pcoord->x3v(k);
      for (int j = pmb->js; j <= pmb->je; j++) {
        x2 = pmb->pcoord->x2v(j);
        for (int i = pmb->is; i <= pmb->ie; i++) {
          x1 = pmb->pcoord->x1v(i);
          rsecn = std::sqrt(x1 * x1 + R0 * R0 - 2 * R0 * x1 * std::cos(x2));
          if (rsecn < rSink) {
            pmb->phydro->w(IVX, k, j, i) = 0;
            pmb->phydro->w(IVY, k, j, i) = 0;
            pmb->phydro->w(IVZ, k, j, i) = 0;
            // USES EOS EQUATION OF STATE
            pmb->phydro->w(IDN, k, j, i) = dfloor;
            pmb->phydro->w(IPR, k, j, i) = dfloor * CS02;

            pmb->phydro->u(IDN, k, j, i) = dfloor;
            pmb->phydro->u(IM1, k, j, i) = 0;
            pmb->phydro->u(IM2, k, j, i) = 0;
            pmb->phydro->u(IM3, k, j, i) = 0;
            pmb->phydro->u(IDN, k, j, i) = dfloor;
          }
        }
      }
    }
    //}
  }
  return;
}

void DiskSourceFunction(MeshBlock *pmb, const Real time, const Real dt,
                        const AthenaArray<Real> &prim,
                        const AthenaArray<Real> &prim_scalar,
                        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                        AthenaArray<Real> &cons_scalar) {
  Real vk, vr, vth, Sig, Sig0, rprim, rsecn;
  Real x1, x2, x3;
  Real Fpr, Cs;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    x3 = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      x2 = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        x1 = pmb->pcoord->x1v(i);

        rprim = x1;
        rsecn =
            std::sqrt(rprim * rprim + R0 * R0 - 2 * R0 * rprim * std::cos(x2));
        vk = VelProf(rprim);
        Sig0 = DenProf(rprim);
        Sig = prim(IDN, k, j, i);

        // primary gravity
        Fpr = -1. / rprim / rprim;
        cons(IM1, k, j, i) += dt * Sig * Fpr;
        cons(IM2, k, j, i) += 0;

        // centrifugal
        cons(IM1, k, j, i) += dt * Sig * Omega0 * Omega0 * rprim;

        // coriolis
        vr = prim(IM1, k, j, i);
        vth = prim(IM2, k, j, i);
        cons(IM1, k, j, i) += 2 * dt * Sig * Omega0 * vth;
        cons(IM2, k, j, i) += -2 * dt * Sig * Omega0 * vr;

        // satellite gravity
        Cs = gm1 / (rsecn * rsecn + soft_sat * soft_sat) /
             std::sqrt(rsecn * rsecn + soft_sat * soft_sat);
        cons(IM1, k, j, i) += Sig * dt * Cs * (std::cos(x2) * R0 - x1);
        cons(IM2, k, j, i) += -Sig * dt * Cs * R0 * std::sin(x2);

        // indirect
        // cons(IM1, k, j, i) += gm1*std::cos(x2)/(R0*R0+soft_sat*soft_sat);
        // cons(IM2, k, j, i) += -gm1*std::sin(x2)/(R0*R0+soft_sat*soft_sat);

        // wave damping regions
        if ((rprim <= innerbdy) and (innerbdy != x1min)) {
          Real x = (rprim - x1min) / (innerbdy - x1min);
          Real factor = splineKernel(x);
          cons(IDN, k, j, i) +=
              -factor * dt * (cons(IDN, k, j, i) - Sig0) / T_damp_in;
          cons(IM1, k, j, i) += -factor * dt * (cons(IM1, k, j, i)) / T_damp_in;
          cons(IM2, k, j, i) +=
              -factor * dt * (cons(IM2, k, j, i) - Sig0 * vk) / T_damp_in;
        }
        if ((rprim >= WDL1) and (WDL2 != WDL1)) {
          Real x = 1 - (rprim - WDL1) / (WDL2 - WDL1);
          Real factor = splineKernel(x);
          cons(IDN, k, j, i) +=
              -factor * dt * (cons(IDN, k, j, i) - Sig0) / T_damp_bdy;
          cons(IM1, k, j, i) +=
              -factor * dt * (cons(IM1, k, j, i)) / T_damp_bdy;
          cons(IM2, k, j, i) +=
              -factor * dt * (cons(IM2, k, j, i) - Sig0 * vk) / T_damp_bdy;
        }
      }
    }
  }
}

int RefinementCondition(MeshBlock *pmb) {
  Real x1min = pmb->pcoord->x1v(pmb->is);
  Real x1max = pmb->pcoord->x1v(pmb->ie);
  Real x2min = pmb->pcoord->x2v(pmb->js);
  Real x2max = pmb->pcoord->x2v(pmb->je);
  Real cond = 0;
  if ((x1min < (R0 + l_refine)) and (x1max > (R0 - l_refine))) {
    // radial band
    if ((x2min < 2 * PI * l_refine / R0) and
        (x2max > -2 * PI * l_refine / R0)) {
      // azimuthal band
      cond = 1;
    }
  }
  return (cond);
}

void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                 FaceField &b, Real time, Real dt, int il, int iu, int jl,
                 int ju, int kl, int ku, int ngh) {
  Real Sig, vk, rprim, x1, x2, x3, Cx, Cy;
  for (int k = kl; k <= ku; ++k) {
    x3 = pco->x3v(k);
    for (int j = jl; j <= ju; ++j) {
      x2 = pco->x2v(j);
      for (int i = 1; i <= ngh; ++i) {
        x1 = pco->x1v(il - i);

        rprim = x1;
        Sig = DenProf(rprim);
        vk = VelProf(rprim);

        prim(IDN, k, j, il - i) = Sig;
        prim(IM1, k, j, il - i) = 0.;
        prim(IM2, k, j, il - i) = vk;
        prim(IM3, k, j, il - i) = 0.;
      }
    }
  }
}

void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                 FaceField &b, Real time, Real dt, int il, int iu, int jl,
                 int ju, int kl, int ku, int ngh) {
  Real Sig, vk, rprim, x1, x2, x3, Cx, Cy;
  for (int k = kl; k <= ku; ++k) {
    x3 = pco->x3v(k);
    for (int j = jl; j <= ju; ++j) {
      x2 = pco->x2v(j);
      for (int i = 1; i <= ngh; ++i) {
        x1 = pco->x1v(iu + i);

        rprim = x1;
        Sig = DenProf(rprim);
        vk = VelProf(rprim);

        prim(IDN, k, j, iu + i) = Sig;
        prim(IM1, k, j, iu + i) = 0.;
        prim(IM2, k, j, iu + i) = vk;
        prim(IM3, k, j, iu + i) = 0.;
      }
    }
  }
}
