// C++ headers
#include <cmath>     // sqrt()
#include <fstream>   // ofstream
#include <iomanip>   // setprecision
#include <iostream>  // cout, endl
#include <sstream>   // stringstream
#include <stdexcept> // runtime_error
#include <string>    // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../outputs/outputs.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

namespace
{
    Real A, m, qshear, Omega0, cs0;
}

Real DenProf(const Real rad);
Real AzimVelProf(const Real rad);
Real RadVelProf(const Real rad);
void DiskSourceFunction(MeshBlock *pmb, const Real time, const Real dt,
                        const AthenaArray<Real> &prim,
                        const AthenaArray<Real> &prim_scalar,
                        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                        AthenaArray<Real> &cons_scalar);

void OuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int il, int iu, int jl,
             int ju, int kl, int ku, int ngh);
void InnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int il, int iu, int jl,
             int ju, int kl, int ku, int ngh);
void OuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int il, int iu, int jl,
             int ju, int kl, int ku, int ngh);
void InnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int il, int iu, int jl,
             int ju, int kl, int ku, int ngh);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
    // Source function with perturbation + rotating frame terms.
    EnrollUserExplicitSourceFunction(DiskSourceFunction);

    // m-wave component parameters
    m = pin->GetReal("problem", "m");
    A = pin->GetReal("problem", "A");
    // shearing box parameters
    Omega0 = pin->GetReal("problem", "Omega0");
    qshear = pin->GetReal("problem", "qshear");
    // hydro parameters
    cs0 = pin->GetReal("hydro", "iso_sound_speed");

    // Boundary conditions.
    if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user"))
    {
        EnrollUserBoundaryFunction(BoundaryFace::inner_x1, InnerX1);
    }
    if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user"))
    {
        EnrollUserBoundaryFunction(BoundaryFace::outer_x1, OuterX1);
    }
    if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user"))
    {
        EnrollUserBoundaryFunction(BoundaryFace::inner_x2, InnerX2);
    }
    if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user"))
    {
        EnrollUserBoundaryFunction(BoundaryFace::outer_x2, OuterX2);
    }
}

Real DenProf(const Real x1)
{
    return 1.0;
}
Real AzimVelProf(const Real x1)
{
    return -qshear * Omega0 * x1;
}
Real RadVelProf(const Real x1)
{
    return 0.0;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
    Real x1, x2, x3;
    Real den, vy, vx;
    for (int k = ks; k <= ke; ++k)
    {
        x3 = pcoord->x3v(k);
        for (int j = js; j <= je; ++j)
        {
            x2 = pcoord->x2v(j);
            for (int i = is; i <= ie; ++i)
            {
                x1 = pcoord->x1v(i);

                den = DenProf(x1);
                vy = AzimVelProf(x1);
                vx = RadVelProf(x1);

                phydro->u(IDN, k, j, i) = den;
                phydro->u(IM1, k, j, i) = den * vx;
                phydro->u(IM2, k, j, i) = den * vy;
                phydro->u(IM3, k, j, i) = 0.;
            }
        }
    }
}

void DiskSourceFunction(MeshBlock *pmb, const Real time, const Real dt,
                        const AthenaArray<Real> &prim,
                        const AthenaArray<Real> &prim_scalar,
                        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                        AthenaArray<Real> &cons_scalar)
{
    Real x1, x2, x3, den, vx, vy;
    for (int k = pmb->ks; k <= pmb->ke; ++k)
    {
        x3 = pmb->pcoord->x3v(k);
        for (int j = pmb->js; j <= pmb->je; ++j)
        {
            x2 = pmb->pcoord->x2v(j);
            for (int i = pmb->is; i <= pmb->ie; ++i)
            {
                x1 = pmb->pcoord->x1v(i);
                den = prim(IDN, k, j, i);
                vx = prim(IVX, k, j, i);
                vy = prim(IVY, k, j, i);

                // Centrifugal force.
                cons(IM1, k, j, i) += dt * den * 2 * qshear * x1 * Omega0 * Omega0;

                // Coriolis force.
                cons(IM1, k, j, i) += dt * den * 2 * Omega0 * vy;
                cons(IM2, k, j, i) += -dt * den * 2 * Omega0 * vx;

                // m-component wave perturbation.
                cons(IM2, k, j, i) += dt * den * A * m * std::sin(m * x2);
            }
        }
    }
}

void OuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int il, int iu, int jl,
             int ju, int kl, int ku, int ngh)
{
    Real x1, x2, x3;
    Real den, vx, vy;
    for (int k = kl; k <= ku; ++k)
    {
        x3 = pco->x3v(k);
        for (int j = jl; j <= ju; ++j)
        {
            x2 = pco->x2v(j);
            for (int i = 1; i <= ngh; ++i)
            {
                x1 = pco->x1v(iu + i);

                den = DenProf(x1);
                vx = RadVelProf(x1);
                vy = AzimVelProf(x1);

                prim(IDN, k, j, iu + i) = den;
                prim(IVX, k, j, iu + i) = vx;
                prim(IVY, k, j, iu + i) = vy;
                prim(IVZ, k, j, iu + i) = 0.;
            }
        }
    }
}
void InnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int il, int iu, int jl,
             int ju, int kl, int ku, int ngh)
{
    Real x1, x2, x3;
    Real den, vx, vy;
    for (int k = kl; k <= ku; ++k)
    {
        x3 = pco->x3v(k);
        for (int j = jl; j <= ju; ++j)
        {
            x2 = pco->x2v(j);
            for (int i = 1; i <= ngh; ++i)
            {
                x1 = pco->x1v(il - i);

                den = DenProf(x1);
                vx = RadVelProf(x1);
                vy = AzimVelProf(x1);

                prim(IDN, k, j, il - i) = den;
                prim(IVX, k, j, il - i) = vx;
                prim(IVY, k, j, il - i) = vy;
                prim(IVZ, k, j, il - i) = 0.;
            }
        }
    }
}
void OuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int il, int iu, int jl,
             int ju, int kl, int ku, int ngh)
{
    Real x1, x2, x3;
    Real den, vx, vy;
    for (int k = kl; k <= ku; ++k)
    {
        x3 = pco->x3v(k);
        for (int j = 1; j <= ngh; ++j)
        {
            x2 = pco->x2v(ju + j);
            for (int i = il; i <= iu; ++i)
            {
                x1 = pco->x1v(i);

                den = DenProf(x1);
                vx = RadVelProf(x1);
                vy = AzimVelProf(x1);

                // x2 (y) > 0
                if (x1 > 0)
                {
                    // inflow fixed to initial
                    prim(IDN, k, ju + j, i) = den;
                    prim(IVX, k, ju + j, i) = vx;
                    prim(IVY, k, ju + j, i) = vy;
                    prim(IVZ, k, ju + j, i) = 0.;
                }
                else
                {
                    // outflow copied from last active
                    prim(IDN, k, ju + j, i) = prim(IDN, k, ju, i);
                    prim(IVX, k, ju + j, i) = prim(IVX, k, ju, i);
                    prim(IVY, k, ju + j, i) = prim(IVY, k, ju, i);
                    prim(IVZ, k, ju + j, i) = prim(IVZ, k, ju, i);
                }
            }
        }
    }
}
void InnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int il, int iu, int jl,
             int ju, int kl, int ku, int ngh)
{
    Real x1, x2, x3;
    Real den, vx, vy;
    for (int k = kl; k <= ku; ++k)
    {
        x3 = pco->x3v(k);
        for (int j = 1; j <= ngh; ++j)
        {
            x2 = pco->x2v(jl - j);
            for (int i = il; i <= iu; ++i)
            {
                x1 = pco->x1v(i);

                den = DenProf(x1);
                vx = RadVelProf(x1);
                vy = AzimVelProf(x1);

                // x2 (y) < 0
                if (x1 > 0)
                {
                    // outflow copied from last active
                    prim(IDN, k, jl - j, i) = prim(IDN, k, jl, i);
                    prim(IVX, k, jl - j, i) = prim(IVX, k, jl, i);
                    prim(IVY, k, jl - j, i) = prim(IVY, k, jl, i);
                    prim(IVZ, k, jl - j, i) = prim(IVZ, k, jl, i);
                }
                else
                {
                    // inflow fixed to initial
                    prim(IDN, k, jl - j, i) = den;
                    prim(IVX, k, jl - j, i) = vx;
                    prim(IVY, k, jl - j, i) = vy;
                    prim(IVZ, k, jl - j, i) = 0.;
                }
            }
        }
    }
}