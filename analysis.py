import sys
import numpy as np
sys.path.insert(0, "/home/astrosun/jtlaune/athena/")
import athplot

def reduceProfiles(filename):
    """Takes an hdf5 athdf filename, restricts it to the base grid, and calculates Jdots, Mdot, the j profile.

    Returns: dictionary with keys
    "JdotAdv": The advective component of Jdot.
    "JdotVisc": The viscous component of Jdot.
    "dTgravdr": Differential torque from the planet's gravity.
    "Tgrav_gtrr": The gravitational component of Jdot.
    "Mdot": The Mdot profile.
    "jProf": The specific AM (j) profile.
    "r": The radial coordinates.
    """
    data = athplot.rawDataRestricted(filename)
    cccs, lds = data.get_level(0)
    RR, PHI = np.meshgrid(cccs[0], cccs[1])
    Nphi, Nr = RR.shape
    DENS = lds["dens"]
    VR = lds["mom1"] / DENS
    VPHI = lds["mom2"] / DENS + RR

    dVRdPHI = np.gradient(VR, PHI[:, 0], axis=0, edge_order=2)
    dVPHI_RRdRR = np.gradient(VPHI / RR, RR[0, :], axis=1, edge_order=2)

    JdotAdv = np.zeros(RR.shape[1])
    JdotVisc = np.zeros(RR.shape[1])
    dTgravdr = np.zeros(RR.shape[1])
    Tgrav_gtrr = np.zeros(RR.shape[1])
    Mdot = np.zeros(RR.shape[1])
    jProf = np.zeros(RR.shape[1])

    dRR = np.gradient(RR, axis=1, edge_order=2)
    dPHI = np.gradient(PHI, axis=0, edge_order=2)

    q = 1e-4
    GPOT = -q / np.sqrt(1 + RR**2 - 2 * RR * np.cos(PHI)) - 1 / RR
    dGPOTdPHI = np.gradient(GPOT, PHI[:, 0], axis=0, edge_order=2)

    nu = 1e-4
    for i in range(Nr):
        for j in range(Nphi):
            r = RR[j, i]
            phi = PHI[j, i]
            Sig = DENS[j, i]
            vr = VR[j, i]
            vphi = VPHI[j, i]
            dr = dRR[j, i]
            dphi = dPHI[j, i]

            dvrdphi = dVRdPHI[j, i]
            dvphi_rdr = dVPHI_RRdRR[j, i]
            dpotdphi = dGPOTdPHI[j, i]

            JdotAdv[i] += -(r**2) * Sig * vr * vphi * dphi
            JdotVisc[i] += (
                -(r**3) * nu * Sig * (dvphi_rdr + 1 / r**2 * dvrdphi) * dphi
            )
            dTgravdr[i] += -r * Sig * dpotdphi * dphi
            Mdot[i] += -r * vr * Sig * dphi
            jProf[i] += (r * vphi * dphi) / (2 * np.pi)

    for i in range(Nr):
        for j in range(Nr - i):
            dr = dRR[0, i + j]
            Tgrav_gtrr[i] += dTgravdr[i + j] * dr

    return {
        "JdotAdv": JdotAdv,
        "JdotVisc": JdotVisc,
        "dTgravdr": dTgravdr,
        "Tgrav_gtrr": Tgrav_gtrr,
        "Mdot": Mdot,
        "jProf": jProf,
        "r":RR[0,:],
    }