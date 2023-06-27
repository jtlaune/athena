from matplotlib import ticker
from matplotlib import colors
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_styles import analytic
import yt

sys.path.insert(0, "/home/astrosun/jtlaune/athena/vis/python/")
import athena_read as athr


def plotAdjustKwargs(fig, ax, **kwargs):
    if "ybds" in kwargs.keys():
        ax.set_ylim(kwargs["ybds"][0], kwargs["ybds"][1])
    if "symlog" in kwargs.keys():
        if (type(kwargs["symlog"]) is bool) and (kwargs["symlog"]):
            if "linthresh" in kwargs.keys():
                ax.set_yscale("symlog", linthresh=kwargs["linthresh"])
            else:
                ax.set_yscale("symlog")
    if "yMajorMaxN" in kwargs.keys():
        ax.yaxis.set_major_locator(ticker.MaxNLocator(kwargs["yMajorMaxN"]))
    if "yMinorMaxN" in kwargs.keys():
        ax.yaxis.set_minor_locator(ticker.MaxNLocator(kwargs["yMinorMaxN"]))
    if "yTickLabelSize" in kwargs.keys():
        ax.yaxis.set_tick_params(labelsize=kwargs["yTickLabelSize"])


class rawDataRestricted:
    """
    Author: Rixin Li, 16 May 2023
    """

    def __init__(self, filename, quantities=None, **kwargs):
        ds_raw = athr.athdf(
            filename,
            raw=True,
            quantities=quantities,
        )
        self.ds_raw = ds_raw
        self.filename = filename

        if (
            quantities is None
        ):  # todo: check if input quantities are included in VariableNames
            self.quantities = [
                x.decode("ascii", "replace") for x in ds_raw["VariableNames"][:]
            ]
        elif isinstance(quantities, str):
            self.quantities = [
                quantities,
            ]
        else:
            self.quantities = quantities

        # mostly copied from athena_read.py
        self.block_size = ds_raw["MeshBlockSize"]
        self.root_grid_size = ds_raw["RootGridSize"]
        self.levels = ds_raw["Levels"][:]
        self.logical_locations = ds_raw["LogicalLocations"][:]
        self.max_level = ds_raw["MaxLevel"]

        # for now, if ndim = 2, we assume only the first two dimensions are meaningful
        self.ndim = self.get_ndim()

        # w/o knowing this, we'll running into the following error:
        # "Block boundaries at finest level must be cell boundaries at desired level for subsampling or fast restriction to work"
        block_refine_limit = np.log2(self.block_size)
        self.block_refine_limit = block_refine_limit[block_refine_limit > 0].min()

        self.restrict_data(min_lev=kwargs.get("min_lev_to_restrict", 0))

    def get_level(
        self,
        lev,
    ):
        # construct mesh data from blocks

        sel_mb_lev = np.where(self.levels == lev)[0]
        logi_locs = self.logical_locations[sel_mb_lev]
        anchor = logi_locs.min(axis=0)
        logi_locs -= anchor
        Nx_mb = self.block_size
        Nx_lev = Nx_mb * (logi_locs.max(axis=0) + 1)  # b/c locs starts from 0

        # reconstruct cell center coordinates
        ccx1, ccx2, ccx3 = (
            np.zeros(Nx_lev[0], dtype=np.float32),
            np.zeros(Nx_lev[1], dtype=np.float32),
            np.zeros(Nx_lev[2], dtype=np.float32),
        )

        if self.ndim == 2:
            Nx_lev = Nx_lev[:2]

        level_data = {}
        for _q in self.quantities:
            level_data[_q] = np.zeros(Nx_lev[::-1], dtype=np.float32)

        for idx_sel_mb, idx_mb in enumerate(sel_mb_lev):
            # print(idx_sel_mb, idx_mb)
            _ccx1, _ccx2, _ccx3 = (
                self.ds_raw["x1v"][idx_mb],
                self.ds_raw["x2v"][idx_mb],
                self.ds_raw["x3v"][idx_mb],
            )
            ccx1[
                Nx_mb[0]
                * logi_locs[idx_sel_mb][0] : Nx_mb[0]
                * (logi_locs[idx_sel_mb][0] + 1)
            ] = _ccx1
            ccx2[
                Nx_mb[1]
                * logi_locs[idx_sel_mb][1] : Nx_mb[1]
                * (logi_locs[idx_sel_mb][1] + 1)
            ] = _ccx2
            ccx3[
                Nx_mb[2]
                * logi_locs[idx_sel_mb][2] : Nx_mb[2]
                * (logi_locs[idx_sel_mb][2] + 1)
            ] = _ccx3

            for _q in self.quantities:
                if self.ndim == 2:
                    (level_data[_q])[
                        Nx_mb[1]
                        * logi_locs[idx_sel_mb][1] : Nx_mb[1]
                        * (logi_locs[idx_sel_mb][1] + 1),
                        Nx_mb[0]
                        * logi_locs[idx_sel_mb][0] : Nx_mb[0]
                        * (logi_locs[idx_sel_mb][0] + 1),
                    ] = self.ds_raw[_q][idx_mb]
                if self.ndim == 3:
                    (level_data[_q])[
                        Nx_mb[2]
                        * logi_locs[idx_sel_mb][2] : Nx_mb[2]
                        * (logi_locs[idx_sel_mb][2] + 1),
                        Nx_mb[1]
                        * logi_locs[idx_sel_mb][1] : Nx_mb[1]
                        * (logi_locs[idx_sel_mb][1] + 1),
                        Nx_mb[0]
                        * logi_locs[idx_sel_mb][0] : Nx_mb[0]
                        * (logi_locs[idx_sel_mb][0] + 1),
                    ] = self.ds_raw[_q][idx_mb]

        return [ccx1, ccx2, ccx3], level_data

    def restrict_data(self, min_lev=0):
        # restrict data level by level, from the finest level to root level

        for lev in range(self.max_level, min_lev, -1):
            logi_locs = self.logical_locations[self.levels == lev]
            logi_locs_parent = logi_locs // 2

            # to find and group fine mesh blocks that can be merged into one coarse mesh blocks
            unq, count = np.unique(logi_locs_parent, axis=0, return_counts=True)
            repeated_groups = unq[count > 1]

            re_levels = []
            re_logi_locs = []
            re_data = {
                "x1f": [],
                "x1v": [],
                "x2f": [],
                "x2v": [],
                "x3f": [],
                "x3v": [],
            }
            for _q in self.quantities:
                re_data[_q] = []

            for repeated_group in repeated_groups:
                repeated_idx = np.argwhere(
                    np.all(logi_locs_parent == repeated_group, axis=1)
                )
                # print(repeated_idx.ravel()) # one can check this so we know it is 2D or 3D

                # hard-coded for 3D (but seems to also work in 2D so far)
                idx_to_merge = np.argwhere(self.levels == lev)[
                    repeated_idx.ravel()
                ].ravel()

                # athr.athdf uses face coordinates to find the enclosure boundaries, so center-coordiantes are fine to capture the mesh blocks
                bounding_box = np.array(
                    [
                        [
                            self.ds_raw["x1v"][idx_to_merge].min(),
                            self.ds_raw["x1v"][idx_to_merge].max(),
                        ],
                        [
                            self.ds_raw["x2v"][idx_to_merge].min(),
                            self.ds_raw["x2v"][idx_to_merge].max(),
                        ],
                        [
                            self.ds_raw["x3v"][idx_to_merge].min(),
                            self.ds_raw["x3v"][idx_to_merge].max(),
                        ],
                    ]
                )

                _ds = athr.athdf(
                    self.filename,
                    level=lev - 1,
                    fast_restrict=True,
                    quantities=self.quantities,
                    max_level=min(self.max_level, lev - 1 + self.block_refine_limit),
                    x1_min=bounding_box[0][0],
                    x1_max=bounding_box[0][1],
                    x2_min=bounding_box[1][0],
                    x2_max=bounding_box[1][1],
                    x3_min=bounding_box[2][0],
                    x3_max=bounding_box[2][1],
                )

                re_levels.append(lev - 1)
                re_logi_locs.append(repeated_group)
                for _coord in [
                    "x1f",
                    "x1v",
                    "x2f",
                    "x2v",
                    "x3f",
                    "x3v",
                ]:
                    re_data[_coord].append(_ds[_coord])
                for _q in self.quantities:
                    re_data[_q].append(_ds[_q])

            self.levels = np.hstack([self.levels, np.atleast_1d(re_levels)])
            self.logical_locations = np.vstack(
                [self.logical_locations, np.array(re_logi_locs)]
            )
            for _coord in [
                "x1f",
                "x1v",
                "x2f",
                "x2v",
                "x3f",
                "x3v",
            ]:
                self.ds_raw[_coord] = np.vstack(
                    [self.ds_raw[_coord], np.array(re_data[_coord])]
                )
            for _q in self.quantities:
                self.ds_raw[_q] = np.vstack([self.ds_raw[_q], np.array(re_data[_q])])

    def get_ndim(self, num_ghost=0):
        # mostly copied from athena_read.py
        nx_vals = []
        for d in range(3):
            if self.block_size[d] == 1 and self.root_grid_size[d] > 1:  # sum or slice
                other_locations = [
                    location
                    for location in zip(
                        self.levels,
                        self.logical_locations[:, (d + 1) % 3],
                        self.logical_locations[:, (d + 2) % 3],
                    )
                ]
                if len(set(other_locations)) == len(other_locations):  # effective slice
                    nx_vals.append(1)
                else:  # nontrivial sum
                    num_blocks_this_dim = 0
                    for level_this_dim, loc_this_dim in zip(
                        self.levels, self.logical_locations[:, d]
                    ):
                        if level_this_dim <= level:
                            possible_max = (loc_this_dim + 1) * 2 ** (
                                level - level_this_dim
                            )
                            num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                        else:
                            possible_max = (loc_this_dim + 1) // 2 ** (
                                level_this_dim - level
                            )
                            num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                    nx_vals.append(num_blocks_this_dim)
            elif self.block_size[d] == 1:  # singleton dimension
                nx_vals.append(1)
            else:  # normal case
                nx_vals.append(
                    self.root_grid_size[d] * 2**self.max_level + 2 * num_ghost
                )
        self.nx1 = nx_vals[0]
        self.nx2 = nx_vals[1]
        self.nx3 = nx_vals[2]
        lx1 = self.nx1 // self.block_size[0]
        lx2 = self.nx2 // self.block_size[1]
        lx3 = self.nx3 // self.block_size[2]
        num_extended_dims = 0
        for nx in nx_vals:
            if nx > 1:
                num_extended_dims += 1

        return num_extended_dims

    def plot2d(self, var, fig, ax, nlev, vmin, vmax):
        for lev in range(nlev):
            cccs, lds = self.get_level(
                lev
            )  # meaning: cell center coordinates, level data set
            meshr, meshphi = np.meshgrid(cccs[0], cccs[1])
            im = ax.pcolormesh(
                meshr,
                meshphi,
                lds[var],
                shading="nearest",
                cmap="inferno",
                norm=colors.LogNorm(vmin,vmax),
            )


class PpdCylAthhdf5(object):
    def __init__(self, q, rootGrid, filename, outdir, quantities=None, **kwargs):
        """
        q =     planet/star mass ratio
        rootGrid =      (Nx1,Nx2) of the rood grid
        """
        self.q = q
        self.filename = filename
        self.ds = yt.load(self.filename)
        self.nx1 = rootGrid[0]
        self.nx2 = rootGrid[1]
        self.data_dict = None

        _, outname = os.path.split(self.filename)
        outname, _ = os.path.splitext(outname)
        outname = outname + ".npz"
        self.outname = os.path.join(outdir, outname)

    def load(self):
        if os.path.exists(self.outname):
            self.data_dict = np.load(self.outname)
        else:
            self.reduce()
        return self.data_dict

    def reduce(self):
        """
        Calculates "Sig", "Fx_g", "Fy_g", "Fx_g_rHexcl", "Fy_g_rHexcl", "Mdot", radial profiles.
        TODO: "Edot", "Jdot" profiles
        """
        self.radProfs = {}

        q = self.q
        rH = (q / 3) ** (1.0 / 3)

        dd = self.ds.all_data()

        dens_data = np.array(dd["athena_pp", "dens"])
        vr_data = np.array(dd["athena_pp", "mom1"]) / dens_data
        vphi_data = np.array(dd["athena_pp", "mom2"]) / dens_data
        coords = np.array(dd.fcoords)
        fwidths = np.array(dd.fwidth)

        Nr = self.nx1 # root grid
        istart = 0
        iend = Nr + 1

        edge_coords = np.linspace(0.4, 1.6, iend, endpoint=True)
        midpts = (edge_coords[:-1] + edge_coords[1:]) / 2

        mRing = np.zeros(Nr)
        mDotRing = np.zeros(Nr)
        GFxProf = np.zeros(Nr)
        GFyProf = np.zeros(Nr)
        GFx_rHexclProf = np.zeros(Nr)
        GFy_rHexclProf = np.zeros(Nr)
        vrProf = np.zeros(Nr)
        vphiProf = np.zeros(Nr) 

        for i in range(len(dens_data)):
            Sig = dens_data[i]
            vr = vr_data[i]
            vphi = vphi_data[i]
            r = coords[i, 0]
            phi = coords[i, 1]
            for j in range(istart, iend):
                if edge_coords[j] < r < edge_coords[j + 1]:
                    rsecn = np.sqrt(r**2 + 1 - 2 * r * np.cos(phi))
                    x = r * np.cos(phi)
                    y = r * np.sin(phi)

                    dr = fwidths[i, 0]
                    dphi = fwidths[i, 1]
                    dA = r * dr * dphi

                    GdFx = dA * Sig * q * (1 - x) / rsecn**3
                    GdFy = -dA * Sig * q * y / rsecn**3
                    GFxProf[j] += GdFx
                    GFyProf[j] += GdFy

                    vrProf[j] += dA*vr
                    vphiProf[j] += dA*vphi

                    dm = dA * Sig
                    dmDot = r * dphi * Sig * vr
                    mRing[j] += dm
                    mDotRing[j] += dmDot

                    if np.sqrt((x - 1.0) ** 2 + (y) ** 2) > rH:
                        GFx_rHexclProf[j] += GdFx
                        GFy_rHexclProf[j] += GdFy

        SigProf = mRing / (2 * np.pi * midpts[:] * np.diff(edge_coords)[:])
        vrProf = vrProf / (2 * np.pi * midpts[:] * np.diff(edge_coords)[:])
        vphiProf = vphiProf / (2 * np.pi * midpts[:] * np.diff(edge_coords)[:])

        np.savez(
            self.outname,
            Sig=SigProf,
            vr=vrProf,
            vphi=vphiProf,
            Fx_g=GFxProf,
            Fy_g=GFyProf,
            Fx_g_rHexcl=GFx_rHexclProf,
            Fy_g_rHexcl=GFy_rHexclProf,
            Mdot=mDotRing,
            rEdgeCoords=edge_coords,
            rMidpts=midpts,
        )
        self.data_dict = np.load(self.outname)


class athhst(object):
    def __init__(self, q, filename):
        self.q = q
        self.data = np.loadtxt(filename)

        t = self.data[:, 0] / np.pi / 2
        it = -1

        self.t = t[:it]
        self.FP_x = self.data[:it, -7]
        self.FP_y = self.data[:it, -6]
        self.Fsgrav_x = self.data[:it, -5]
        self.Fsgrav_y = self.data[:it, -4]
        self.accrate = self.data[:it, -3]
        self.momx_accrate = self.data[:it, -2]
        self.momy_accrate = self.data[:it, -1]

    @mpl.rc_context(analytic)
    def plotGravTorq(self, label, **kwargs):
        """
        kwargs
        figax=      (fig, ax)
        ybds=       (lower y, upper y)
        """

        if "figax" not in kwargs.keys():
            fig, ax = plt.subplots(figsize=(8, 5))
        else:
            fig = kwargs["figax"][0]
            ax = kwargs["figax"][1]

        if "c" in kwargs.keys():
            cl = kwargs["c"]
        else:
            cl = "k"

        totTorq = self.Fsgrav_y + self.momy_accrate + self.FP_y
        ax.plot(
            self.t,
            totTorq,
            c=cl,
            lw=2,
            ls="-",
            label=f"{label}" + r", $\mathtt{tot}$",
            zorder=-5,
        )
        ax.plot(
            self.t,
            self.Fsgrav_y,
            c=cl,
            lw=2,
            ls="--",
            label=f"{label}" + r", $\mathtt{grav}$",
            zorder=-5,
        )
        ax.plot(
            self.t,
            self.momy_accrate,
            c=cl,
            lw=2,
            ls="-.",
            label=f"{label}" + r" $\mathtt{accr}$",
            zorder=-5,
        )
        ax.plot(
            self.t,
            self.FP_y,
            c=cl,
            lw=2,
            ls=":",
            label=f"{label}" + r", $\mathtt{pres}$",
            zorder=-5,
        )

        ax.set_xlim((self.t[0], self.t[-1]))
        plotAdjustKwargs(fig, ax, **kwargs)
        ax.grid(True)
        ax.set_ylabel(r"$\Gamma/(M_pR_0^2\Omega^2)$")
        ax.set_xlabel(r"$T/(2\pi\Omega^{-1})$")

        ax.legend(bbox_to_anchor=[1.0, 1.0, 0, 0])

        fig.tight_layout()
        return (fig, ax)
