import sys
import numpy as np
sys.path.insert(0, '/Users/jtlaune/athena/vis/python/')
import athena_read as athr

class Raw_Data_Restricted:
    """
        Author: Rixin Li, 16 May 2023
    """    
    def __init__(self, filename, q=None, **kwargs):
    
        ds_raw = athr.athdf(filename, raw=True, quantities=q, )
        self.ds_raw = ds_raw
        self.filename = filename
        
        if q is None:  # todo: check if input quantities are included in VariableNames
            self.quantities = [x.decode('ascii', 'replace') for x in ds_raw['VariableNames'][:]]
        elif isinstance(q, str):
            self.quantities = [q, ]
        else:
            self.quantities = q
    
        # mostly copied from athena_read.py
        self.block_size = ds_raw['MeshBlockSize']
        self.root_grid_size = ds_raw['RootGridSize']
        self.levels = ds_raw['Levels'][:]
        self.logical_locations = ds_raw['LogicalLocations'][:]
        self.max_level = ds_raw['MaxLevel']
        
        # for now, if ndim = 2, we assume only the first two dimensions are meaningful
        self.ndim = self.get_ndim()
        
        # w/o knowing this, we'll running into the following error:
        # "Block boundaries at finest level must be cell boundaries at desired level for subsampling or fast restriction to work" 
        block_refine_limit = np.log2(self.block_size)
        self.block_refine_limit = block_refine_limit[block_refine_limit>0].min()
        
        self.restrict_data(min_lev=kwargs.get("min_lev_to_restrict", 0))
        
    def get_level(self, lev, ):
        # construct mesh data from blocks
        
        sel_mb_lev = np.where(self.levels == lev)[0]
        logi_locs = self.logical_locations[sel_mb_lev]
        anchor = logi_locs.min(axis=0)
        logi_locs -= anchor
        Nx_mb = self.block_size
        Nx_lev = Nx_mb * (logi_locs.max(axis=0) + 1)  # b/c locs starts from 0
        
        # reconstruct cell center coordinates
        ccx1, ccx2, ccx3 = np.zeros(Nx_lev[0], dtype=np.float32), np.zeros(Nx_lev[1], dtype=np.float32), np.zeros(Nx_lev[2], dtype=np.float32)
        
        if self.ndim == 2:
            Nx_lev = (Nx_lev[:2])
            
        level_data = {}
        for _q in self.quantities:
            level_data[_q] = np.zeros(Nx_lev[::-1], dtype=np.float32)
        
        for idx_sel_mb, idx_mb in enumerate(sel_mb_lev):
            #print(idx_sel_mb, idx_mb)
            _ccx1, _ccx2, _ccx3 = self.ds_raw['x1v'][idx_mb], self.ds_raw['x2v'][idx_mb], self.ds_raw['x3v'][idx_mb]
            ccx1[Nx_mb[0]*logi_locs[idx_sel_mb][0]:Nx_mb[0]*(logi_locs[idx_sel_mb][0]+1)] = _ccx1
            ccx2[Nx_mb[1]*logi_locs[idx_sel_mb][1]:Nx_mb[1]*(logi_locs[idx_sel_mb][1]+1)] = _ccx2
            ccx3[Nx_mb[2]*logi_locs[idx_sel_mb][2]:Nx_mb[2]*(logi_locs[idx_sel_mb][2]+1)] = _ccx3
            
            for _q in self.quantities:
                if self.ndim == 2:
                    (level_data[_q])[Nx_mb[1]*logi_locs[idx_sel_mb][1]:Nx_mb[1]*(logi_locs[idx_sel_mb][1]+1),
                                     Nx_mb[0]*logi_locs[idx_sel_mb][0]:Nx_mb[0]*(logi_locs[idx_sel_mb][0]+1)] = self.ds_raw[_q][idx_mb]
                if self.ndim == 3:
                    (level_data[_q])[Nx_mb[2]*logi_locs[idx_sel_mb][2]:Nx_mb[2]*(logi_locs[idx_sel_mb][2]+1), 
                                     Nx_mb[1]*logi_locs[idx_sel_mb][1]:Nx_mb[1]*(logi_locs[idx_sel_mb][1]+1),
                                     Nx_mb[0]*logi_locs[idx_sel_mb][0]:Nx_mb[0]*(logi_locs[idx_sel_mb][0]+1)] = self.ds_raw[_q][idx_mb]
        
        return [ccx1, ccx2, ccx3], level_data
    
        
    def restrict_data(self, min_lev = 0):
        # restrict data level by level, from the finest level to root level
        
        for lev in range(self.max_level, min_lev, -1):
            
            logi_locs = self.logical_locations[self.levels==lev]
            logi_locs_parent = logi_locs // 2

            # to find and group fine mesh blocks that can be merged into one coarse mesh blocks
            unq, count = np.unique(logi_locs_parent, axis=0, return_counts=True)
            repeated_groups = unq[count>1]

            re_levels = []
            re_logi_locs = []
            re_data = {'x1f': [], 'x1v': [], 'x2f': [], 'x2v': [], 'x3f': [], 'x3v': [], }
            for _q in self.quantities:
                re_data[_q] = []

            for repeated_group in repeated_groups:
                repeated_idx = np.argwhere(np.all(logi_locs_parent == repeated_group, axis=1))
                #print(repeated_idx.ravel()) # one can check this so we know it is 2D or 3D

                # hard-coded for 3D (but seems to also work in 2D so far)
                idx_to_merge = np.argwhere(self.levels==lev)[repeated_idx.ravel()].ravel()

                # athr.athdf uses face coordinates to find the enclosure boundaries, so center-coordiantes are fine to capture the mesh blocks
                bounding_box = np.array([[self.ds_raw['x1v'][idx_to_merge].min(), self.ds_raw['x1v'][idx_to_merge].max()], 
                                         [self.ds_raw['x2v'][idx_to_merge].min(), self.ds_raw['x2v'][idx_to_merge].max()], 
                                         [self.ds_raw['x3v'][idx_to_merge].min(), self.ds_raw['x3v'][idx_to_merge].max()], 
                                        ])
                
                _ds = athr.athdf(self.filename, level=lev-1, fast_restrict=True, quantities=self.quantities, 
                                 max_level=min(self.max_level, lev-1 + self.block_refine_limit),
                                 x1_min=bounding_box[0][0], x1_max=bounding_box[0][1], 
                                 x2_min=bounding_box[1][0], x2_max=bounding_box[1][1], 
                                 x3_min=bounding_box[2][0], x3_max=bounding_box[2][1], )

                re_levels.append(lev-1)
                re_logi_locs.append(repeated_group)
                for _coord in ['x1f', 'x1v', 'x2f', 'x2v', 'x3f', 'x3v', ]:
                    re_data[_coord].append(_ds[_coord])
                for _q in self.quantities:
                    re_data[_q].append(_ds[_q])

            self.levels = np.hstack([self.levels, np.atleast_1d(re_levels)])
            self.logical_locations = np.vstack([self.logical_locations, np.array(re_logi_locs)])
            for _coord in ['x1f', 'x1v', 'x2f', 'x2v', 'x3f', 'x3v', ]:
                self.ds_raw[_coord] = np.vstack([self.ds_raw[_coord], np.array(re_data[_coord])])
            for _q in self.quantities:
                self.ds_raw[_q] = np.vstack([self.ds_raw[_q], np.array(re_data[_q])])
        
    def get_ndim(self, num_ghost = 0):
        
        # mostly copied from athena_read.py
        nx_vals = []
        for d in range(3):
            if self.block_size[d] == 1 and self.root_grid_size[d] > 1:  # sum or slice
                other_locations = [location
                                   for location in zip(self.levels,
                                                       self.logical_locations[:, (d+1) % 3],
                                                       self.logical_locations[:, (d+2) % 3])]
                if len(set(other_locations)) == len(other_locations):  # effective slice
                    nx_vals.append(1)
                else:  # nontrivial sum
                    num_blocks_this_dim = 0
                    for level_this_dim, loc_this_dim in zip(self.levels,
                                                            self.logical_locations[:, d]):
                        if level_this_dim <= level:
                            possible_max = (loc_this_dim+1) * 2**(level-level_this_dim)
                            num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                        else:
                            possible_max = (loc_this_dim+1) // 2**(level_this_dim-level)
                            num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                    nx_vals.append(num_blocks_this_dim)
            elif self.block_size[d] == 1:  # singleton dimension
                nx_vals.append(1)
            else:  # normal case
                nx_vals.append(self.root_grid_size[d] * 2**self.max_level + 2 * num_ghost)
        nx1 = nx_vals[0]
        nx2 = nx_vals[1]
        nx3 = nx_vals[2]
        lx1 = nx1 // self.block_size[0]
        lx2 = nx2 // self.block_size[1]
        lx3 = nx3 // self.block_size[2]
        num_extended_dims = 0
        for nx in nx_vals:
            if nx > 1:
                num_extended_dims += 1
                
        return num_extended_dims