#Parameters specific to the EDML ice core
self.age_top = 70.
self.depth = np.arange(18., 2564.+0.01, 1.)
self.corr_deporate_age = np.arange(self.age_top, 300000+self.age_top+0.01, self.age_step)
self.corr_lid_age = np.arange(self.age_top, 300000+self.age_top+0.01, self.age_step)
self.corr_thinning_depth = np.arange(self.depth[0], self.depth[-1]+0.01,\
                                (self.depth[-1]-self.depth[0])/(self.corr_thinning_nodes-1))
