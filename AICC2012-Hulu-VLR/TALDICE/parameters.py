#Parameters specific to the TALDICE ice core.
self.age_top = -54.
self.depth = np.arange(0., 1486.+0.01, 1.)
self.corr_deporate_age = np.arange(self.age_top, 500000+self.age_top+0.01, self.age_step)
self.corr_lid_age = np.arange(self.age_top, 500000+self.age_top+0.01, self.age_step)
self.corr_thinning_depth = np.arange(self.depth[0], self.depth[-1]+0.01,\
                              (self.depth[-1]-self.depth[0])/(self.corr_thinning_nodes-1))
