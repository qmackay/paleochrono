#Parameters specific to the NGRIP ice core
self.udepth_top = 3.313966267738
self.age_top = -30.
self.depth = np.arange(8., 3084.+0.01, 1.)
self.corr_deporate_age = np.arange(self.age_top, 150000+self.age_top+0.01, self.age_step)
self.corr_lid_age = np.arange(self.age_top, 150000+self.age_top+0.01, self.age_step)
self.corr_thinning_depth = np.arange(self.depth[0], self.depth[-1]+0.01,\
                                (self.depth[-1]-self.depth[0])/(self.corr_tau_nodes-1))