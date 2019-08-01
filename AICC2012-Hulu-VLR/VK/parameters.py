#Parameters specific to the Vostok ice core
self.age_top_prior = -50.
self.age_top_sigma = 10. #Error on age_top
self.depth = np.arange(0., 3501+0.01, 1.)
self.corr_deporate_age = np.arange(0., 450000, self.age_step)
self.corr_lid_age = np.arange(0., 450000, self.age_step)
self.corr_thinning_depth = np.arange(self.depth[0], self.depth[-1]+0.01,\
                              (self.depth[-1]-self.depth[0])/(self.corr_thinning_nodes-1))
