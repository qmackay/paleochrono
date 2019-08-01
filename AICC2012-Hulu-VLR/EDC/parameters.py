#Parameters specific to the EDC ice core
self.age_top_prior = -55. #age at the top of the core
self.age_top_sigma = 10. #Error on age_top
self.depth = np.arange(0., 3259.3+0.01, 0.55) #Define the depth grid for the age calculation
#Age grid for the accu correction function
self.corr_deporate_age = np.arange(0., 1000000, self.age_step)
#Age grid for the LID correction function
self.corr_lid_age = np.arange(0., 1000000, self.age_step)
#Depth grid for the thinning correction function
self.corr_thinning_depth = np.arange(self.depth[0], self.depth[-1]+0.01,\
                                (self.depth[-1]-self.depth[0])/(self.corr_thinning_nodes-1))
