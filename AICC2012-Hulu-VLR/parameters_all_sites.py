self.calc_deporate = False           #Use False for now.
self.calc_thinning = False         #Use False for now.
self.calc_lid = False         #Use False for now.
self.corr_thinning_nodes = 51  #Define the number of nodes of the thinning function.
self.age_step = 10000.	#Define the age step for the LID and accu correction functions.
self.dens_firn = 0.7            #Average density of the firn
self.start = 'prior'  #prior, restart or random

#Parameters needed to define the covariance matrices as in AICC2012.
self.lambda_thinning = 70
self.lambda_deporate = 4000
self.lambda_lid = 4000