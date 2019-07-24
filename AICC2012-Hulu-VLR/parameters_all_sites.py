self.calc_deporate = False           #Use False for now.
self.calc_tau = False         #Use False for now.
self.calc_lid = False         #Use False for now.
self.corr_tau_nodes = 26  #Define the number of nodes of the thinning function.
self.age_step = 20000.	#Define the age step for the LID and accu correction functions.
self.dens_firn = 0.7            #Average density of the firn
self.start = 'default'  #default, restart or random

#Parameters needed to define the covariance matrices as in AICC2012.
self.lambda_tau = 70
self.lambda_deporate = 4000
self.lambda_lid = 4000