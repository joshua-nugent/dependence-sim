
############## EDITED for SEARCH-IPT by Josh Nugent
# - get.weights function removed
# - some changes to get.IC.variance()
# - coded for one-sided p-values
# Automatically re-scales Y values to not exceed 1


##############
# Stage2_Functions.R 
# R code to implement all Stage2 analyses to compare intervention effects between arms
#
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Lead Statistician for SEARCH
#
# Input: Stage1 estimated cluster-level outcome & cluster-level covariates/exposure 
# Outputs: an estimate of the intervention effect
# Goal option specifies if arithmetic risk ratio (aRR) or risk difference (RD)
# Weighting option specifies to weight clusters or individuals equally when computing
#		the intervention effect. 
##---------------------------
#################



#-----------------------------------------------------#-----------------------------------------------------
# Stage2: function to estimate the intervention effect
#
# input: 
#	  goal;  aRR= arithmetic risk ratio; otherwise risk difference (RD)
#	  weighting; "clust" weights clusters equally; otherwise weights indv equally,
#	  observed data (data.input),
#	  name of column corresponding to the outcome (outcome),
#	  survival: indicator prob survival inputted but outcome is prob of badness
#		  currently applicable to TB/mortality
#	  set of candidate adjustment variables; must be at cluster-level (clust.adj)
#	  prespecified adjustment variables for estimation of the conditional mean outcome with main terms (e.g. Y~A+W1+W2) (QAdj), 
#	  prespecified adjustment variables for estimation of the propensity score with main terms (e.g. A~W1+W2) (gAdj),
#	  indicator to do adaptive prespecification (do.data.adapt),
#	  indicator to break the pair-match (break.match)
#	  indicator to print updates (verbose)
# output: point estimate and inference
#
# note the data.input must have
#		1) a column 'id' indicating cluster membership 
#		2) a column 'U' as dummy adjustment variable (=1 for all)
#		both of these requirements are taken care with preprocess function
#
# note2: if you are running adaptive prespecification to select the candidate C-TMLE 
#   minimizing the variance (which nearly all Stage 2 analyses do), 
#   please set do.data.adapt=T and specify the candidate adjustment variables in clust.adj... 
#   Please do NOT use (QAdj, gAdj) which are reserved for running non-adaptive glms only. 
#
# Note3: when specifying clust.adj, please include the dummy variable U 
#   to make sure that the unadjusted estimator is considered as a candidate. 
#   e.g. clust.adj <- c(‘U’, ‘hiv_prev_0’, ‘mc_cover_0’)
#
# note4:to run the unadjusted estimator than set clust.adj=NULL and do.data.adapt=F.
#
# EXAMPLES
## Stage2 estimation of the arithmetic risk ratio (ARR), weighting clusters equally, & preserving the matches
#
## First with TMLE with adaptive pre-specification
# clust.adj <- c('U', 'chc_prev', 'bmi')
# Yc1.adj <- Stage2(goal=‘aRR’, weighting='clust’, data.input=pri, 
#                   outcome='Yc', clust.adj=clust.adj,  do.data.adapt=T,  break.match=F) 
#                   
## Then with the unadjusted estimator
#  Yc1.unadj <- Stage2(goal=‘aRR’, weighting=‘cluster-level’, data.input=pri, 
#               outcome= 'YcU', clust.adj=NULL, do.data.adapt=F,  break.match=F)
#-----------------------------------------------------#-----------------------------------------------------
Stage2_cluster <- function(goal='aRR', weighting,# 'Yc_ind' or 'Yc_clust'
                   data.input, outcome = NULL,
                   clust.adj, QAdj = NULL, gAdj = NULL, verbose = F){	

  ###### TRANSFORM data.input here........................ scale by max value and min value
	outcome.string <- make.full.word(outcome)
	names(data.input)[grep(outcome.string, names(data.input) ) ] <- 'Y'
	
	if(max(data.input[,'Y']) > 1){
	  scale_value <- max(data.input[,'Y'])
	} else {
	  scale_value <- 1
	}
	if(min(data.input[,'Y']) < 0){
	  scale_value_min <- min(data.input[,'Y'])
	} else {
	  scale_value_min <- 0
	}
	data.input[,'Y'] <- (data.input[,'Y'] - scale_value_min) / (scale_value - scale_value_min)
	#print(b)
	#print(scale_value)
	#print(scale_value_min)
	
	nIndv.string <- make.full.word( paste('nIndv', outcome, sep='_') )
	names(data.input)[grep(nIndv.string, names(data.input) ) ] <- 'nIndv'

	
	# Run full TMLE algorithm with prespecified adjustment set
	# Output 'est' has been unscaled within the 'do.TMLE' function
	est <- do.TMLE(goal = goal, weighting = weighting, train = data.input,  QAdj = QAdj,  
			 gAdj = gAdj, verbose = verbose,
			 scale_value = scale_value, scale_value_min = scale_value_min)
		
	# Get point estimates and inference
	R1 = est$R1
	R0 = est$R0
	
	n.clust <- length(unique(data.input$id)) 
	n.pair <- length(unique(data.input$pair))
	
	df <- n.clust - 2
	#print("Getting inference...")
	inference <- get.inference(goal = goal, R1 = R1, R0 = R0, var.IC = est$var.break, df = df)
	
	Txt <- get.CI(psi.hat = R1, var.IC = est$var.R1, df = (n.clust - 2))
	Con <- get.CI(psi.hat = R0, var.IC = est$var.R0, df = (n.clust - 2))

	return(data.frame(Txt=Txt, Con=Con, inference,  QAdj=est$QAdj, gAdj=est$gAdj))
}


#-----------------------------------------------------#-----------------------------------------------------
# do.TMLE: master function to run TMLE 
#
# input: 
#		goal - aRR: arithmetic risk ratio; O/W risk difference,
#		training data (train),
#		weighting; "clust" weights clusters equally; O/W weights indv equally,
#		prespecified adjustment variables for the conditional mean outcome (QAdj), 
#		prespecified adjustment variables for the propensity score (gAdj),
#		initial estimator of the conditional mean outcome (Q.out),
#		estimator of the propensity score (p.out),
#		indicator to print updates (verbose)
#
# output: list 
#		training data augmented with estimates,
#		prespecified adjustment variables for the conditional mean outcome (QAdj), 
#		prespecified adjustment variables for the propensity score (gAdj),
#		initial estimator of the conditional mean outcome (Q.out),
#		estimator of the propensity score (p.out),
#		estimted fluctutation coef (epsilon),
#		estimated risk under txt and control (R1, R0),
#		estimated variance preserving the match (var.pair),
# 		estimated variance breaking the match (var.break)
#-----------------------------------------------------#-----------------------------------------------------

do.TMLE <- function(goal, train, weighting, 
                    QAdj, gAdj=NULL, Q.out=NULL, p.out=NULL,
                    verbose=F, scale_value = 1, scale_value_min = 0){	
	
	J <- length(unique(train$id))
	n <- nrow(train)

	################ adding weights
	#print("adding weights... ")
	if(weighting == 'Yc_mean'){
	  #print("... Yc_mean: Cluster-level values are weighted by each cluster mean")
	  train$alpha  <- 1 / train$N_j
	} else if (weighting == 'Yc_ind'){
	  #print("... Yc_ind: Cluster-level values are weighted to make each individual have equal weight regardless of cluster")
	  train$alpha  <- J/n #train$N_j * J / n
	  #print(sum(train$alpha))
	} else {
	  print("please specify a weighting scheme for the Yc values")
	  break
	}
	# do they sum to J? Yes.
	#print(sum(train$alpha))
	#print(J)
	
	#=====================================================
	# Step1 - initial estimation of E(Y|A,W)= Qbar(A,W)
	#=====================================================
	# run glm on the adjustment set
	Q <- do.Init.Qbar(train=train, QAdj = QAdj, glm.out = Q.out, verbose = verbose)
	train <- Q$train

	#==========================================================
	# Step2: Calculate the clever covariate
	#==========================================================	
  #print("Calculating the clever covariate...")
	G <- get.clever.cov(train=train, gAdj = gAdj, p.out = p.out, verbose = verbose)
  #print(G)
	train <- G$train
	
	#==========================================================
	# Step3: Targeting
	#==========================================================
	#print("Getting epsilon...")
	eps <- get.epsilon(train=train, goal=goal, verbose=verbose)
	#print("doing targeting; updating the predictions...")
	train <- do.targeting(train=train, eps=eps, goal=goal)
	
	#==========================================================
	# Step4: Parameter and variance estimation
	#==========================================================
	#print("getting the parameter estimates and variance...")
	variance.out <- get.IC.variance(goal=goal, Vdata=train, #R1=R1, R0=R0,
	                                weighting = weighting,
	                                scale_value = scale_value, scale_value_min = scale_value_min)
	
	RETURN <- list(train=train, 	
		QAdj = Q$QAdj, Q.out = Q$glm.out,
		gAdj = G$gAdj, p.out = G$p.out, 
		eps = eps, 
		R1 = variance.out$R1,
		R0 = variance.out$R0,
		var.R1 = variance.out$var.R1, 
		var.R0 = variance.out$var.R0,
		var.break = variance.out$var.break)	
	RETURN
}

#-----------------------------------------------------#-----------------------------------------------------
# do.Init.Qbar - function to do initial estimation of E[Y|A,W] = Qbar(A,W)
# 	input: 
#		data set, adjustment variable(s), outcome regression fit on training set, verbose
# 	output: 
#		adjustment variable(s)
#		outcome regression fit on training set
#		data set augmented w/ initial predictions: Qbar(A,W), Qbar(1,W) and Qbar(0,W)
#-----------------------------------------------------#-----------------------------------------------------
do.Init.Qbar <- function(train, QAdj, glm.out = NULL, verbose = F){
	
	train.temp <- train[, c(QAdj, 'A', 'Y')]
	X1 <- X0 <- train.temp
	X1$A <-1; X0$A <- 0	
	if(is.null(glm.out)){
		# fit using the training data. Confirmed with LB that weights are needed, though that is surprising to me
	  # why binomial? If we didn't scale could we use LM?
		glm.out<- suppressWarnings(glm(Y~ ., family='binomial', data=train.temp, weights = train$alpha))
		if(verbose) print(glm.out)
	}	
	
	# get initial predictions
	QbarAW <- predict(glm.out, newdata=train.temp, type='response')
	Qbar1W <- predict(glm.out, newdata=X1, type='response')
	Qbar0W <- predict(glm.out, newdata=X0, type='response')
	
	list(QAdj=QAdj, glm.out=glm.out, train=data.frame(train, QbarAW ,Qbar1W ,Qbar0W))
}

#-----------------------------------------------------#-----------------------------------------------------
# get.clever.cov - function calculate the clever covariate
# 	input: 
#		data set, adjustment variable(s), pscore regression fit on training set, verbose
# 	output: 
#		adjustment variable(s), pscore regression fit on training set, 
#		data set augmented with pscore & clever covariate (H.AW, H.1W, H.0W)
#-----------------------------------------------------#-----------------------------------------------------
get.clever.cov<- function(train, gAdj, p.out = NULL, verbose = F){
	if(is.null(gAdj)){
		gAdj <- 'U'
	}
	train.temp <- train[, c(gAdj, 'A')]  
	if( is.null(p.out) ){
		# fit pscore on training set 	
		p.out <- suppressWarnings(glm(A ~ ., family='binomial', data = train.temp, weights = train$alpha))
		if(verbose){print(p.out)}
	}
	# now use p.out to get estimated pscores
	pscore <- predict(p.out, newdata= train.temp,  type="response")
	# Note if gAdj=U & train$alpha !=1, then this will differ from 0.5

	# bound g - should not apply for a randomized trial
	pscore [pscore < 0.025] <- 0.025
	pscore [pscore > 0.975] <- 0.975

	A.train <- train$A
	# Clever covariate is two-dimensional; 
	H.1W <- A.train/pscore 
	H.0W <- (1-A.train)/(1-pscore)
	# via delta method
	H.AW <- H.1W - H.0W

	list(gAdj=gAdj, p.out=p.out,  train=data.frame(train, pscore, H.1W , H.0W , H.AW) ) 
}	
	
	
#-----------------------------------------------------#-----------------------------------------------------
# get.epsilon - function calculate the fluctuation coefficient
# 	input: 
#		data set, goal with 'aRR'=arithmetic RR, verbose
# 	output: 
#		estimated fluctuation coefficient (eps)
#-----------------------------------------------------#-----------------------------------------------------
get.epsilon <- function(train, goal, verbose=F){
	
	A.train<- train$A
	Y.train<- train$Y
		
	# Skip fitting if outcome=0 for all observations in either txt or control group
	Skip.update <- mean(Y.train[A.train==1])==0 | mean(Y.train[A.train==0])==0 |  
		mean(Y.train[A.train==1])==1 | mean(Y.train[A.train==0])==1 
			
	if(goal!='aRR'){
		# if not going after aRR, then use a 1-dim clever covariate
		if(!Skip.update){
			logitUpdate<- suppressWarnings( 
				glm(Y.train ~ -1 +offset(qlogis(train$QbarAW )) + train$H.AW, family="binomial",  weights=train$alpha))
			eps<-logitUpdate$coef
		} else{
			eps<- 0
		}
		names(eps) <- 'H.AW'
	} else {
		# targeting the aRR requires a two-dimensional clever covariate
		if(!Skip.update){
			logitUpdate<- suppressWarnings(
				glm(Y.train ~ -1 +offset(qlogis(train$QbarAW )) + train$H.0W + train$H.1W, family="binomial", weights=train$alpha))
			eps<-logitUpdate$coef
		} else{
			eps <- c(0,0)
		}
		names(eps)<- c('H.0W', 'H.1W')	
	}
	if(verbose) print(eps)
	return(eps)
}

#-----------------------------------------------------#-----------------------------------------------------
# do.targeting - function to update initial estimators of QbarAW
# 	input: 
#		data set (train), fluctuation coefficient (eps), goal (aRR= arithmetic risk ratio; otherwise RD)
# 	output: 
#		data.frame w/ targeted predictions: Qbar*(A,W), Qbar*(1,W), Qbar*(0,W)
#-----------------------------------------------------#-----------------------------------------------------

do.targeting <- function(train, eps, goal){
	g1W <- train$pscore
	g0W <- (1 - g1W)
	if(goal != 'aRR'){
		# updated QbarAW estimates for training set. 
		QbarAW.star <- plogis( qlogis(train$QbarAW ) + eps*train$H.AW)	
		Qbar1W.star <- plogis( qlogis(train$Qbar1W ) + eps/g1W )
		Qbar0W.star <- plogis( qlogis(train$Qbar0W ) - eps/g0W )
	} else {
		# updated QbarAW estimates for training set. 
	  print("Josh: You are going for risk ratio for some reason. this is bad.")
		QbarAW.star <- plogis( qlogis(train$QbarAW) + eps['H.0W']*train$H.0W + eps['H.1W']*train$H.1W)	
		Qbar0W.star <- plogis( qlogis(train$Qbar0W) + eps['H.0W']/g0W )
		Qbar1W.star <- plogis( qlogis(train$Qbar1W) + eps['H.1W']/g1W )
	}
	train <- data.frame(train, QbarAW.star, Qbar1W.star, Qbar0W.star)		
	return(train)
}
	

#-----------------------------------------------------#-----------------------------------------------------
# get.IC.variance - function to do influence curve-based variance estimate 
# 	input: 
#		goal (aRR= arithmetic risk ratio; otherwise RD)
#		dataset
#		risk estimates under txt and control R1 & R0
# 	output: on log scale for if goal='aRR'
#		estimated IC & variance - preserving/breaking the match
#-----------------------------------------------------#-----------------------------------------------------
get.IC.variance <- function(goal, Vdata, R1, R0,
                            scale_value = 1, scale_value_min = 0, weighting){
	
  # aggregated data: the weights are 1 for cluster-effect & n_j*J/nTot for individual-level effect
  # (JN adds: weights are J/nTot, but then they are summed in clusters, making it equiv-ish to above)
  # note weights are normalized to sum to J		
  
  #print('AGGREGATING ALL DATA')
  J <- length(unique(Vdata$id))
  n <- nrow(Vdata)
  
  if(weighting == 'Yc_mean' & n > J){ # need to aggregate PEs and IC separately...
    # IF DATA ARE AT THE INDV-LEVEL, BUT GOAL IS THE CLUSTER-LEVEL EFFECT 
    # get point estimates by aggregating to the cluster level 
    #   (by taking the weighted sum)
    # then take mean of cluster-level endpoints
    #print("Using cluster-level means, unweighted...")
    #Vdata <- aggregate(Vdata, by=list(Vdata$id), mean)[,-1]
    #Vdata$alpha  <- 1
    
    R1 <- mean( aggregate(data.frame(Vdata$alpha*Vdata$Qbar1W.star), by=list(Vdata$id), sum)[,2] )
    R0 <- mean( aggregate(data.frame(Vdata$alpha*Vdata$Qbar0W.star), by=list(Vdata$id), sum)[,2] )
    #print("Getting mean predicted values...")
    #OLD:R1 <- mean(Vdata$alpha * Vdata$Qbar1W.star)
    #OLD:R0 <- mean(Vdata$alpha * Vdata$Qbar0W.star)
    
    #print("Unscaling the results...")
    R1 <- R1 * (scale_value - scale_value_min) + scale_value_min
    R0 <- R0 * (scale_value - scale_value_min) + scale_value_min
    
    #print("Getting IC and unscaling it...")
    DY1 <- Vdata$alpha * Vdata$H.1W *
      (Vdata$Y - Vdata$Qbar1W.star) * (scale_value - scale_value_min) + scale_value_min
    DY0 <- Vdata$alpha * Vdata$H.0W *
      (Vdata$Y - Vdata$Qbar0W.star) * (scale_value - scale_value_min) + scale_value_min
    DY1 <- aggregate(DY1, by=list(Vdata$id), sum)[,-1]
    DY0 <- aggregate(DY0, by=list(Vdata$id), sum)[,-1]
    #OLD:DY1 <- aggregate(x = DY1, by = list(id = Vdata$id), mean)[,2]
    #OLD:DY0 <- aggregate(x = DY0, by = list(id = Vdata$id), mean)[,2]
  } else if (weighting == 'Yc_ind'){
    #print(names(Vdata))
    #print(Vdata)
    #print("Using cluster-level values, weighted by the number of individuals in the cluster...")
    #print("Getting cluster-level weighted predicted values...")
    Rs <- Vdata %>% 
      group_by(id) %>% summarise(R1 = sum(Qbar1W.star*alpha),
                                 R0 = sum(Qbar0W.star*alpha))
    R1 <- sum(Rs$R1) / J # was     #R1 <- mean(Vdata$alpha * Vdata$Qbar1W.star)
    R0 <- sum(Rs$R0) / J # was     #R0 <- mean(Vdata$alpha * Vdata$Qbar0W.star)
    #print("Unscaling the results...")
    R1 <- R1 * (scale_value - scale_value_min) + scale_value_min
    R0 <- R0 * (scale_value - scale_value_min) + scale_value_min
    #print(paste(R1, R0))

    #print("Getting individual-level IC estimates and unscaling them...")
    #print(paste("Sum of weights:",sum(Vdata$alpha)))
    DY1 <-# Vdata$alpha * this happens later at the aggregation level
      Vdata$H.1W *
      (Vdata$Y - Vdata$Qbar1W.star) * (scale_value - scale_value_min) + scale_value_min
    DY0 <-# Vdata$alpha * this happens later at the aggregation level
      Vdata$H.0W *
      (Vdata$Y - Vdata$Qbar0W.star) * (scale_value - scale_value_min) + scale_value_min
    # print(paste("mean DY1: ", mean(DY1)))
    # print(paste("mean DY0: ", mean(DY0)))
    #print("Aggregating ICs accounting for individual-level weighting")
    agg_Ds <- cbind.data.frame(DY1 = DY1, DY0 = DY0, alpha = Vdata$alpha, id = Vdata$id) %>% 
      group_by(id) %>% summarise(D1 = sum(DY1 * alpha),
                                 D0 = sum(DY0 * alpha))
    DY1 <- agg_Ds$D1 # was: #DY1 <- c(ltmle:::HouseholdIC(as.matrix(DY1), id = Vdata$id))
    DY0 <- agg_Ds$D0 # was: #DY0 <- c(ltmle:::HouseholdIC(as.matrix(DY0), id = Vdata$id))
  }
  ########################################## something with nIndv later?
  
	if(goal !='aRR'){
		# going after RD, easy IC
		DY <-  DY1 - DY0
		print(mean(DY))
		print(var(DY))
		#print(paste("mean DY: ", mean(DY)))
	} else { 
		# going after aRR, then get IC estimate on log scale;	i.e. Delta method for log(aRR) = log(R1) - log(R0)
	  print("Josh: Beware: apparently you are going for risk ratio. you don't want this.")
		DY <- 1/R1*DY1 - 1/R0*DY0
	}
	
	# estimated variance for txt specific means
	var.R1 <- var(DY1) / J
	var.R0 <- var(DY0) / J
	var.break <- var(DY) / J

	return(list(R1 = R1, R0 = R0, var.R1=var.R1, var.R0=var.R0, DY=DY, var.break=var.break))
}


#-----------------------------------------------------#-----------------------------------------------------
# get.CI: mini function to calculate two-sided (1-sig.level)% confidence intervals 
#	input: 
# 		pt estimate
# 		var.IC (variance estimate)
#		df (degrees of freedom for Student's t-dist ) 
#		sig.level (significance level)
# output: 
#		data.frame with point est, lower/higher 95% CI
#-----------------------------------------------------#-----------------------------------------------------	


get.CI <- function(psi.hat,  var.IC, df, sig.level=0.05){
	
	# cutoff based on t-dist for testing and CI	
	cutoff <- qt(sig.level/2, df = df, lower.tail = F)
	
	# standard error (square root of the variance)
	se <- sqrt(var.IC)
	# 95% confidence interval 
	CI.lo <- (psi.hat - cutoff*se)
	CI.hi <- (psi.hat + cutoff*se)
  out <- data.frame(est = psi.hat, CI.lo, CI.hi, var.hat = var.IC)
	out
}

#-----------------------------------------------------#-----------------------------------------------------
# get.Inference: function to calculate (1-sig.level)% confidence intervals & test the null hypothesis
#	input: 
#		goal (aRR= arithmetic risk ratio; otherwise RD)
#		risk estimates on txt and control (R1, R0)
# 		var (variance estimate - on log scale for aRR), 
#		df (degrees of freedom if using a Student's t-dist ) 
#		sig.level (significance level)
# output: 
#		variance, test statistic, confidence intervals, pval, indicator reject null
# 		note: if goal=aRR, variance & test stat are on log-scale
#-----------------------------------------------------#-----------------------------------------------------	

get.inference <- function(goal, R1, R0, var.IC, df, sig.level=0.05,
                          one_sided = T, alternate_null_negative = F){
	
	if( goal=='aRR' ){
		psi.hat <- log(R1/R0)
	} else{
		psi.hat <- R1 - R0	
	}
  # standard error (square root of the variance)
  se <- sqrt(var.IC)
  # test statistic (if goal=aRR then on the transformed scale)
  tstat <- psi.hat/se
  
  #print(tstat)
  #print(df)
  # cutoff based on t-dist for testing and CI	
  if(one_sided){
    pval <- pt(tstat, df = df, lower.tail = alternate_null_negative) 
  } else {
    pval <- 2*pt(abs(tstat), df=df, lower.tail = F)
  }
  # reject the null
  reject <- pval < sig.level 
  cutoff_2_sided <- qt(sig.level/2, df = df, lower.tail=F)
  
  CI.lo <- (psi.hat - cutoff_2_sided*se)
  CI.hi <- (psi.hat + cutoff_2_sided*se)
  
 	if(goal=='aRR'){
  			psi.hat<- exp(psi.hat)
  			CI.lo <- exp(CI.lo)
  			CI.hi <- exp(CI.hi)
  	}  
 	
 	out<- data.frame(Effect.est=psi.hat, Effect.CI.lo=CI.lo, Effect.CI.hi=CI.hi, var.hat=var.IC, tstat, pval, reject)
  	out

}




############
make.full.word <- function(string){
	paste("\\<", string, "\\>", sep="")
}







