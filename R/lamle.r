##### Class definitions
setClass("lamleout", representation(Data = "list", Estimates = "list", Model = "list", Optim = "list", Timing = "data.frame"))

###### Functions
optCR <- function(start, tol, maxit, y, J, cats, p, model, modelpars, modeltype, link, estimator, mu, invSigma){
    # SJ: This function estimates the latent variable by solving dh / d(latent variable) = 0.
    res <- optC(start, tol, maxit, y, J, cats, p, model, modelpars, modeltype, link, estimator, mu, invSigma)
    return(res)
}

DGP <- function(a, b, modeltype, z){
	if(is.vector(z)) z <- as.matrix(z)
	m <- unlist(lapply(b, length)) + 1
	Nf <- nrow(z)
	Ni <- nrow(a)
	out <- matrix(NA, nrow(z), length(b))
	linpred <- a %*% t(z)
	for(i in 1:Ni){
		if(modeltype[i] == "GPCM"){
			bi <- c(NA, b[[i]])
			propi <- matrix(NA, nrow = Nf, ncol = m[i])
			denom <- rep(1.0, Nf)
			for(k in 2:m[i]){
				tmp <- rep(0, Nf)
				for(v in 2:k) tmp <- tmp + linpred[i,] + bi[v]
				denom <- denom + exp(tmp)
			}
			for(yi in 1:m[i]){
				if(yi == 1){
					propi[, yi] <- 1 / denom
				} else{
					tmp <- rep(0, Nf)
					for(v in 2:yi) tmp <- tmp + linpred[i,] + bi[v]
					propi[, yi] <- exp(tmp) / denom
				}
			}
			for(f in 1:Nf) out[f, i] <- which(rmultinom(1, 1, propi[f,]) == 1)
		}
		if(modeltype[i] == "GRM"){
			bi <- c(NA, b[[i]])
			propi <- matrix(NA, nrow = Nf, ncol = m[i])
			propis <- matrix(0, nrow = Nf, ncol = m[i] + 1)
			propis[,1] <- 1.0			
			for(k in 2:(m[i])){
				tmp <- linpred[i,] + bi[k]
				tmp <- exp(tmp)
				propis[, k] <- tmp / (1.0 + tmp);
			}
			#print(propis)
			for(k in 1:m[i]){
				propi[, k] <- propis[, k] - propis[, k + 1]
			}
			#print(propi)
			for(f in 1:Nf) out[f, i] <- which(rmultinom(1, 1, propi[f,]) == 1)
		}
		#Should be: size = 1.0 / bi[3]
		#if(modeltype[i] == "negbin"){
		#	bi <- c(NA, b[[i]])
		#	for(f in 1:Nf) out[f, i] <- rnbinom(1, size = bi[3], mu = exp(linpred[i,f] + bi[2]))
		#}
		if(modeltype[i] == "negbin"){
			bi <- c(NA, b[[i]])
			for(f in 1:Nf) out[f, i] <- rnbinom(1, size = 1.0 / bi[3], mu = exp(linpred[i,f] + bi[2]))
		}
		if(modeltype[i] == "normal"){
			bi <- c(NA, b[[i]])
			for(f in 1:Nf) out[f, i] <- rnorm(1, mean = (linpred[i,f] + bi[2]), sd = sqrt(bi[3]))
		}
		if(modeltype[i] == "poisson"){
			bi <- c(NA, b[[i]])
			for(f in 1:Nf) out[f, i] <- rpois(1, exp(linpred[i,f] + bi[2]))
		}
	}
	return(out)
}

lamle.sim <- function(obj, seed = NULL, N = 1, z = NULL){
	p <- length(unique(unlist(obj@Model$model)))
	model <- obj@Model$model
	G <- length(unique(obj@Data$group))
	J <- length(obj@Estimates$modelpars[[1]])
	if(is.null(N)) N <- 1
	ng <- numeric(G)
	for(g in 1:G) ng[g] <- sum(obj@Data$group == (g - 1))
	ng <- ng * N
	set.seed(seed)
	if(is.null(z)){
		z <- vector("list", G)
		for(g in 1:G) z[[g]] <- rmvnorm(ng[g], obj@Estimates$mu[g,], as.matrix(obj@Estimates$covmat[g,,]))
	}
	out <- vector("list", G)
	apar <- vector("list", G)
	bpar <- vector("list", G)
	for(g in 1:G){
		aparg <- matrix(0, nrow = J, ncol = p)
		bparg <- vector("list", J)
		for(i in 1:J){
			aparg[i, model[[i]] + 1] <- obj@Estimates$modelpars[[g]][[i]][1:length(model[[i]]), 1]
			bparg[[i]] <- obj@Estimates$modelpars[[g]][[i]][-(1:length(model[[i]])), 1]
		}
		apar[[g]] <- aparg
		bpar[[g]] <- bparg
		out[[g]] <- DGP(apar[[g]], bpar[[g]], obj@Model$modeltype, z[[g]])
	}
	return(out)
}

lamle.fit <- function(obj, obs = NULL, N = NULL, z = NULL, seed = NULL, use = "complete", residmat = FALSE){
	G <- length(unique(obj@Data$group))
	if(is.null(obs)) obs <- obj@Data$y
	myseed <- set.seed(seed)
	mysim <- lamle.sim(obj, seed = myseed, N = N, z = z)
	groupout <- numeric(G)
	for(g in 1:G){
		simcorg <- cor(mysim[[g]])
		obscorg <- cor(obs[which(obj@Data$group == (g - 1)),], use = use)
		groupout[g] <- (sqrt(mean((simcorg[lower.tri(simcorg)] - obscorg[lower.tri(obscorg)])^2)))
	}
	mydata <- mysim[[1]]
	if(G != 1) for(g in 2:G) rbind(mydata, mysim[[g]])
	simcor <- cor(mydata, use = use)
	obscor <- cor(obs, use = use)
	if(residmat){
		retobj <- list(overall = sqrt(mean((simcor[lower.tri(simcor)] - obscor[lower.tri(obscor)])^2)), group = groupout, residmat = simcor - obscor)
	} else retobj <- list(overall = sqrt(mean((simcor[lower.tri(simcor)] - obscor[lower.tri(obscor)])^2)), group = groupout)
	return(retobj)
}

#Make covstruct from number of dimensions
makecov <- function(ndim, resdim = NULL){
	#resdim indicates uncorrelated latent variables
	if(!is.null(resdim)) newdim <- ndim - length(resdim) else newdim <- ndim
	outmat <- matrix(NA, nrow = newdim * (newdim - 1) / 2, ncol = 2)
	parindex <- 1
	for(j in 1:(ndim - 1)){
		for(k in j:(ndim - 1)){
			if(j %in% resdim || (k + 1) %in% resdim) next
			outmat[parindex, ] <- c(j, k + 1)
			parindex <- parindex + 1
		}
	}
	outmat - 1
}

#This function creates needed objects to accommodate parameter restrictions in the main R and C++ code (between groups, between parameters, and fixed to specified constants)
preparepar <- function(groupequal, parequal, parfix, obscat, modeltype, G, p, model, y){
	#'groupequal' is a list with group equality constraints for each type of parameter
	#'parequal' is a list with equality constraints between parameters of each type
	#'parfix' is a list with entries that indobscate with parameters should be set to specified constants (and hence not be estimated)
	#'obscat' indicates the number of categories for each variable (count and continuous data: 1)
	#'modeltype' indicates the observed variable measurement model (categorical data: "GPCM", "GRM", "NRM"; count data: "negbin", "poisson"; continuous data: "normal")
	#'G' indicates the number of groups
	#'p' indicates the number of dimensions
	#'model' is a list that indicates the dimensions that have non-zero slopes for each variable
	#'y' is the observed data matrix
	
	meany <- apply(y, 2, mean, na.rm = TRUE)
	vary <- apply(y, 2, var, na.rm = TRUE)
	
	dimy <- length(obscat)
	modelstructure <- matrix(NA, dimy, p)
	for(i in 1:dimy) modelstructure[i, model[[i]]] <- 1
	
	dimx <- 0
	npari <- numeric(dimy)
	for(i in 1:dimy){
		if(modeltype[i] %in% c("GPCM", "GRM")){
			npari[i] <- p + obscat[i] + 1
		}
		if(modeltype[i] == "NRM"){
			npari[i] <- obscat[i] * p + obscat[i] + 1
		}
		#Check
		if(modeltype[i] == "normal"){
			npari[i] <- p + obscat[i] + 1
		}
		if(modeltype[i] == "negbin"){
			npari[i] <- p + obscat[i] + 1
		}
		#Check
		if(modeltype[i] == "poisson"){
			npari[i] <- p + obscat[i] + 1
		}
	}
	#print(npari)
	npars <- G * sum(npari) + p * G + p * G + (p * (p - 1) / 2 * G)
	
	mypartable <- data.frame(parname = character(npars), group = character(npars), label = character(npars), toest = logical(npars), fixed = logical(npars), value = numeric(npars), variable = numeric(npars), type = numeric(npars))
	
	#Number of restrictions of each type
	npareq <- unlist(lapply(parequal, length))
	nparfix <- unlist(lapply(parfix, length))
	ngroupeq <- unlist(lapply(groupequal, length))
	
	#Measurement model parameters
	parindex <- 1
	for(i in 1:dimy){
		for(g in 1:G){
			#Initialize slope parameters and assign labels
			tempindex <- parindex
			if(modeltype[i] == "NRM"){
				for(k in 1:obscat[i]){
					for(j in 1:p){
						mypartable[tempindex, 1] <- paste0("G",g,"I",i,"a",j,"c",k)
						mypartable[tempindex, 2] <- paste0(g)
						mypartable[tempindex, 4] <- FALSE
						mypartable[tempindex, 5] <- FALSE
						mypartable[tempindex, 6] <- 0.0
						mypartable[tempindex, 7] <- i
						mypartable[tempindex, 8] <- "slope"
						tempindex <- tempindex + 1
					}
				}
			} else{
				for(j in 1:p){
					mypartable[tempindex, 1] <- paste0("G",g,"I",i,"a",j)
					mypartable[tempindex, 2] <- paste0(g)
					mypartable[tempindex, 4] <- FALSE
					mypartable[tempindex, 5] <- FALSE
					mypartable[tempindex, 6] <- 0.0
					mypartable[tempindex, 7] <- i
					mypartable[tempindex, 8] <- "slope"
					tempindex <- tempindex + 1
				}
			}
			
			#Set default starting values and baseline estimation configuration
			if(modeltype[i] %in% c("GRM", "GPCM")){
				for(j in 1:p){
					if(j %in% model[[i]]){
						mypartable[parindex, 3] <- paste0("G",g,"I",i,"a",j)
						mypartable[parindex, 4] <- TRUE
						mypartable[parindex, 5] <- FALSE
						mypartable[parindex, 6] <- 0.8
					} else{
						mypartable[parindex, 3] <- NA
					}
					parindex <- parindex + 1
				}
			}
			if(modeltype[i] %in% c("NRM")){
				for(k in 1:obscat[i]){
					for(j in 1:p){
						if(k == 1){
							mypartable[parindex, 3] <- NA
						} else{					
							if(j %in% model[[i]]){
								mypartable[parindex, 3] <- paste0("G",g,"I",i,"a",j,"c",k)
								mypartable[parindex, 4] <- TRUE
								mypartable[parindex, 5] <- FALSE
								mypartable[parindex, 6] <- 0.8
							} else{
								mypartable[parindex, 3] <- NA
							}
						}
						parindex <- parindex + 1
					}
				}			
			}
			if(modeltype[i] %in% c("negbin", "normal", "poisson")){
				for(j in 1:p){				
					if(j %in% model[[i]]){
						mypartable[parindex, 3] <- paste0("G",g,"I",i,"a",j)
						mypartable[parindex, 4] <- TRUE
						mypartable[parindex, 5] <- FALSE
						mypartable[parindex, 6] <- 0.5
					} else{
						mypartable[parindex, 3] <- NA
					}
					parindex <- parindex + 1
				}	
			}
			
			#Assign labels for equalities of slopes
			for(m in 1:npareq[1]){
				if(i %in% parequal[[1]][[m]]$variables){
					if(modeltype[i] %in% c("GRM", "GPCM", "negbin", "normal", "poisson")){
						for(n in 1:length(parequal[[1]][[m]]$dimension)){
							if(!is.na(modelstructure[i, parequal[[1]][[m]]$dimension[n]])){
								mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + parequal[[1]][[m]]$dimension[n]
								mypartable[mytempindex, 3] <- paste0("G",g,paste0("I",parequal[[1]][[m]]$variables, collapse = ""),"a",parequal[[1]][[m]]$dimension[1])
							}						
						}				
					}
					#NRM needs to account for different category-specific discrimination pars
					if(modeltype[i] %in% c("NRM")){
						for(n in 1:length(parequal[[1]][[m]]$dimension)){
							if(!is.na(modelstructure[i, parequal[[1]][[m]]$dimension[n]])){
								for(k in 2:obscat[i]){
									mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p * (k - 1) + parequal[[1]][[m]]$dimension[n]
									mypartable[mytempindex, 3] <- paste0("G",g,paste0("I",parequal[[1]][[m]]$variables, collapse = ""),"a",parequal[[1]][[m]]$dimension[1], "c", k)
								}
							}						
						}
					}
				
				}			
			}		
			#Define fixed slope parameters
			for(m in 1:nparfix[1]){
				if(g %in% parfix[[1]][[m]]$groups){
					if(i %in% parfix[[1]][[m]]$variables){
						if(modeltype[i] %in% c("GRM", "GPCM", "negbin", "normal", "poisson")){
							for(n in 1:length(parfix[[1]][[m]]$dimension)){
								if(!is.na(modelstructure[i, parfix[[1]][[m]]$dimension[n]])){
									mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + parfix[[1]][[m]]$dimension[n]
									mypartable[mytempindex, 4] <- FALSE
									mypartable[mytempindex, 5] <- TRUE
									mypartable[mytempindex, 6] <- parfix[[1]][[m]]$value[n]
								}						
							}				
						}
						#NRM needs to account for different category-specific discrimination pars
						if(modeltype[i] %in% c("NRM")){
							for(n in 1:length(parfix[[1]][[m]]$dimension)){
								if(!is.na(modelstructure[i, parfix[[1]][[m]]$dimension[n]])){
									for(k in 2:obscat[i]){
										mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p * (k - 1) + parfix[[1]][[m]]$dimension[n]
										mypartable[mytempindex, 4] <- FALSE
										mypartable[mytempindex, 5] <- TRUE
										mypartable[mytempindex, 6] <- parfix[[1]][[m]]$value[n]
									}
								}						
							}
						}
					}			
				}
			}
			
			#Initialize intercept parameters
			if(modeltype[i] %in% c("GRM", "GPCM", "NRM")){
				for(k in 1:obscat[i]){
					if(k == 1){
						mypartable[parindex, 1] <- paste0("G",g,"I",i,"b",k)
						mypartable[parindex, 2] <- paste0(g)
						mypartable[parindex, 3] <- NA
						mypartable[parindex, 4] <- FALSE
						mypartable[parindex, 5] <- FALSE
						mypartable[parindex, 6] <- 0.0
						mypartable[parindex, 7] <- i
						mypartable[parindex, 8] <- "intercept"
					} else{
						mypartable[parindex, 1] <- paste0("G",g,"I",i,"b",k)
						mypartable[parindex, 2] <- paste0(g)
						mypartable[parindex, 3] <- paste0("G",g,"I",i,"b",k)
						mypartable[parindex, 4] <- TRUE
						mypartable[parindex, 5] <- FALSE
						mypartable[parindex, 6] <- -seq(-(obscat[i] - 2), (obscat[i] - 2), by = 2)[k - 1]
						mypartable[parindex, 7] <- i
						mypartable[parindex, 8] <- "intercept"
					}
					parindex <- parindex + 1
				}
			}
			#Count data and continuous data only have a single intercept
			if(modeltype[i] %in% c("negbin", "normal", "poisson")){
				mypartable[parindex, 1] <- paste0("G",g,"I",i,"b",0)
				mypartable[parindex, 2] <- paste0(g)
				mypartable[parindex, 3] <- paste0("G",g,"I",i,"b",0)
				mypartable[parindex, 4] <- TRUE
				mypartable[parindex, 5] <- FALSE
				if(modeltype[i] %in% c("negbin", "poisson")) mypartable[parindex, 6] <- log(meany[i])
				if(modeltype[i] == "normal") mypartable[parindex, 6] <- meany[i]				
				mypartable[parindex, 7] <- i
				mypartable[parindex, 8] <- "intercept"
				parindex <- parindex + 1
			}
			
			#Assign labels for equalities of intercepts
			for(m in 1:npareq[2]){
				if(i %in% parequal[[2]][[m]]$variables){
					if(modeltype[i] %in% c("negbin", "normal", "poisson")){
						mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p + 1
						mypartable[mytempindex, 3] <- paste0("G",g,paste0("I", parequal[[2]][[m]]$variables, collapse = ""),"b",0)
					} else{
						for(k in 2:obscat[i]){
							if(modeltype[i] == "NRM") mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p * obscat[i] + k else mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p + k
							mypartable[mytempindex, 3] <- paste0("G",g,paste0("I", parequal[[2]][[m]]$variables, collapse = ""),"b",k)
						}
					}
				}
			}
			#Define fixed intercept parameters
			for(m in 1:nparfix[2]){
				if(g %in% parfix[[2]][[m]]$groups){
					if(i %in% parfix[[2]][[m]]$variables){
						if(modeltype[i] %in% c("negbin", "normal", "poisson")){
							mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p + 1
							mypartable[mytempindex, 4] <- FALSE
							mypartable[mytempindex, 5] <- TRUE
							mypartable[mytempindex, 6] <- parfix[[2]][[m]]$value[1]
						} else{
							for(k in 2:obscat[i]){
								if(modeltype[i] == "NRM") mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p * obscat[i] + k else mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p + k
								mypartable[mytempindex, 4] <- FALSE
								mypartable[mytempindex, 5] <- TRUE
								mypartable[mytempindex, 6] <- parfix[[2]][[m]]$value[k-1]
							}
						}
					}
				}
			}
			#Scale parameters
			mypartable[parindex, 1] <- paste0("G",g,"I",i,"phi")
			mypartable[parindex, 2] <- paste0(g)
			if(modeltype[i] %in% c("negbin", "normal")){
				mypartable[parindex, 3] <- paste0("G",g,"I",i,"phi")
				mypartable[parindex, 4] <- TRUE
				mypartable[parindex, 5] <- FALSE
				if(modeltype[i] == "negbin") mypartable[parindex, 6] <- 1.0
				if(modeltype[i] == "normal") mypartable[parindex, 6] <- vary[i]
				mypartable[parindex, 7] <- i
				mypartable[parindex, 8] <- "scale"
			}
			if(modeltype[i] %in% c("GRM", "GPCM", "NRM", "poisson")){
				mypartable[parindex, 3] <- NA
				mypartable[parindex, 4] <- FALSE
				mypartable[parindex, 5] <- FALSE
				mypartable[parindex, 6] <- 0.0
				mypartable[parindex, 7] <- i
				mypartable[parindex, 8] <- "scale"
			}
			parindex <- parindex + 1
			
			for(m in 1:npareq[3]){
				if(i %in% parequal[[3]][[m]]$variables){
					if(modeltype[i] == "NRM") mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p * obscat[i] + obscat[i] + 1 else mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p + obscat[i] + 1
					mypartable[mytempindex, 3] <- paste0("G",g,paste0("I", parequal[[2]][[m]]$variables, collapse = ""),"phi")
				}
			}
			
			for(m in 1:nparfix[3]){
				if(g %in% parfix[[3]][[m]]$groups){
					if(i %in% parfix[[3]][[m]]$variables){
						if(modeltype[i] == "NRM") mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p * obscat[i] + obscat[i] + 1 else mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p + obscat[i] + 1
						mypartable[mytempindex, 4] <- FALSE
						mypartable[mytempindex, 5] <- TRUE
						mypartable[mytempindex, 6] <- parfix[[3]][[m]]$value
					}
				}
			}
			
			#Equality restrictions between groups: Slope parameters
			for(m in 1:length(groupequal[[1]])){
				if(i %in% groupequal[[1]][[m]]$variables){
					if(g %in% groupequal[[1]][[m]]$groups){
						if(modeltype[i] == "NRM"){
							for(k in 2:obscat[i]){
								mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p * (k - 1) + groupequal[[1]][[m]]$dimension
								if(g < 10) topattern <- ".."
								if(g >= 10) topattern <- "..."
								if(g >= 100) topattern <- "...."	
								mypartable[mytempindex, 3] <- paste0(paste0("G", groupequal[[1]][[m]]$groups, collapse = ""),sub(topattern, "",mypartable[mytempindex, 3]))
							}
						} else{
							mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + groupequal[[1]][[m]]$dimension
							if(g < 10) topattern <- ".."
							if(g >= 10) topattern <- "..."
							if(g >= 100) topattern <- "...."	
							mypartable[mytempindex, 3] <- paste0(paste0("G", groupequal[[1]][[m]]$groups, collapse = ""),sub(topattern, "",mypartable[mytempindex, 3]))
						}
					}
				}
			}
			#Equality restrictions between groups: Intercept parameters
			for(m in 1:length(groupequal[[2]])){
				if(i %in% groupequal[[2]][[m]]$variables){
					if(g %in% groupequal[[2]][[m]]$groups){
						if(modeltype[i] %in% c("normal", "negbin", "poisson")){
							mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p + 1
							if(g < 10) topattern <- ".."
							if(g >= 10) topattern <- "..."
							if(g >= 100) topattern <- "...."	
							mypartable[mytempindex, 3] <- paste0(paste0("G", groupequal[[2]][[m]]$groups, collapse = ""),sub(topattern, "",mypartable[mytempindex, 3]))
						} else{
							for(k in 2:obscat[i]){
								if(modeltype[i] == "NRM") mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p * obscat[i] + k else mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p + k
								if(g < 10) topattern <- ".."
								if(g >= 10) topattern <- "..."
								if(g >= 100) topattern <- "...."	
								mypartable[mytempindex, 3] <- paste0(paste0("G", groupequal[[2]][[m]]$groups, collapse = ""),sub(topattern, "",mypartable[mytempindex, 3]))
							}
						}
					}
				}
			}
			#Equality restrictions between groups: Scale parameters
			for(m in 1:length(groupequal[[3]])){
				if(i %in% groupequal[[3]][[m]]$variables){
					if(g %in% groupequal[[3]][[m]]$groups){
						if(modeltype[i] == "NRM") mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p * obscat[i] + obscat[i] + 1 else mytempindex <- sum(npari[0:(i-1)]) * G + (g - 1) * npari[i] + p  + obscat[i] + 1
						if(g < 10) topattern <- ".."
						if(g >= 10) topattern <- "..."
						if(g >= 100) topattern <- "...."	
						mypartable[mytempindex, 3] <- paste0(paste0("G", groupequal[[3]][[m]]$groups, collapse = ""),sub(topattern, "",mypartable[mytempindex, 3]))
					}
				}
			}
		}
	}
	#print("Measurement model parameters OK")
	#Initialize mean parameters
	for(g in 1:G){
		for(j in 1:p){
			mypartable[parindex, 1] <- paste0("G",g,"mu",j)
			mypartable[parindex, 2] <- paste0(g)
			mypartable[parindex, 3] <- paste0("G",g,"mu",j)
			mypartable[parindex, 4] <- TRUE
			mypartable[parindex, 5] <- FALSE
			mypartable[parindex, 6] <- 0.0
			mypartable[parindex, 8] <- "mean"
			parindex <- parindex + 1
		}
	}
	#Initialize variance parameters
	for(g in 1:G){
		for(j in 1:p){
			mypartable[parindex, 1] <- paste0("G",g,"sigma","d",j,"d",j)
			mypartable[parindex, 2] <- paste0(g)
			mypartable[parindex, 3] <- paste0("G",g,"sigma","d",j,"d",j)
			mypartable[parindex, 4] <- TRUE
			mypartable[parindex, 5] <- FALSE
			mypartable[parindex, 6] <- 1.0
			mypartable[parindex, 8] <- "variance"
			parindex <- parindex + 1
		}
	}
	#Initialize covariance parameters
	for(g in 1:G){
		if(p > 1){
			for(j in 1:(p - 1)){
				for(k in (j + 1):p){
					mypartable[parindex, 1] <- paste0("G",g,"sigma","d",j,"d",k)
					mypartable[parindex, 2] <- paste0(g)
					mypartable[parindex, 3] <- paste0("G",g,"sigma","d",j,"d",k)
					mypartable[parindex, 4] <- TRUE
					mypartable[parindex, 5] <- FALSE
					mypartable[parindex, 6] <- 0.25
					mypartable[parindex, 8] <- "covariance"
					parindex <- parindex + 1
				}
			}
		}
	}
	#print("Initialize distribution parameters OK")
	#Assign labels for equalities of means
	for(g in 1:G){
		for(m in 1:npareq[5]){
			if(g %in% parequal[[5]][[m]]$groups){
				for(j in 1:length(parequal[[5]][[m]]$mean)){
					mytempindex <- G * sum(npari) + (g - 1) * p + parequal[[5]][[m]]$mean[j]
					mypartable[mytempindex, 3] <- paste0("G",g,"mu",parequal[[5]][[m]]$mean[1])
				}
			}
		}
	}
	#Assign fixed means
	for(g in 1:G){
		for(m in 1:nparfix[5]){
			if(g %in% parfix[[5]][[m]]$groups){
				for(j in 1:length(parfix[[5]][[m]]$mean)){
					mytempindex <- G * sum(npari) + (g - 1) * p + parfix[[5]][[m]]$mean[j]
					mypartable[mytempindex, 4] <- FALSE
					mypartable[mytempindex, 5] <- TRUE
					mypartable[mytempindex, 6] <- parfix[[5]][[m]]$value[j]
				}
			}
		}
	}
	#Set group equality constraints for means
	for(g in 1:G){
		for(m in 1:length(groupequal[[5]])){
			if(g %in% groupequal[[5]][[m]]$groups){
				for(j in 1:length(groupequal[[5]][[m]]$mean)){
					mytempindex <- G * sum(npari) + (g - 1) * p + groupequal[[5]][[m]]$mean[j]
					if(g < 10) topattern <- ".."
					if(g >= 10) topattern <- "..."
					if(g > 100) topattern <- "...."	
					mypartable[mytempindex, 3] <- paste0(paste0("G", groupequal[[5]][[m]]$groups, collapse = ""),sub(topattern, "",mypartable[mytempindex, 3]))
				}
			}
		}
	}
	#print("Mean parameters OK")
	#Assign labels for equalities of variances
	for(g in 1:G){
		for(m in 1:npareq[6]){
			if(g %in% parequal[[6]][[m]]$groups){
				for(j in 1:length(parequal[[6]][[m]]$variance)){
					mytempindex <- G * sum(npari) + G * p + (g - 1) * p + parequal[[6]][[m]]$variance[j]
					mylabelindex <- G * sum(npari) + G * p + (g - 1) * p + parequal[[6]][[m]]$variance[1]
					mypartable[mytempindex, 3] <- mypartable[mylabelindex, 3]
				}
			}
		}
	}
	#Assign fixed variances
	for(g in 1:G){
		for(m in 1:nparfix[6]){
			if(g %in% parfix[[6]][[m]]$groups){
				for(j in 1:length(parfix[[6]][[m]]$variance)){
					mytempindex <- G * sum(npari) + G * p + (g - 1) * p + parfix[[6]][[m]]$variance[j]
					mypartable[mytempindex, 4] <- FALSE
					mypartable[mytempindex, 5] <- TRUE
					mypartable[mytempindex, 6] <- parfix[[6]][[m]]$value[j]
				}
			}
		}
	}
	#Set group equality constraints for variances
	for(g in 1:G){
		for(m in 1:length(groupequal[[6]])){
			if(g %in% groupequal[[6]][[m]]$groups){
				for(j in 1:length(groupequal[[6]][[m]]$variance)){
					mytempindex <- G * sum(npari) + G * p + (g - 1) * p + groupequal[[6]][[m]]$variance[j]
					if(g < 10) topattern <- ".."
					if(g >= 10) topattern <- "..."
					if(g >= 100) topattern <- "...."	
					mypartable[mytempindex, 3] <- paste0(paste0("G", groupequal[[6]][[m]]$groups, collapse = ""),sub(topattern, "",mypartable[mytempindex, 3]))
				}
			}
		}
	}
	#print("Variance parameters OK")
	#Assign labels for equalities of covariances
	for(g in 1:G){
		for(m in 1:npareq[7]){
			if(g %in% parequal[[7]][[m]]$groups){
				for(j in 1:length(parequal[[7]][[m]]$covariance)){
					#Need to find this index location..
					mytempindex <- G * sum(npari) + G * p + G * p + (g - 1) * (p * (p - 1) / 2) + parequal[[7]][[m]]$covariance[j]
					mylabelindex <- G * sum(npari) + G * p + G * p + (g - 1) * (p * (p - 1) / 2) + parequal[[7]][[m]]$covariance[1]
					mypartable[mytempindex, 3] <- mypartable[mylabelindex, 3]
				}
			}
		}
	}
	#Assign fixed covariances
	for(g in 1:G){
		for(m in 1:nparfix[7]){
			if(g %in% parfix[[7]][[m]]$groups){
				for(j in 1:length(parfix[[7]][[m]]$covariance)){
					#Need to find this index location..
					mytempindex <- G * sum(npari) + G * p + G * p + (g - 1) * (p * (p - 1) / 2) + parfix[[7]][[m]]$covariance[j]
					mypartable[mytempindex, 4] <- FALSE
					mypartable[mytempindex, 5] <- TRUE
					mypartable[mytempindex, 6] <- parfix[[7]][[m]]$value[j]
				}
			}
		}
	}
	#Set group equality constraints for covariances
	for(g in 1:G){
		for(m in 1:length(groupequal[[7]])){
			if(g %in% groupequal[[7]][[m]]$groups){
				for(j in 1:length(groupequal[[7]][[m]]$covariance)){
					#Need to find this index location..
					mytempindex <- G * sum(npari) + G * p + G * p + (g - 1) * (p * (p - 1) / 2) + groupequal[[7]][[m]]$covariance[j]
					print(mytempindex)
					if(g < 10) topattern <- ".."
					if(g >= 10) topattern <- "..."
					if(g >= 100) topattern <- "...."	
					mypartable[mytempindex, 3] <- paste0(paste0("G", groupequal[[7]][[m]]$groups, collapse = ""),sub(topattern, "",mypartable[mytempindex, 3]))
				}
			}
		}
	}
	#print("Covariance parameters OK")
	#print("Distribution parameters OK")
	#Output is like previous lamle code but with additional parameters and with object indobscating fixed parameters
	estparameters <- (mypartable[mypartable$toest,])
	fixedparameters <- data.frame(label = mypartable[mypartable$fixed,]$label)
	estandfixedparameters <- mypartable[mypartable$fixed | mypartable$toest,]
	
	fixedparameters$variable <- numeric(length(fixedparameters$label))
	fixedparameters$group <- numeric(length(fixedparameters$label))
	fixedparameters$value <- numeric(length(fixedparameters$label))
	fixedparameters$fixedtoall <- numeric(length(fixedparameters$label))
	fixedparameters$type <- numeric(length(fixedparameters$label))
	uniqueparameters <- data.frame(label = unique(estparameters$label))
	uniqueparameters$value <- numeric(length(uniqueparameters$label))
	uniqueparameters$uniquetoall <- numeric(length(uniqueparameters$label))
	uniqueparameters$uniquetoest <- numeric(length(uniqueparameters$label))
	uniqueparameters$uniquetoestandfix <- numeric(length(uniqueparameters$label))
	uniqueparameters$type <- numeric(length(uniqueparameters$label))
	for(i in 1:nrow(uniqueparameters)){
		myindobscator1 <- which(mypartable$label == uniqueparameters$label[i])
		myindobscator2 <- which(estparameters$label == uniqueparameters$label[i])
		myindobscator3 <- which(estandfixedparameters$label == uniqueparameters$label[i])
		uniqueparameters$uniquetoall[i] <- toString(myindobscator1)
		uniqueparameters$value[i] <- mypartable[myindobscator1[1],]$value
		uniqueparameters$type[i] <- mypartable[myindobscator1[1],]$type
		uniqueparameters$uniquetoest[i] <- toString(myindobscator2)
		uniqueparameters$uniquetoestandfix[i] <- toString(myindobscator3)
	}
	estparameters$esttounique <- (match(estparameters$label, uniqueparameters$label))
	estandfixedparameters$estandfixtounique <- (match(estandfixedparameters$label, uniqueparameters$label))
	
	for(i in 1:nrow(fixedparameters)){
		myindobscator1 <- which(mypartable$label == fixedparameters$label[i])
		fixedparameters$fixedtoall[i] <- toString(myindobscator1)
		fixedparameters$value[i] <- mypartable[myindobscator1[1],]$value
		fixedparameters$variable[i] <- mypartable[myindobscator1[1],]$variable
		fixedparameters$group[i] <- mypartable[myindobscator1[1],]$group
		fixedparameters$type[i] <- mypartable[myindobscator1[1],]$type
	}

	mypartable$esttounique <- rep(NA, nrow(mypartable))
	mypartable$esttounique[match(estparameters$parname, mypartable$parname)] <- estparameters$esttounique
	mypartable$estandfixtounique <- rep(NA, nrow(mypartable))
	mypartable$estandfixtounique[match(estandfixedparameters$parname, mypartable$parname)] <- estandfixedparameters$estandfixtounique
	
	estandfixedparameters[is.na(estandfixedparameters[,9]),9] <- 0
	estandfixedparameters[,4] <- as.numeric(estandfixedparameters[,4])	
	
	return(list(allpars = mypartable, estpars = estparameters, uniquepars = uniqueparameters, fixpars = fixedparameters, estandfixpars = estandfixedparameters))
}

#Computes the gradient by calling C++ functions
#We require knowing which parameters are estimated and which are set to (potentially non-zero) constants 
#'pars' are the unique parameters
#'estfixpars' includes the estimated parameters and the fixed parameters, both of which which need to be passed to the C++ code
#'estfixpars' includes the correspondence between the unique parameters and the estimated parameters (esttounique), this is needed because parameters can be set equal and we need to identify these
#'estfixpars' entries are ordered by item-by group and then by parameter type-by group
#Generally for item parameters: slopes / intercepts / scale
#For distribution parameters: means, variances, covariances
computeGrad <- function(pars, estfixpars, y, z, J, mi, p, model, modeltype, link, N, covstruct, G, group, z.tol = 1e-8, filters = NULL, X = NULL, modeupdate = TRUE, quadpt = 1, adapt = TRUE, fullexp = TRUE, method = "ghq", npartype){
    if(method == "lap"){
        lap <- quadpt
        if(lap > 2) lap <- 2
    }
    #Make modelpars a list of lists
    modelpars <- vector("list", G)
    for(g in 1:G) modelpars[[g]] <- vector("list", J)
    parindex <- 1
	for(i in 1:J){
		npi <- length(model[[i]])
		pivec <- numeric(npi)
        for(g in 1:G){
            modelpars[[g]][[i]] <- vector("list", 2)
            if(modeltype[i] %in% c("GPCM", "GRM")){
				#Number of non-zero slopes: npi
				#Number of non-zero intercepts: mj[j] - 1
				npari <- npi + mi[i] - 1
                gindex <- parindex
                for(j in 1:npi){
					pivec[j] <- gindex
                    gindex <- gindex + 1
                }
				modelpars[[g]][[i]][[1]] <- estfixpars[pivec, 6]
                modelpars[[g]][[i]][[2]] <- c(0, estfixpars[gindex:(gindex + mi[i] - 2), 6])
            }
			if(modeltype[i] == "NRM"){
				#Number of non-zero slopes: npi * (mj[j] - 1)
				#Number of non-zero intercepts: mj[j] - 1
				npari <- npi * (mi[i] - 1) + mi[i] - 1
                gindex <- parindex
                for(j in 1:(npi * (mi[i] - 1))){
					pivec[j] <- gindex
                    gindex <- gindex + 1
                }
				modelpars[[g]][[i]][[1]] <- estfixpars[pivec, 6]
				modelpars[[g]][[i]][[2]] <- c(0, estfixpars[gindex:(gindex + mi[i] - 2), 6])
			}
			if(modeltype[i] %in% c("negbin", "normal")){
				#Number of non-zero slopes: npi
				#Number of non-zero intercepts: 1
				#Number of scale parameters: 1
				npari <- npi + 2
				gindex <- parindex
                for(j in 1:npi){
					pivec[j] <- gindex
                    gindex <- gindex + 1
                }
				modelpars[[g]][[i]][[1]] <- estfixpars[pivec, 6]
				modelpars[[g]][[i]][[2]] <- c(0, estfixpars[gindex:(gindex + 1), 6])
			}
			if(modeltype[i] %in% c("poisson")){
				#Number of non-zero slopes: npi
				#Number of non-zero intercepts: 1
				#Number of scale parameters: 0
				npari <- npi + 1
				gindex <- parindex
                for(j in 1:npi){
					pivec[j] <- gindex
                    gindex <- gindex + 1
                }
				modelpars[[g]][[i]][[1]] <- estfixpars[pivec, 6]
				modelpars[[g]][[i]][[2]] <- c(0, estfixpars[gindex, 6])
			}
			parindex <- parindex + npari
		}
    }
	
    #print(modelpars)
    mu <- array(0, dim = c(G, p))
    Sigma <- array(0, dim = c(G, p, p))
    invSigma <- array(0, dim = c(G, p, p))
    #Order is the same
    #print(modelpars)
    #print(estfixpars)
	
    for(g in 1:G){
		mu[g, ] <- estfixpars[parindex:(parindex + p - 1), 6]
		parindex <- parindex + p
    }
	#print(mu)
	#Regression model only works with single group for now
	if(!is.null(X)){
		ntotpar <- length(pars)
		nbetapar <- (ncol(X) - 1) * p
		betamat <- matrix(0, nrow = p, ncol = ncol(X))
		betaparind <- ntotpar - nbetapar
		for(tt in 1:p){
			betamat[tt, 2:ncol(X)] <- pars[(betaparind + 1):(betaparind + ncol(X) - 1)]
			betaparind <- betaparind + ncol(X) - 1
		}
		for(ss in 1:p){
			temp1 <- X[1:N, 2:ncol(X)]
			temp2 <- betamat[ss, 2:(ncol(X))]
			betamat[ss, 1] <- -mean(temp2 %*% t(temp1))
		}
	}
    
    for(g in 1:G){
		Sigma[g,,] <- diag(estfixpars[parindex:(parindex + p - 1), 6], p)
        parindex <- parindex + p
	}
	for(g in 1:G){
		if(p > 1){
			for(j in 1:(p - 1)){
				for(k in (j + 1):p){
					Sigma[g, j, k] <- estfixpars[parindex, 6]
					Sigma[g, k, j] <- estfixpars[parindex, 6]
					parindex <- parindex + 1
				}
			}
		}
    }
	
	eigencheck <- logical(G)
    for(g in 1:G) eigencheck[g] <- any(eigen(Sigma[g,,])$values < 0)
	if(method == "ghq"){
		if(any(eigencheck)) stop("The covariance matrix for the latent variables became indefinite and the algorithm had to stop. Lower the step length, try new starting values, or use a different estimation method.")
	}
	
    for(g in 1:G) invSigma[g,,] <- solve(Sigma[g,,])
    zopt <- matrix(0, N, p)
    #print(X)
	#BA 2022-07-01: No update needed here? (But update needed for optCR_String)
    if(modeupdate){
        ####  SJ: Sometimes (error occurs), optim is used to estimate latent variables.
		for(g in 1:G){
			#Can easily parallelize here
			for(n in 1:N){
				if((group[n] + 1) != g) next
				if(is.null(X)) ztemp <- try(optCR(start = z[n,], tol = z.tol, maxit = 100, y = y[n,], J = J, cats = mi, p = p, model = model, modelpars = modelpars[[g]], modeltype = modeltype, link = link, estimator = "MAP", mu = mu[g,], invSigma = as.matrix(invSigma[g,,])))
				if(!is.null(X)) ztemp <- try(optCR(start = z[n,], tol = z.tol, maxit = 100, y = y[n,], J = J, cats = mi, p = p, model = model, modelpars = modelpars[[g]], modeltype = modeltype, link = link, estimator = "MAP", mu = as.numeric(X[n,] %*% t(betamat)), invSigma = as.matrix(invSigma[g,,])))
				# SJ: I have not programmed the following...
				if("try-error" %in% class(ztemp)) ztemp <- 99
				if(ztemp[1] == 99){
					if(is.null(X)) ztemp <- optim(par = z[n,], fn = hoptC, gr = dhoptC, y = y[n,], J = J, cats = mi, p = p, model = model, modelpars = modelpars[[g]], modeltype = modeltype, link = link, mu = mu[g,], Sigma = as.matrix(Sigma[g,,]), method = "BFGS", control = list(reltol = z.tol))$par
					if(!is.null(X)) ztemp <- optim(par = z[n,], fn = hoptC, gr = dhoptC, y = y[n,], J = J, cats = mi, p = p, model = model, modelpars = modelpars[[g]], modeltype = modeltype, link = link, mu = as.numeric(X[n,] %*% t(betamat)), Sigma = as.matrix(Sigma[g,,]), method = "BFGS", control = list(reltol = z.tol))$par
				}
				zopt[n,] <- ztemp
			}
		}
        z <- zopt
    }
    #Need to pass parameter structure object to C++
	#This includes the estimated and fixed parameters, and the connection between the estimated parameters and the unique parameters
    passestfixpars <- as.matrix(estfixpars[,c(4, 6, 9)])
	#print(passestfixpars)
	#First column: Estimated (1) or not estimated (0)
	#Second column: Parameter value (estimate or constant)
	#Third column: Correspondence between unique parameters and estimated parameters
	#This is the same as previous code, except that we now also have fixed parameters to consider

    #Theta includes mu
    if(is.null(X)) X <- matrix(0, 0, 1)
    ####  SJ: Computes the log-likelihood and gradient, both for each individual.
	
	#Can parallelize here with multiple calls and combination of output
    if(method == "lap"){
        out <- mglogLGrad(pars = pars, estfixpars = passestfixpars, y = y, theta = z, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, filters = filters, X = X, approx = "Laplace", accuracy = lap, npartype = npartype)
    } else if(method == "ghq"){
        gh.xw <- fastGHQuad::gaussHermiteData(quadpt)
        templistx <- templistw <- vector("list", p)
        for(j in 1:p){
            templistx[[j]] <- gh.xw$x
            templistw[[j]] <- gh.xw$w
        }
        gh.point <- expand.grid( templistx )
        gh.weight <- expand.grid( templistw )
        gh.weight <- apply(cbind(gh.weight, exp(gh.point^2)), 1, function(x) prod(x))
        gh.point <- t(gh.point)
        
        out <- try(mglogLGrad_quad(pars = pars, estfixpars = passestfixpars, y = y, theta = z, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, filters = filters, X = X, adapt = adapt, fullexp = fullexp, quadp = gh.point, quadw = gh.weight, method = method, npartype = npartype))
		if(inherits(out, "try-error")) return(out)
    }

    out$zopt <- z
    out$modelpars <- modelpars
    out$mu <- mu
    out$sigma <- Sigma
    return(out)
}


##### lamle - Latent variable model maximum likelihood estimation using adaptive quadrature and Laplace approximations
#### Item-level measurement models supported
### Categorical data: "GPCM", "GRM"
### Count data: "negbin", "poisson"
### Continuous data: "normal"
### Change to not refer to variables as "items", since we have many types of observed variables supported now. ()
#Instead make control argument like:
#control = list(adapt = TRUE, fullexp = TRUE, z.tol = 1e-8, maxit = 500, step.length = c(0.5, 1), first.step = 25, inithess = "NCW", startval = NULL, startmode = NULL, starteval = FALSE, checkgrad = FALSE, modeupdate = TRUE, filtering = "advanced")

lamle <- function(y, model, modeltype = rep("GRM", ncol(y)), link = rep("logit", ncol(y)), mi = NULL, X = NULL, group = NULL, parequal = NULL, parfix = NULL, groupequal = NULL, resdim = NULL, sw = NULL, method = "lap", accuracy = 2, tol = 1e-4, maxit = 500, optimizer = "BFGS", obsinfo = FALSE, adapt = TRUE, fullexp = TRUE, z.tol = 1e-8, maxdiff = 0.25, step.length = c(0.5, 1), first.step = 25, inithess = "NCW", startval = NULL, startmode = NULL, starteval = FALSE, checkgrad = FALSE, modeupdate = TRUE, verbose = FALSE, filtering = "advanced"){

	#### BA: Define number of categories automatically
	#### BA: (Also add some check for coding from (1, ..., max)?)
	starttime <- proc.time()[3]
	if(is.null(mi)) mi <- apply(y, 2, function(x) length(na.omit(unique(x))))
	if(any(modeltype %in% c("normal", "negbin", "poisson"))) mi[which(modeltype %in% c("normal", "negbin", "poisson"))] <- 1
    if(method == "lap"){
        lap <- accuracy
        if(lap > 2){
            lap <- 2
            accuracy <- 2
        }
    }

    success <- FALSE
    p <- ncol(model)
    J <- nrow(model)
    N <- nrow(y)
	if(is.null(group)){
        group <- rep(0, N)
        G <- 1
    } else G <- length(unique(group))
	
	##### Check input for issues
	##### a) Consistency data matrix and model definition (x)
	##### b) Consistency theta matrix and model definition (x)
	##### c) Consistency data matrix and number of item categories (x)
	##### d) Consistency modeltype and data matrix (x)
	##### e) Consistency link and data matrix (x)
	##### f) Consistency group and data matrix (x)
	if(is.data.frame(y)) y <- as.matrix(y)
	if(J >= N) stop("There are too few cases in relation to the number of variables, which makes the model not identified. Please check input data.")
	if(J != ncol(y)) stop("Data/argument mismatch. Check consistency between number of observed variables in model definition and number of observed variables in data matrix.")
	if(!(is.null(startmode))){
		if(p != ncol(startmode)) stop("Data/argument mismatch. Argument 'startmode' must be a matrix with number of columns equal to the number of latent variables.")
		if(N != nrow(startmode)) stop("Data/argument mismatch. Argument 'startmode' must be a matrix with number of rows equal to the number of cases in the observed data.")
	}
	if(!is.null(mi)){
		compto <- mi[mi != 1]
		ncats <- apply(y[, mi != 1], 2, function(x) length(na.omit(unique(x))))
		if(any(ncats > compto)) stop("Data/argument mismatch. Check consistency between number of categories specified for each observed categorical variable and the observed data.")		
	}
	
	if(length(modeltype) != ncol(y)) stop("Data/argument mismatch. Check consistency between length of 'modeltype' and the number of observed variables in the data.")
	if(length(link) != ncol(y)) stop("Data/argument mismatch. Check consistency between length of 'link' and the number of observed variables in the data.")
	if(length(group) != N) stop("Data/argument mismatch. Check consistency between length of 'group' and the number of cases in the data.")
	
	##### End check input for issues
		
	#Fix sw (weight each gradient/log-likelihood entry by survey weight)
	if(is.null(sw)) sw <- rep(1.0, N)
	
	#Make covstruct (default if missing is all correlated)
	#Question is if different correlation structures between groups makes sense?
	if(p != 1) covstruct <- makecov(ndim = p, resdim = resdim) else covstruct <- matrix(0, 0, 0)
	
	#Make default parameter equality restrictions (nothing)
	if(is.null(parequal)){
		parequal <- vector("list", 7)
		for(i in 1:7) parequal[[i]] <- vector("list", 1)
	}
	#Make default fixed parameters (mean parameters set to 0 in group 1, variance parameters set to 1 in group 1)
	if(is.null(parfix)){
		parfix <- vector("list", 7)
		for(i in c(1, 2, 3, 4, 5, 6, 7)) parfix[[i]] <- vector("list", 1)
		parfix[[5]][[1]] <- list(mean = 1:p, value = rep(0, p), groups = 1)
		parfix[[6]][[1]] <- list(variance = 1:p, value = rep(1, p), groups = 1)
	}
    newmodel <- vector("list", J)
    for(i in 1:J) newmodel[[i]] <- which(!is.na(model[i,]))
	
	#Make default group equality restrictions (all items invariant, none for distributions)
	if(is.null(groupequal)){
		groupequal <- vector("list", 7)
		groupequal[[1]] <- vector("list", p)
		for(i in 1:p){
			groupequal[[1]][[i]] <- list(variables = 1:J, dimension = i, groups = 1:G)
		}		
		groupequal[[2]][[1]] <- list(variables = 1:J, groups = 1:G)
		for(i in 3:7) groupequal[[i]] <- vector("list", 1)
	}
	
	if(verbose) print(groupequal)
    if(!is.null(X)){
        #Starting values for regression coefficients
        betastart <- runif(p * ncol(X))
        #Add constant to matrix of covariates
        X <- cbind(1, X)
        nbetapar <- p * (ncol(X) - 1)
    } else{
        betastart <- NULL
        nbetapar <- 0
    }
	
	tospec <- preparepar(groupequal = groupequal, parequal = parequal, parfix = parfix, obscat = mi, modeltype = modeltype, G = G, p = p, model = newmodel, y = y)
	
	#Estimated and fixed parameters
	estandfixedpars <- tospec$estandfixpars
	#Unique parameters
	myuniquepars <- tospec$uniquepars
	#Parameters of each type
	npartype <- numeric(5)
	npartype[1] <- sum(myuniquepars[,6] %in%  c("slope", "intercept", "scale")); npartype[2] <- sum(myuniquepars[,6] == "regression"); npartype[3] <- sum(myuniquepars[,6] == "mean"); npartype[4] <- sum(myuniquepars[,6] == "variance"); npartype[5] <- sum(myuniquepars[,6] == "covariance")

    if(verbose) print(tospec)
    myfilters <- vector("list", 7)
    if(method == "lap" && accuracy == 2){
        message('Filtering zero and repeated entries in the likelihood function.')
        #Define filters
        if(filtering == "advanced"){
            model0 <- model
            model0[is.na(model0)] <- 0
            d3hfilter <- FindUniqComb_jlmrst_GLLVM(loadmat = model0, check_zero = TRUE)
            d4hfilter <- FindUniqComb_jlmr_GLLVM(loadmat = model0, check_zero = TRUE)
            uniqi3 <- uniqi4 <- vector("list", J)
            for(i in 1:J){
                d3hifilter <- FindUniqComb_jlmrst_GLLVM(loadmat = matrix(model0[i,], nrow = 1), check_zero = TRUE)
                d4hifilter <- FindUniqComb_jlmr_GLLVM(loadmat = matrix(model0[i,], nrow = 1), check_zero = TRUE)
                uniqi3[[i]] <- unique(c(d3hifilter$UniqComb_3rd_4[,1], d3hifilter$UniqComb_3rd_4[,2], d3hifilter$UniqComb_3rd_6[,1], d3hifilter$UniqComb_3rd_6[,2]))
                uniqi4[[i]] <- unique(d4hifilter$UniqComb[,1])
            }
            myfilters[[1]] <- uniqi3
            myfilters[[2]] <- uniqi4
            myfilters[[3]] <- d3hfilter$Uniq3
            myfilters[[4]] <- d4hfilter$Uniq4
            myfilters[[5]] <- d3hfilter$UniqComb_3rd_4
            myfilters[[6]] <- d3hfilter$UniqComb_3rd_6
            myfilters[[7]] <- d4hfilter$UniqComb
        }
    } else{
		myfilters[[1]] <- myfilters[[2]] <- vector("list", J)
		myfilters[[3]] <- myfilters[[4]] <- myfilters[[5]] <- myfilters[[6]] <- myfilters[[7]] <- matrix(0, 1, 1)
	}
    model <- lapply(newmodel, function(x) x - 1)
    
	#2021-10-06: Fixed startval input used correctly.
	#BA 2022-07-01: Make sure setup correctly.
    if(is.null(startval)){
		startval <- c(tospec$uniquepars[, 2], betastart)
	} else{
		startval[is.na(startval)] <- c(tospec$uniquepars[, 2][is.na(startval)], betastart)
        for(j in 1:(length(startval) - nbetapar)){
            estandfixedpars[as.numeric(unlist(strsplit(tospec$uniquepars[j, 5], ", "))), 6] <- startval[j]
        }
	}
    if(is.null(startmode)) startmode <- matrix(0, nrow = N, ncol = p)
    #Recode NA to 9999
    newy <- y
    for(i in 1:J) newy[is.na(y[,i]),i] <- 9999
    y <- newy
    
	if(starteval){		
        for(j in 1:(length(startval) - nbetapar)){
            estandfixedpars[as.numeric(unlist(strsplit(tospec$uniquepars[j, 5], ", "))), 6] <- startval[j]
        }
        #### SJ: I added if-else here
		#### BA (2021-10-07): fixed to 'modeltype' and 'link'
		#BA 2022-07-01: covstruct input adjust or remove?
        if( method == "lap" ){
            if(filtering == "advanced") tempobj <- computeGrad(pars = startval, estfixpars = estandfixedpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = modeupdate, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "lap", npartype = npartype)
        } else if( method == "ghq" ){
            tempobj <- try(computeGrad(pars = startval, estfixpars = estandfixedpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = modeupdate, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "ghq", npartype = npartype))			
			if(inherits(tempobj, "try-error")) return(tempobj)
        }
		return(tempobj)
	}
	
    #Check analytical gradient against numerical gradient. checkgrad = TRUE/FALSE
    if(checkgrad){
		#BA 2022-07-01: change based on uniquepars def?
        for(j in 1:(length(startval) - nbetapar)){
            estandfixedpars[as.numeric(unlist(strsplit(tospec$uniquepars[j, 5], ", "))), 6] <- startval[j]
        }
        #### SJ: I added if-else here
		#### BA (2021-10-07): fixed to 'modeltype' and 'link'
		#BA 2022-07-01: covstruct input adjust or remove?
        if( method == "lap" ){
            if(filtering == "advanced") tempobj <- computeGrad(pars = startval, estfixpars = estandfixedpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = modeupdate, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "lap", npartype = npartype)
        } else if( method == "ghq" ){
            tempobj <- try(computeGrad(pars = startval, estfixpars = estandfixedpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = modeupdate, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "ghq", npartype = npartype))
			if(inherits(tempobj, "try-error")) return(tempobj)
        }
        
        tograd <- function(x, estfixpars, y, startmode, J, mi, p, model, modeltype, N, covstruct, G, group, z.tol, filtering, filters, X = X, modeupdate, lap, npartype){
            for(j in 1:(length(x)- nbetapar)){
                estfixpars[as.numeric(unlist(strsplit(tospec$uniquepars[j,5], ", "))),6] <- x[j]
            }
            if( method == "lap" ){
                if(filtering == "advanced") tempobj <- computeGrad(pars = startval, estfixpars = estfixpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = modeupdate, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "lap", npartype = npartype)
            } else if( method == "ghq" ){
                tempobj <- computeGrad(pars = startval, estfixpars = estfixpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = modeupdate, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "ghq", npartype = npartype)
            }
            return(sum(tempobj$loglik))
        }
        res <- numDeriv::grad(tograd, startval, method = "Richardson", estfixpars = estandfixedpars, y = y, startmode = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filtering = filtering, filters = myfilters, X = X, modeupdate = modeupdate, lap = lap, npartype = npartype)	
        return(list(analytical = rowSums(tempobj$gradient), numerical = res, loglik = sum(tempobj$loglik)))
    }
    
	preptime <- proc.time()[3] - starttime
	
	message('Starting quasi-Newton method, using ', optimizer, ' with tolerance ', tol, '.')
    #### SJ: I only changed the BFGS algorithm, not the BHHH.
	#### BA: I updated BHHH.
    if(optimizer == "BFGS"){
        if( method == "lap" ){
            if(filtering == "advanced") tempobj <- computeGrad(pars = startval, estfixpars = estandfixedpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = TRUE, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "lap", npartype = npartype)
        } else if( method == "ghq" ){
            tempobj <- try(computeGrad(pars = startval, estfixpars = estandfixedpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = TRUE, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "ghq", npartype = npartype))
			if(inherits(tempobj, "try-error")) return(tempobj)
        }
        
        #print(tempobj$loglik)
        #print(rowSums(tempobj$gradient))
        #print(tempobj$thetaopt)
        x <- startval
        alpha <- step.length[1]
        dfxk0 <- newgrad <- rowSums(tempobj$gradient)
        dfxk <- dfxk0
		Id <- diag(length(x))
		#This sets the initial Hessian to the empirical cross-product matrix
		if(inithess == "crossprod"){
			Hk <- solve(tcrossprod(tempobj$gradient))
		}
		#This sets the initial Hessian to the identity matrix
		if(inithess == "identity"){
			Hk <- Id / N
		}
		#This sets the initial Hessian like p. 142 in Nocedal and Wright (2006)
		if(inithess == "NCW"){
			Hk <- Id / sqrt(sum(dfxk0^2))
		}
		        
        lltrace <- numeric(maxit)
        glltrace <- matrix(NA, maxit, length(startval))
        partrace <- matrix(NA, maxit, length(startval))
        for(iter in 1:maxit){
            xk <- x
            pk <- -Hk %*% dfxk
			if(max(abs(pk)) > maxdiff){
                pk <- maxdiff * pk / max(abs(pk))
            }
            sk <- as.numeric(alpha * pk)
            x <- xk - sk
            lltrace[iter] <- sum(tempobj$loglik)
            glltrace[iter,] <- newgrad
            partrace[iter,] <- xk
            for(j in 1:(length(x)- nbetapar)){
                estandfixedpars[as.numeric(unlist(strsplit(tospec$uniquepars[j, 5], ", "))), 6] <- x[j]
            }
            if(max(abs(x - xk)) < tol){
				success <- TRUE
				message('Iteration ',iter, ': Maximum difference is ',round(max(abs(x - xk)), 5),' and the log-likelihood is ', round(lltrace[iter], 5),'.')
				message('Algorithm stopped within tolerance ', tol, ' after ', iter, ' iterations, computing standard errors.')
				break
			}
			message('Iteration ',iter, ': Maximum difference is ',round(max(abs(x - xk)), 5),' and the log-likelihood is ', round(lltrace[iter], 5),'.')
			if(verbose) print(estandfixedpars)
            #print(rowSums(tempobj$gradient), digits = 2)
            if( method == "lap" ){
                if(filtering == "advanced") tempobj <- computeGrad(pars = x, estfixpars = estandfixedpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = TRUE, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "lap", npartype = npartype)
            } else if( method == "ghq" ){
                tempobj <- try(computeGrad(pars = x, estfixpars = estandfixedpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = TRUE, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "ghq", npartype = npartype))
				if(inherits(tempobj, "try-error")) return(tempobj)
            }
            
            #output here is named ok
			#Need to update yk based on actual direction used?
            newgrad <- rowSums(tempobj$gradient)
            dfxk <- newgrad
            yk <- as.numeric(dfxk - dfxk0)
            rhok <- as.numeric(1 / (t(yk) %*% sk))
            Hk <- (Id - rhok * sk %*% t(yk)) %*% Hk %*% (Id - rhok * yk %*% t(sk)) + rhok * sk %*% t(sk)
            dfxk0 <- dfxk
            startmode <- tempobj$zopt
            if(iter == first.step) alpha <- step.length[2]
        }
    }
    
    if(optimizer == "BHHH"){
        x <- xk <- startval
        alpha <- step.length[1]
        lltrace <- numeric(maxit)
        glltrace <- matrix(NA, maxit, length(startval))
        partrace <- matrix(NA, maxit, length(startval))
        for(iter in 1:maxit){
			if( method == "lap" ){
				if(filtering == "advanced") tempobj <- computeGrad(pars = x, estfixpars = estandfixedpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = TRUE, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "lap", npartype = npartype)
            } else if( method == "ghq" ){
                tempobj <- try(computeGrad(pars = x, estfixpars = estandfixedpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = TRUE, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "ghq", npartype = npartype))
				if(inherits(tempobj, "try-error")) return(tempobj)
            }
			
			
			# tempobj$gradient is the gradient of the approximated likelihood for each parameter (rows) and individual (columns)
			#if(!is.null(sw)){
			#
			#	f.grad <- sw %*% t(tempobj$gradient)
			#	f.grad <- rowSums(tempobj$gradient %*% diag(sw))
            #	f.hess <- tcrossprod(tempobj$gradient %*% diag(sw))
			#
			#}
			
            f.grad <- rowSums(tempobj$gradient)
            f.hess <- tcrossprod(tempobj$gradient)
			pk <- alpha * as.numeric(solve(f.hess, f.grad))
			if(max(abs(pk)) > maxdiff){
                pk <- maxdiff * pk / max(abs(pk))
            }
            x <- xk + pk
			if(!is.null(sw)) lltrace[iter] <- sum(sw * tempobj$loglik) else lltrace[iter] <- sum(tempobj$loglik)
            glltrace[iter,] <- f.grad
            partrace[iter,] <- xk
			
            for(j in 1:(length(x) - nbetapar)){
                estandfixedpars[as.numeric(unlist(strsplit(tospec$uniquepars[j, 5], ", "))), 6] <- x[j]
            }
            if(max(abs(x - xk)) < tol){
				success <- TRUE
				message('Iteration ',iter, ': Maximum difference is ',round(max(abs(x - xk)), 5),' and the log-likelihood is ', round(lltrace[iter], 5),'.')
				message('Algorithm stopped within tolerance ', tol, ' after ', iter, ' iterations, computing standard errors.')
				break
			}
			message('Iteration ',iter, ': Maximum difference is ',round(max(abs(x - xk)), 5),' and the log-likelihood is ', round(lltrace[iter], 5),'.')

            # print(rowSums(tempobj$gradient), digits = 2)
            xk <- x
            startmode <- tempobj$zopt
            if(iter == first.step) alpha <- step.length[2]
        }
    }
    esttime <- proc.time()[3] - preptime - starttime
		
	#Point of these duplicate calls? (same function called)
    if(method == "lap"){
        if(filtering == "advanced") tempobj <- computeGrad(pars = x, estfixpars = estandfixedpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = TRUE, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "lap", npartype = npartype)
    } else if(method == "ghq"){
        tempobj <- try(computeGrad(pars = x, estfixpars = estandfixedpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = TRUE, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "ghq", npartype = npartype))
		if(inherits(tempobj, "try-error")) return(tempobj)
    }
    
    #Approximated observed information matrix
    #With and without accounting for theta dependence in numerical derivatives (modeupdate = TRUE/FALSE)
    if(obsinfo && success == TRUE){
        tohess <- function(x, estfixpars, y, startmode, J, mi, p, model, modeltype, link, N, covstruct, G, group, z.tol, filtering, filters, X = X, modeupdate, accuracy, npartype){
            for(j in 1:(length(x)- nbetapar)){
                estfixpars[as.numeric(unlist(strsplit(tospec$uniquepars[j, 5], ", "))), 6] <- x[j]
            }
            if(method == "lap"){
                if(filtering == "advanced") tempobj <- computeGrad(pars = x, estfixpars = estfixpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = modeupdate, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "lap", npartype = npartype)
            } else if(method == "ghq"){
                tempobj <- computeGrad(pars = x, estfixpars = estfixpars, y = y, z = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filters = myfilters, X = X, modeupdate = modeupdate, quadpt = accuracy, adapt = adapt, fullexp = fullexp, method = "ghq", npartype = npartype)
            }
            
            return(rowSums(tempobj$gradient))
        }
        Amat <- jacobian(tohess, x, method = "simple", estfixpars = estandfixedpars, y = y, startmode = startmode, J = J, mi = mi, p = p, model = model, modeltype = modeltype, link = link, N = N, covstruct = covstruct, G = G, group = group, z.tol = z.tol, filtering = filtering, filters = myfilters, X = X, modeupdate = modeupdate, accuracy = accuracy, npartype = npartype)
    }
    Bmat <- tcrossprod(tempobj$gradient)
    if(obsinfo && success == TRUE){
        acov <- solve(-Amat)
        separ <- sqrt(diag(acov)) 
    } else{
        Amat <- NULL
        acov <- solve(Bmat)
    }
    separ <- sqrt(diag(acov))	
	setime <- proc.time()[3] - preptime - starttime - esttime
    #Put item parameters from estpars in output object
    #estpars[k,4] denotes which pars[] / separs[] entry for k-th row in estpars
    #Also give label as names
	estandfixedpars$se <- rep(NA, nrow(estandfixedpars))
	for(j in 1:(length(x) - nbetapar)){
		estandfixedpars[as.numeric(unlist(strsplit(tospec$uniquepars[j, 5], ", "))), 10] <- separ[j]
	}
	
	modelparsoptC <- tempobj$modelpars
    modelpars <- vector("list", G)
    for(g in 1:G) modelpars[[g]] <- vector("list", J)
		
    k <- 1
    for(i in 1:J){
		for(g in 1:G){
            npi <- length(model[[i]])
			if(modeltype[i] %in% c("NRM")){
				npari <- npi * (mi[i] - 1) + mi[i] - 1
			}
            if(modeltype[i] %in% c("GPCM","GRM")){
				npari <- npi + mi[i] - 1
			}
			if(modeltype[i] %in% c("negbin", "normal")){
				npari <- npi + mi[i] + 1
			}
			if(modeltype[i] %in% c("poisson")){
				npari <- npi + mi[i]
			}
			modelpars[[g]][[i]] <- data.frame(par = rep(NA, npari), se = rep(NA, npari), label = rep(NA, npari))
			modelpars[[g]][[i]]$par <- estandfixedpars[k:(k + npari - 1), 6]
			modelpars[[g]][[i]]$se <- estandfixedpars[k:(k + npari - 1), 10]
			modelpars[[g]][[i]]$label <- estandfixedpars[k:(k + npari - 1), 3]
			k <- k + npari
        }
    }

	##### Output to-do-list:
	##### 1. Hessian matrices for latent variables
	##### 2. Function call
	##### 3. Fit measures
	##### 4. Settings used
	##### 5. Reliability measures
	
	loadg <- vector("list", G)
	for(g in 1:G){
		tmpload <- matrix(0, nrow = J, ncol = p)
		for(i in 1:J){
			npi <- length(model[[i]])
			if(modeltype[i] != "NRM") tmpload[i,model[[i]] + 1] <- modelpars[[g]][[i]]$par[1:npi]
		}
		loadg[[g]] <- tmpload
	}

	Data <- list(y = y,
		group = group, 
		N = N, 
		m = mi)

	Estimates <- list(par = x, 
		separ = separ, 
		acov = acov, 
		Amat = Amat, 
		Bmat = Bmat, 
		loadings = loadg, 
		covmat = tempobj$sigma, 
		mu = tempobj$mu, 
		partable = estandfixedpars, 
		modelpars = modelpars, 
		modelparsoptC = modelparsoptC, 
		map = tempobj$zopt)

	Model <- list(model = model, 
		covstruct = covstruct, 
		groupequal = groupequal, 
		parequal = parequal, 
		parfix = parfix, 
		modeltype = modeltype,
		filters = myfilters, 
		method = method, 
		accuracy = accuracy)
		
	Optim <- list(loglik = sum(tempobj$loglik), 
		grad = rowSums(tempobj$gradient), 
		logliki = tempobj$loglik, 
		gradi = tempobj$gradient, 
		AIC = 2 * length(x) - 2 * sum(tempobj$loglik), 
		BIC = log(N) * length(x) - 2 * sum(tempobj$loglik), 
		lltrace = lltrace[1:iter], 
		glltrace = glltrace[1:iter,], 
		partrace = partrace[1:iter,], 
		iter = iter,
		call = match.call())

	Timing <- data.frame(prep = as.numeric(preptime), 
		est = as.numeric(esttime), 
		se = as.numeric(setime))
		
    #res <- list(par = x, separ = separ, acov = acov, Amat = Amat, Bmat = Bmat, loglik = sum(tempobj$loglik), grad = rowSums(tempobj$gradient), logliki = tempobj$loglik, gradi = tempobj$gradient, modelpars = modelpars, modelparsoptC = modelparsoptC, loadings = loadg, covmat = tempobj$sigma, mu = tempobj$mu, partable = estandfixedpars, map = tempobj$thetaopt, group = group, N = N, model = model, covstruct = covstruct, groupequal = groupequal, parequal = parequal, parfix = parfix, modeltype = modeltype, m = mi, AIC = 2 * length(x) - 2 * sum(tempobj$loglik), BIC = log(N) * length(x) - 2 * sum(tempobj$loglik), lltrace = lltrace[1:iter], glltrace = glltrace[1:iter,], partrace = partrace[1:iter,], iter = iter, filters = myfilters, method = method, accuracy = accuracy, timing = data.frame(prep = as.numeric(preptime), est = as.numeric(esttime), se = as.numeric(setime)))
	res <- new("lamleout", Data = Data, Estimates = Estimates, Model = Model, Optim = Optim, Timing = Timing)
    return(res)
}

#Conditional probabilities, supports 1-PL, 2-PL, 3-PL, GPCM, GRM
Pi <- function(x, cats, model, z){
	JX <- length(cats)
	probs <- vector("list", JX)
	for(i in 1:JX){
		if(model[i] %in% c("1-PL", "2-PL", "3-PL", "1PL", "2PL", "3PL")){
			out <- matrix(0, nrow = cats[i], ncol = length(z))
			apar <- x[[i]][1]
			bpar <- x[[i]][2]
			if(model[i] %in% c("1-PL", "2-PL", "1PL", "2PL")) cpar <- 0 else cpar <- exp(x[[i]][3]) / (1.0 + exp(x[[i]][3]))
			out <- matrix(0, nrow = cats[i], ncol = length(z))
			out[1, ] <- 1.0 - (cpar + (1.0 - cpar) / (1 + exp(-apar * z - bpar)))
			out[2, ]  <- (cpar + (1.0 - cpar) / (1 + exp(-apar * z - bpar)))
			probs[[i]] <- out
		}
		if(model[i] == "GPCM"){
			out <- matrix(0, nrow = cats[i], ncol = length(z))
			apar <- x[[i]][1]
			bpar <- x[[i]][-1]
			denom <- 0
			for(j in 1:(cats[i] - 1)){
				temp <- 0
				for(l in 1:j) temp <- apar * z + bpar[l] + temp 
				denom <- exp(temp) + denom
			}
			out[1, ] <- 1 / (1 + denom)
			for(j in 1:(cats[i] - 1)){
				numer <- exp(j * apar * z + sum(bpar[1:j]))
				out[j + 1, ] <- numer / (1 + denom)
			}
			probs[[i]] <- out			
		}
		if(model[i] == "GRM"){
			out <- matrix(0, nrow = cats[i], ncol = length(z))
			apar <- x[[i]][1]
			bpar <- x[[i]][-1]
			out[1,] <- 1 - 1 / (1 + exp(-apar * z - bpar[1]))
			out[cats[i], ] <- 1 / (1 + exp(-apar * z - bpar[cats[i] - 1]))
			if(cats[i] > 2) for(j in 2:(cats[i] - 1)) out[j, ] <- 1 / (1 + exp(-apar * z - bpar[j-1])) - 1 / (1 + exp(-apar * z - bpar[j]))
			probs[[i]] <- out
		}
		if(model[i] == "NRM") stop("NRM not implemented yet.")
	}
	return(probs)
}

#Derivatives of the conditional probabilities wrt theta
dPidz <- function(x, cats, model, z){ 
  JX <- length(cats)
  dPidz <- vector("list", JX)
  probs <- vector("list", length(cats))
  probs <- Pi(x = x, cats = cats, model = model, z = z)
  for (i in 1:JX){
    if(model[i] %in% c("1-PL", "2-PL", "3-PL", "1PL", "2PL", "3PL")){
      apar <- x[[i]][1]
      bpar <- x[[i]][2]
      if(model[i] %in% c("1-PL", "2-PL", "1PL", "2PL")) cpar <- 0 else cpar <- exp(x[[i]][3]) / (1.0 + exp(x[[i]][3]))
      res <-  matrix(NA, nrow = cats[i], ncol = length(z))
      res[1, ] <- - apar[i] * (1 - probs[[i]][2, ]) * (probs[[i]][2, ] - cpar[i]) / (1 - cpar[i])
      res[2, ] <- apar[i] * (1 - probs[[i]][2, ]) * (probs[[i]][2, ] - cpar[i]) / (1 - cpar[i])
      dPidz[[i]] <- res
    }

    sumkprobs <- numeric(length(z))
    if (model[i] == "GPCM"){
      apar <- x[[i]][1]
      bpar <- x[[i]][-1]
      res <- matrix(NA, nrow = cats[i], ncol = length(z))
      sumkprobs <- 0
      for (k in 1:cats[i]){ # for k in number of categories per item j
        sumkprobs <- k * probs[[i]][k, ] + sumkprobs 
      }
      for (k in 1:cats[i]){
        res[k, ] <- probs[[i]][k, ] * apar[i] * (k - sumkprobs)
      }
      dPidz[[i]] <- res
    }
    if (model[i] == "GRM"){
	  Pstar <- matrix(0, nrow = cats[i] + 1, ncol = length(z))
	  Pstar[1, ] <- 1.0
	  for(k in 1:cats[i]){
		Pstar[k + 1, ] <- Pstar[k, ] - probs[[i]][k, ]
	  }
	  apar <- x[[i]][1]
      bpar <- x[[i]][-1]
      res <- matrix(NA, nrow = cats[i], ncol = length(z))
      for (k in 1:(cats[i])) res[k, ] <- (apar[i] * Pstar[k, ] * (1.0 - Pstar[k, ])) - (apar[i] * Pstar[k + 1, ] * (1.0 - Pstar[k + 1, ])) 
      dPidz[[i]] <- res
    }
  }
  return(dPidz)
} 

#Item information function
Infoi <- function(alphai, z, model, catsi){
	info <- matrix(1e-15, nrow = length(z), ncol = catsi, dimnames = list(z, 1:catsi))
	if(model %in% c("1-PL", "2-PL", "GPCM")){
		P1 <- Pi(list(alphai), catsi, model, z)
		tempRF1 <- numeric(length(z))
		for(k in 1:catsi) tempRF1 <- k * P1[[1]][k,] + tempRF1
		for(k in 1:catsi) info[,k] <- P1[[1]][k,] * alphai[1]^2 * (k - tempRF1)^2
	}
	if(model == "GRM"){
		Pz <- Pi(list(alphai), catsi, model, z)[[1]]
		Pstar <- matrix(0, nrow = catsi + 1, ncol = length(z))
		Pstar[1, ] <- 1.0
		for(k in 1:catsi) Pstar[k + 1, ] <- Pstar[k, ] - Pz[k, ]
		dPdz <- dPidz(list(alphai), catsi, model, z)[[1]]
		d2Pdz2 <- matrix(0, nrow = catsi, ncol = length(z))
		for(k in 1:catsi) d2Pdz2[k, ] <- alphai[1]^2 * Pstar[k, ] * (1.0 - Pstar[k, ]) * (1.0 - 2.0 * Pstar[k, ]) - alphai[1]^2 * Pstar[k + 1, ] * (1.0 - Pstar[k + 1, ]) * (1.0 - 2.0 * Pstar[k + 1, ])
		for(k in 1:catsi) info[,k] <- (dPdz[k, ]^2 / Pz[k, ] - d2Pdz2[k, ])
	}
	if(model == "normal"){
		info[,1] <- alphai[1]^2 / alphai[3]
	}
	if(model == "poisson"){
		linpred <- alphai[2] + alphai[1] * z
		info[,1] <- alphai[1]^2 * exp(linpred)
	}
	if(model == "negbin"){
		linpred <- alphai[2] + alphai[1] * z
		explinpred <- exp(linpred)
		info[,1] <- alphai[1]^2 * (alphai[3] * explinpred + 1.0) * explinpred / (1.0 + explinpred)^2
	}
	return(info)
}

lamle.info <- function(x, variables = NULL, group = NULL, z = NULL, eachcat = FALSE){
	if(is.null(z)) z <- seq(-6, 6, by = 0.01)
	if(is.null(group)) group <- 1
	if(any(apply(x@Estimates$loadings[[1]], 1, function(x) sum(x != 0)) > 1)) stop("Function only supports unidimensional or independent-clusters multidimensional models.")
	if(is.null(variables)) variables <- 1:length(x@Data$m)
	out <- vector("list", length(variables))
	names(out) <- variables
	k <- 1
	for(i in variables){
		#if(x$modeltype[i] %in% c("poisson", "negbin", "normal")) stop("Function only supports models for ordinal data (not 'normal', 'poisson', or 'negbin'.")
		out[[k]] <- Infoi(x@Estimates$modelpars[[group]][[i]][,1], z, x@Model$modeltype[i], x@Data$m[i])
		k <- k + 1
	}
	if(eachcat) return(out) else return(sapply(out, function(x) rowSums(x)))
}

#Output density/probability mass function for specified points w/ continuous and count data? (Question is if someone wants this?)
lamle.prob <- function(x, variables = NULL, group = NULL, z = NULL){
	if(is.null(z)) z <- seq(-6, 6, by = 0.01)
	if(is.null(group)) group <- 1
	if(is.null(variables)) variables <- 1:length(x@Data$m)
	if(any(apply(x@Estimates$loadings[[1]], 1, function(x) sum(x != 0)) > 1)) stop("Function supports only unidimensional or independent-clusters multidimensional models.")
	out <- vector("list", length(variables))
	names(out) <- variables
	k <- 1
	for(i in variables){
		if(x@Model$modeltype[i] %in% c("poisson", "negbin", "normal")) stop("Function only supports models for ordinal data (not 'normal', 'poisson', or 'negbin'.")
		tempmat <- t(Pi(list(x@Estimates$modelpars[[group]][[i]][,1]), x@Data$m[i], x@Model$modeltype[i], z)[[1]])
		dimnames(tempmat) <- list(z, 1:x@Data$m[i])
		out[[k]] <- tempmat
		k <- k + 1
	}
	out
}

#Visible function that should be documented
lamle.compute <- function(x, compute = "prob", ...){
	#Add possible alternative scoring vector (new input to lamle)
	#Expected scores for Poisson and Negative-binomial: exp(linpred)
	#Expected scores for Normal: linpred
	if(compute %in% c("prob", "probabilities", "probtrace", "trace", "icc", "ICC", "irf", "IRF", "icrf", "ICRF")){
		out <- lamle.prob(x, ...)
		return(out)
	} else if(compute %in% c("itemscore", "expectedscore", "testscore", "tcc", "TCC")){
		probs <- lamle.prob(x, ...)
		variables <- as.numeric(names(probs))
		itemscores <- matrix(0, nrow = nrow(probs[[1]]), ncol = length(variables))
		for(i in 1:length(variables)){
			scorei <- 1:x@Data$m[variables[i]]
			for(j in scorei){
				itemscores[,i] <- probs[[i]][,j] * j + itemscores[,i]
			}
		}
		if(compute %in% c("itemscore", "expectedscore")){
			colnames(itemscores) <- variables
			rownames(itemscores) <- rownames(probs[[1]])
			return(itemscores)
		} else if(compute %in% c("testscore", "tcc", "TCC")){
			testscores <- rowSums(itemscores)
			names(testscores) <- rownames(probs[[1]])
			return(testscores)
		}
	} 
	else if(compute %in% c("info", "information", "iteminfo", "iif", "IIF")){
		out <- lamle.info(x, ...)
		return(out)
	} else if(compute %in% c("testinfo", "testinformation", "tif", "TIF")){
		out <- lamle.info(x, ...)
		return(rowSums(out))
	} else if(compute %in% c("srmsr", "SRMSR")){
		return(lamle.fit(obj = x, ...))		
	} else return("Requested output not recognized.")
}

#Visible function that should be documented
lamle.plot <- function(x, toplot = "prob", palette = "Dark 3", ...){
	#Restore graphical settings upon exit.
	oldpar <- par(no.readonly = TRUE) 
	on.exit(par(oldpar))
	if(toplot %in% c("prob", "probabilities", "probtrace", "trace", "icc", "ICC", "irf", "IRF", "icrf", "ICRF")){
		#Add support for single item only
		out <- lamle.compute(x, compute = toplot, ...)
		xplot <- as.numeric(rownames(out[[1]]))
		colour <- hcl.colors(ncol(out[[1]]), palette = palette)
		plot(x = xplot, y = out[[1]][,1], type = "l", lty = 1, main = paste0("Variable ", names(out)[1]), xlab = "Latent variable value", ylab = "Probability", xlim = range(xplot), ylim = c(0, 1), lwd = 2.5, col = colour[1])
		for(k in 2:ncol(out[[1]])){
			par(new = TRUE)
			plot(x = xplot, y = out[[1]][,k], type = "l", lty = k, main = "", xlab = "", ylab = "", xlim = range(xplot), ylim = c(0, 1), lwd = 2.5, col = colour[k], axes = FALSE)
		}
	} else if(toplot %in% c("itemscore", "expectedscore", "testscore", "tcc", "TCC")){
		out <- lamle.compute(x, compute = toplot, ...)
		if(toplot %in% c("itemscore", "expectedscore")){
			xplot <- as.numeric(rownames(out))
			colour <- hcl.colors(ncol(out), palette = palette)
			plot(x = xplot, y = out[,1], type = "l", lty = 1, main = "", xlab = "Latent variable value", ylab = "Expected score", xlim = range(xplot), ylim = c(range(out)), lwd = 2.5, col = colour[1])
			if(ncol(out) > 1){
				for(k in 2:ncol(out)){
					par(new = TRUE)
					plot(x = xplot, y = out[,k], type = "l", lty = k, main = "", xlab = "", ylab = "", xlim = range(xplot), ylim = c(range(out)), lwd = 2.5, col = colour[k], axes = FALSE)
				}
			}
		} else if(toplot %in% c("testscore", "tcc", "TCC")){
			xplot <- as.numeric(names(out))
			colour <- hcl.colors(1, palette = palette)
			plot(x = xplot, y = out, type = "l", lty = 1, main = "", xlab = "Latent variable value", ylab = "Test characteristic curve", xlim = range(xplot), ylim = c(range(out)), lwd = 2.5, col = colour[1])		
		}
	} else if(toplot %in% c("info", "information", "iteminfo", "iif", "IIF")){
		#Add support for single or multiple items
		out <- lamle.compute(x, compute = toplot, ...)
		xplot <- as.numeric(rownames(out))
		colour <- hcl.colors(ncol(out), palette = palette)
		plot(x = xplot, y = out[,1], type = "l", lty = 1, main = "", xlab = "Latent variable value", ylab = "Information", xlim = range(xplot), ylim = c(0, max(out)), lwd = 2.5, col = colour[1])
		if(ncol(out) > 1){
			for(k in 2:ncol(out)){
				par(new = TRUE)
				plot(x = xplot, y = out[,k], type = "l", lty = k, main = "", xlab = "", ylab = "", xlim = range(xplot), ylim = c(0, max(out)), lwd = 2.5, col = colour[k], axes = FALSE)
			}
		}
	} else if(toplot %in% c("testinfo", "testinformation", "tif", "TIF")){
		#Add support for single or multiple groups
		out <- lamle.compute(x, compute = toplot, ...)
		xplot <- as.numeric(names(out))
		colour <- hcl.colors(1, palette = palette)
		plot(x = xplot, y = out, type = "l", lty = 1, main = "", xlab = "Latent variable value", ylab = "Test information", xlim = range(xplot), ylim = c(0, max(out)), lwd = 2.5, col = colour[1])		
	} else return("Requested output not recognized.")
}

#Visible function that should be documented
#Instead add to lamle.compute()? (Need additional input.)
lamle.predict <- function(x, obs, estimator = "MAP", information = NULL, z.tol = 1e-8, group = NULL){
	#Take parameter estimates from lamle object x and estimate the latent variable based on the observed data in obs.
	#Support MAP and MLE (can add quickly)
	#Add WLE and EAP?
	y <- obs
	J <- ncol(obs)
	mi <- x@Data$m
	p <- ncol(x@Estimates$loadings[[1]])
	N <- nrow(obs)
	if(is.null(group)) group <- rep(0, N)
	G <- length(unique(x@Data$group))
	model <- x@Model$model
	modelparsoptC <- x@Estimates$modelparsoptC
	modeltype <- x@Model$modeltype
	link <- rep("logit", ncol(obs))
	mu <- x@Estimates$mu
	Sigma <- x@Estimates$covmat
	if(nrow(x@Estimates$map) == nrow(obs)) z <- x@Estimates$map else z <- matrix(0, N, p)
	invSigma <- array(0, dim = c(G, p, p))
	for(g in 1:G) invSigma[g,,] <- solve(Sigma[g,,])
	zopt <- matrix(0, N, p)
	for(g in 1:G){
		for(n in 1:N){
			if((group[n] + 1) != g) next
			ztemp <- try(optCR(start = z[n,], tol = z.tol, maxit = 100, y = y[n,], J = J, cats = mi, p = p, model = model, modelpars = modelparsoptC[[g]], modeltype = modeltype, link = link, estimator = estimator, mu = mu[g,], invSigma = as.matrix(invSigma[g,,])))
			if("try-error" %in% class(ztemp)) ztemp <- 99
			if(ztemp[1] == 99){
				if(estimator == "MAP") ztemp <- optim(par = z[n,], fn = hoptC, gr = dhoptC, y = y[n,], J = J, cats = mi, p = p, model = model, modelpars = modelparsoptC[[g]], modeltype = modeltype, link = link, mu = mu[g,], Sigma = as.matrix(Sigma[g,,]), method = "BFGS", control = list(reltol = z.tol))$par
				if(estimator == "MLE") ztemp <- optim(par = z[n,], fn = mleoptC, gr = dmleoptC, y = y[n,], J = J, cats = mi, p = p, model = model, modelpars = modelparsoptC[[g]], modeltype = modeltype, link = link, method = "BFGS", control = list(reltol = z.tol))$par
			}
			zopt[n,] <- ztemp
		}
	}
	zvcov <- NULL
	if(!is.null(information)){
		zvcov <- vector("list", N)
		if(!(information %in% c("expected", "observed"))) stop("Argument information must be either 'expected' or 'observed'.")
		for(g in 1:G){
			for(n in 1:N){
				if((group[n] + 1) != g) next
				zvcovtemp <- try(infoC(start = zopt[n,], tol = z.tol, maxit = 100, y = y[n,], J = J, cats = mi, p = p, model = model, modelpars = modelparsoptC[[g]], modeltype = modeltype, link = link, estimator = estimator, information = information, mu = mu[g,], invSigma = as.matrix(invSigma[g,,])))
				zvcov[[n]] <- zvcovtemp				
			}
		}
	}
	return(list(zopt, zvcov))
}



