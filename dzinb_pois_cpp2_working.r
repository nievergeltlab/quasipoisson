#Adam: Pcor is now a list of matrices


library(Rcpp)
library(RcppArmadillo)
library(extraDistr)
#library(poismodel)
library(quasipoissonprobs)
#Rcpp.package.skeleton("quasipoissonprobs2",code_files=c("poismodel.cpp"))
#sourceCpp("poismodel.cpp")

quasipois <- function(p,y,xx,PCLzeromatlist,...) #take categorical covariate names to split 
 {
    ncoefp <- ncol(xx)
    coefp <- p[1:ncoefp] #get the coefficients for the count term           
    maxy=max(y)      
    L <- exp(xx %*% coefp)  #predicted count values      

    #exp of negative lambda part of likelihood
    enegL1 <- exp(-L)
    
    #This is the computation of E^-lambda * lambda^y/factorial(y)
    
   l1 <- lambdakfactorial_alt3(L,enegL1,maxy) #fastest to already multply by exp-L #this doesnt have to include maxy, can just derive from data.
  
    
   #PCLzeromat2 <- diag(41)  this is the matrix to use if trying to match a default model

   negllik <- l1*PCLzeromatlist[y + 1,]   
   #each cell is now probability of (y true ) * (y obs | y true) 

   unitvector <- matrix(1,ncol=1,nrow=(maxy+1)) #using a unit vector instead of apply provides a large speedup
   -sum(log(negllik %*% unitvector)) #the sum and log operations are only a marginal task
   #-sum(log(matprod(negllik, unitvector))) #no speed gains by sending this out to a program, this operation is efficient!
 }

   quasipois_gr <- function(p, y, xx, ww,PCLzeromatlist)
     {
        maxy=max(y) # max value of y       
        ## gradient
        ncoefp <- ncol(xx)
        coefp <- p[1:ncoefp]
        Sb <- xx %*% coefp
        #lambdas
        L <- exp(Sb)
        
       #The likelihood forms the denominator of the derivative (because d/dB ln(lik) is 1/lik * d/dB lik)
      
       enegL1 <- exp(-L)
       l1 <- lambdakfactorial_alt3(L,enegL1,maxy) 
       negllik <- l1*PCLzeromatlist[y + 1,]   
       unitvector <- matrix(1,ncol=1,nrow=(maxy+1)) 
       lik <- negllik %*% unitvector
        
       #The derivative is just the probabilities of y multiplied by lambda - y... so I just need to calculate this and then take the basic product

       #This is expected minus obs
       negllik2 <- obs_min_expected(L,maxy)

       #top part of derivative formula, not yet multipled by Xi  
       likprod <- (negllik * negllik2 ) 

       #this is the numerator of the derivative for each Bi, then divided by the denominator (total likelihood is denominator)
       nnliktot <- list (length(ncoefp))
       for (i in 1:ncoefp)
       {
        nnliktot[[i]] <- sum( (xx[,i] * likprod) %*% unitvector  /lik)
       }
               
       grad <- c(unlist(nnliktot))        
     }              

quasipoislm  <- function(parm=NULL,y = y,xx = xx, ntries = 3,method = "BFGS")
{    
   size <- nrow(xx)
   reg0 <- glm.fit(xx, y, family = poisson(),intercept=F) #initial values for count distribution estimation #intercept is turned off! we already made it

   par0 <- c(reg0$coef)  # ,
   print(par0)
   pnames <- c(paste(dimnames(xx)[[2]],"x",sep="")) # "fitval",

   distr.name <- "Quasipoisson corrected for misclassification probabilities"
   
   ######generate P matrices    
   maxy=max(y)
   PCLzeromat <- matrix(NA, nrow=maxy+1,ncol=maxy+1)
   row.names(PCLzeromat) <- as.character(c(0:maxy)) 
   colnames(PCLzeromat)  <- as.character(c(0:maxy))

   for (truescore in c(0:maxy) ) #this is score + 1, to make indexing easier
   {
    underprediction=.5 # + .1*(max(ww) - traumacount) #.3
    #print(underprediction)
    underprediction <- max(underprediction,0) #if trauma max out, no underprediction!!
   #  print(underprediction)
    #this is optimized under ml, since i have nothing informative, in a future test it will be a function of the data itself..
    #missclass_p= 0.5  +.5*(obsscore/max_ptsd_score ) 
    #I can probably do ML to do parameter choosing, but thats going to be complicated
    adj_p <- .9 #+ #.3*(traumacount/hight)
    adj_score = min(max(truescore - underprediction,0),maxy) #this is what we are likely to observe
    ev=round(adj_score/adj_p)
    #I still think these probabilities should be calibrated with bayes theorem and population prevalences
    
    #By default i am using a truncated binomial distribution.. others may be more appropriate. This is essentially normal-ish
    #this is the loss points, this is why the prediction sucked.. 
    sim_p <- dtbinom(c(0:maxy),ev,adj_p)
    #what is then the best estimate for p? and how do I make it so P varies with trauma?
    #I think that the offset is actually more varaible given low trauma, in spite of p increasing the likelihood of an offset. So it is counter 
    names(sim_p) <- as.character(c(0:maxy))
    #becuase not all possible values will  be assigned probs, add them here to make a conformable vector
    addx1 <- names(PCLzeromat)[-(which(names(PCLzeromat) %in% names(sim_p)))]
    #since they are 0 probs, just add 0 to the table 
    sim_p[addx1] <- 0
    #make sure that the order is right..
    sim_p2 <-   sim_p[order(as.numeric(names(sim_p)))]       
    #these probabilities are p(y = x | y = true score) , where x is a value in the range of true scores.
    PCLzeromat[, truescore+1] <- sim_p2           
  }

    # if no vector of initial values is given then use par0
    if(is.null(parm)) par0 <- par0 else par0 <- parm

    z0 <- optim(par0,fn=quasipois,  hessian = T, method = method, gr=quasipois_gr,
    control = list(trace=1,maxit = 200),  y=y,xx=xx,PCLzeromatlist=PCLzeromat) # gr=gr.qp, # grad=, 

    t <- 1 
    #print(z0)
    cov <- solve(z0$hessian)
    se <- sqrt(diag(cov))
    names(se) <- pnames
    estimate <- array(z0$par)
    names(estimate) <- pnames
    print(estimate)
    print(se)
    corr <- cov/(se %o% se)
    colnames(corr) <- pnames
    rownames(corr) <- pnames


    #this is for estimating the dispersion parameter for quasipoisson
    expected <- exp(xx %*% estimate)
    resid <- (y - expected) 
    dfparm <- (length(y) - length(estimate))
    resid_pearson <- resid/sqrt(expected)   
    resid1 <- sqrt(sum(resid_pearson^2)/dfparm)

    #se = se,  poisson se is useless, so i only put in the quasipoisson SE
    z1 <- list(distribution = distr.name, coefficients = estimate, se= se*resid1, dispersion = resid1^2, maxlike = z0$
        value, corr = corr, aic = round(2 * (z0$value + length(estimate))), df = 
        size - length(estimate), counts = z0$counts, convergence = z0$convergence,
         ntries = t, message = z0$message)

    class(z1) <- c("quasipoisson")
    return(z1)
}

#NOte: range is truncated to the observed maximum of Y, which is not 100% correct because there is no adjustment for truncation
#even thiough there should be, but much faster than going over the total range of values.
 
