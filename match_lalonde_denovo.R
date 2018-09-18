############
# Code to create a matched dataset for use with the denovo package using the lalonde dataset
# FOr illustrative purposes, not enough data for meaningful analysis
# Authors: Kwonsang Lee and Ben Sabath
############

unmatched=read.csv("lalonde.csv")
out_filename <- "lalonde_matched.csv"

## data description
# treatment: an indicator variable for treatment status (whether a subject received job training programs or not).
# age: age in years.
# education: years of schooling.
# black: indicator variable for blacks.
# hispanic: indicator variable for Hispanics.
# married: indicator variable for martial status.
# nodegree: indicator variable for high school diploma.
# re74: real earnings in 1974.
# re75: real earnings in 1975.
# re78: real earnings in 1978. -> outcome 


#### propensity score estimation
outcome_var <- "re78"
exact_vars<- c("black","hispanic","married","nodegree") # Variables that observations should be exactly matched on
matching_vars <- c("age","education") # Variables used to match within the sets of exactly matched variables
treatment_var <- "treatment"
match_formula <- as.formula(paste0(treatment_var, " ~ ",paste(exact_vars, collapse="+"),"+", 
                                   paste(matching_vars, collapse="+")))

propscore.model=glm(match_formula, family=binomial, data=unmatched)
logit.ps=predict(propscore.model, newdata=unmatched)

#### check overlap
hist(logit.ps[unmatched[[treatment_var]]==0], col=rgb(1,0,0,0.2))
hist(logit.ps[unmatched[[treatment_var]]==1], col=rgb(0,0,1,0.2), add=T)

unmatched$logit.ps=logit.ps ## add a variable of estimated propensity scores


#####
# Function for computing 
# rank based Mahalanobis distance.  Prevents an outlier from
# inflating the variance for a variable, thereby decreasing its importance.
# Also, the variances are not permitted to decrease as ties 
# become more common, so that, for example, it is not more important
# to match on a rare binary variable than on a common binary variable
# z is a vector, length(z)=n, with z=1 for treated, z=0 for control
# X is a matrix with n rows containing variables in the distance

smahal=
  function(z,X){
    X<-as.matrix(X)
    n<-dim(X)[1]
    rownames(X)<-1:n
    k<-dim(X)[2]
    m<-sum(z, na.rm=T)
    for (j in 1:k) X[,j]<-rank(X[,j])
    cv<-cov(X)
    vuntied<-var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    message("here")
    out<-matrix(NA,m,n-m)
    message("here")
    Xc<-X[z==0,]
    Xt<-X[z==1,]
    message(m)
    message(dim(out))
    rownames(out)<-rownames(X)[z==1]
    colnames(out)<-rownames(X)[z==0]
    library(MASS)
    icov<-ginv(cv)
    for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
    out
  }

# Function for adding a propensity score caliper to a distance matrix dmat
# calipersd is the caliper in terms of standard deviation of the logit propensity scoe
addcaliper=function(dmat,z,logitp,calipersd=.5,penalty=1000){
  sd.logitp=sd(logitp)
  adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat=dmat+adif*penalty
  dmat
}

exact_vals <- list()
## Create data.frame of all permutations of unique exact match values
all_cases <- data.frame(init = 0)
for (exact in exact_vars) {
  message(exact)
  unique_vals <- data.frame(unique(unmatched[[exact]]))
  names(unique_vals) <- exact
  all_cases <- merge(all_cases, unique_vals , all = T) ## data.frame created to allow for all possible case merge
}
all_cases$init <- NULL

print(names(all_cases))

for (i in 1:nrow(all_cases)) {
  print(i)
  ## boolean vector holding all cases where the exact matches are present
  subset_index <- unmatched[[exact_vars[1]]] == all_cases[[exact_vars[1]]][i]
  for(j in 2:length(exact_vars)) {
    subset_index <- subset_index & (unmatched[[exact_vars[j]]] == all_cases[[exact_vars[j]]][i])
  }
  ## Get row numbers where all variables being used for exact matcing fit the permutation
  subset_index <- which(subset_index)
  print(all_cases[i,])
  if (length(subset_index) == 0) {
    next
  }
  subset <- unmatched[subset_index,]
  
  if ((sum(subset[[treatment_var]] == 1, na.rm=T) <= 1) || (sum(subset[[treatment_var]] == 0, na.rm=T) <= 1 )) {
    next
  }

  ###################################################
  ### Matching
  Xmat=subset[, matching_vars]
             
  # Rank based Mahalanobis distance
  treatment=subset[[treatment_var]]
  
  ## remove observations where treatment is unknown
  Xmat <- Xmat[!is.na(treatment),]
  subset <- subset[!is.na(treatment),]
  treatment <- treatment[!is.na(treatment)]
  distmat=smahal(treatment,Xmat)
  message("Distance!")

  # Add caliper (give a penalty if two propensity scores are too different)
  distmat2=addcaliper(distmat,treatment, subset$logit.ps)
  message("Distance 2")

  ### Create a subject index and name the rows and columns of distance matrix by ### this subject index
  subject.index=seq(1,length(treatment),1)
  rownames(distmat2)=subject.index[treatment==1]
  colnames(distmat2)=subject.index[treatment==0]

  treated.indices=subject.index[treatment==1]
  control.indices=subject.index[treatment==0]


  # Pair Matching
  library(optmatch)
  options("optmatch_max_problem_size" = Inf)
  matchvec=pairmatch(distmat2)
  message("We matched fam!")

  # Create vectors of the subject indices of the treatment units ordered by
  # their matched set and corresponding control unit
  subject.indices=as.numeric(names(matchvec))
  treated.subject.index=rep(0,sum(treatment==1))
  matched.control.subject.index=rep(0,length(treated.subject.index))
  included.already=rep(0,length(matchvec))
  count=1
  for(i in 1:length(matchvec)){
    if(included.already[i]==0 & !is.na(matchvec[i])){
      temp.indices=which(matchvec==matchvec[i])
      if(treatment[subject.indices[temp.indices[1]]]==1){
        treated.subject.index[count]=subject.indices[temp.indices[1]]
        matched.control.subject.index[count]=subject.indices[temp.indices[2]]
      }
      if(treatment[subject.indices[temp.indices[1]]]==0){
        treated.subject.index[count]=subject.indices[temp.indices[2]]
        matched.control.subject.index[count]=subject.indices[temp.indices[1]]
      }
      included.already[temp.indices]=1
      count=count+1
    }
  }
  matched.pair.treated=subset[treated.subject.index,]
  matched.pair.control=subset[matched.control.subject.index,]
  
  if (!exists("out.treated")) {
    out.treated <- matched.pair.treated
  } else {
    out.treated <- rbind(out.treated, matched.pair.treated)
  }
  
  if (!exists("out.control")) {
    out.control <- matched.pair.control
  } else {
    out.control <- rbind(out.control, matched.pair.control)
  }
  

### matching balance check
treated.mat=Xmat[treatment==1,]
control.mat.before=Xmat[treatment==0,] ## containing all control units before matching
control.mat.after=Xmat[matched.control.subject.index,] ## containing matched controls

before.matching.control.mean=apply(control.mat.before, 2, mean)
new.sd=sqrt((apply(treated.mat, 2, var)+apply(control.mat.before, 2, var))/2) ## compute sds for variables 

matching.treated.mean=apply(treated.mat, 2, mean)
after.matching.control.mean=apply(control.mat.after, 2, mean)

std.diff.mat=round(cbind((matching.treated.mean-before.matching.control.mean)/new.sd, (matching.treated.mean-after.matching.control.mean)/new.sd), 3)
colnames(std.diff.mat)=c("before", "after")
std.diff.mat
}
rm(unmatched)
gc()

out <- data.frame(init = rep_len(NA, nrow(out.treated)))
out[[paste0(outcome_var, "_t")]] <- out.treated[[outcome_var]]
out[[paste0(outcome_var, "_c")]] <- out.control[[outcome_var]]
out$init <- NULL
for (exact in exact_vars) {
  out[[exact]] <- out.treated[[exact]]
}

for (match in matching_vars) {
  out[[match]] <- (out.control[[match]] + out.treated[[match]])/2
}

write.csv(out, out_filename, row.names = F)

### outcome comparison
##mean(unmatched[unmatched[[treatment_var]]==1, "re78"]) ## average earnings for treated subjects
##mean(unmatched[unmatched[[treatment_var]]==0, "re78"]) ## before matching
##mean(matched.pair.control$re78) ## only for matched controls
