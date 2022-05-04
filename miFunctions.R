# ------------------------------------------------------------------------------
# Program: miFunctions.R  
#  Author: Hermine Maes
#    Date: 02 21 2018 
#
# Set of my options & functions used in basic twin methodology scripts
#   Email: hmaes@vcu.edu
# -------|---------|---------|---------|---------|---------|---------|---------|

# Options
#mxOption( NULL, "Default optimizer", "NPSOL" )
#mxOption(NULL, "Default optimizer", "SLSQP")
#mxOption(NULL, "Default optimizer", "CSOLNP")
mxOption(NULL, 'Number of Threads', parallel::detectCores()) #now
#mxOption( NULL, "Checkpoint Prefix", filename )
mxOption( NULL, "Checkpoint Units", "iterations" )
mxOption( NULL, "Checkpoint Count", 1 )
options(width=120)
options(digits=8)
mxVersion()

# Functions to assign labels
labLower  <- function(lab,nv) { paste(lab,rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="") }
labSdiag  <- function(lab,nv) { paste(lab,rev(nv+1-sequence(1:(nv-1))),rep(1:(nv-1),(nv-1):1),sep="") }
labOdiag  <- function(lab,nv) { paste(lab,c(rev(nv+1-sequence(1:(nv-1))),rep(1:(nv-1),(nv-1):1)),c(rep(1:(nv-1),(nv-1):1),rev(nv+1-sequence(1:(nv-1)))),sep="") }
labFullSq <- function(lab,nv) { paste(lab,1:nv,rep(1:nv,each=nv),sep="") }
labDiag   <- function(lab,nv) { paste(lab,1:nv,1:nv,sep="") } 
labSymm   <- function(lab,nv) { paste(lab,rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="") }
labFull   <- function(lab,nr,nc) { paste(lab,1:nr,rep(1:nc,each=nr),sep="") }
labFullR  <- function(lab,nr,nc) { paste(lab,rep(1:nr,each=nc),1:nc,sep="") }
labVect   <- function(lab,nv) { paste(lab,1:nv, sep="") }
labVars   <- function(lab,vars) { paste(lab,vars,sep="") }
labTh     <- function(lab,vars,nth) { paste(paste("t",1:nth,lab,sep=""),rep(vars,each=nth),sep="") }


# Functions to assign values
valDiag   <- function(valD,dim) {
valF      <- diag(valD,dim,dim)        # values for diagonal of covariance matrix
valF
}
valDiagO  <- function(valD,valOD,dim) {
valF      <- diag(valD,dim,dim)        # values for diagonal of covariance matrix
valF[lower.tri(valF)] <- valOD         # values for below diagonal elements 
valF[upper.tri(valF)] <- valOD         # values for above diagonal elements
valF
}
valDiagLU <- function(valD,valLD,valUD,dim) {
valF      <- diag(valD,dim,dim)        # values for diagonal of covariance matrix
valF[lower.tri(valF)] <- valLD         # values for below diagonal elements 
valF[upper.tri(valF)] <- valUD         # values for above diagonal elements
valF
}

myMean    <- function(x) { mean(as.numeric(x),na.rm=TRUE) }
myCov     <- function(x) { cov(as.numeric(x),use="complete.obs") }
myCor     <- function(x) { cor(as.numeric(x),use="everything") }

# Functions to generate output

fitGofs   <- function(fit) {
          summ <- summary(fit)
          cat(paste("Mx:", fit$name,"  os=", summ$ob,"  ns=", summ$nu,"   ep=", summ$es,
                    "   co=", sum(summ$cons),"  df=", summ$de, "  ll=", round(summ$Mi,4), 
                    "  cpu=", round(summ$cpu,4),"  opt=", summ$op,"  ver=", summ$mx,
                    "  stc=", fit$output$status$code, "\n",sep=""))
}

fitGofS   <- function(fit) {
          summ <- summary(fit)
          cat(paste("Mx:", fit$name,"  #statistics=", summ$ob,"  #records=", summ$nu,"   #parameters=", summ$es,
                    "   #constraints=", sum(summ$cons),"  df=", summ$de, "  -2LL=", round(summ$Mi,4), 
                    "  cpu=", round(summ$cpu,4),"  optim=", summ$op,"  version=", summ$mx,
                    "  code=", fit$output$status$code, "\n",sep=""))
}

fitGofT   <- function(fit) {
          summ <- summary(fit)
          cat(paste("Mx:", fit$name," ", summ$ob," ", summ$nu," ", summ$es," ", sum(summ$cons)," ", summ$de," ", round(summ$Mi,4)," ", 
                    round(summ$Mi-2*summ$de,4)," ", fit$output$status$code, "\n",sep=""))
}


fitEsts   <- function(fit) {
	 print(round(fit$output$estimate,4)) 
}

fitEstCis   <- function(fit) {
	 print(round(fit$output$estimate,4)) 
	 print(round(fit$output$confidenceIntervals,4))
}

fitEstVC   <- function(fit) {
	 print(round(fit$output$estimate,4)) 
	 print(round(fit$output$confidenceIntervals,4))
	 print(round(fit$VC$result,4))
}

fitEstVCfm   <- function(fit) {
	 print(round(fit$output$estimate,4)) 
	 print(round(fit$output$confidenceIntervals,4))
	 print(round(fit$VCf$result,4))
	 print(round(fit$VCm$result,4))
}

lrtSAT    <- function(fit,llSAT,dfSAT) {
 	 chi <- summary( fit, SaturatedLikelihood = llSAT, SaturatedDoF = dfSAT)$Chi
 	 df  <- summary( fit, SaturatedLikelihood = llSAT, SaturatedDoF = dfSAT)$ChiDoF
 	 p   <- 1-pchisq(chi,df) 
 	 cat(paste( "lrtX2(df=",df,")=",round(chi,4),", p=",round(p,4),sep="")) 
 	 }
		 
fitExpc   <- function(fit) {
print(round(fit$MZ$meanMZ$values,4))
print(round(fit$DZ$meanDZ$values,4))
print(round(fit$MZ$covMZ$values,4))
print(round(fit$DZ$covDZ$values,4)) }

fitExpc2   <- function(fit) {
print(round(fit$MZ$meanMZ$values,4))
print(round(fit$DZ$meanDZ$values,4))
print(round(fit$MZ$covMZ$result,4))
print(round(fit$DZ$covDZ$result,4)) }

fitExpb   <- function(fit) {
print(round(fit$MZ$threMZ$values,4))
print(round(fit$DZ$threDZ$values,4))
print(round(fit$MZ$corMZ$values,4))
print(round(fit$DZ$corDZ$values,4)) }

fitExpo   <- function(fit) {
print(round(fit$MZ$thinMZ$values,4))
print(round(fit$DZ$thinDZ$values,4))
print(round(fit$MZ$threMZ$values,4))
print(round(fit$DZ$threDZ$values,4))
print(round(fit$MZ$corMZ$values,4))
print(round(fit$DZ$corDZ$values,4)) }

fitExpm   <- function(fit) {
print(round(fit$MZ$thinMZ$values,4))
print(round(fit$DZ$thinDZ$values,4))
print(round(fit$MZ$threMZ$values,4))
print(round(fit$DZ$threDZ$values,4))
print(round(fit$MZ$corMZ$values,4))
print(round(fit$DZ$corDZ$values,4))
print(round(fit$MZ$meanMZ$values,4))
print(round(fit$DZ$meanDZ$values,4))
print(round(fit$MZ$covMZ$values,4))
print(round(fit$DZ$covDZ$values,4)) }

# Function "parameterSpecifations()" prints labels of a MxMatrix with
# square brackets surrounding free parameters; returns a matrix of strings
# -----------------------------------------------------------------------
parameterSpecifications <- function(model) {
	resultsList <- .collectParameterSpecifications(model)
	if(length(resultsList) > 0) {
		resultsNames <- names(resultsList)
		for(i in 1:length(resultsList)) {
			cat(resultsNames[[i]],'\n')
			print(resultsList[[i]], quote=FALSE)
			cat('\n')
		}
	}
}

.collectParameterSpecifications <- function(model) {
	listReturn <- list()
	if(length(model@matrices) > 0) {
		for(i in 1:length(model@matrices)) {
			current <- model@matrices[[i]]
			extract <- is(current, "FullMatrix") ||
				is(current, "LowerMatrix") ||
				is(current, "DiagMatrix") ||
				is(current, "SymmMatrix") ||
				is(current, "StandMatrix")
			if(extract) {
				retval <- mapply(.parameterSpecificationsHelper, 
					current@labels, current@free, current@values)
				retval <- matrix(retval, nrow(current), ncol(current))
				dimnames(retval) <- dimnames(current)
				storeName <- paste('model:', model@name,', matrix:', current@name, sep='')
				listReturn[[storeName]] <- retval
			}
		}
	}
	names(model@submodels) <- NULL
	matrices <- lapply(model@submodels, .collectParameterSpecifications)
	listReturn <- append(listReturn, unlist(matrices, FALSE))
	return(listReturn)
}

.parameterSpecificationsHelper <- function(label, free, value) {
	if(free) return(paste('[', label, ']', sep = ''))
	else return(value)
}

# Function "formatOutputMatrices()" prints matrix with specified labels and
# number of decimals
# -----------------------------------------------------------------------
#parse(text=matricesList[k]) == matricesList[[k]]
formatOutputMatrices <- function(fittedModel,matricesList,labelsList,vars,digits) {
	if(length(matricesList) > 0) {
	for(k in 1:length(matricesList)) {
		print(paste("Matrix",matricesList[[k]]))
		print(formatOutputMatrix(
			evalQuote(matricesList[[k]], fittedModel),
			labelsList[[k]],vars,digits), quote=FALSE)
			cat('\n')
		}
	}
}

formatOutputMatrix <- function(matrix,label,vars,digits) {
	#table <- round(eval(substitute(mxEval(Matrix,Model))),ND)
	matrix <- apply(matrix, c(1,2), round, digits = digits)
	retval <- apply(matrix, c(1,2), format, scientific=FALSE, nsmall = digits)

	cols <- character(ncol(retval))
	for(i in 1:ncol(retval)) {paste(label,i,sep="")} -> cols[i]
	colnames(retval) <- cols
	if (nrow(retval) == length(vars)) {
	rownames(retval) <- vars
	} else {
	rows <- character(nrow(retval))
	for(j in 1:nrow(retval)) {paste("LP",j,sep="")} -> rows[j]
	rownames(retval) <- rows
	}
	return(retval)
}



# Function "formatMatrix()" returns a matrix with specified dimnames and # of decimal places
# -----------------------------------------------------------------------
formatMatrix <- function(matrix, dimnames, digits) {
	retval <- apply(matrix, c(1,2), round, digits)
	dimnames(retval) <- dimnames
	return(retval)
}

evalQuote <- function(expstring, model, compute = FALSE, show = FALSE) {
	return(eval(substitute(mxEval(x, model, compute, show),
			list(x = parse(text=expstring)[[1]]))))
}

