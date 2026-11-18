## function for computing phylogenetic signal by the lambda (Pagel 1999)
## or K (Blomberg et al. 2003) methods
## written by Liam J. Revell 2011/2012, 2019, 2020, 2021, 2023, ...

# function matches data to tree
# written by Liam J. Revell 2011, 2015
matchDatatoTree<-function(tree,x,name){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.matrix(x)) x<-x[,1]
	if(is.null(names(x))){
		if(length(x)==Ntip(tree)){
			print(paste(name,"has no names; assuming x is in the same order as tree$tip.label"))
			names(x)<-tree$tip.label
		} else

			stop(paste(name,"has no names and is a different length than tree$tip.label"))

	}
	if(any(is.na(match(names(x),tree$tip.label)))){
		print(paste("some species in",name,"are missing from tree, dropping missing taxa from",name))
		x<-x[intersect(tree$tip.label,names(x))]
	}
	return(x)
}

# function matches tree to data
# written by Liam J. Revell 2011, 2015
matchTreetoData<-function(tree,x,name){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(any(is.na(match(tree$tip.label,names(x))))){
		print(paste("some species in tree are missing from",name,", dropping missing taxa from the tree"))
		tree<-drop.tip(tree,setdiff(tree$tip.label,names(x)))
	}
	if(any(is.na(x))){
		print(paste("some data in",name,"given as 'NA', dropping corresponding species from tree"))
		tree<-drop.tip(tree,names(which(is.na(x))))
	}
	return(tree)
}

get_pagel_cov <- function(lam, tau) {
  
  n <- length(tau$tip.label)
  sig_vcv <- ape::vcv.phylo(tau)
  diag_vcv <- diag(sig_vcv)
  
  sig_lam <- sig_vcv - diag(diag_vcv, nrow = n)
  out_vcv <- diag(diag_vcv, nrow = n) + lam*sig_lam
  
  out_vcv
}

maxLambda<-function(tree){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.ultrametric(tree)){
		H<-nodeHeights(tree)
		return(max(H[,2])/max(H[,1]))
	} else return(1)
}

phylosig<-function(tree,x,method="K",test=FALSE,nsim=1000,
	se=NULL,start=NULL,control=list(),niter=10){
	# some minor error checking
	if(!inherits(tree,"phylo")) 
		stop("tree should be an object of class \"phylo\".")
	x<-matchDatatoTree(tree,x,"x")
	tree<-matchTreetoData(tree,x,"x")
	if(!is.null(se)){
		se<-matchDatatoTree(tree,se,"se")
		tree<-matchTreetoData(tree,se,"se")
		me=TRUE
		M<-diag(se^2)
		rownames(M)<-colnames(M)<-names(se)
	} else me=FALSE
	if(!is.null(start)&&!is.null(se)){
		if(start[1]<=0||start[2]<0||start[2]>maxLambda(tree)){
			message("some of the elements of 'start' are invalid, resetting to random")
			start<-NULL
		}
	}
	# done error handling
	if(method=="K"){
		C<-vcv.phylo(tree)
		x<-x[rownames(C)]
		n<-nrow(C)
		if(!me){
			invC<-solve(C)
			a<-sum(invC%*%x)/sum(invC)
			K<-(t(x-a)%*%(x-a)/(t(x-a)%*%invC%*%(x-a)))/((sum(diag(C))-
				n/sum(invC))/(n-1)) # calculate K
			if(!test){
				object<-as.numeric(K)
			} else {
				P=0.0
				simX<-x
				simK<-vector()
				for(i in 1:nsim){
					a<-sum(invC%*%simX)/sum(invC)
					simK[i]<-(t(simX-a)%*%(simX-a)/(t(simX-a)%*%invC%*%
						(simX-a)))/((sum(diag(C))-n/sum(invC))/(n-1))
					## calculate P-value for randomization test
					if(simK[i]>=K) P<-P+1/nsim 
					simX<-sample(simX) # randomize x
				}
				object<-list(K=as.numeric(K),P=P,sim.K=simK)
			}
		} else {
			likelihoodK<-function(theta,C,M,y){
				Ce<-theta*C+M
				invCe<-solve(Ce)
				a<-as.numeric(sum(invCe%*%y)/sum(invCe))
				logL<--t(y-a)%*%invCe%*%(y-a)/2-n*log(2*pi)/2-
					determinant(Ce,logarithm=TRUE)$modulus/2
				logL[1,1]
			}
			M<-M[rownames(C),colnames(C)]
			invC<-solve(C)
			maxSig2<-as.numeric(t(x-as.numeric(sum(invC%*%x)/
				sum(invC)))%*%invC%*%(x-as.numeric(sum(invC%*%x)/
				sum(invC)))/n)
			res<-optimize(f=likelihoodK,interval=c(0,maxSig2),y=x,
				C=C,M=M,maximum=TRUE) # optimize sig2
			sig2<-res$maximum*n/(n-1)
			Ce<-sig2*C+M
			invCe<-solve(Ce)
			a<-as.numeric(sum(invCe%*%x)/sum(invCe))
			K<-(t(x-a)%*%(x-a)/(t(x-a)%*%invCe%*%(x-a)))/((sum(diag(Ce))-
				n/sum(invCe))/(n-1)) # calculate K
			if(!test){
				object<-list(K=as.numeric(K),sig2=as.numeric(sig2),
					logL=res$objective,
					lik=function(sig2) likelihoodK(sig2,C=C,M=M,y=x))
			} else {
				P=0.0
				simX<-x
				simK<-vector()
				for(i in 1:nsim){
					maxSig2<-as.numeric(t(simX-as.numeric(sum(invC%*%
						simX)/sum(invC)))%*%invC%*%(simX-as.numeric(sum(invC%*%
						simX)/sum(invC)))/n)
					simRes<-optimize(f=likelihoodK,interval=c(0,maxSig2),y=simX,
						C=C,M=M,maximum=TRUE) # optimize sig2
					simSig2<-simRes$maximum*n/(n-1)
					Ce<-simSig2*C+M
					invCe<-solve(Ce)
					a<-as.numeric(sum(invCe%*%simX)/sum(invCe))
					# calculate K
					simK[i]<-(t(simX-a)%*%(simX-a)/(t(simX-a)%*%invCe%*%
						(simX-a)))/((sum(diag(Ce))-n/sum(invCe))/(n-1)) 
					# calculate P-value for randomization test
					if(simK[i]>=K) P<-P+1/nsim 
					o<-sample(1:n)
					simX<-x[o]
					M<-diag(diag(M)[o]^2) # randomize x & errors
				}
				object<-list(K=as.numeric(K),P=P,sim.K=simK,
					sig2=as.numeric(sig2),logL=res$objective,
					lik=function(sig2) likelihoodK(sig2,C=C,M=M,y=x))
			}
		}
	} else if(method=="lambda"){
		# function to compute C with lambda
		lambda.transform<-function(C,lambda){
			dC<-diag(diag(C))
			C<-lambda*(C-dC)+dC
			C
		}
		# fast likelihood functions, accepting eigen decomposition of C
		fast_lhoodLambda<-function(theta,y, eig_values, eig_vectors, z_1, z_x, tree_age){
			n<-nrow(eig_vectors)
			f_eig <- 1/(1+theta*as.numeric(eig_values))
			a<-as.numeric(sum(z_1*f_eig*z_x)/sum(f_eig*z_1^2))
			sig_2<-as.numeric((1/tree_age)*sum(f_eig*(t(eig_vectors)%*%(y-a))^2)/n)
			logL <- -n/2 -n*log(2*pi)/2 - n*log(sig_2)/2 + 0.5*sum(log((1/tree_age)*(1/(1+theta*eig_values))))
			logL
		}

		# standard likelihood
		likelihoodLambda<-function(theta,C,y){
			Cl<-lambda.transform(C,theta)
			invCl<-solve(Cl)
			n<-nrow(Cl)
			y<-y[rownames(Cl)]
			a<-as.numeric(sum(invCl%*%y)/sum(invCl))
			sig2<-as.numeric(t(y-a)%*%invCl%*%(y-a)/n)
			#logL<--t(y-a)%*%(1/sig2*invCl)%*%(y-a)/2-n*log(2*pi)/2-
			#	determinant(sig2*Cl,logarithm=TRUE)$modulus/2
			logL<- -n/2-n*log(2*pi)/2-
				determinant(sig2*Cl,logarithm=TRUE)$modulus/2
			logL
		}
		# likelihood function with error
		likelihoodLambda.me<-function(theta,C,y,M){
			Cl<-theta[1]*lambda.transform(C,theta[2])
			V<-Cl+M
			invV<-solve(V)
			n<-nrow(Cl)
			y<-y[rownames(Cl)]
			a<-as.numeric(sum(invV%*%y)/sum(invV))
			logL<--t(y-a)%*%invV%*%(y-a)/2-n*log(2*pi)/2-
				determinant(V,logarithm=TRUE)$modulus/2
			logL[1,1]
		}
		C<-vcv.phylo(tree)
		x<-x[rownames(C)]

		# check C has constant diagonal, required for fast lhood
		d <- diag(C)
		ultraMetric <- diff(range(d)) < 1e-10 * max(d)

		if (ultraMetric & !me) {
			likelihoodLambda <- fast_lhoodLambda

			# precompute eigen decomposition and associated products
			tree_age <- C[1, 1]
			off_diag <- (1/tree_age) * (C - diag(diag(C), nrow = length(tree$tip.label)))
			eigen_C <- eigen(off_diag)

			# precompute products required for likelihood
			z_1 <- t(eigen_C$vectors) %*% rep(1, length(tree$tip.label))
            z_x <- t(eigen_C$vectors) %*% x
		} else {
			if (!ultraMetric) {
				warning("Tree not ultrametric, using slow likelihood function for lambda")
			}
		}

		maxlam<-maxLambda(tree)
		if(!me){
			INTERVAL<-cbind(
				seq(0,maxlam-maxlam/niter,length.out=niter),
				seq(maxlam/niter,maxlam,length.out=niter))
			if (ultraMetric){
				FITS<-apply(INTERVAL,1,optimize,f=likelihoodLambda,
					y=x, eig_values=eigen_C$values, eig_vectors=eigen_C$vectors,
					z_1=z_1, z_x=z_x, tree_age=tree_age, maximum=TRUE)
			} else {
				FITS<-apply(INTERVAL,1,optimize,f=likelihoodLambda,
					y=x,C=C,maximum=TRUE)
			}
			LIKS<-sapply(FITS,function(x) x$objective)
			res<-FITS[[which(LIKS==max(LIKS))[1]]]
			if (ultraMetric) {
				   lik <- function(theta) fast_lhoodLambda(theta, y=x, eig_values=eigen_C$values, eig_vectors=eigen_C$vectors, z_1=z_1, z_x=z_x, tree_age=tree_age)
			} else {
				   lik <- function(theta) likelihoodLambda(theta, C=C, y=x)
			}
			if(!test){
				object<-list(lambda=res$maximum,logL=res$objective,
					lik=lik)
			} else {
				# compute likelihood of lambda=0
				
				logL0<-lik(0)
				P<-as.numeric(pchisq(2*(res$objective-logL0),df=1,
					lower.tail=FALSE)) # P-value
				object<-list(lambda=res$maximum,logL=res$objective,
					logL0=logL0,P=P,lik=lik)
			}
		} else {
			FITS<-list()
			control$fnscale=-1
			M<-M[rownames(C),colnames(C)]
			for(i in 1:niter){
				if(is.null(start)) s<-c(0.5*runif(n=1)*mean(pic(x,
					multi2di(tree,random=FALSE))^2),runif(n=1))
				else s<-start
				FITS[[i]]<-optim(s,likelihoodLambda.me,C=C,y=x,M=M,
					method="L-BFGS-B",lower=c(0,0),upper=c(Inf,maxlam),
					control=control)
			}
			LIKS<-sapply(FITS,function(x) x$value)
			res<-FITS[[which(LIKS==max(LIKS))[1]]]
			if(!test){
				object<-list(lambda=res$par[2],sig2=res$par[1],
					logL=res$value,convergence=res$convergence,
					message=res$message,lik=function(lambda,
					sig2=res$par[1]) likelihoodLambda.me(c(sig2,lambda),
					C=C,M=M,y=x))
			} else {
				for(i in 1:niter){
					s<-0.5*runif(n=1)*mean(pic(x,multi2di(tree,
						random=FALSE))^2)
					FITS[[i]]<-optim(c(s,0),likelihoodLambda.me,C=C,
						y=x,M=M,method="L-BFGS-B",lower=c(0,0),
						upper=c(Inf,1e-10),control=control)
				}
				LIKS<-sapply(FITS,function(x) x$value)
				res0<-FITS[[which(LIKS==max(LIKS))[1]]]
				P<-as.numeric(pchisq(2*(res$value-res0$value),df=1,
					lower.tail=FALSE))
				object<-list(lambda=res$par[2],sig2=res$par[1],
					logL=res$value,convergence=res$convergence,
					message=res$message,logL0=res0$value,P=P,
					lik=function(lambda,sig2=res$par[1]) likelihoodLambda.me(c(sig2,
					lambda),C=C,M=M,y=x))
			}
		}
	} else
		stop(paste("do not recognize method = \"",method,
			"\"; methods are \"K\" and \"lambda\"",sep=""))
	attr(object,"class")<-"phylosig"
	attr(object,"method")<-method
	attr(object,"test")<-test
	attr(object,"se")<-!is.null(se)
	object
}

print.phylosig<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-6
	cat("\n")
	if(attr(x,"method")=="K"){
		if(attr(x,"test")||attr(x,"se")){
			cat(paste("Phylogenetic signal K :",signif(x$K,
				digits),"\n"))
			if(attr(x,"se")){
				cat(paste("MLE(sig2) :",signif(x$sig2,digits),
					"\n"))
				cat(paste("logL(sig2) :",signif(x$logL,digits),
					"\n"))
			}
			if(attr(x,"test"))
				cat(paste("P-value (based on",length(x$sim.K),
					"randomizations) :",signif(x$P,digits),"\n"))
		} else cat(paste("Phylogenetic signal K :",signif(x[1],
			digits),"\n"))
	} else if(attr(x,"method")=="lambda"){
		cat(paste("Phylogenetic signal lambda :",signif(x$lambda,
			digits),"\n"))
		cat(paste("logL(lambda) :",signif(x$logL,digits),"\n"))
		if(attr(x,"se")) cat(paste("MLE(sig2) :",signif(x$sig2,
			digits),"\n"))
		if(attr(x,"test")){
			cat(paste("LR(lambda=0) :",signif(2*(x$logL-x$logL0),
				digits),"\n"))
			cat(paste("P-value (based on LR test) :",signif(x$P,
				digits),"\n"))
		}
	}
	cat("\n")			
}

plot.phylosig<-function(x,...){
	if(hasArg(what)) what<-list(...)$what
	else what<-if(attr(x,"method")=="lambda") "lambda" else 
		if(attr(x,"method")=="K"&&attr(x,"test")) "K" else "sig2"
	if(hasArg(res)) res<-list(...)$res
	else res<-100
	if(hasArg(las)) las<-list(...)$las
	else las<-par()$las
	if(hasArg(cex.lab)) cex.lab<-list(...)$cex.lab*par()$cex
	else cex.lab<-par()$cex.lab
	if(hasArg(cex.axis)) cex.axis<-list(...)$cex.axis*par()$cex
	else cex.axis<-par()$cex.axis
	if(hasArg(bty)) bty<-list(...)$bty
	else bty<-par()$bty
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-NULL
	if(attr(x,"method")=="lambda"){
		lambda<-seq(0,max(c(1,x$lambda)),length.out=res)
		logL<-sapply(lambda,x$lik)
		plot(lambda,logL,xlab=expression(lambda),ylab="log(L)",
			type="l",bty=bty,las=las,cex.lab=cex.lab,cex.axis=cex.axis,
			xlim=xlim)
		lines(rep(x$lambda,2),c(par()$usr[3],x$logL),lty="dotted")
		text(x=x$lambda+0.01*diff(par()$usr[1:2]),
			par()$usr[3]+0.5*diff(par()$usr[3:4]),
			expression(paste("MLE(",lambda,")")),srt=90,adj=c(0.5,1))
		if(attr(x,"test")){
			lines(c(0,x$lambda),rep(x$logL0,2),lty="dotted")
			text(x=0.5*x$lambda,
				y=x$logL0+
				if(x$logL0>(par()$usr[3]+0.5*(diff(par()$usr[3:4]))))
					-0.01*diff(par()$usr[3:4])
				else 0.01*diff(par()$usr[3:4]),
				expression(paste("logL(",lambda,"=0)")),
				adj=if(x$logL0>(par()$usr[3]+0.5*(diff(par()$usr[3:4]))))
					c(0.5,1)
				else c(0.5,0))
		}
	} else if(attr(x,"method")=="K"){
		if(what=="K"){
			if(attr(x,"test")==FALSE)
				cat("Sorry. This is not a valid plotting option for your object.\n\n")
			else {
				hist(x$sim.K,breaks=min(c(max(12,round(length(x$sim.K)/10)),
					20)),bty=bty,col="lightgrey",border="lightgrey",
					main="",xlab="K",ylab="null distribution of K",
					las=las,cex.lab=cex.lab,cex.axis=cex.axis)
				arrows(x0=x$K,y0=par()$usr[4],y1=0,length=0.12,
					col=make.transparent("blue",0.5),lwd=2)
				text(x$K,0.95*par()$usr[4],"observed value of K",
					pos=if(x$K>mean(range(x$sim.K))) 2 else 4)
			}
		} else if(what=="sig2"){
			if(attr(x,"se")==FALSE)
				cat("Sorry. This is not a valid plotting option for your object.\n\n")
			else {
				sig2<-seq(0.25*x$sig2,1.75*x$sig2,length.out=res)
				logL<-sapply(sig2,x$lik)
				plot(sig2,logL,xlab=expression(sigma^2),ylab="log(L)",
					type="l",bty=bty,las=las,cex.lab=cex.lab,
					cex.axis=cex.axis,xlim=xlim)
				lines(rep(x$sig2,2),c(par()$usr[3],x$logL),lty="dotted")
				text(x=x$sig2+0.01*diff(par()$usr[1:2]),
					par()$usr[3]+0.5*diff(par()$usr[3:4]),
					expression(paste("MLE(",sigma^2,")")),srt=90,
					adj=c(0.5,1))
			}
		}	
	}
}