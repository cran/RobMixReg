
#library(glmnet)

#' Simulate high dimension data for RBSL algorithm validation.
#'
#' @param n Patient number.
#' @param bet The coefficient matrix.
#' @param pr A vector of probability threshold which simulate the sampling based on uniform distribution.
#' @param sigma A vector of noise level. The length should be equal to the component number.
#' @return A list object consist of x, y, true cluster label.
simu_data_sparse<-function(n,bet, pr, sigma ){
  nc=length(pr)
  x=t(replicate(ncol(bet)-1,rnorm(n,0,1)))
  rownames(x)=paste("feature",1:(ncol(bet)-1),sep="")
  colnames(x)=paste("Sample",1:n,sep="")
  u=runif(n,0,1);
  clusters=rep(0,n)
  y=rep(0,n)
  pr_tmp=c(0,cumsum(pr))
  for(j in 1:nc){
    ee=rnorm(n,0,sigma[j])
    yy=bet[j,]%*%rbind(1,x)+ee
    y=y+(u<=pr_tmp[j+1] & u>pr_tmp[j])*yy
    clusters=clusters+(u<=pr_tmp[j+1] & u>pr_tmp[j])*j
  }
  x=t(x)
  names(y)=paste("S",1:n,sep="")
  sss=order(clusters)
  x=x[sss,]
  y=y[sss]
  clusters=clusters[sss]
  tmp_list=list(x=x,y=y,clusters=clusters);names(tmp_list)=c("x","y","clusters")
  return(tmp_list)
}


#' Adaptive lasso.
#'
#' @param XX The independent variable.
#' @param yy The dependent variable.
#' @return A list object consist of index of selected variable and coefficient for all variables.
Rec_Lm<-function(XX,yy){
  ok<-complete.cases(XX,yy)
  XX<-XX[ok,]                            # get rid of na's
  yy<-yy[ok]                             # since regsubsets can't handle na's
  m<-ncol(XX)
  n<-nrow(XX)
  XX=as.matrix(XX)

  pp=ncol(XX)+1
  nn=nrow(XX)
  coeff_all = rep(0, ncol(XX))
  names(coeff_all) = colnames(XX)
  if(pp >=nn -2){
    aaa=lars(XX,yy)
    lars_coff = coef(aaa,s=nn-2,mode="step")
    ZZ=XX[,which(lars_coff!=0),drop=FALSE]
    zero_ind = which(lars_coff==0)

    # cvfit <- cv.glmnet(XX, yy,alpha=1)
    # # # lars_coff = as.vector(coef(cvfit, s = "lambda.1se"))[-1]
    # lars_coff = as.vector(coef(cvfit, cvfit$lambda.min))[-1]
    # ZZ=XX[,which(lars_coff!=0),drop=FALSE]
    # zero_ind = which(lars_coff==0)
    # if(length(zero_ind)==0) print("mywarning!!!")
  }else{
    ZZ=XX
  }

  out <- lm(yy~ZZ) #20200604 remove global variable
  # print("####")
  # print(out)
  out.lsa = lsa(out)
  #	print("####")
  coeff<-out.lsa$beta.bic
  coeff2<-coeff[-1]               # get rid of intercept
  if(pp >=nn -2){
    coeff_all[zero_ind] = 0
    coeff_all[-zero_ind] = coeff2
  }else{
    coeff_all=coeff2
  }
  XX.ind=which(coeff_all!=0)
  return(list(XX.ind=XX.ind,coeff_all=coeff_all))
}


#' Perform the RBSL algorithm one times.
#'
#' @param x The matrix
#' @param y The external supervised variable.
#' @param nit xxx?
#' @param nc The component number in the mixture model.
#' @param max_iter The maximum iteration number.
#' @return A list object consist of coefficient, clustering membership, data x, external variable y, predicted y based on regression model.
CSMR_one<-function(x,y,nit=1,nc,max_iter){
  data=data.frame((x),y)
  mycall = match.call();
  nx=ncol(data)-1; nobs = nrow(data)            # number of observations
  flag=0; ccc=0
  x_small=x
  f_small = nobs/(nc) - 3
  if(f_small < nx){
    sel_IND = order(abs(cor(x,y)),decreasing=TRUE)[1:f_small]
    x_small=x[,sel_IND]
  }
  f_list=sapply(1:nc, function(ii)as.formula(paste("y~",paste(colnames(x_small),collapse="+"),sep="")))
  ppp = mixtureReg(regData = data.frame(x_small,y),formulaList = f_list,mixingProb = "Constant",silently=TRUE,max_iter = 100)
  #www=cbind(ppp$posterior[[1]], ppp$posterior[[2]])
  www=do.call(cbind,  ppp$posterior)
  coffs_old = matrix(rep(Inf,nc*(nx)),ncol=nc)
  while(flag==0 ){
    #print(ccc)
    coffs=matrix(0,nrow=ncol(x),ncol=nc)
    rownames(coffs)=colnames(x)
    f_list=vector("list",nc)
    for(ii in 1:nc){
      inds_in_tmp=which(apply(www, 1, which.max)==ii & apply(www, 1, max)>0)
      x_tmp=x[inds_in_tmp,,drop=FALSE]
      y_tmp=y[inds_in_tmp]
      mod_cv = Rec_Lm(x_tmp,y_tmp)
      coffs[,ii]=mod_cv$coeff_all
      f_list[[ii]]=as.formula(paste("y~",paste(colnames(x)[coffs[,ii]!=0],collapse="+"),sep=""))
    }
    ppp = mixtureReg(regData = data.frame(x,y),formulaList = f_list, mixingProb = "loess",silently=TRUE,max_iter = 100)
    #www=cbind(ppp$posterior[[1]], ppp$posterior[[2]])
    www=do.call(cbind,  ppp$posterior)
    #uu1=coefficients(ppp$lmList[[1]])[-1]
    #uu2=coefficients(ppp$lmList[[2]])[-1]
    #	coffs=coffs*0
    #	coffs[match(names(uu1),rownames(coffs)),1]=uu1
    #	coffs[match(names(uu2),rownames(coffs)),2]=uu2
    ccc=ccc+1
    if(ccc>max_iter  ){flag=1}
    if( sum((coffs_old-coffs)^2)<0.00001 & ccc>1){
      flag=1;
    }else{
      coffs_old=coffs
    }
  }
  predit = sapply(ppp$lmList, function(ii)predict(ii))
  clus = apply(www, MARGIN=1, which.max)
  yhat = unlist(sapply(1:length(clus), function(gg)predit[gg,clus[gg]]))
  mx.model = ppp   #also return model
  return(list(coffs=coffs, clus=clus, yhat=yhat, mx.model=mx.model))
}


#' The main function of the RBSL algorithm.
#'
#' @param x The matrix
#' @param y The external supervised variable.
#' @param nit xxx?
#' @param nc The component number in the mixture model.
#' @param max_iter The maximum iteration number.
#' @return A list object consist of coefficient, clustering membership, data x, external variable y, predicted y based on regression model.
CSMR<-function(x,y,nit,nc,max_iter){
  rrr = CSMR_one(x,y,nit,nc,max_iter); #print(rrr$coffs)
  x1=x
  while(length(which(rowSums(abs(rrr$coffs))==0))>0){
    x1=x1[,which(rowSums(abs(rrr$coffs))>0),drop=FALSE]
    rrr = CSMR_one(x1,y,nit,nc,max_iter); #print(rrr$coffs)
  }
  coffs_final = matrix(0,nrow=ncol(x),ncol=nc)
  rownames(coffs_final)=colnames(x)
  coffs_final[match(rownames(rrr$coffs),rownames(coffs_final)),]=rrr$coffs
  return(list(coffs=coffs_final,clus=rrr$clus,x=x,y=y,yhat=rrr$yhat))
}

#' The train function of the CSMR algorithm.
#'
#' @param x The matrix
#' @param y The external supervised variable.
#' @param nit xxx
#' @param nc The component number in the mixture model.
#' @param max_iter The maximum iteration number.
#' @return A list object consist of coefficient, clustering membership, data x, external variable y, predicted y based on regression model.
CSMR_train<-function(x,y,nit,nc,max_iter){
  if(nc==1){
    print("Only have 1 component.[CSMR train].")
    ccc=Rec_Lm(x,y)
    return(ccc)   #note:return is different to situation of nc>1.
  }
  rrr = CSMR_one(x,y,nit,nc,max_iter);# print(rrr$coffs)
  x1=x
  while(length(which(rowSums(abs(rrr$coffs))==0))>0){
    x1=x1[,which(rowSums(abs(rrr$coffs))>0),drop=FALSE]
    rrr = CSMR_one(x1,y,nit,nc,max_iter); #print(rrr$coffs)
  }
  coffs_final = matrix(0,nrow=ncol(x),ncol=nc)
  rownames(coffs_final)=colnames(x)
  coffs_final[match(rownames(rrr$coffs),rownames(coffs_final)),]=rrr$coffs
  return(list(coffs=coffs_final,clus=rrr$clus,x=x,y=y,yhat=rrr$yhat,CSMR.model=rrr$mx.model,infoAll=rrr))
}

#' The predict function of the CSMR algorithm.
#'
#' @param CSMR_coffs The coefficient matrix.
#' @param CSMR.model The trained model.
#' @param xnew x variable.
#' @param ynew y variable.
#' @param singleMode A parameter to set the component to one.
#' @return A list object consist of coefficient, clustering membership, data x, external variable y, predicted y based on regression model.
CSMR_predict <- function(CSMR_coffs, CSMR.model, xnew, ynew,singleMode=F){
  if(singleMode==T){  #nc=1 component.
    yhat= t(CSMR_coffs) %*% t(xnew)
    yhat=as.vector(yhat)
    return(yhat)
  }

  conditionalP <- function(res, lambda, sigma) {
    # res: vector of residuals from predictions
    # lambda: vector or scalar of prior probability
    # sigma: variance estimate
    # return: vector of posterior probability
    lambda * dnorm(x = res, mean = 0, sd = sigma)
  }

  #lambda0=do.call(cbind,  CSMR.model$posterior)
  #lambda=apply(lambda0,2,sum)
  lambdaList=lapply(CSMR.model$posterior, sum)  #new prior
  sigmaList = lapply(X =CSMR.model$lmList, FUN = function(m) summary(m)$sigma)
  Bet=CSMR_coffs
  intercept0=lapply(CSMR.model$lmList, FUN=function(li) li$coefficients[1])
  intercept=do.call(cbind, intercept0)
  k=length(CSMR.model$lmList)
  yhatM=matrix(0,length(ynew),k)
  for(i in 1:nrow(xnew)){
    yhatM[i,] = xnew[i,] %*% Bet + intercept
  }
  resList=list()
  for(j in 1:ncol(yhatM)){
    resList[[j]]=ynew-yhatM[,j]
  }
  PList = mapply(conditionalP, resList, lambdaList, sigmaList, SIMPLIFY = FALSE)
  sumsP = rowSums(do.call(cbind, PList))
  ClassList = lapply(PList, function(p) p/sumsP)
  Class.id0=do.call(cbind, ClassList)
  Class.id0[is.nan(Class.id0)] <- 0  #avoid NaN as bug below.
  Class.id=apply(Class.id0, 1, which.max)
  yhat = unlist(sapply(1:length(Class.id), function(gg) yhatM[gg,Class.id[gg]])) #bug solved:return diff length(NA),thus it is a list instead of vector.

  return(yhat)
}


#' Plot the coefficient matrix.
#'
#' @param rrr The result from CSMR function
blockMap <- function(rrr){
  # x = rrr$x
  # mat = matrix(0, nrow(x), ncol(x))
  # rownames(mat) = rownames(x); colnames(mat) = colnames(x)
  #
  # color_bar = c("pink","lightblue","yellow","green")
  #
  # nc = length(table(rrr$clus))
  # for(i in 1:nc){
  #   gene_set = which(rrr$coffs[,i] != 0)
  #   sample_set = which(rrr$clus == i)
  #   mat[sample_set, gene_set] = i
  # }
  rowNum = nrow(rrr$coff)
  rSet = c(rowNum-1,rowNum)
  mat = rrr$coff[-rSet,,drop=F] # new return list

  # if(nc==2){
  #   eg_color = c("aliceblue","blue","dodgerblue")
  # }else{
  #   eg_color = colorRampPalette(c('blue', 'dodgerblue'))(nc)
  # }
  #eg_color = brewer.pal(n=8, name = "PiYG")
  eg_color = colorRampPalette(c('blue', 'aliceblue', 'dodgerblue'))
  #dev.new(width=2, height=4)
  # heatmap.2(mat, trace='none', Rowv=F,Colv=F)
  heatmap.2(mat, trace='none', Rowv = F, Colv = F,
            #col = colorRampPalette(c('blue', 'yellow'))(3),
            col = eg_color,
            margins=c(3,1.5), # ("margin.Y", "margin.X") #marginX can solve the problem overlapping the 'Patient'.
            symkey=FALSE,
            symbreaks=FALSE,
            dendrogram='none',
            density.info='none',
            #denscol="black",
            srtCol=0,
            key=T,
            keysize=1,
            key.title = "Subset",
            labRow=NA,
            #labCol=NA,
            #xlab = 'Gene',
            ylab = 'Feature',
            cexCol=1.2,
            #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
            #key.par=list(mar=c(bottom, left, top, right)),
            key.par=list(mar=c(6,5,3,5)),
            # lmat -- added 2 lattice sections (5 and 6) for padding
            lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(5, 8, 5)
            #lwid control the left/right margin of the heatmap
            )
  # heatmap.2(mat, Rowv=NULL,Colv=NULL,
  #           col = eg_color,
  #           scale="none",
  #           margins=c(3,0), # ("margin.Y", "margin.X")
  #           trace='none',
  #           symkey=FALSE,
  #           symbreaks=FALSE,
  #           dendrogram='none',
  #           density.info='histogram',
  #           denscol="black",
  #           keysize=1,
  #           #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
  #           key.par=list(mar=c(3.5,0,3,0)),
  #           # lmat -- added 2 lattice sections (5 and 6) for padding
  #           lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1)
  # )


}

#-------------------CSMR example----------------

# n=400####need to loop through 100, 200 400
# bet1=bet2=rep(0,101)
# bet1[2:21]=sign(runif(20,-1,1))*runif(20,2,5)
# bet2[22:41]=sign(runif(20,-1,1))*runif(20,2,5)
# bet=rbind(bet1,bet2)
# pr=c(1,1)*0.5
# sigs=c(1,1)####need to loop through 0.5, 1, 2
# tmp_list = simu_data_sparse(n=n,bet=bet,pr=pr,sigma=sigs)
# nit=1
# nc=2
# max_iter=50
# x=tmp_list$x
# y=tmp_list$y
# rrr=CSMR(x,y,nit,nc,max_iter)
