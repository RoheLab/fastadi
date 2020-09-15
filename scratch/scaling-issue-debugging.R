require(RSpectra)
require(Matrix)

##############################
## subfunctions for fastAdi ##
##############################

#########################################
#  eigInit and its subfunction Mx
#  are to make the svd in the initalizer.
#  given sparse M compute eigen of
#  t(M)%*%M/p^2 - (1-p)*diag(t(M)%*%M)
#  where p is the mean nnz of M.

eigInit = function(M, r){
  #  given sparse M compute eigen of
  #  t(M)%*%M - (1-p)*diag(t(M)%*%M)
  #  where p is the mean nnz of M.

  return(
    eigs_sym(
      Mx, r, n = ncol(M),
      args = list(M=M,
                  p = nnzero(M)/prod(dim(M))
      )
    )
  )
}

Mx = function(x, args)
{
  ## Has the effect of  [ t(M)%*%M - (1-p)*diag(t(M)%*%M) ] *x
  return(drop(crossprod(args$M, args$M%*%x)/args$p^2 - (1-args$p)*Diagonal(ncol(args$M),colSums(args$M^2))%*%x/args$p^2))
}

#######################################################################################
# svdSUV and its subfunctions Ax and Atx
# are source files for taking the svd of
# M = S + U%*%V^T
# where S is sparse with nonzeros that match the Data, it is the residuals:
# S = Data - U%*%V^T  (where Data is sparse and this operation is only on the nonzeros)
# Thank you to Yixuan Qiu for help.

svdSUV = function(U,V, Data){
  # U and V are skinny matrices, probably the parameters.
  # Data is a sparse data matrix.
  # returns svd of
  # M = S + U%*%V^T
  # where S is sparse with nonzeros that match the Data
  # S = Data - U%*%V^T  (where Data is sparse and this operation is only on the nonzeros)

  S = graphMult(U, V, Data, resid = T)
  args = list(S = S, U=U, V=V)
  ei = eigs_sym(f, ncol(args$U), n = ncol(args$S), args = args)
  u = as.matrix(Ax(ei$vec, args))
  u = apply(u,2,function(x) return(x/sqrt(sum(x^2))))
  return(list(u =u , v = ei$vec, d = sqrt(ei$val)))

}

Ax = function(x, args)
{
  ## Has the effect of A * x,
  return(drop(args$S %*% x + args$U%*%crossprod(args$V,x)))
}

Atx = function(x, args)
{
  ## Has the effect of A' * x,
  return(drop(t(args$S) %*% x + args$V%*%crossprod(args$U,x)))
}

f = function(x, args)
{
  Atx(Ax(x, args), args)
}

##########################################################################
# graphMult is a subfunction in several pieces for computing the residuals
# S = G - U%*%t(V)
# but only computing it on the nonzero elements of G.

graphMult = function(U,V,G, resids = F){
  # U and V are matrices with U %*% V^T defined.
  # G is a sparse matrix
  # if resids = F, this function computes the elements of U%*%V^T on *only* non-zero elements of G.
  # if resids = T, this function computes the elements of G - U%*%V^T on *only* non-zero elements of G.

  mt = as(G, "dgTMatrix")
  # the indices for which we want to compute the matrix multiplication:
  i = mt@i+1
  j = mt@j+1

  # the rows of U and the columns of V for which we want to compute the multiplication:
  left= U[i,]
  right =  V[j,]  #  Note:  there is most certainly a more memory efficient way to do this.

  # the inner products to compute the elements of the U%*%V:
  uv = rowSums(left*right)

  if(!resids) return(sparseMatrix(i=i, j=j, x = uv))
  if(resids) return(sparseMatrix(i=i, j=j, x = mt@x - uv))

}

######################################################################
# The initializer currently uses a different "sign" matching function,
#  which is faster.
# If the old matching functions need to be adapted, then
# they functions need to be adapted to take an svds and return an svds.
#   Also, they cannot use the (dense) ind matrix.

positive = function(x){
  x[x<0] = .001
  return(x)
}


fastInitial=function(M.p,r){
  ei = eigInit(M.p, r)

  p = nnzero(M.p)/prod(dim(M.p))
  tau = (sum(M.p@x^2)/p-sum(ei$val))/(ncol(M.p)-r)
  s = list(u=c(),d= c(),v=c())
  s$d=sqrt(positive(ei$val-tau))
  s$v=ei$vectors
  s$u=eigInit(t(M.p), r)$vec

  B = svd(t(s$u)%*%(M.p%*%s$v))
  s$u = s$u %*%B$u
  s$v = s$v %*% B$v
  return(s)
  return(M.hat)}


####################################################
# threshSVD is the key internal function for fastAdi

threshSVD= function(so, M.p, frobMp){

  # uses sparse+lowrank svd
  # computes thresholding parameters
  # returns the thresholded SVD
  # M.p is the data
  # frobMp only needs to be computed once.
  # It is the squared frobeinus norm of M.p (frobMp = sum(M.p@x^2))
  d = ncol(M.p)
  s = svdSUV(so$u, so$v%*%diag(so$d), M.p)
  MTM = frobMp + sum(so$d^2) - sum(graphMult(so$u,so$v%*%diag(so$d),M.p)@x^2)
  tau = (sum(MTM) - sum(s$d^2))/(d-r)
  s$d=sqrt(positive(s$d^2 - tau))
  return(s)}

#############################
## fastAdi (main function) ##
#############################

fastAdi = function(M.p,r,tol=1e-07,itmax=200){

  n = nrow(M.p); d = ncol(M.p)
  M.p = as(M.p, "dgCMatrix")
  p.hat = nnzero(M.p)/prod(dim(M.p))
  s = fastInitial(M.p, r)

  frobMp = sum(M.p@x^2)
  itr=0; error=Inf
  oldResids = M.p@x
  resids = graphMult(s$u,s$v%*%diag(s$d),M.p, resids = T)@x

  while(error> tol){
    oldResids = resids
    s=threshSVD(s,M.p,frobMp)
    resids = graphMult(s$u,s$v%*%diag(s$d),M.p, resids = T)@x
    print(sd(resids))
    error=sum((oldResids-resids)^2)/sum(resids^2)

    itr=itr+1
    if(itr%%10==0)cat(".")
    if(itr>itmax){break}
  }
  return(s)
}

n = 500
d = 100
r = 5

A = matrix(runif(n*r,-5,5),n,r)
B = matrix(runif(d*r,-5,5),d,r)
M0 = A %*% t(B)

err = matrix(rnorm(n*d),n,d)
Mf = M0+err

p = 0.1
y = matrix(rbinom(n*d,1,p),n,d)
dat = Mf*y

out = fastAdi(M.p=dat,r=r)
str(out)

MM <- Matrix(dat)

error <- fastadi:::masked_approximation(out, MM) - MM

sum(error^2)

out2 <- fastadi::sparse_adaptive_impute(MM, r)


error2 <- fastadi:::masked_approximation(out2, MM) - MM

# okay, so sparse adaptive impute appears to be working reasonably well

sum(error2^2)

Matrix(dat)

A

sd <- function(U, V) {
  norm(tcrossprod(U) - tcrossprod(V), type = "F")^2
}

sd(A, out$u)
sd(B, out$v)

A
out$u


true <- svd(tcrossprod(A, B))

true$d[1:r]
out2$d
