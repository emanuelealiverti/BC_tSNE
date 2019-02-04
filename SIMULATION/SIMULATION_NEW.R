###############################################################
# Batch - corrected t-SNE 
###############################################################

#rm(list=ls())
set.seed(2731)

k.rid = 20
#n = 500
#p = 1000
p = 3000
n = 200

W = matrix(rnorm(p*k.rid), k.rid)
S = matrix(rnorm(n*k.rid), n)
z= sample(1:3, rep = T,size = n)
Z = model.matrix(~as.factor(z))[,-1]
#lambda = rnorm(k.rid)
lambda = matrix(rnorm(k.rid*NCOL(Z)), k.rid)


#lambda = sample(c(-1:2), k.rid, rep=T)
#lambda
#A = matrix(rnorm(n*p,sd=2),n) + S%*%W + 4*z

#A = matrix(rnorm(n*p,sd=2),n) 
#A = (S - lambda * z ) %*% W
A = (S - Z %*% t(lambda) ) %*% W
#A = matrix(rnorm(n*p,sd=2),n) 
#A = (S + lambda*z)%*%W
X = scale(A)
#lattice::levelplot(cor(X))
require(Rtsne)
sn = Rtsne(X = X, initial_dims = 30)
str(sn)
res = list()
res$'Unadjusted' = data.frame(sn$Y)
plot(sn$Y, col = z+1)


library(sOG)
X_tilde = OG(X = X, z = Z, K = 10)
sn_tilde = Rtsne(X = X_tilde$S, pca = FALSE)
res$'Preprocessing Adjustment' = data.frame(sn_tilde$Y)
plot(sn_tilde$Y, col = z+1)


## GRADIENT implementation
s = svd(X, nu = 10, nv = 10)
X_red = scale(s$u %*% diag(s$d[1:10]))
write.table(X_red, file = 'X.txt',row.names = F,col.names = F)
write.table(Z, file = 'z.txt',row.names = F,col.names = F)
system('python ../tsne_python/tsneGRADIENT.py')
#system('python ../tsne_python/tsne.py')
snPY = read.table("X_tilde.txt",col.names = c('X1','X2'))
snPY = data.frame(scale(snPY))
res$'Gradient Adjustment' = snPY
plot(snPY, col = z+1)
str(res)
df = reshape2::melt(res, id.var = c('X1','X2'))
df$z = rep(z, 3)
df$z = factor(df$z)
df$L1 = factor(df$L1, ordered = T, levels = names(res))
require(ggplot2)
ggplot(df) + geom_point(aes(X1,X2,col=z,shape=z),size=5) + 
	facet_wrap(~L1, scales = 'free')+
	theme_bw()+
	theme(legend.position = 'none',
	      axis.title = element_blank(),
	      strip.text.x = element_text(face = 'bold', size = 35))
#ggsave(filename = "../../PAPER/sim.pdf", height = 9, width = 27)
graphics.off()

###############################################################
# Batch - corrected t-SNE - high order 
###############################################################

set.seed(2731)

k.rid = 20
#n = 500
#p = 1000
p = 3000
n = 200

W = matrix(rnorm(p*k.rid), k.rid)
S = matrix(rnorm(n*k.rid), n)
z = rnorm(n)
#lambda = rnorm(k.rid)
#lambda = matrix(rnorm(k.rid*NCOL(Z)), k.rid)


#lambda = sample(c(-1:2), k.rid, rep=T)
#lambda
#A = matrix(rnorm(n*p,sd=2),n) + S%*%W + 4*z

#A = matrix(rnorm(n*p,sd=2),n) 
A = (S - exp(z)*1.5 ) %*% W
#A = (S - exp(Z) %*% t(lambda)) %*% W
#A = matrix(rnorm(n*p,sd=2),n) 
#A = (S + lambda*z)%*%W
X = scale(A)
#lattice::levelplot(cor(X))
require(Rtsne)
sn = Rtsne(X = X, initial_dims = 30)
str(sn)
res = list()
res$'Unadjusted' = data.frame(sn$Y)
plot(sn$Y, col = (z > 0) +1)


library(sOG)
X_tilde = OG(X = X, z = Z, K = 10)
sn_tilde = Rtsne(X = X_tilde$S, pca = FALSE)
res$'Preprocessing Adjustment' = data.frame(sn_tilde$Y)
plot(sn_tilde$Y, col = (z > 0) + 1)


## GRADIENT implementation
s = svd(X, nu = 10, nv = 10)
X_red = scale(s$u %*% diag(s$d[1:10]))
write.table(X_red, file = 'X.txt',row.names = F,col.names = F)
write.table(z, file = 'z.txt',row.names = F,col.names = F)

dir()
system('python ../tsne_python/tsneGRADIENT.py')

#system('python ../tsne_python/tsne.py')
snPY = read.table("X_tilde.txt",col.names = c('X1','X2'))
snPY = data.frame(scale(snPY))
res$'Gradient Adjustment' = snPY
plot(snPY, col = (z>0)+1)
str(res)
df = reshape2::melt(res, id.var = c('X1','X2'))
df$z = rep(z>0, 3)
df$z = factor(df$z)
df$L1 = factor(df$L1, ordered = T, levels = names(res))

require(ggplot2)
ggplot(df) + geom_point(aes(X1,X2,col=z,shape=z),size=5) + 
	facet_wrap(~L1, scales = 'free')+
	theme_bw()+
	theme(legend.position = 'none',
	      axis.title = element_blank(),
	      strip.text.x = element_text(face = 'bold', size = 35))
plot3D::scatter3D(z, res$Unadjusted$X1, res$Unadjusted$X2, col = 1, xlab = 'Systmetatic')
plot3D::scatter3D(z, res$`Gradient Adjustment`$X1, res$`Gradient Adjustment`$X2, col = 1, xlab = 'Systmetatic')
plot3D::scatter3D(z, res$Unadjusted$X1, res$Unadjusted$X2, col = 1, xlab = 'Systmetatic')
install.packages('plot3D')
#ggsave(filename = "../../PAPER/sim2.pdf", height = 9, width = 27)
