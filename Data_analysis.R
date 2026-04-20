###### libraries ######
library(Rfast)
library(Matrix)
library(highmean)
library(readxl)

######### set up data ###########
load("FFPE_ovarian_data.RData")
load("FFPE_ovarian_SR.RData")
load("pathways.RData")

DataX <- t(data1[rownames(data1) %in% union(go$HALLMARK_INTERFERON_GAMMA_RESPONSE, go$HALLMARK_INTERFERON_ALPHA_RESPONSE), senref==1])
DataY <- t(data1[rownames(data1) %in% union(go$HALLMARK_INTERFERON_GAMMA_RESPONSE, go$HALLMARK_INTERFERON_ALPHA_RESPONSE), senref==0])

ptrc<- read_xlsx("Supplementary_table_4.xlsx", sheet=2)
sig<- ptrc[ptrc$fdr < 0.1,]
sig.up<- sig$pathway[sig$updown=="up"]
sig.dn<- sig$pathway[sig$updown=="down"]

go.up<- go1[names(go1) %in% sig.up]
go.dn<- go1[names(go1) %in% sig.dn]

alpha<- 0.05

count.dn<- NULL;
pvals.dn<- NULL
for(i in 1:length(go.dn)){
  DataY <- t(data1[rownames(data1) %in% go.dn[[i]],  senref==1])
  DataX <- t(data1[rownames(data1) %in% go.dn[[i]], senref==0])
  
  X<-DataX  ##### Data matrix of order n1X p
  Y<-DataY #### ##### Data matrix of order n2X p
  p<-ncol(DataX)
  n1<-nrow(DataX) ### sample size for popu1
  n2<-nrow(DataY) ### sample size for popu2
  
  #######################################################
  l <-p^{1/2}  # window size   we choose epsilon=1/2   Based on theorem 2.4
  
  fun_rearrange <- function(R){
    p <- ncol(R)
    RL<-R
    
    RL[lower.tri(RL)]<-0  #### Upper triangular matrix
    
    for (i in 1:(p-2)){
      sl<-seq(1,i,1)
      rl<-cbind(RL[i,][-sl],RL[i,][-sl]<0)
      rl[,2]<-replace(rl[,2],rl[,2]==1,-1)
      rl[,2]<-replace(rl[,2],rl[,2]==0,1)
      rl<-cbind(abs(rl[,1]),rl[,2])
      rl<-rl[order(rl[,1],decreasing=TRUE),]
      RL[i,(i+1):p]<-as.vector(rl[,1]*rl[,2])
    }
    R.new <- RL + t(RL)
    diag(R.new)<- diag(R) #### rearranged matrix
    return(R.new)
  }
  
  f<- function(u){
    return(u*(1/4 + 1/(2*pi)*asin(u)) + (sqrt(1 - u^2) - 1)/(2*pi))
  }
  
  W1<-apply(data.frame(cbind((colMeans(X) - colMeans(Y))/sqrt(apply(X, 2, var)/n1+ apply(Y, 2, var)/n2), rep(0,p))), 1, max)
  W <- mean(W1)
  
  ### mu Hat of W under the null hypothesis
  WMu<- sqrt(1/(2*pi))
  
  ################### variance estimate of the W ##################################
  
  eta.matrix<- mat.mult(mat.mult(diag(sqrt(1/diag(cov(X)/n1+cov(Y)/n2)),nrow=p),(cov(X)/n1+cov(Y)/n2)),diag(sqrt(1/diag(cov(X)/n1+cov(Y)/n2)),nrow=p))
  
  diag(eta.matrix)<-0
  
  gamma.hat<-apply(eta.matrix,2,f) ####### computes gamma.hat matrix with diagonal elemnst 0
  
  diag(gamma.hat)<-1/2*(1-1/pi)  ####### computes gamma.hat matrix  
  
  
  M<-fun_rearrange(gamma.hat)
  
  cvv<-sum(band(gamma.hat, 1, l))
  
  WVar <-(sum(diag(M))+2*cvv)/p       # estimate variance of p^{1/2} W
  
  Tobs<- sqrt(p)*(W - WMu) / sqrt(WVar)
  
  pvals.dn <- c(pvals.dn, (1 - pnorm(Tobs, 0, 1)))

  ########## other methods ###########  
  SD<-apval_Sri2008(X, Y)                               
  CQ<-apval_Chen2010(X,Y)
  CXL<-apval_Cai2014(X,Y)
  
  if((SD$pval< 2*alpha) & (sum(colMeans(X) - colMeans(Y))>0)) count1<-1 else count1<-0    ### count1<-1 means rejection of H_0 by the method SD
  if((CQ$pval< 2*alpha)& (sum(colMeans(X) - colMeans(Y))>0)) count2<-1 else count2<-0       ### count2<-1 means rejection of H_0  by the method CQ
  if((CXL$pval< 2*alpha) & (sum(colMeans(X) - colMeans(Y))>0)) count3<-1 else count3<-0      ### count3<-1 means rejection of H_0 by the method CXL
  
  count.dn<- rbind(count.dn, c(count1, count2, count3))
}
df.dn<- cbind.data.frame(matrix(round(pvals.dn, 3),length(go.dn),1), count.dn)
rownames(df.dn)<- names(go.dn)
colnames(df.dn)<- c("SMC", "SD", "CQ", "Cai")

########### overlap plot ############
df<- rbind.data.frame(df.up, df.dn)
path.list<- list("SMC"= rownames(df)[df$SMC < 0.05], "SD"= rownames(df)[df$SD ==1], "CQ"=rownames(df)[df$CQ==1], "Cai"=rownames(df)[df$Cai==1])

upset(fromList(path.list), order.by = "freq", main.bar.color = "black", sets.bar.color = "black", text.scale = 2)

mat<- df[rownames(df) %in% unique(unlist(path.list)),]
mat2<- mat
mat2$SMC<- ifelse(mat$SMC < 0.05, 1, 0)
mat3<- mat2[c(4, 12, 28, 31, 33, 39, 88, 109, 113, 114, 117, 120, 135, 137, 138, 142, 144, 170:171, 216, 219, 252),]
Heatmap(as.matrix(mat3), col = colorRamp2(c(0,1), c("white","red")), show_column_dend = F, rect_gp = gpar(col = "white", lwd = 1), border_gp = gpar(col = "black", lty = 1), height = nrow(mat3)*unit(7, "mm"), width = ncol(mat3)*unit(7, "mm"),show_row_dend = F, cluster_columns = F, cluster_rows = T, row_order = rownames(mat3), column_order = colnames(mat3), show_row_names = T, show_column_names = T, show_heatmap_legend=F)

########### heatmap of immune pathway genes ############
mat.ord<- cbind(t(DataX), t(DataY))
mat.ord4<- t(apply(mat.ord, 1, function(x) (x-mean(x))/sd(x)))
senref<- as.factor(senref)

ha<- HeatmapAnnotation(response = senref[c(which(senref==0), which(senref==1))])

Heatmap(mat.ord4, show_column_dend = F, show_row_dend = F, cluster_columns = F, cluster_rows = T, row_order = rownames(mat.ord4), column_order = colnames(mat.ord4), top_annotation = ha, show_row_names = F, show_column_names = F, show_heatmap_legend=F)

