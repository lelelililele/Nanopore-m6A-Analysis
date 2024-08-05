###lasso回归
library(data.table)
irlnc_sig<-read.table("UC_diffgenesig.txt",header=T)
irlnc_sig<-as.data.frame(irlnc_sig)


library(glmnet)
library(GGally)
library(survival)
rownames(irlnc_sig)<-irlnc_sig$sample
irlnc_sig<-irlnc_sig[,-1]
data<-irlnc_sig[,-c(2)]

x=as.matrix(data[,-c(1,2)])
y<-data[,c('OS.time','OS')]
y$time<-as.double(data$OS.time)
y$status<-as.double(data$OS)
y<-as.matrix(survival::Surv(y$time,y$status))

###training set
set.seed(123)
ind<- sample(2,nrow(data),replace=TRUE,prob=c(0.9,0.1))
df.train<-data[ind==1, ]
df.validate<- data[ind==2, ]

Y = ifelse(data$OS=="1",1,0)
trainY<-Y[ind==1]
validateY<-Y[ind==2]

lasso_fit <- cv.glmnet(x, y, family="cox", type.measure='default', alpha=1,nfolds=10) ##type.measure='default'
plot(lasso_fit, xvar = "norm", label = TRUE)
coefficient <- coef(lasso_fit, s=lasso_fit$lambda.min)
Active.Index <- which(as.numeric(coefficient) != 0)
sig_gene_multi_cox <- rownames(coefficient)[Active.Index]
cox<-as.matrix(coefficient)

#candidates build COX model
formula_for_multivariate <- as.formula(paste0('Surv(OS.time, OS)~', paste(sig_gene_multi_cox, sep = '', collapse = '+')))
multi_variate_cox <- coxph(formula_for_multivariate, data = df.train)
#check if variances are supported by PH hypothesis.
ph_hypo_multi <- cox.zph(multi_variate_cox)
#ggcoxzph(ph_hypo_multi)
#The last row of the table records the test results on the GLOBAL model. Delete it.
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
#Remove variances not supported by ph hypothesis and perform the 2nd regression.
formula_for_multivariate <- as.formula(paste0('Surv(OS.time, OS)~', paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05], sep = '', collapse = '+')))
multi_variate_cox <- coxph(formula_for_multivariate, data = df.train)
library('survminer')
ggforest(multi_variate_cox)
length(multi_variate_cox$coefficients)


multi_variate_cox_2<-step(multi_variate_cox,direction = "both")

ggforest(multi_variate_cox_2)

C_index <- multi_variate_cox_2$concordance['concordance']
if(C_index >= 0.9){
  print('High accuracy')
}else{ 
  if(C_index < 0.9 & C_index >= 0.7){
    print('Medium accuracy')
  }else{
    print('Low accuracy')
  }
}
#calculate the risk score of each sample.



  coxscore<-predict(multi_variate_cox_2,type = "risk",newdata = df.train)
  
  coxrisk<-as.vector(ifelse(coxscore>median(coxscore),"high","low"))
  
  riskdata<-cbind(df.train[,c(1,2)],coxscore,coxrisk)
  
  colnames(riskdata)<-c("OS","OS.time","riskscore","risk")
  
  coxtestscore<-predict(multi_variate_cox_2,type='risk',newdata=df.validate)
  
  coxtestrisk<-as.vector(ifelse(coxtestscore>median(coxscore),"1","0"))
  
  testriskdata<-cbind(df.validate[,c(1,2)],coxtestscore,coxtestrisk)
  
  colnames(testriskdata)<-c("OS","OS.time","riskscore","risk")


library(pROC)

  r_train<-roc(trainY,riskdata$riskscore)
  plot(r_train,main = "Training dataset",col="#0072B5FF",print.auc=TRUE,print.thres=TRUE)
  r_train$auc
  train_result <- coords(r_train, "best")
  train_result$threshold
  
  train_risk<-as.numeric(ifelse(riskdata$riskscore>train_result$threshold,"1","0"))
  
  riskdata<-cbind(riskdata,train_risk)
  
  colnames(riskdata)<-c("OS","OS.time","riskscore","risk","risk_roc")
  
  #
  test_risk<-as.numeric(ifelse(testriskdata$riskscore>train_result$threshold,"1","0"))
  
  testriskdata<-cbind(testriskdata,test_risk)
  
  colnames(testriskdata)<-c("OS","OS.time","riskscore","risk","risk_roc")
  
  #testriskdata$coxtestrisk<-as.numeric(testriskdata$coxtestrisk)
  r_test<-roc(validateY,testriskdata$riskscore)
  
  plot(r_test,main = "Testing dataset",col="#0072B5FF",print.auc=TRUE,print.thres=TRUE)
  
  r_test$auc




  train_surv.cut <- surv_cutpoint(riskdata,time = "OS.time",
                                  event = "OS",
                                  variables = "riskscore")
  summary(train_surv.cut)
  plot(train_surv.cut, "riskscore", palette = "npg")
  surv.cat <- surv_categorize(train_surv.cut) 
  
  fit_curve_test<-survfit(Surv(OS.time,OS)~risk_roc,data = riskdata)
  ##
 #ggsurvplot(fit_curve_test, data = riskdata, legend.title = "survival curve of Pseudogene", legend.labs = c("Low","High"),pval = TRUE,palette =c("#a00020","#4575b4"))
  ggsurvplot(fit_curve_test, data = riskdata, legend.title = "survival curve of Pseudogene", legend.labs = c("Low","High"),pval = TRUE,palette =c("#a00020","#4575b4"), conf.int = TRUE,ncensor.plot = TRUE)
  res_cox<-coxph(Surv(OS.time,OS)~risk_roc,data = riskdata)
  summary(res_cox)
 
ggcoxfunctional(multi_variate_cox_2,data=df.train,point.col = "#4575b4", point.alpha = 1)
ggcoxdiagnostics(multi_variate_cox_2)
ggsurvevents(fit=fit_curve_test)
###ggcoxzph()


 
  ####
  diff_test<-survdiff(Surv(OS.time,OS)~risk,data = testriskdata)
  
  p_valuetest<-signif((1 - pchisq(diff_test$chisq,df=1)),3)
  
  p_valuetest
  
  fit_curve_test<-survfit(Surv(OS.time,OS)~risk_roc,data = testriskdata)##
  
  ggsurvplot(fit_curve_test, data = testriskdata,
             legend.title = "survival curve of validation dataset",
             legend.labs = c("Low","High"),
             pval = TRUE,#Pֵ
             ylab = "Overall survival", 
             xlab = "Time(days)", 
             # Add risk table
             risk.table = TRUE,
             tables.height = 0.2,
             tables.theme = theme_cleantable(),
             palette =c("jco" ),#"jco","jama","nejm",npg","hue","grey","aaas","lancet","jco", "ucscgb","uchicago","simpsons"??"rickandmorty"
             ggtheme = theme(panel.grid =element_blank(),
                             panel.spacing.y = unit(0.5, "mm"),
                             panel.spacing=unit(0.5,"cm"),
                             panel.background = element_rect(fill = 'white',colour = "black")) # Change ggplot2 theme
  )
  
  res_cox<-coxph(Surv(OS.time,OS)~risk_roc,data = testriskdata)
  summary(res_cox)