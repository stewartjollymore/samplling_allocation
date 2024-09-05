## Creator: Stewart Jollymore
## based on work by Jeff Schinder, Jill Devers PhD.


# This program reads in the following files and
# determines the overall integer sample allocation:
# 1) costs, eligibility, and completion
# 2) domains
# 3) precisions
# 4) population data
# 
# -----------------------------------------------------

library(dplyr)
library(reshape2)
library(ggplot2)

#-----------------------------------------------------
#  Step 1: Read all of the Four data sets into R:
#
#  1. Costs
#  2. Domains
#  3. Precisions
#  4. Data
#-----------------------------------------------------

## Synthaize Data
{
strata_domain_count_synth<-cbind(strata_domain_count[,c(1,2)],
                                 floor(rexp(2070,
                                            1/mean(strata_domain_count$STRDOMSIZE))))

names(strata_domain_count_synth)<-c("DOMAIN","STRATA","STRDOMSIZE")

strata_information_synth<- strata_domain_count_synth %>% 
                            group_by(STRATA) %>% 
                            summarise(STRATSIZE=max(STRDOMSIZE))

domain_info_synth<-strata_domain_count_synth %>% 
                    group_by(DOMAIN) %>% 
                    summarise(DOMSIZE=sum(STRDOMSIZE))

domain_info_synth$precision<-sample(c(0.0007455,0.000031351, 0.0064321),
                                    size = nrow(domain_info_synth),
                                    prob = c(0.2,0.3,0.5),
                                    replace = TRUE)

domain_info_synth$prevalence<-0.5
names(domain_info_synth)<-c("DOMAIN","DOMSIZE", "precision", "prevalence")

costs_synth<-as.data.frame(cbind(seq(1,length(unique(strata_domain_count$STRATA)),1),
                   runif(90,min=(min(costs$COST)),max=(max(costs$COST))),
                   runif(90,min=(min(costs$RESP.RATE)),max=(max(costs$RESP.RATE))),
                   runif(90,min=(min(costs$ELIG.RATE)),max=(max(costs$ELIG.RATE)))))

names(costs_synth)<-c("STRATA","COST","RESP.RATE","ELIG.RATE")


strata_domain_count_synth<-merge(strata_domain_count_synth,
                                 costs_synth,
                                 by="STRATA",
                                 all.x=T)

#Merge on Strata Information
strata_domain_count_synth<-merge(strata_domain_count_synth,
                                 strata_information_synth,
                                 by="STRATA",
                                 all.x=T)

#Merge on Domain Infomration
strata_domain_count_synth<-merge(strata_domain_count_synth,
                                 domain_info_synth,
                                 by="DOMAIN",
                                 all.x=T)
head(strata_domain_count)

strata_domain_count_synth$precision<-(strata_domain_count_synth$precision/1.96)^2

head(strata_domain_count_synth)
}


## Create components and inital lambda
{
## Renames values to coincide with literature making reading 
## equations easier
p=strata_domain_count_synth$prevalence
N_d=strata_domain_count_synth$DOMSIZE
N_h=strata_domain_count_synth$STRATSIZE
N_dh=strata_domain_count_synth$STRDOMSIZE
e=strata_domain_count_synth$ELIG.RATE


## This is varaince of a population proportion by strata component of 
## the Lagrange Multiplier 
popvar<-(p*(N_dh/N_h)*(1-(p*(N_dh/N_h))))

## VARCOMP is the full variance component which comes from solving objective 
## function for Lagrange Multiplier
strata_domain_count_synth$VARCOMP<-((N_h/N_d)^2)*N_h/(N_h-1))*popvar

strata_domain_count_synth<- strata_domain_count_synth %>% arrange(STRATA,DOMAIN)
head(strata_domain_count_synth)

## Calulate Initial Lambda (see Chromy paper - TOP Left of Page 196)
## Initial lambda is derivated in the readme file but comes from 
## maximizing with respect to allocation sample size using partial derivatives
## of the objective function and solving for lambda
i_lambda<- strata_domain_count_synth %>% group_by(DOMAIN) %>%
  mutate(num = (sqrt(VARCOMP*COST)),
         den=precision) %>%
  summarise(num = sum(num),
            den=mean(den)) %>%  # mean is taken here since precision is domain specific
                                # and mean of a constant is its self
  mutate(lambda = (num/(den)),
         initial_lambda=lambda^2)

domain_info<-merge(domain_info,i_lambda,by="DOMAIN",all=T)
domain_info$lambda_i<-domain_info$lambda^2
head(domain_info)

#  Set up for the iteration

store_ss<-c()
store_allocation<-c()
store_lambda<-data.frame(x=c(1:nrow(domain_info)))
head(store_lambda)
}

## Minimization Allgorithm
{
for (i in 1:30) {


  cat(paste0("current iteration: ",i,"\n"))

  strata_domain_count_synth$lambda_i<-NULL
  strata_domain_count_synth$n_h<-NULL
  strata_domain_count_synth$ALLOCATION<-NULL
  strata_domain_count_i<-merge(strata_domain_count_synth,
                               domain_info[,c("DOMAIN","lambda_i")],
                               by="DOMAIN",
                               all.x=T)

  strata_domain_count_i$n<-strata_domain_count_i$VARCOMP*strata_domain_count_i$lambda_i

  store_lambda<-cbind(store_lambda,domain_info$lambda_i)

  hm<- strata_domain_count_i %>%
    group_by(STRATA) %>% 
    summarise(STRATSIZE = mean(STRATSIZE),
              n=sum(n),
              COST=mean(COST),
              RESP.RATE=mean(RESP.RATE))


  hm$n_h<-sqrt(hm$n/hm$COST) #KKT ERROR HERE See Mason 1995 page 20


  hm$noKKTflag<-ifelse(hm$n_h>=hm$STRATSIZE,
                       hm$STRATA,
                       NA)

  strata_domain_count_i<-merge(strata_domain_count_i,hm[,c("STRATA","n_h")],by="STRATA")

  #compute domain level variances
  strata_domain_count_i$STRDOMVAR<-(strata_domain_count_i$VARCOMP/strata_domain_count_i$N2)

  dom_var<- strata_domain_count_i %>% 
    select(DOMAIN,STRDOMVAR,precision) %>% 
    group_by(DOMAIN) %>% 
    summarise(DOM_VAR=sum(STRDOMVAR),
              precision=mean(precision))

  domain_info$DOMVAR<-dom_var$DOM_VAR

  domain_info$SS<-sqrt(domain_info$lambda_i)*(domain_info$DOMVAR-dom_var$precision)^2
  conv.crit =  0.000000001

  domain_info$FIN = ifelse(domain_info$SS < conv.crit,
                           TRUE,
                           FALSE)
  store_ss<-rbind(store_ss,sum(domain_info$SS))

  #TERMINATION CRITERION
  domain_info$lambda_i<-(domain_info$lambda_i)*(domain_info$DOMVAR/dom_var$precision)^2

  domain_info$lambda_i<-ifelse(domain_info$lambda_i < 0.0000001,
                                  0,
                                  domain_info$lambda_i)


  sum_squares<-sum(domain_info$SS)

  #BUILD ALLOCATION AFTER ALLOCATION COMPUTED MHM *NODS
  max_iters<-i


  allocation_interal<- hm %>%
    select(STRATA, n_h, STRATSIZE, RESP.RATE) %>%
    group_by(STRATA) %>%
    summarise(N_sum = sum(ceiling(n_h)), 
              STRATASIZE = mean(STRATSIZE),
              RESP=mean(RESP.RATE))

  allocation_interal$N_sum<-ifelse(allocation_interal$N_sum < 2,
                                   2,
                                   allocation_interal$N_sum)

  allocation_interal$sample_size<-round(allocation_interal$N_sum/allocation_interal$RESP)

  allocation_interal$sample_size<-ifelse(allocation_interal$sample_size >= allocation_interal$STRATASIZE,
                                         allocation_interal$STRATASIZE,
                                         allocation_interal$sample_size)

  store_allocation<-rbind(store_allocation,cbind(sum(allocation_interal$sample_size),sum(allocation_interal$N_sum)))

  if(sum(domain_info$SS) < conv.crit){
    break }
}
}



allocation<- hm %>%
  group_by(STRATA) %>%
  summarise(N_sum = sum(ceiling(n_h)), 
            STRATASIZE = mean(STRATSIZE),
            RESP=mean(RESP.RATE))

allocation$N_sum<-ifelse(allocation$N_sum < 2,
                         2,
                         allocation$N_sum)

allocation$sample_size<-round(allocation$N_sum/allocation$RESP)

allocation$sample_size<-ifelse(allocation$sample_size >= allocation$STRATASIZE,
                               allocation$STRATASIZE,
                               allocation$sample_size)

iter<-paste0("iter",c(1:max_iters))
names(store_lambda)<-c("Domain",iter)
store_lambda<-melt(store_lambda,id="Domain")
store_lambda$value<-store_lambda$value^.5
store_lambda<- store_lambda %>% 
  group_by(Domain) %>% 
  mutate(Iterations=c(1:n()))


store_ss<-data.frame(store_ss)
store_ss$Iterations<-c(1:max_iters)

store_allocation<-data.frame(store_allocation)

store_allocation$Iterations<-c(1:max_iters)


df_lambda<-store_lambda
df_ss<-store_ss
df_allocation<-allocation
max_iters<-max_iters                              

domain_info$RATIO<-(domain_info$iter_lambda/domain_info$init_lambda)*100

domain_info<-merge(domain_info,domains[,c("DOMAIN", "DOMAIN.LABEL")], by="DOMAIN")

domain_info<-merge(domain_info,tool_precisions, by="DOMAIN")

dominfo_out<-domain_info[,c("DOMAIN","DOMAIN.LABEL","DOMSIZE","precision","DOMVAR","RATIO")]

dominfo_out$DOMVAR<-sqrt(dominfo_out$DOMVAR)*1.96

dominfo_out$RATIO<-round(dominfo_out$RATIO,digits=1)
dominfo_out$DOMVAR<-round(dominfo_out$DOMVAR,digits=3)

names(dominfo_out)<-c("Domain","Domain Label","Domain Size","Set Precision","Obtained Precision","Final Lambda")
dominfo<-dominfo_out
strdomcnt<-strata_domain_count


df.s<-dominfo
df.s$DOMAIN.LABEL<-factor(df.s$`Domain Label`,levels=df.s$`Domain Label`[order((df.s$RATIO))])

ggplot(data=df.s)+geom_bar(aes(x=`Domain Label`,y=`Final Lambda`),stat="identity",fill="steelblue",alpha=0.7)+
  coord_flip()+xlab("")+theme_minimal()+ggtitle("Allocation Drivers")+ylab("Initial to Final Lambda Ratio")




ggplot(data=df_ss,aes(x=factor(Iterations),y=store_ss,group=1))+geom_line(alpha=1,size=1)+theme_minimal()+
  xlab("Iteration")+ylab("Sum of Squares")+ggtitle("Sum of Squares vs Iteration")+
  geom_hline(yintercept = conv.crit,linetype=2,color="red",size=1.25,alpha=0.5) #+scale_y_continuous(expand=c(0,0))



df_lambda<-df_lambda[1:max_iters*nrow(domains),]
ggplot(data=store_lambda,aes(x=factor(Iterations),y=value,group=Domain))+geom_path(alpha=0.5)+theme_minimal()+
  xlab("Iteration")+ylab("Lambda")+ggtitle("Lambda vs Iteration") #+scale_y_continuous(expand=c(0,0))
