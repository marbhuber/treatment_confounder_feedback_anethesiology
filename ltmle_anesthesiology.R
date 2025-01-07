##############################################
##############################################
# ltmle_anesthesiology.R
#
# Simulation study of longitudinal
# anesthetic care
# - Only binary random variables
# - No unmeasured confounding
#
# Date: 23 Dec 2024
# Author: Markus Huber (markus.huber@insel.ch)
##############################################
##############################################

# clear workspace and load packages
rm(list=ls())
library(tidyverse)
library(lmtp)
library(ltmle)
library(compareGroups)
library(dagitty)
library(ggdag)

# helper function for binary outcome
rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))

#########################################
# Cohort
# L0: MAP below baseline
# A0: Norepinephrine administration
# L0: MAP below baseline intraoperatively
# A1: Norepinephrine administration
# Y:  MAP below baseline postoperatively

set.seed(1234)
n    <- 10000
L0   <- rbinom(n,1,0.35)
A0   <- rexpit(0.8*L0) 
L1   <- rexpit(0.3*L0 - 1.2*A0)
A1   <- rexpit(-3+0.3+L0 + 0.8*A0 + 0.4*L0 + 0.8+L1)
Y    <- rexpit(0+0.2*L0 - 0.7*A0 - 1.2*A1 + 0.2*L0 + 0.2*L1)
df   <- data.frame(L0, A0, L1, A1, L2 = Y, Y)


########################
# True underlying effect
  
# Treat all
set.seed(1234)
n      <- 1e7
L0     <- rbinom(n,1,0.35)
A0     <- 1
L1     <- rexpit(0.3*L0 - 1.2*A0)
A1     <- 1
Y      <- rexpit(0+0.2*L0 - 0.7*A0 - 1.2*A1 + 0.2*L0 + 0.2*L1)
df.all <- data.frame(L0, A0, L1, A1, L2 = Y, Y)
 
# Treat none
set.seed(1234)
n       <- 1e7
L0      <- rbinom(n,1,0.35)
A0      <- 0
L1      <- rexpit(0.3*L0 - 1.2*A0)
A1      <- 0
Y       <- rexpit(0+0.2*L0 - 0.7*A0 - 1.2*A1 + 0.2*L0 + 0.2*L1)
df.none <- data.frame(L0, A0, L1, A1, L2 = Y, Y)
  
df.true <- mean(df.all$Y)-mean(df.none$Y)

########
########
# Figure
# DAG

L0.stat <- paste0("[",format(round(100*sum(df$L0)/nrow(df),1),nsmall=1),"%]")
L1.stat <- paste0("[",format(round(100*sum(df$L1)/nrow(df),1),nsmall=1),"%]")
Y.stat  <- paste0("[",format(round(100*sum(df$Y)/nrow(df),1),nsmall=1),"%]")
A0.stat <- paste0("[",format(round(100*sum(df$A0)/nrow(df),1),nsmall=1),"%]")
A1.stat <- paste0("[",format(round(100*sum(df$A1)/nrow(df),1),nsmall=1),"%]")

mydag <- dagitty("dag {
L0 [pos=\"0,0.01\"]
A0 [pos=\"1,0.03\"]
L1 [pos=\"2,0.01\"]
A1 [pos=\"3,0.03\"]
Y  [pos=\"4,0.01\"]
L0 -> A1
L0 -> L1
L0 -> Y
L0 -> A0
L1 -> A1
L1 -> Y
A0 -> L1
A0 -> A1
A0 -> Y
A1 -> Y
}") %>% 
  tidy_dagitty() %>%
  dag_label(
    labels = c(
      "L0" = paste0("MAP\nbelow preoperative\n(after induction)\n",L0.stat),
      "L1" = paste0("MAP\nbelow preoperative\n(intraoperative)\n",L1.stat),
      "A0" = paste0("Norepinephrine\n(intraoperative)\n",A0.stat),
      "A1" = paste0("Norepinephrine\n(intraoperative)\n",A1.stat),
      "Y"  = paste0("MAP\nbelow preoperative\n(postoperative)\n",Y.stat)))
      
fig1A <- mydag$data %>% 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_dag_text(col = "black") +
  theme_dag()+
  geom_dag_edges_arc(data = function(x) filter(x, (name == "L0")),curvature = 0.002)+
  geom_dag_edges_arc(data = function(x) filter(x, (name == "A0")),curvature = 0.001)+
  geom_dag_edges_arc(data = function(x) filter(x, (name == "A1")),curvature = 0.001)+
  geom_dag_edges_arc(data = function(x) filter(x, (name == "L1")),curvature = 0.001)+
  geom_text(data = function(x) filter(x, (name == "A0")), aes(label = label),check_overlap = T,nudge_y=0.005)+
  geom_text(data = function(x) filter(x, (name == "A1")), aes(label = label),check_overlap = T,nudge_y=0.005)+
  geom_text(data = function(x) filter(x, (name == "L0")), aes(label = label),check_overlap = T,nudge_y=-0.006,nudge_x=0.0)+
  geom_text(data = function(x) filter(x, (name == "L1")), aes(label = label),check_overlap = T,nudge_y=-0.006,nudge_x=0.0)+
  geom_text(data = function(x) filter(x, (name == "Y")), aes(label = label),check_overlap = T,nudge_y=-0.006,nudge_x=0.0)+
  scale_y_continuous(limits = c(-0.002,.04))+
  scale_x_continuous(limits = c(-0.3,4.3))

fig1A
ggsave("figure1.jpg",dpi=1200,width=9,height=6)
ggsave("figure1.pdf",dpi=1200,width=9,height=6)

###################
###################
# Figure
# Treatment effects

set.seed(1234)
df.analysis <- c()

n.boot = 1000

for (myboot in 1:n.boot){

  print(myboot)
  
  # helper data frame
  rm(dummy)
  dummy = df %>%
    sample_n(size=nrow(df),replace=T) %>% 
    mutate(A.sum = factor(A0+A1)) %>% 
    mutate(L.sum = (L0+L1))

  # only treatments as covariates
  rm(mymod,myeff)
  mymod = glm(Y~A0+A1,dummy %>% mutate(A0=factor(A0),A1=factor(A1)),family=binomial())
  myeff = mean(
    predict(mymod,newdata=mutate(dummy,A0=factor(1),A1=factor(1)),type="response")-
    predict(mymod,newdata=mutate(dummy,A0=factor(0),A1=factor(0)),type="response")
  )
  df.analysis <- rbind(df.analysis,data.frame(
    Type = "Traditional Regression (unadjusted)",
    Analysis = "Y~A0+A1",
    estimate = myeff
  ))
  
  # only sum of treatments as covariates
  rm(mymod,myeff)
  mymod = glm(Y~A.sum,dummy,family=binomial())
  myeff = mean(
    predict(mymod,newdata=mutate(dummy,A.sum=factor(2)),type="response")-
    predict(mymod,newdata=mutate(dummy,A.sum=factor(0)),type="response")
  )
  df.analysis <- rbind(df.analysis,data.frame(
    Type = "Traditional Regression (unadjusted)",
    Analysis = "Y~A(1+2)",
    estimate = myeff
  ))
  
  # confounding adjustment (as factors)
  rm(mymod,myeff)
  mymod = glm(Y~A0+A1+L0+L1,dummy %>% mutate(A0=factor(A0),A1=factor(A1)),family=binomial())
  myeff = mean(
    predict(mymod,newdata=mutate(dummy,A0=factor(1),A1=factor(1)),type="response")-
    predict(mymod,newdata=mutate(dummy,A0=factor(0),A1=factor(0)),type="response")
  )
  df.analysis <- rbind(df.analysis,data.frame(
    Type = "Traditional Regression (adjusted)",
    Analysis = "Y~A0+A1+L0+L1",
    estimate = myeff
  ))
  
  # confounding adjustment (as factors)
  rm(mymod,myeff)
  mymod = glm(Y~A.sum+L0+L1,dummy,family=binomial())
  myeff = mean(
    predict(mymod,newdata=mutate(dummy,A.sum=factor(2)),type="response")-
      predict(mymod,newdata=mutate(dummy,A.sum=factor(0)),type="response")
  )
  df.analysis <- rbind(df.analysis,data.frame(
    Type = "Traditional Regression (adjusted)",
    Analysis = "Y~A(1+2)+L0+L1",
    estimate = myeff
  ))
  
  # confounding adjustment (cumulative)
  rm(mymod,myeff)
  mymod = glm(Y~A0+A1+L.sum,dummy %>% mutate(A0=factor(A0),A1=factor(A1)),family=binomial())
  myeff = mean(
    predict(mymod,newdata=mutate(dummy,A0=factor(1),A1=factor(1)),type="response")-
      predict(mymod,newdata=mutate(dummy,A0=factor(0),A1=factor(0)),type="response")
  )
  df.analysis <- rbind(df.analysis,data.frame(
    Type = "Traditional Regression (adjusted)",
    Analysis = "Y~A0+A1+L(0+1)",
    estimate = myeff
  ))
  
  # confounding adjustment (cumulative)
  rm(mymod,myeff)
  mymod = glm(Y~A.sum+L.sum,dummy,family=binomial())
  myeff = mean(
    predict(mymod,newdata=mutate(dummy,A.sum=factor(2)),type="response")-
      predict(mymod,newdata=mutate(dummy,A.sum=factor(0)),type="response")
  )
  df.analysis <- rbind(df.analysis,data.frame(
    Type = "Traditional Regression (adjusted)",
    Analysis = "Y~A(1+2)+L(0+1)",
    estimate = myeff
  ))
  
  # TMLE
  dummy <- dummy %>% select(-A.sum,-L.sum,-L2)
  rm(a)
  a = ltmle(dummy, Anodes=c("A0", "A1"), Lnodes=c("L0","L1"), Ynodes="Y", abar=list(c(1,1),c(0,0)), SL.library="glm") %>% summary()
  df.analysis <- rbind(df.analysis,data.frame(
    Type = "Accounting for treatment-confounder feedback",
    Analysis = "Targeted Maximum\nLikelihood Estimation\n(TMLE)",
    estimate = a$effect.measures$ATE$estimate
  ))
  
  rm(a)
  a = ltmle(dummy, Anodes=c("A0", "A1"), Lnodes=c("L0","L1"), Ynodes="Y", abar=list(c(1,1),c(0,0)), SL.library="glm") %>% summary(estimator="iptw")
  df.analysis <- rbind(df.analysis,data.frame(
    Type = "Accounting for treatment-confounder feedback",
    Analysis = "G-methods\n(IPTW)",
    estimate = a$effect.measures$ATE$estimate
  ))
  
  rm(a)
  a = ltmle(dummy, Anodes=c("A0", "A1"), Lnodes=c("L0","L1"), Ynodes="Y", abar=list(c(1,1),c(0,0)), SL.library="glm",gcomp=TRUE) %>% summary(estimator="gcomp")
  df.analysis <- rbind(df.analysis,data.frame(
    Type = "Accounting for treatment-confounder feedback",
    Analysis = "G-methods\n(G-computation)",
    estimate = a$effect.measures$ATE$estimate
  ))
  
}

df.analysis %>% 
  group_by(Analysis) %>% 
  mutate(
    mymean  = mean(estimate),
    mylower = quantile(estimate,c(0.025)),
    myupper = quantile(estimate,c(1-0.025))
  ) %>% 
  select(Type,contains("my"),Analysis) %>% 
  unique() %>% 
  mutate(
    mylabel = paste0(
      format(round(100*mymean,0),nsmall=0),
      "%, 95%-CI:\n(",
      format(round(100*mylower,0),nsmall=0),
      "% to ",
      format(round(100*myupper,0),nsmall=0),
      "%)"
    )
  ) %>% 
  mutate(
    Analysis = factor(Analysis,levels = c(
      "Y~A0+A1",
      "Y~A(1+2)",
      "Y~A0+A1+L0+L1",            
      "Y~A(1+2)+L0+L1" ,              
      "Y~A0+A1+L(0+1)",
      "Y~A(1+2)+L(0+1)",
      "G-methods\n(IPTW)",
      "G-methods\n(G-computation)",
      "Targeted Maximum\nLikelihood Estimation\n(TMLE)" 
    ))
  ) %>% 
  ggplot(aes(x=mymean,xmin=mylower,xmax=myupper,y=Analysis,color=Type,label=mylabel))+
  geom_vline(xintercept = df.true,color="black",linetype="dashed")+
  geom_point(size=3)+
  geom_errorbarh(height=0.2)+
  theme_classic()+
  scale_x_continuous(labels = scales::percent,n.breaks=10)+
  theme(
    plot.margin = margin(0.5, 2, , , "cm"),
    legend.position = "bottom",
    text = element_text(size=14)
  )+
  ylab("Method")+
  xlab("Treatment effect\n(treat always versus treat none)")+
  ggsci::scale_color_jama()+
  coord_flip()+
  scale_y_discrete(guide = guide_axis(n.dodge=2))+
  annotate(geom="text", y=1.6, x=df.true-0.005, label="True, simulated value (-40.5%)",
           color="black")+
  geom_text(aes(x=-0.31),check_overlap = F,size=3.5)

ggsave("figure2.jpg",dpi=600,width=10.5,height=7)
ggsave("figure2.pdf",dpi=600,width=10.5,height=7)


