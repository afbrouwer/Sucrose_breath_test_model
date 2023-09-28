library(deSolve)
library(ggplot2)
library(ggpubr)

model = function(t,par){
  rho = par[1]
  pi = par[2]/par[1]
  kappa = exp(par[3])/(1+exp(par[3]))
  
  if (par[1]==par[2]){ #if pi == 1, the model simplifies
    y= 50 * kappa * rho^3 * t^2 * exp(-rho *t)
  } else{ #if pi =/= 1, use the model simplifies
    y= 100 * kappa*pi*rho/((pi-1)^2)*(exp(-pi*rho*t)+((pi*rho-rho)*t-1)*exp(-rho*t))
  }
  
  return(y)
}

sim_model = function(times,data,params){
  
  # Simulate the model
  sim = cbind(times, model(times,params))

  #Find where simulation matches the data times
  index=which(is.element(round(sim[,1],digits=2),round(data$Time,digits=2)))

  n = nrow(data)
  sd = 0.555
  
  NLL = ((n/2)*log(2*pi) + n*log(sd) + 
           (1/(2*sd^2))*sum((data$PDRr-sim[index,2])^2))
  
  return(NLL)
  
}

fitfun = function(par,data){
  
  NLL = sim_model(seq(0,8,0.01),data,par)

  return(NLL)
}

data_all = read.csv("Data_SFG.csv")
IDs = paste0("SFG_",1:19)

results = matrix(0, nrow = 19, ncol = 10)
colnames(results) = c("Participant","n","NLL", "SIC","par[1]","par[2]","par[3]","rho","pi","kappa")


#Set initial parameter guess. 
#The model can be sensitive to the initial guess, so multiple guesses may be needed.
param_init = c(2,0.5,1)

for (j in 1:19) {
  id = IDs[j]

  print(id)
  
  data = data_all[data_all$Participant==id & data_all$Label == "S",]
  
  mle=optim(c(2,0.5,1),fitfun,data = data,control = list(maxit=10000))
  while (mle$convergence!=0){
    mle=optim(mle$par,fitfun,data = data,control = list(maxit=1000))
  }

  results[j,] = c(id, nrow(data), mle$value, 
                  2*mle$value + log(nrow(data))*length(mle$par),
                  mle$par[1], mle$par[2], mle$par[3],
                  mle$par[1], mle$par[2]/mle$par[1], exp(mle$par[3])/(1+exp(mle$par[3])))

}


write.csv(results,"Results_S.csv")

################################################################################
## Plotting output

data_all = read.csv("Data_SFG.csv")
IDs = paste0("SFG_",1:19)
params = read.csv("Results_S.csv")[,6:8]

model = function(t,par){
  rho = par[1]
  pi = par[2]/par[1]
  kappa = exp(par[3])/(1+exp(par[3]))
  
  if (par[1]==par[2]){ #if pi == 1, the model simplifies
    y= 50 * kappa * rho^3 * t^2 * exp(-rho *t)
  } else{ #if pi =/= 1, use the model simplifies
    y= 100 * kappa*pi*rho/((pi-1)^2)*(exp(-pi*rho*t)+((pi*rho-rho)*t-1)*exp(-rho*t))
  }
  
  return(y)
}

sim_plot = function(data, sim){
  p = ggplot() +
    geom_point(aes(x=data$Time,y=data$PDRr,bg="S"),pch=21,size=1,stroke=0.5) +
    geom_line(aes(x=sim[,1],y=sim[,2],col = "S"),alpha=0.8) +
    scale_color_manual(guide=guide_legend("Model"),values = c("#9c0046")) +
    scale_fill_manual(guide="none",values = c("#9c0046"))+
    xlab("Time (hr)")+ylab("PDRr (1/hr)        ")+
    theme_classic()
}


plot_function=function(j){
  id = IDs[j]
  
  data = data_all[data_all$Participant==id & data_all$Label == "S",]
  
  par= unlist(params[j,])

  times=seq(0,8,0.01)
  sim = cbind(times,model(times,par))
  
  p  =  sim_plot(data, sim)
  
  return(p)
  
}

plot_list=lapply(1:19,plot_function)


p_out = ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
                  plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],
                  plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
                  plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],
                  plot_list[[17]],plot_list[[18]],plot_list[[19]],nrow=5, ncol=4,
                  common.legend = TRUE,labels="auto")
ggsave("Figure_S.pdf",width=6.5,height=7)
