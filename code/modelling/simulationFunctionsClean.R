


# Simulation functions ----------------------------------------------------


### Sample mutation times ####

#after sampling Wpre, weak/strong arrive at rate bpre*muwd*We^(lpre*t)
#prob weak has not occurred is exp(-bpre*muwd*Wpre/lpre*( exp(lpre*t)-1 ) )
#prob mut has occured is 1-exp(-bpre*muwd*Wpre/lpre*( exp(lpre*t)-1 ) )
#invert this for inverse quantile sampling
sample_time_driverclone = function(bp,dp,Wp,lp,mud){
  nud = bp*mud
  unifval = runif(1)
  sampletime = 1/lp * log(1+ (Wp*nud/lp)^(-1)*log( (1-unifval)^(-1)   )   )
  return(sampletime)
}


#### clone sizes given timings and random amps ####
approx_clonesize <- function(t,l,W,t0){
  #for clone initiated as time t0, growing as ~W*exp(l(t-t0))
  if (t<t0){
    return(0)
  } else {
    
    return(W*exp(l* ( t - t0) ) )
  }
}



#sizes of precursor, weak, strong as function of time
approx_sim_all = function(t,
                          lp,
                          Wp,
                          lwd,
                          Wwd,
                          twd,
                          lsd,
                          Wsd,
                          tsd){
  
  return( approx_clonesize(t,lp,Wp,0)+
            approx_clonesize (t,lwd,Wwd,twd)+
            approx_clonesize (t,lsd,Wsd,tsd))
}

#time population hits size N
approx_sim_find_time_hitN =  function(N,
                                      lp,
                                      Wp,
                                      lwd,
                                      Wwd,
                                      twd,
                                      lsd,
                                      Wsd,
                                      tsd){
  
  
  rooteqn = function(t){
    mall = approx_sim_all(t,
                          lp,
                          Wp,
                          lwd,
                          Wwd,
                          twd,
                          lsd,
                          Wsd,
                          tsd)
    if (is.na(mall) == F){
      return( mall - N)
    } else {
      return (10^300 - N)
    }
  }
  
  
  #if this fails it is likely to be because the effective birth rate for precursors
  #is so close to death rate-precur, that dp/lp>>1 and so population reaches 
  #N at roughly 0 time or even negative time
  return( suppressWarnings(uniroot(rooteqn,lower = .001, upper = 100*1/lp*log(N))$root))
}

#population sizes over size grid until
#total population hits N
simApproxSizegrid = function(bp,
                             dp,
                             muwd,
                             musd,
                             swd,
                             ssd,
                             sizegrid){
  
  # sizegrid = 10^seq(3,7)
  # bp = 1
  # dp = 0
  # muwd = .1
  # musd = .01
  # swd = 1
  # ssd = 3
  # N = 10^5
  
  #derived parameters
  bpeff = bp *(1 - muwd - musd)
  lp = bpeff - dp
  
  bwd = lp*(1+swd)+dp
  lwd = bwd-dp
  muwdSurv = muwd*(1-dp/bwd) #only care about surviving mutations.
  
  bsd = lp*(1+ssd)+dp
  lsd = bsd-dp
  musdSurv = musd*(1-dp/bsd) #only care about surviving mutations.
  
  #sample random amplitudes
  Wp = rexp(1,1-dp/bp)
  Wwd = rexp(1,1-dp/bwd)
  Wsd = rexp(1,1-dp/bsd)
  
  #sample mutation times
  twd = sample_time_driverclone(bp,dp,Wp,lp, muwdSurv)
  tsd = sample_time_driverclone(bp,dp,Wp,lp, musdSurv)
  
  #clonesizes
  clonesizes_hitN = sapply(sizegrid, function(N) {
    t = approx_sim_find_time_hitN(N,
                                  lp,
                                  Wp,
                                  lwd,
                                  Wwd,
                                  twd,
                                  lsd,
                                  Wsd,
                                  tsd)
    return(c(N,approx_clonesize(t,lp,Wp,0),
             approx_clonesize(t,lwd,Wwd,twd),
             approx_clonesize(t,lsd,Wsd,tsd)))})
  
  return(clonesizes_hitN)
}

#wrap simApproxSizegrid over multiple runs
simApproxSizegridRuns <- function(bp,
                                  dp,
                                  muwd,
                                  musd,
                                  swd,
                                  ssd,
                                  runs,
                                  sizegrid,
                                  seednum){
  
  
  set.seed(seednum)
  
  simoutput = lapply(1:runs, function(r) {
    # print(r)
    simout_df = simApproxSizegrid(bp,
                                  dp,
                                  muwd,
                                  musd,
                                  swd,
                                  ssd,
                                  sizegrid) %>% t() %>%  data.frame()
    colnames(simout_df)  = c("N","pre","weak","strong")
    simout_df$run = rep(r, nrow(simout_df))
    
    
    return(simout_df)
    
  }) %>% do.call(rbind,.)
  return(simoutput)
}

# Tertiary transition event --------------------------------------------------------


#time until third event occurs with prob mu3 per 
#cell with a 2nd driver dividing
simTimesThirdEvent = function(bp,
                              dp,
                              muwd,
                              musd,
                              swd,
                              ssd,
                              mu3){
  
  
  
  #derived parameters
  bpeff = bp *(1 - muwd - musd)
  lp = bpeff - dp
  
  bwd = lp*(1+swd)+dp
  lwd = bwd-dp
  muwdSurv = muwd*(1-dp/bwd) #only care about surviving mutations.
  
  bsd = lp*(1+ssd)+dp
  lsd = bsd-dp
  musdSurv = musd*(1-dp/bsd) #only care about surviving mutations.
  
  #sample random amplitudes
  Wp = rexp(1,1-dp/bp)
  Wwd = rexp(1,1-dp/bwd)
  Wsd = rexp(1,1-dp/bsd)
  
  #sample mutation times
  twd = sample_time_driverclone(bp,dp,Wp,lp, muwdSurv)
  tsd = sample_time_driverclone(bp,dp,Wp,lp, musdSurv)
  
  t3wd = sample_time_driverclone(bwd,dp,Wwd,lwd,mu3)+twd
  t3sd = sample_time_driverclone(bsd,dp,Wsd,lsd,mu3)+tsd
  
  t3rd = min(c(t3wd,t3sd))
  
  
  size_at_3rd <-  approx_sim_all(t3rd,
                                 lp,
                                 Wp,
                                 lwd,
                                 Wwd,
                                 twd,
                                 lsd,
                                 Wsd,
                                 tsd) %>% sum
  
  #return all times 
  return(c(t3wd,t3sd, size_at_3rd,twd,tsd))
}

#probability 3rd event occurs via 
#weak driver clone
prob3rdEventWeak <- function( bp,
                              dp,
                              muwd,
                              musd,
                              swd,
                              ssd,
                              mu3,
                              runs,
                              seednum){
  set.seed(seednum)
  simout <-  sapply(1:runs, function(x){
    simTimesThirdEvent( bp,
                        dp,
                        muwd,
                        musd,
                        swd,
                        ssd,
                        mu3) %>% return()
  })
  prob3rdviaWeak <- sum(simout[1,]<simout[2,])/runs
  return(prob3rdviaWeak)
}






# Plotting ----------------------------------------------------------------

if (!exists("clone_colors")){
  clone_colors = c("#619CFF",
                   "darkgoldenrod1",
                   "#F8766D" )
}

#### Trajectories ####


#clone sizes over tumour size , multiple trajectories +median.
get_trajectory_plot = function(bp,dp,
                               muwd,musd,
                               swd,ssd,
                               seednum=1,runs=50,
                               sizegrid=10^(seq(3,12,.2)),
                               alp =.05){
  
  # muwd = 0.005
  # musd = 0.00005
  # swd = 1.5
  # ssd = 5
  # seednum = 1
  # runs = 100
  
  
  simoutput <-  simApproxSizegridRuns(bp,
                                      dp,
                                      muwd,
                                      musd,
                                      swd,
                                      ssd,
                                      runs,
                                      sizegrid,
                                      seednum)
  
  simoutput_long = melt(simoutput,id = c("N","run"))
  
  plot_trajectories = ggplot( )+
    geom_line(data = subset(simoutput_long, run ==1),
              aes(x = N,
                  y = value/N,
                  color = variable),
              alpha = alp)+
    scale_x_continuous(trans = "log10",
                       limits =c(min(sizegrid),max(sizegrid)))+
    theme_classic()+
    scale_color_manual(values = clone_colors,
                       labels = c("Precursor", 
                                  "Mini driver",
                                  "Major driver"))+
    xlab("Tumour size (cell number)")+
    ylab("Clone proportion")+
    theme(legend.position = "none")
  # theme(legend.title=element_blank(),
  #       legend.text=element_text(size=10))
  # plot_trajectories1 = plot_trajectories
  # plot_trajectories1
  for (k in 2:runs){
    plot_trajectories = plot_trajectories +   geom_line(data = subset(simoutput_long, run ==k),
                                                        aes(x = N,
                                                            y = value/N,
                                                            color = variable),
                                                        alpha = alp)
  }
  
  
  medians_props_df = lapply(sizegrid,function(n){
    pre = (subset(simoutput_long,N == n & variable == "pre")$value/n) %>% median
    weak = (subset(simoutput_long,N == n & variable == "weak")$value/n) %>% median
    strong = (subset(simoutput_long,N == n & variable == "strong")$value/n) %>% median
    return(c(n,pre,weak,strong))
  }) %>% do.call(rbind,.) %>% data.frame()
  
  colnames(medians_props_df )  = c("N","pre","weak","strong")
  medians_props_df_long = melt(medians_props_df ,id = "N")
  
  
  plot_trajectories = plot_trajectories +  
    geom_line(data = medians_props_df_long,
              aes(x = N,
                  y = value,
                  color = variable),
              linetype = "dashed",
              size = 1)
  
  return(plot_trajectories)
  
}




#### Contour plot probability 3rd event via weak clone ####

plotProb3rdEventWeakGrid <-  function(bp,
                                      dp,
                                      muwd,
                                      mutBiasGrid,
                                      swd,
                                      ssd,
                                      mu3grid,
                                      runs ,
                                      seednum){
  
  #example params
  # bp = .2
  # dp = 0
  # muwd = 5*10^(-3)
  # # musd = .0001
  # swd = 1
  # ssd = 3
  # runs = 100
  # seednum = 1
  # mutBiasGrid = 10^seq(0,log10(500),by = .2)
  # mu3grid= 10^(-(seq(8,1,-.25)))
  # 
  # 
  
  pgrid = expand.grid(rmu = mutBiasGrid,
                      mu3 =  mu3grid)
  pb <- txtProgressBar(min = 0, max = nrow(pgrid), style = 3)
  
  prob3rdEventWeakVals = sapply(1:nrow(pgrid),function(k){
    setTxtProgressBar(pb, k)
    prob3rdEventWeak( bp,
                      dp,
                      muwd,
                      muwd/pgrid[k,"rmu"],
                      swd,
                      ssd,
                      pgrid[k,"mu3"],
                      runs,
                      seednum) %>% return()
  }) 
  
  
  prob3rdEventWeakParams <- cbind("prob3rdWeak" =   
                                    prob3rdEventWeakVals ,pgrid) %>% 
    as.data.frame
  
  mut3rdString <- expression("Tertiary driver rate"~(mu["3"]))
  
  mutBiasString <- expression(Mutational~bias~to~weak~driver~(beta))
  
  mutBiasNeeded3rd = sapply(1:nrow(pgrid), function(k){
    rs = ssd/swd
    u =    pgrid[k,"mu3"]
    return( u^(1/(1+rs*swd)-1/(1+swd)))
  }) 
  
  prob3rdEventWeakParams$mutBiasNeeded <-   mutBiasNeeded3rd
  prob3rdEventWeakParams$probThresh <- (prob3rdEventWeakParams$prob3rdWeak>0.5) %>% 
    as.numeric %>% factor
  plot_prob3rdweak <-  ggplot(data=prob3rdEventWeakParams)+
    geom_tile( aes( mu3,rmu,fill=prob3rdWeak),
               color = "white")+
    scale_x_continuous(trans = "log10",
                       limits = c(.5*min(prob3rdEventWeakParams$mu3),
                                  2*max(prob3rdEventWeakParams$mu3)))+
    scale_y_continuous(trans = "log10",
                       limits = c(min(prob3rdEventWeakParams$rmu),
                                  max(prob3rdEventWeakParams$rmu)))+
    scale_fill_gradient(low = "black", high = clone_colors[2],
                        name="Probability third driver\n via weak clone\n(simulations)",
                        breaks = c(0,.25,.5,.75,1),
                        labels = c(0,.25,.5,.75,1),
                        limits = c(0,1))+
    theme_bw()+
    ylab(mutBiasString)+
    xlab(mut3rdString) +
    geom_line( aes(x = mu3,
                   y=mutBiasNeeded),
               linetype = "dashed",
               color = "black",
               size = .5)
  
  
  return(plot_prob3rdweak)
  
}
