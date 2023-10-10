library(ERP)

dataInpeps = read.table("/Users/bolger/PycharmProjects/BrainIHM_python/AllData_Rformat.csv", header = FALSE, sep = ",", dec = ".")
time_pt<-seq(from=-250, to= 801, by=1.9531)
dim(dataInpeps)

# Condition effect for Human vs. Agent ------------------------------------

dataInpeps1 <- dataInpeps[dataInpeps$V5 != "InCongru-Congru",]  # Exclude factor level corresponding to Incongru-Congru
dataInpeps2 <- dataInpeps1[dataInpeps1$V5 != "InCongru-InCongru",] # Exclude factor level corresponding to Incongru-Incongru

# Define the erp-data and the covariates for this effect of intere --------

covariates.CongAllchans <- dataInpeps2[, 2:4] 
erpdata.CongAllchans <- as.numeric(unlist(dataInpeps2[, -(1:5)])) 
erpdata.CongAllchans <- dataInpeps2[, -(1:5)]
pgroup = factor(dataInpeps2$V4)
condcong = factor(dataInpeps2$V5)
chans = factor(dataInpeps2$V3)
sujets = factor(dataInpeps2$V2)


# Set up design matrix ----------------------------------------------------

#design <- model.matrix(~ C(sujets, sum) + relevel(pgroup, ref="Agent") + condcong   + chans + chans:pgroup +
#                        chans:condcong + relevel(pgroup, ref="Agent"):condcong,data = covariates.CongAllchans)   # Non-Null design
#design0 <-model.matrix(~ C(sujets, sum) + relevel(pgroup, ref = "Agent") + condcong + chans + chans:pgroup + 
#                        chans:condcong, data = covariates.CongAllchans)   # NULL design

design <- model.matrix(~ C(sujets, sum) + pgroup + chans + pgroup:chans, data = covariates.CongAllchans)
design0 <- model.matrix(~ C(sujets, sum)+pgroup+chans, data = covariates.CongAllchans)

# Carry out ftest using current design and design0 ------------------------

F <- erpFtest(erpdata.CongAllchans, design, design0, nbf = NULL, pvalue = "none", 
              wantplot = TRUE)   # p-value not calculated and number of factors estimated. 

# Monte-Carlo p-value calculated with a Satterthwaite approximation. Number of monte-carlo samples used = 200. 
F2 <-erpFtest(erpdata.CongAllchans, design, design0, nbf = F$nbf, pvalue = "Satterthwaite") 

# Average the curves of 20 samples ----------------------------------------

avetest <- erpavetest(erpdata.CongAllchans, design, design0, nintervals = 20, method = "bonferroni", alpha=.01)

erpplot(erpdata.CongAllchans, design, effect = ncol(design), lwd = 2, interval = "simultaneous",
        frames = time_pt, y = rev(c(-4, 8)), xlab = "Time", ylab = "Condition Effect")
title(main = " Group (Agent - Human) - Condition (congruent-Congruent) \n pvalue = .01 (All channels)", font.main = 2)
abline(v=250, col="red")
points(time_pt[avetest$significant], rep(0, length(avetest$significant)), pch = 20, col = "goldenrod")
abline(v = time_pt[avetest$breaks], lty = 2, col = "darkgray")

time_pt = seq(0,800,2)     # sequence of time points (1 time point every 2ms in [0,1000])
nbs = 20                    # Number of B-splines for the plot of the effect curve
effect=which(colnames(design)=="Agent")
erpplot(erpdata.CongAllchans,design=design,frames=time_pt,effect=ncol(design),xlab="Time (ms)",
        ylab=expression(Effect~curve~(mu~V)),bty="l",ylim=rev(c(-3*(10^-6),3*(10^-6))),nbs=nbs,
        cex.axis=1.25,cex.lab=1.25,interval="simultaneous")
# with interval="simultaneous", both the pointwise and the simultaneous confidence bands
# are plotted
abline(v=time_pt[avetest$breaks],lty=2,col="darkgray")
# Add a grid to show breakpoints
points(time_pt[avetest$significant],rep(0,length(avetest$significant)),pch=16,col="blue")


# Apply the adaptive Factor-Adjustment Method -----------------------------

fabh1 <- erpfatest(erpdata.CongAllchans, design, design0, nbf = NULL, wantplot = TRUE)
fabh2 <- erpfatest(erpdata.CongAllchans, design, design0, nbf = fabh1$nbf)


