library(ERP)

dataInpeps = read.table("/Users/bolger/PycharmProjects/BrainIHM_python/AllData_Rformat.csv", header = FALSE, sep = ",", dec = ".")
time_pt<-seq(from=-250, to= 801, by=1.9531)

# Drop two very noisy participants (32 and 37). 
new_data1 <- droplevels(dataInpeps[dataInpeps$V2 != "S37", ])
new_data <- droplevels(new_data1[new_data1$V2 != "S32", ])

T <- length(time_pt)   # Number of time samples
n <- nrow(new_data)  # Number of erp curves

covariates <- new_data[, 2:5] 
erpdata <- as.numeric(unlist(new_data[, -(1:5)]))

sujets = factor(new_data$V2)
chans = factor(new_data$V3)
pgroup = factor(new_data$V4)
condition= factor(new_data$V5)

##==========================================================================================
# Only keep 15 out of the 64 electrodes.

dataInpeps_chans <- droplevels(new_data[new_data$V3 %in% c("F3", "Fz", "F4", "FC3", "FCz", "FC4", "C3", "Cz", "C4",
                                  "CP3", "CPz", "CP4", "P3", "Pz", "P4"), ])
dataInpeps_group <- droplevels(dataInpeps_chans[dataInpeps_chans$V4 %in% c("Human"),])


covariates_subset <- dataInpeps_chans[, 2:5]
erpdata_subset <- dataInpeps_chans[, -(1:5)] # Retain only the mV
Group = factor(dataInpeps_chans$V4)
Conds = factor(dataInpeps_chans$V5)
Channels = factor(dataInpeps_chans$V3)
Subjects = factor(dataInpeps_chans$V2)

with(covariates_subset, table(V3, V5, V4))

## ======================================================================================================
## Plot the difference curve for Agent - Human (Channel C4)
#  Set up a linear model to account for Group and Condition effects.
#  Note that the same participants viewed Agent and Human groups; this is accounted for in the model.

design <- model.matrix(~ C(Subjects, sum)+Group+Conds+Channels + Channels:Conds + Channels:Group + Conds:Group + Channels:Conds:Group,
                       data=covariates_subset)

erpplot(erpdata_subset, design, effect = 38, interval = "simultaneous", 
        nbs = 20, lwd = 2, frames = time_pt, 
        xlab = "Time (ms)", ylab = "Group effect", ylim=rev(c(-3*(10^-6),2*(10^-6))))
title("Human-Agent difference curve")

##=======================================================================================================
# Plot the ERPs for each Agent vs. Human for each condition over the 15 channels of interest.

chanlevels <- c("F3", "Fz", "F4", "FC3", "FCz", "FC4", "C3", "Cz", "C4",
                "CP3", "CPz", "CP4", "P3", "Pz", "P4")
chanmat <- matrix(chanlevels, nrow=5, byrow=TRUE)
conds <- levels(covariates_subset$V5)
sgroups <- levels(covariates_subset$V4)
sujets <- levels(covariates_subset$V2)

#adjust plot margins
par(mar = c(1, 1, 1, 1))
par(mfrow = c(5,3))

for (i in 1:5) {
  for (j in 1:3) {
    
    chansel <- (covariates_subset$V3 == chanmat[i,j]) &
      (covariates_subset$V5 == "InCongru-InCongru") 
    
    erpplot(erpdata_subset[chansel, ],frames = time_pt, col = c("blue", "green") ,
            xlab = "Time (ms)", ylab = "ERP (mV)", ylim=rev(c(-10*(10^-6),10*(10^-6))),
            axes="FALSE",
            main = paste("Channel: ", chanmat[i,j], "- Condition: Congru-Congru", sep = ""))
    axis(1, pos=0); axis(2, pos=0)
    grid()
    
    
    
    legend("topright", bty = "n", lwd = 3, col = c("blue", "green"),
           legend = c("Agent", "Human"))
  }
}

##=========================================================================================
# Plot the difference curve for Agent vs. Human
# The effect is calculated for each condition separately and presents Human-Agent.

DataNew <- droplevels(dataInpeps_chans[dataInpeps_chans$V5 %in% c("InCongru-InCongru"), ])
covariates_new <- DataNew[, 2:5]
erpdata_new <- DataNew[, -(1:5)]*(10^6)
ChannelsL =  c("F3", "Fz", "F4", "FC3", "FCz", "FC4", "C3", "Cz", "C4",
               "CP3", "CPz", "CP4", "P3", "Pz", "P4")


par(mar = c(1, 1, 1, 1))
par(mfrow = c(5,3))
count = 1


for (i1 in 1:5) {
  for (j1 in 1:3) {
    
    covariates_new$V3 <- relevel(covariates_new$V3, ref = ChannelsL[count])
    
    designN <- model.matrix(~ C(V2, sum) + V4 + V3 + V3:V4,
                            data =covariates_new)
    effect = which(colnames(designN)=='V4Human')
    
    erpplot(erpdata_new, designN, effect=effect, interval='simultaneous', alpha = 0.05, nbs=NULL, lwd = 2,
            frames = time_pt, ylim=rev(c(-6,6)), xlab = 'Time (ms)',
            ylab = expression(mu~V), bty="l", axes="FALSE")
    axis(1, pos=0); axis(2, pos=0)
    grid()
    title(paste("Human - Agent effect curve: Channel ", ChannelsL[count], sep= ""), sub="Condition 'Cong-Cong' - 37 subjects", cex.main=1.25, font.main=4,
         cex.sub=0.75, font.sub=3)
    count = count+1
    
  }
}

