library(ERP)

PEPS_data = read.table("/Users/bolger/PycharmProjects/BrainIHM_python/AllData_Rformat.csv", header = FALSE, sep = ",", dec = ".")
time_pt<-seq(from=-250, to= 801, by=1.9531)
dim(PEPS_data)

# Drop two very noisy participants (32 and 37). 
new_data2 <- droplevels(PEPS_data[PEPS_data$V2 != "S37", ])
new_data1 <- droplevels(new_data2[new_data2$V2 != "S32", ])

# Retain only the channels of interest and the conditions of interest. 
PEPS_chansubset <- droplevels(new_data1[new_data1$V3 %in% c("Fz", "FCz", "Cz", "CPz","Pz"), ])
PEPS_congvsincong <- droplevels(PEPS_chansubset[PEPS_chansubset$V5 %in% c("Congru-Congru", "InCongru-Congru"), ])
PEPS_group <- droplevels(PEPS_congvsincong[PEPS_congvsincong$V4 %in% c("Agent"), ])

covariates <- PEPS_congvsincong[,1:5]
mvdata <- PEPS_congvsincong[, -(1:5)]*(10^6)

covariates_group <- PEPS_group[,1:5]
mvdata_group <- PEPS_group[,-(1:5)]*(10^6)

# Define the channels of interest.
ChannelsL =  c("Fz", "FCz", "Cz", "CPz","Pz")
count = 1

par(mfrow = c(5,1)) 

for (i in 1:5) {
    
    covariates$V3 <- relevel(covariates$V3, ref = ChannelsL[count])
    design = model.matrix(~C(V2, sum)+V3+V5+V3:V5, data = covariates_group)
    design0 = model.matrix(~C(V2, sum)+V3, data = covariates_group)
    colnames(design)
    
    #avetest <- erpavetest(mvdata, design, design0, alpha=.001, method="fdr", nintervals = 100)
    gb <- gbtest(mvdata_group, design, design0, graphthresh = .01)
    
    effect = which(colnames(design)=="V5InCongru-Congru")
    erpplot(mvdata_group, design, effect=effect, interval='simultaneous', alpha = 0.01, lwd = 2,
            frames = time_pt, ylim=rev(c(-5,5)), xlab = 'Time (ms)',
            ylab = expression(mu~V), bty="l", axes="FALSE")
    axis(1, pos=0); axis(2, pos=0)
    grid()
    points(time_pt[gb$significant], rep(0, length(gb$significant)), pch=20, col="goldenrod")
    #abline(v=time_pt[avetest$breaks], col="darkgray")
    title(paste("Agent: Congru-Congru vs. InCongru-Congru effect curve.\n Channel: ", ChannelsL[count], sep= ""), cex.main=1.25, font.main=4)
    count = count+1
}
