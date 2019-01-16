# Set plotting parameters
library(grDevices)
library(ggthemes)
library(ggalt)

#change color of CI fill here
fillcolor <- "#230FCB"
transparency <- .6
coverageLevels <- seq(.1,.9, by = .1)

# figure 3 (incidence)
# this feeds in from "singleCoverageRisk"

fig3Data <- as.data.frame(rbind(fig3_10,fig3_20,fig3_30,fig3_40,fig3_50,fig3_60,fig3_70,fig3_80, fig3_90))
fig3Data <- cbind(coverageLevels, fig3Data)
colnames(fig3Data)<-c("coverage","estimate","lowerCI","upperCI")

fig4Data <- as.data.frame(rbind(fig4_10,fig4_20,fig4_30,fig4_40,fig4_50,fig4_60,fig4_70,fig4_80, fig4_90))
fig4Data <- cbind(coverageLevels, fig4Data)
colnames(fig4Data)<-c("level","estimate","lowerCI","upperCI")

#ggplot() + geom_ribbon(aes(ymin = ribspllower$y, ymax = ribsplupper$y, x = ribspllower$x), fill = fillcolor, alpha = transparency) + geom_line(aes(y=fig3Data$estimate, x=coverageLevels))
#plot figures 3 and 4
spl <- as.data.frame(spline(coverageLevels, fig3Data$estimate))
ribspllower <- as.data.frame(spline(coverageLevels,fig3Data$lowerCI))
ribsplupper <- as.data.frame(spline(coverageLevels,fig3Data$upperCI))
fig3plot <- ggplot(fig3Data, aes(x = coverageLevels, y = estimate)) + 
  geom_ribbon(aes(ymin=lowerCI, ymax=upperCI, x=coverageLevels),alpha = transparency, fill = fillcolor)+
  geom_line(size = 1) +
  xlab(expression(paste(alpha,"'",sep = "")))+
  ylab(expression("Estimated Disseminated Effect")) +
  theme_classic()

ggsave("Figure_3.pdf",fig3plot)

fig4plot <- ggplot(fig4Data, aes(x = coverageLevels, y= estimate)) + 
  geom_ribbon(aes(ymin = lowerCI,ymax = upperCI, x = coverageLevels), fill = fillcolor, alpha = transparency) +
  geom_line(size = 1) +
  xlab(expression(paste(alpha,"'",sep = "")))+ scale_x_continuous(breaks = seq(.1:.9, by = .1)) +
  ylab(expression("Estimated Individual Effect")) +
  theme_classic() 
ggsave("Figure_4.pdf",fig4plot)


testsmooth <- ggplot() + geom_ribbon(aes(ymin = ribspllower$y, ymax = ribsplupper$y, x = ribspllower$x)) +
  geom_line(aes(y=fig3Data$estimate, x=fig3Data$coverageLevels))+
  xlab(expression(paste(alpha,"'",sep = "")))+
  ylab(expression(widehat(DE)(alpha)))
  
