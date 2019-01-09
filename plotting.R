# Set plotting parameters
library(grDevices)

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


#plot figures 3 and 4
fig3plot <- ggplot(fig3Data) + geom_ribbon(aes(ymin=lowerCI, ymax=upperCI, x=coverageLevels),alpha = transparency, fill = fillcolor)+
  geom_line(aes(y=estimate, x=coverageLevels))+
  xlab(expression(paste(alpha,"'",sep = "")))+
  ylab(expression(widehat(DE)(alpha)))

ggsave("Figure_3.pdf",fig3plot)

fig4plot <- ggplot(fig4Data, aes(x = coverageLevels, y= estimate)) + 
  geom_ribbon(aes(ymin = lowerCI,ymax = upperCI, x = coverageLevels), fill = fillcolor, alpha = transparency) +
  geom_line() +
  xlab(expression(paste(alpha,"'",sep = "")))+
  ylab(expression(widehat(IE)(alpha)))
ggsave("Figure_4.eps",fig4plot)
