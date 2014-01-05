demog <- read.csv('DemogWithASLVals.csv')
require(reshape2)
require(ggplot2)
RawCBF    <- demog$RawCBFVals.Label5[-127]
FuncCBF   <- demog$FuncCBFVals.Label5[-127]
StructCBF <- demog$StructCBFVals.Label5[-127]
Thick     <- demog$Thickness_AAL.AAL5[-127]
Age       <- demog$AgeAtScan[-127]
Sex       <- demog$Sex[-127]
mydat.thick <- data.frame(RawCBF=RawCBF, ResCBF=FuncCBF, 
  AnatCBF=StructCBF, Thick=Thick)
mydat.thick.m <- melt(mydat, id="Thick")
ggplot(mydat.thick.m, aes(Thick, value, colour=variable)) + geom_point() + geom_smooth(method='lm') + xlab('Thickness (mm)') + 
  ylab('CBF (ml/100g/min') + ggtitle('Raw, Structurally Predicted, and Residual CBF vs. Cortical Thickness')
ggsave('../figure/ThickVsCBF.pdf', width=10, height=7)
mydat.age <- data.frame(RawCBF=RawCBF, ResCBF=FuncCBF, 
  AnatCBF=StructCBF, Age=Age)
mydat.age.m <- melt(mydat.age, id='Age')
ggplot(mydat.age.m, aes(Age, value, colour=variable)) + geom_point() + geom_smooth(method='lm') + xlab('Age (years)') + 
  ylab('CBF (ml/100g/min') + ggtitle('Raw, Structurally Predicted, and Residual CBF vs. Age')
ggsave('../figure/AgeVsCBF.pdf', width=10, height=7)
print(summary(lm(RawCBF ~ Thick)))
print(summary(lm(StructCBF ~ Thick)))
print(summary(lm(FuncCBF ~ Thick)))
