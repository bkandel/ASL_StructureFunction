demog <- read.csv('DemogWithASLVals.csv')
require(reshape2)
require(ggplot2)
raw <- demog$FuncCBFVals.Label87[-127] + demog$StructCBFVals.Label87[-127]
func <- demog$FuncCBFVals.Label87[-127]
struct <- demog$StructCBFVals.Label87[-127]
age <- demog$AgeAtScan[-127]
mydat <- data.frame(raw=raw, age=age, func=func, struct=struct)
mydat.m <- melt(mydat, id="age")
ggplot(mydat.m, aes(age, value, colour=variable)) + geom_point() + geom_smooth(method='lm')

print(summary(lm(raw ~ age)))
print(summary(lm(struct ~ age)))
print(summary(lm(func ~ age)))