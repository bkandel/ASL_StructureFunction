demog <- read.csv('DemogWithASLVals.csv')
structcols <- grepl('Struct', names(demog))
thickcols <- grepl('Thickness_AAL.AAL', names(demog))
funccols <- grepl('Func', names(demog))
rawcols <- grepl('RawCBF', names(demog))
badsubj <- 127

subjcors <- NULL
for( mysubj in 1:dim(demog)[1]){
  cor.struct <- cor.test(as.numeric(demog[mysubj, thickcols]), 
           as.numeric(demog[mysubj, structcols][, 2:91]))$estimate
  #plot(as.numeric(demog[mysubj, thickcols]), 
  #     as.numeric(demog[mysubj, structcols][, 2:91]))
  cor.func <- cor.test(as.numeric(demog[mysubj, thickcols]), 
           as.numeric(demog[mysubj, funccols][, 2:91]))$estimate
  cor.raw <- cor.test(as.numeric(demog[mysubj, thickcols]), 
                      as.numeric(demog[mysubj, rawcols][, 2:91]))$estimate
  #plot(as.numeric(demog[mysubj, thickcols]), 
  #     as.numeric(demog[mysubj, funccols][, 2:91]))
  subjcors <- rbind(subjcors, data.frame(StructCor=cor.struct, FuncCor=cor.func, 
                                         RawCor=cor.raw))
}

roicors.thick <- NULL
for(myroi in 2:91){
  cor.struct <- cor.test(as.numeric(demog[, thickcols][, (myroi-1)][-badsubj]), 
                         as.numeric(demog[, structcols][, myroi][-badsubj]))$estimate
  cor.raw <- cor.test(as.numeric(demog[, thickcols][, myroi-1][-badsubj]), 
                      as.numeric(demog[, rawcols][, myroi][-badsubj]))$estimate
  cor.func <- cor.test(as.numeric(demog[, thickcols][, myroi-1][-badsubj]), 
                       as.numeric(demog[, funccols][, myroi][-badsubj]))$estimate
  mylab <- strsplit(colnames(demog[, structcols])[myroi], '.', T)[[1]][2]
  roicors.thick <- rbind(roicors.thick, data.frame(StructCor=cor.struct, FuncCor=cor.func,
                                       RawCor=cor.raw, row.names=mylab))
}
roicors.age <- NULL
for(myroi in 2:91){
  cor.struct <- cor.test(as.numeric(demog$AgeAtScan[-badsubj]), 
                         as.numeric(demog[, structcols][, myroi][-badsubj]))$estimate
  cor.raw <- cor.test(as.numeric(demog$AgeAtScan[-badsubj]), 
                      as.numeric(demog[, rawcols][, myroi][-badsubj]))$estimate
  cor.func <- cor.test(as.numeric(demog[, thickcols][, myroi-1][-badsubj]), 
                       as.numeric(demog[, funccols][, myroi][-badsubj]))$estimate
  cor.thick <- cor.test(as.numeric(demog[, thickcols][, (myroi-1)][-badsubj]), 
             as.numeric(demog$AgeAtScan[-badsubj]))$estimate
  mylab <- strsplit(colnames(demog[, structcols])[myroi], '.', T)[[1]][2]
  roicors.age <- rbind(roicors.age, data.frame(StructCor=cor.struct, FuncCor=cor.func,
                    RawCor=cor.raw, ThickCor=cor.thick, row.names=mylab))
}
