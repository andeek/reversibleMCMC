#### Organizing Soccer Data ####
s0506 <- read.csv("data/soccer0506.csv")
s0607 <- read.csv("data/soccer0607.csv")
s0708 <- read.csv("data/soccer0708.csv")

fullsoccer <- do.call(rbind, list(s0506, s0607, s0708))

fullsoccer$TotalGoals <- fullsoccer$FTHG + fullsoccer$FTAG

soccer.sub <- fullsoccer[,c(2,3,8)]

write.csv(soccer.sub, file="data/soccer.csv", row.names=F)
