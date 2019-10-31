pep3 <- readRDS("data/pepthree.rds")
tree <- read.tree("data/Vascular_Plants_rooted.dated.tre")

rownames(pep3) <- pep3$nrow
tdforpic <- treedata(phy = tree, data = pep3)

vcmax_pic <- ape::pic(tdforpic$data[,2], tdforpic$phy)
jmax_pic <- ape::pic(tdforpic$data[,3], tdforpic$phy)
larea_pic <- ape::pic(tdforpic$data[,5], tdforpic$phy)
lmass_pic <- ape::pic(tdforpic$data[,6], tdforpic$phy)
nitro_pic <- ape::pic(tdforpic$data[,7], tdforpic$phy)
carb_pic <- ape::pic(tdforpic$data[,8], tdforpic$phy)
## remember these should be logged

## test for a relationship between
vjfit <- lm(log(abs(vcmax_pic)) ~ log(abs(jmax_pic)))
amfit <- lm(log(abs(larea_pic)) ~ log(abs(lmass_pic)))
vafit <- lm(log(abs(vcmax_pic)) ~ log(abs(larea_pic)))
vmfit <- lm(log(abs(vcmax_pic)) ~ log(abs(lmass_pic)))
ajfit <- lm(log(abs(larea_pic)) ~ log(abs(jmax_pic)))
mjfit <- lm(log(abs(lmass_pic)) ~ log(abs(jmax_pic)))

summary(vjfit)
summary(amfit)
summary(vafit)
summary(vmfit)
summary(ajfit)
summary(mjfit)

plot(x = vcmax_pic, y = jmax_pic)
plot(x = larea_pic, y = lmass_pic)
plot(x = larea_pic, y = vcmax_pic)
plot(x = lmass_pic, y = vcmax_pic)
plot(x = lmass_pic, y = jmax_pic)
plot(x = larea_pic, y = jmax_pic)
plot(x = nitro_pic, y = carb_pic)
plot(x = lmass_pic, y = nitro_pic)
