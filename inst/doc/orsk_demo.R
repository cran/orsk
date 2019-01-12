### R code from vignette source 'orsk_demo.Rnw'

###################################################
### code chunk number 1: orsk_demo.Rnw:154-155
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: orsk_demo.Rnw:157-158 (eval = FALSE)
###################################################
## install.packages("orsk")


###################################################
### code chunk number 3: orsk_demo.Rnw:161-163 (eval = FALSE)
###################################################
## library("orsk")
## vignette("orsk_demo",package = "orsk")


###################################################
### code chunk number 4: orsk_demo.Rnw:166-167 (eval = FALSE)
###################################################
## edit(vignette("orsk_demo",package = "orsk"))


###################################################
### code chunk number 5: orsk_demo.Rnw:177-178
###################################################
library("orsk")


###################################################
### code chunk number 6: orsk_demo.Rnw:180-183
###################################################
### analysis of Table 2 with grid method
res1 <- orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au= 3.03, method="grid")
summary(res1)


###################################################
### code chunk number 7: orsk_demo.Rnw:190-191
###################################################
plot(res1, type="OR")


###################################################
### code chunk number 8: orsk_demo.Rnw:199-200
###################################################
plot(res1, type="RR")


###################################################
### code chunk number 9: orsk_demo.Rnw:219-223
###################################################
### analysis of Table 2 with optim method
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
seed=579))


###################################################
### code chunk number 10: orsk_demo.Rnw:225-226 (eval = FALSE)
###################################################
## res2 <- orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au= 3.03, method="optim")


###################################################
### code chunk number 11: orsk_demo.Rnw:228-229
###################################################
res2 <- orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au= 3.03, method="optim")


###################################################
### code chunk number 12: orsk_demo.Rnw:231-233
###################################################
summary(res2)
summary(res2$res$RR)


###################################################
### code chunk number 13: orsk_demo.Rnw:235-238
###################################################
#compare the computing speed between the two methods of estimation. Time to finish the modeling
#time.optim <- system.time(orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au=3.03, method="optim"))[1]
time.grid <- system.time(orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au=3.03, method="grid"))[1]


