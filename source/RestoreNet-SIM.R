cat("\nInstall/load packages")
inst.pkgs <- installed.packages()

## required packages:
l.pkgs <- c("devtools",
            "xtable",
            "scatterpie",
            "RColorBrewer",
            "scales",
            "Matrix",
            "stringr",
            "ggplot2",
            "egg")
## check if packages are installed
lapply(l.pkgs, function(pkg){
  if(!(pkg %in% rownames(inst.pkgs))){
    install.packages(pkg)
  }
})

lapply(l.pkgs, function(pkg){library(pkg,character.only=TRUE)})

install_github("delcore-luca/RestoreNet",
               ref = "master")
library(RestoreNet)

rm(list = ls())

# s20 <- 1
s20 <- as.numeric(commandArgs(trailingOnly = TRUE)) ## true noise variance

##############
## currTime ##
##############

currTime <- Sys.time()
currTime <- gsub(" ", "_", currTime)
currTime <- gsub("-", "", currTime)
currTime <- gsub(":", "", currTime)
currTime <- substr(currTime, start = 1, stop = 13)

#############
## folders ##
#############

currDir <- getwd()
setwd(currDir)

#############
## folders ##
#############

resPath <- "./results/"
ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
resPath <- paste(resPath, currTime, sep = "/")
figPath <- paste(resPath, "figures/", sep = "/")

ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
ifelse(!dir.exists(file.path(figPath)), dir.create(file.path(figPath)), FALSE)


#################
## tau-leaping ##
#################
##

rcts <- c("A->1", "B->1", "C->1", "D->1",
          "A->0", "B->0", "C->0", "D->0",
          "A->B", "A->C", "C->D")
ctps <- head(LETTERS, 4)
theta <- c(.2,.15,.17,.09,
           .001, .007, .004, .002,
           .13, .15, .08) ## set theta
names(theta) <- rcts
nC <- 3 ## set number of clones
S <- 100 ## set the number of simulations
tau <- 1

u_1 <- c(.2,.15,.17,.09*5,
         .001, .007, .004, .002,
         .13, .15, .08)
u_2 <- c(.2,.15,.17,.09,
         .001, .007, .004, .002,
         .13, .15, .08)
u_3 <- c(.2,.15,.17*3,.09,
         .001, .007, .004, .002,
         .13, .15, .08)
theta_allcls <- cbind(u_1, u_2, u_3)
rownames(theta_allcls) <- rcts

Y <- array(data = NA,
           dim = c(S + 1, length(ctps), nC),
           dimnames = list(seq(from = 0, to = S*tau, by = tau),
                           ctps,
                           1:nC)) ## empty array to store simulations
Y0 <- c(100,0,0,0) ## initial state
names(Y0) <- ctps

#########################
## matrix for IC results
#########################
nSim <- 100
nModels <- 2
resIC <- array(data = NA,
               dim = c(nModels, 11, nSim),
               dimnames = list(c("null", "re"),
                               c("nPar", "cll", "mll", "cAIC", "mAIC", "Chi2", "p-value", "KLdiv", "KLdiv/N", "BhattDist_nullCond", "BhattDist_nullCond/N"),
                               1:nSim))
res.base.par <- matrix(data = NA, nrow = length(rcts) + 1, ncol = nSim)
rownames(res.base.par) <- c(rcts, "s2")
res.re.par <- matrix(data = NA, nrow = length(rcts)*2 + 1, ncol = nSim)
rownames(res.re.par) <- c(rcts, paste(rcts, "d2", sep = "_"), "s2")
res.re.Euy <- matrix(data = NA, nrow = length(rcts)*nC, ncol = nSim)
rownames(res.re.Euy) <- c(sapply(1:nC, function(cl){paste(rcts, cl, sep = "_")}))


for (sim in 1:nSim) {
  Y[] <- NA
  cat("Simulation n. ", sim, "...", sep = "")
  cat("tau-leaping simulation...", sep = "")
  for (cl in 1:nC) {
    Y[,,cl] <- get.sim.tl(Yt = Y0,
                          theta = theta_allcls[,cl],
                          S = S,
                          s2 = s20,
                          tau = tau,
                          rct.lst = rcts)
  }
  cat("DONE\n", sep = "")

  Y <- Y[which(apply(Y != 0, 1, sum) > 0),,]
  dimnames(Y)[[3]] <- 1:dim(Y)[3]
  Y_NA <- Y
  Y_NA[Y_NA == 0] <- NA


  pdf(file = paste(figPath, "/Ysim", sim,".pdf", sep = ''), height=10, width=5)
  par(mar = c(5,6,2,2), mfrow = c(3,1))
  matplot(Y[,,1], ylim = range(as.numeric(unlist(Y))),
          lty = 1, lwd = 3, type = 'l',
          cex.axis = 2, cex.lab = 2, main = "clone 1", cex.main = 2,
          xlab = "t", ylab = expression(Y[t])) ## plot simulated trajectories
  matplot(Y[,,2], ylim = range(as.numeric(unlist(Y))),
          lty = 1, lwd = 3, type = 'l',
          cex.axis = 2, cex.lab = 2,  main = "clone 2", cex.main = 2,
          xlab = "t", ylab = expression(Y[t])) ## plot simulated trajectories
  matplot(Y[,,3], ylim = range(as.numeric(unlist(Y))),
          lty = 1, lwd = 3, type = 'l',
          cex.axis = 2, cex.lab = 2,  main = "clone 3", cex.main = 2,
          xlab = "t", ylab = expression(Y[t])) ## plot simulated trajectories
  dev.off()

  pdf(file = paste(figPath, "/Ysim_all", sim,".pdf", sep = ''), height=5, width=7)
  par(mar = c(5,6,2,2), mfrow = c(1,1))
  matplot(Y[,,1], ylim = range(as.numeric(unlist(Y))),
          lty = 1, lwd = 3, type = 'l',
          cex.axis = 2, cex.lab = 2,
          xlab = "t", ylab = expression(Y[t])) ## plot simulated trajectories
  matplot(Y[,,2], ylim = range(as.numeric(unlist(Y))),
          lty = 2, lwd = 5, type = 'l',
          cex.axis = 2, cex.lab = 2,
          xlab = "t", ylab = expression(Y[t]), add = T) ## plot simulated trajectories
  matplot(Y[,,3], ylim = range(as.numeric(unlist(Y))),
          lty = 3, lwd = 5, type = 'l',
          cex.axis = 2, cex.lab = 2,
          xlab = "t", ylab = expression(Y[t]), add = T) ## plot simulated trajectories
  text(x = 30,
       y = 210,
       labels = "cell type", cex = 1.5)
  legend(x = 45,
         y = 230,
         legend = ctps,
         col = 1:4, pch = 19, cex = 1.5, horiz = T, bty = "n")
  text(x = 30,
       y = 190,
       labels = "clone", cex = 1.5)
  legend(x = 40,
         y = 210,
         legend = 1:3,
         col = 1, lwd = c(3,3,3), lty = 1:3, cex = 1.5, horiz = T, bty = "n")
  dev.off()


  ######################
  ## null model fitting:
  ######################
  null.res <- fit.null(Y = Y, rct.lst = rcts)

  resIC["null", names(null.res$stats), sim] <- null.res$stats
  res.base.par[,sim] <- null.res$fit$par

  #################################
  ## random-effects model fitting:
  #################################

  re.res <- fit.re(theta_0 = null.res$fit$par,
                   Y = Y,
                   rct.lst = rcts,
                   maxemit = 100)

  resIC["re", , sim] <- re.res$stats
  res.re.par[,sim] <- re.res$fit$par
  res.re.Euy[,sim] <- as.numeric(re.res$fit$VEuy$euy)

  pdf(file = paste(figPath, "/thetaSim", sim,".pdf", sep = ''), height=6, width=6)
  par(mar = c(5,6,2,2))
  plot(log(c(u_1, u_2, u_3)), log(re.res$fit$VEuy$euy),
       # pch = 20,
       pch = c(rep(1,4), rep(0,4), rep(2,3)),
       col = alpha(rep(1:nC, each = length(u_1)), 1),
       lwd = 2,
       cex.axis = 2, cex.lab = 2, cex = 2,
       xlab = expression(log~theta[true]), ylab = expression(log~E[hat(theta)](u*"|"*y)))
  lines(x = log(range(c(u_1, u_2, u_3), as.numeric(re.res$fit$VEuy$euy))),
        y = log(range(c(u_1, u_2, u_3), as.numeric(re.res$fit$VEuy$euy))), type = 'l', lty = 1, lwd = 2, col = "red")
  legend(x = log(min(c(u_1, u_2, u_3, as.numeric(re.res$fit$VEuy$euy))))*1,
         y = log(max(c(u_1, u_2, u_3, as.numeric(re.res$fit$VEuy$euy))))*.8,
         legend = c(expression(alpha),
                    expression(delta),
                    expression(lambda)),
         pch = c(1, 0, 2),
         cex = 2,
         horiz = TRUE,
         bty = "n")
  text(x = log(min(c(u_1, u_2, u_3, as.numeric(re.res$fit$VEuy$euy))))*.9,
       y = log(max(c(u_1, u_2, u_3, as.numeric(re.res$fit$VEuy$euy))))*1.2,
       cex = 1.5,
       labels = "rates")
  legend(x = log(max(c(u_1, u_2, u_3, as.numeric(re.res$fit$VEuy$euy))))*6.5,
         y = log(min(c(u_1, u_2, u_3, as.numeric(re.res$fit$VEuy$euy))))*.9,
         legend = c("1",
                    "2",
                    "3"),
         pch = 20,
         col = 1:3,
         cex = 2,
         horiz = TRUE,
         bty = "n")
  text(x = log(max(c(u_1, u_2, u_3, as.numeric(re.res$fit$VEuy$euy))))*2,
       y = log(min(c(u_1, u_2, u_3, as.numeric(re.res$fit$VEuy$euy))))*.93,
       cex = 1.5,
       labels = "clones")
  dev.off()

  cbind(re.res$fit$VEuy$euy, c(u_1, u_2, u_3))
  relErr(re.res$fit$VEuy$euy, c(u_1, u_2, u_3))

  cat("DONE\n")

  ########### BOXPLOTS ###########
  pdf(file = paste(figPath, "s2_", s20, "sim_", sim, "_boxplots.pdf", sep = ''), height=7, width=7)
  get.boxplots(re.res)
  dev.off()

  ########### SCATTERPIE PLOT ###########
  pdf(file = paste(figPath, "s2_", s20, "sim_", sim, "_scatterpie.pdf", sep = ''), height=7, width=7)
  get.scatterpie(re.res, txt = TRUE)
  dev.off()

  ###### rank the random effects with their variance:
  euy_allIS <- matrix(data = re.res$fit$VEuy[["euy"]],
                      nrow = ncol(re.res$design$M_bdiag)/ncol(re.res$design$M),
                      byrow = TRUE)
  rownames(euy_allIS) <- unique(rownames(re.res$design$M))
  colnames(euy_allIS) <- colnames(re.res$design$V)
  print(xtable(t(euy_allIS), digits = 3), include.rownames=FALSE, file = paste(resPath, "/s2_", s20, "sim_", sim, "Euy_ForLatex.txt", sep = ""))

}

save.image(paste(resPath, "/simulations_s2", s20,".RData", sep = ""))
gc()

## ADDITIONAL PLOTS

pdf(file = paste(figPath, "/AICsBoxplot_s2", s20,".pdf", sep = ''), height=5, width=5)
boxplot(t(resIC[,"mAIC",]), cex.axis = 2, cex.lab = 2, pch = 20, cex = 2, lwd = 2)
dev.off()

pdf(file = paste(figPath, "/ParametersBoxplot_s2", s20,".pdf", sep = ''), height=5, width=14)
boxplot(t(res.re.Euy/c(u_1, u_2, u_3)), xaxt = 'n', cex.axis = 2, lwd = 2)
abline(h = 1, col = "red", lwd = 2)
# axis(1, at=seq(1, nrow(res.re.Euy), by=1), labels = FALSE)
# text(seq(1, nrow(res.re.Euy), by=1) -.8,
#      par("usr")[3] - 1.5,
#      labels = rownames(res.re.Euy),
#      font = 2, cex = 1.5, srt = 45, pos = 1, xpd = TRUE)
dev.off()
