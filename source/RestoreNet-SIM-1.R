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
            "egg",
            "parallel",
            "RestoreNet")
## check if packages are installed
lapply(l.pkgs, function(pkg){
  if(!(pkg %in% rownames(inst.pkgs))){
    install.packages(pkg)
  }
})

lapply(l.pkgs, function(pkg){library(pkg,character.only=TRUE)})

rm(list = ls())

s20_all <- c(.1,1,10)

# set.seed(123)

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

currDir_time <- paste(currDir, currTime, sep = "/")
ifelse(!dir.exists(file.path(currDir_time)), dir.create(file.path(currDir_time)), FALSE)

#########################
## parallel processors ##
#########################

nCores <- 3
cl <- makeCluster(nCores, type = "PSOCK")
clusterExport(cl = cl,
              varlist = ls())
clusterEvalQ(cl, library("Matrix"))
clusterEvalQ(cl, library("RestoreNet"))
clusterEvalQ(cl, sink(paste0(currDir_time, "/", Sys.getpid(), ".txt")))

res_alls20 <- parLapply(cl = cl,
                        X = s20_all,
                        fun = function(s20){

                          set.seed(123)
                          fit.null <- function (Y, rct.lst, maxit = 10000, factr = 1e+07, pgtol = 1e-08,
                                                lmm = 100, trace = TRUE, verbose = TRUE)
                          {
                            check.input(Y = Y, rct.lst = rct.lst)
                            compile.h(rct.lst = rct.lst, envir = environment())
                            V <- get.V(rct.lst = rct.lst)
                            Y <- Y[, rownames(V), , drop = FALSE]
                            Y_rbind <- Reduce(rbind, lapply(1:dim(Y)[3], function(cl) {
                              return(head(Y[, , cl], -1))
                            }))
                            dT <- rep(rep(diff(as.numeric(rownames(Y))), each = ncol(Y)),
                                      times = dim(Y)[3])
                            if (verbose) {
                              cat("creating design matrix...")
                            }
                            M <- eval(parse(text = "get.M(y = Y_rbind, V = V, dT = dT, get.h = get.h)"),
                                      envir = environment())
                            dx <- unlist(lapply(1:dim(Y)[3], function(cl) {
                              return(get.dx(Y[, , cl]))
                            }))
                            VCNs <- rep(1, length(dx))
                            clones <- as.vector(unlist(dimnames(Y)[3]))
                            y <- dx
                            y_all <- y
                            M_all <- M
                            clones_all <- clones <- rep(clones, each = (nrow(Y) - 1) *
                                                          ncol(Y))
                            VCNs_all <- VCNs
                            yfull <- y
                            Vfull <- V
                            yM <- cbind(y, M)
                            yM <- yM[rownames(M) == "T" | names(y) == "T", ]
                            clones <- clones[rownames(M) == "T" | names(y) == "T"]
                            VCNs <- VCNs[rownames(M) == "T" | names(y) == "T"]
                            dT <- dT[rownames(M) == "T" | names(y) == "T"]
                            y <- dx <- yM[, 1]
                            Mfull <- M <- yM[, -1]
                            rownames(M) <- names(y) <- names(VCNs) <- rownames(yM) <- clones
                            nrow(M) == length(dx) & length(dx) == length(VCNs) & length(VCNs) ==
                              length(dT)
                            if (verbose) {
                              cat(" DONE\n")
                            }
                            nObs <- apply(apply(Y != 0, 3, rowSums) != 0, 2, sum)/max(apply(apply(Y !=
                                                                                                    0, 3, rowSums) != 0, 2, sum))
                            if (verbose) {
                              cat(paste("Fitting null model...\n", sep = ""))
                            }
                            p <- ncol(V)
                            epsLB <- 1e-07
                            resNullFullV <- nullModelFitting(theta_start = c(rep(1, p), 1e-3),
                                                             M = M,
                                                             y = y,
                                                             V = V,
                                                             psiLB = rep(epsLB, p + 1),
                                                             psiUB = rep(+Inf, p + 1),
                                                             maxit = maxit,
                                                             factr = factr,
                                                             pgtol = pgtol,
                                                             lmm = lmm,
                                                             VCNs = VCNs,
                                                             nObs = nObs,
                                                             trace = trace)
                            if (verbose) {
                              cat(" DONE\n")
                            }
                            return(resNullFullV)
                          }

                          assignInNamespace("fit.null",fit.null,ns="RestoreNet")
                          rm(fit.null)

                          #############
                          ## folders ##
                          #############

                          resPath <- paste("./results-sim1-s2", s20, "/", sep = "")
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
                                                         c("nPar",
                                                           "cll",
                                                           "mll",
                                                           "cAIC",
                                                           "mAIC",
                                                           "Chi2",
                                                           "p-value",
                                                           "KLdiv",
                                                           "KLdiv/N",
                                                           "BhattDist_nullCond",
                                                           "BhattDist_nullCond/N"),
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

                            saveRDS(Y, paste(resPath, "/sim", sim, "_S", S, "_s20", s20, ".rds", sep = ""))

                            Y <- Y[which(apply(Y != 0, 1, sum) > 0),,]
                            dimnames(Y)[[3]] <- 1:dim(Y)[3]
                            Y_NA <- Y
                            Y_NA[Y_NA == 0] <- NA

                            ######################
                            ## null model fitting:
                            ######################
                            null.res <- fit.null(Y = Y,
                                                 rct.lst = rcts,
                                                 maxit = 1000)

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
                            cat("DONE\n")
                          }

                          save.image(paste(resPath, "/simulations_s2", s20,".RData", sep = ""))
                          gc()

                          res <- list()
                          res$res.base.par <- res.base.par
                          res$res.re.Euy <- res.re.Euy


                          return(res)
                        })

stopCluster(cl)

##########
## plots:
##########

norm2 <- function(x){
  return(sqrt(sum((x)^2)))
}

relErr <- function(x, y){
  norm2(x - y)/norm2(y)
}

nSim <- 100

relErr.all <- sapply(1:length(s20_all), function(j){
  rcts <- c("A->1", "B->1", "C->1", "D->1",
            "A->0", "B->0", "C->0", "D->0",
            "A->B", "A->C", "C->D")

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

  s20 <- s20_all[j]
  currResPath <- currDir

  GLS.file <- read.table(paste(currResPath, paste0("/results-sim1-s2",s20,"/sim1-s2",s20,"_GLSResults.tsv"), sep = ""), sep = "\t", header = F)
  GLS.res <- Reduce(rbind, lapply(1:3, function(j){t(GLS.file[,1:length(u_1)])}))
  GLS.res <- cbind(GLS.res, GLS.res[,1])

  res.re.Euy <- res_alls20[[j]]$res.re.Euy
  RestoreNet.res <- res.re.Euy

  relErr.GLS <- sapply(1:nSim, function(j){
    relErr(GLS.res[,j], c(u_1, u_2, u_3))
  })

  relErr.RestoreNet <-sapply(1:nSim, function(j){
    relErr(RestoreNet.res[,j], c(u_1, u_2, u_3))
  })

  relErr.all_s2 <- cbind(GLS = relErr.GLS,
                         RestoreNet = relErr.RestoreNet)

  return(relErr.all_s2)
}, simplify = "array")

dimnames(relErr.all)[[2]] <- c("GLS", "RestoreNet")
dimnames(relErr.all)[[3]] <- c(0.1, 1, 10)

pdf(file = paste(currDir_time, "/allBoxplots_sim1.pdf", sep = ""), width = 10, height = 3)
par(mar = c(3,3,2,1), mfrow = c(1,3))
boxplot(relErr.all[,,"0.1"],
        outline = FALSE,
        range = 0,
        lwd = 1.5,
        cex.axis = 1.5,
        main = bquote(sigma^2~"="~.1),
        cex.main = 2,
        ylab = "",
        cex.lab = 2,
        ylim = range(relErr.all[,,1:3]))
boxplot(relErr.all[,,"1"],
        outline = FALSE,
        range = 0,
        lwd = 1.5,
        cex.axis = 1.5,
        main = bquote(sigma^2~"="~1),
        cex.main = 2,
        ylim = range(relErr.all[,,1:3]))
boxplot(relErr.all[,,"10"],
        outline = FALSE,
        range = 0,
        lwd = 1.5,
        cex.axis = 1.5,
        main = bquote(sigma^2~"="~10),
        cex.main = 2,
        ylim = range(relErr.all[,,1:3]))
dev.off()

save.image(paste(currDir_time, "/sim1_alls20.RData", sep = ""))
gc()







