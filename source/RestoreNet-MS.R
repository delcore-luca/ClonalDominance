cat("\nInstall/load packages")

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
            "RestoreNet")

lapply(l.pkgs, function(pkg){library(pkg,character.only=TRUE)})

rm(list = ls())

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

load(file = "./data/Y_MS.rda")

#############
## vectorID
#############
# vectorID <- "PGK"
vectorID <- commandArgs(trailingOnly = TRUE)
Y <- Y_MS[[vectorID]]

topClones <- 1000
Y <- Y[,,names(head(sort(apply(Y!=0, 3, sum), decreasing = T), topClones)),drop=FALSE]
cat(paste("n. of clones: ", dim(Y)[3], "\n", sep = ""))

tps <- as.numeric(rownames(Y))
tps <- (tps - min(tps))/(max(tps) - min(tps))
rownames(Y) <- tps

Y <- Y[which(apply(Y != 0, 1, sum) > 0),,]
Y_NA <- Y
Y_NA[Y_NA == 0] <- NA

####################################
## Initialize matrix for IC results
####################################
nModels <- 2
resIC <- matrix(data = NA, nrow = nModels, ncol = 11)
colnames(resIC) <- c("nPar",
                     "cll",
                     "mll",
                     "cAIC",
                     "mAIC",
                     "Chi2",
                     "p-value",
                     "KLdiv",
                     "KLdiv/N",
                     "BhattDist_nullCond",
                     "BhattDist_nullCond/N")
rownames(resIC) <- c("null", "re")
modsPars <- list()
resPAR <- list()
P <- list()


######################
## null model fitting:
######################
null.res <- fit.null(Y = Y,
                     rct.lst = c(
                       "M->1",
                       "B->1",
                       "T->1",
                       "M->0",
                       "B->0",
                       "T->0"))

resIC["null", names(null.res$stats)] <- null.res$stats
resPAR[["null"]] <- null.res$fit$par
P[["null"]] <- null.res$stats["nPar"]

#################################
## random-effects model fitting:
#################################
re.res <- fit.re(theta_0 = null.res$fit$par,
                 Y = Y,
                 rct.lst = c(
                   "M->1",
                   "B->1",
                   "T->1",
                   "M->0",
                   "B->0",
                   "T->0"),
                 maxemit = 1000)

resIC["re", ] <- re.res$stats
resPAR[["re"]] <- re.res$fit$par
P[["re"]] <- length(re.res$fit$par)

########### STORE RESULTS ###########
write.table(resIC, file = paste(resPath, "/", vectorID, "_AICs.tsv", sep = ""), sep = "\t", row.names = TRUE, col.names = TRUE)
print(xtable(resIC[, c("nPar", "mAIC", "KLdiv", "KLdiv/N")]), file = paste(resPath, "/", vectorID, "_AICsforLatex.txt", sep = ""))


########### STORE ENVIRONMENT ###########
save.image(paste(resPath, "/", vectorID, "_LLAwithAIC.RData", sep = ""))
gc()


########### CLONAL TRACKING PLOT ###########
linColMap <- brewer.pal(nrow(null.res$design$V), "Set2")
names(linColMap) <- rownames(null.res$design$V)

pdf(file = paste(figPath, "vectorID_", vectorID, "_clonalTrackingData.pdf", sep = ''), height=5, width=7)
par(mar = c(5,5,2,2))
matplot(as.numeric(rownames(Y_NA)), log(Y_NA[,,1]), lty = 1, pch = 20, col = alpha(linColMap, alpha = .2), cex = 3, lwd = 3,
        ylim = range(c(log(Y_NA)), na.rm = T), type = 'b', cex.axis = 2, cex.lab = 2, ylab = expression(logY[t]), xlab = "time (months)", xaxt = "n")
lapply(1:dim(Y)[3], function(cl){
  matplot(as.numeric(rownames(Y_NA)), log(Y_NA[,,cl]), lty = 1, pch = 20, add = T, col = alpha(linColMap, alpha = .2), cex = 2, type = 'b', lwd = 3)
})
axis(1, at=as.numeric(rownames(Y_NA)),
     labels=rownames(Y_RM[[vectorID]]),
     cex.axis = 2, cex.lab = 2)
matplot(as.numeric(rownames(Y_NA)), apply(log(Y_NA), c(1,2), mean, na.rm = T), type = 'l', lty = 1, col = linColMap, add = T, lwd = 10)
# legend(x = "topleft", legend = colnames(Y), col = linColMap, pch = 20, lwd = 5, cex = 1.5)
dev.off()


linColMap <- brewer.pal(nrow(null.res$design$V), "Set2")
names(linColMap) <- rownames(null.res$design$V)

rownames(Y) <- rownames(Y_RM[[vectorID]])

########### STORE PARAMETERS ###########
xsol <- resPAR$re
psol <- (length(xsol) - 1)/2
rndEffParsMat <- matrix(xsol[1:(2*psol)], nrow = psol, ncol = 2)
rownames(rndEffParsMat) <- colnames(null.res$design$V)

colnames(rndEffParsMat) <- c("theta", "sd")
rndEffParsMat[,"sd"] <- sqrt(rndEffParsMat[,"sd"])
print(xtable(rndEffParsMat, digits = 3), include.rownames=FALSE, file = paste(resPath, "/", vectorID, "REparForLatex.txt", sep = ""))

########### BOXPLOTS ###########
pdf(file = paste(figPath, "vectorID_", vectorID, "_boxplots.pdf", sep = ''), height=10, width=10)
get.boxplots(re.res)
dev.off()

########### SCATTERPIE PLOT ###########
pdf(file = paste(figPath, "vectorID_", vectorID, "_scatterpie.pdf", sep = ''), height=10, width=10)
get.scatterpie(re.res, legend = TRUE)
dev.off()

####### stacked bar plots #######
sc_CTs <- lapply(colnames(Y), function(CT){
  pal <- colorRampPalette(c(linColMap[CT], "black"))
  COLS = pal(length(as.vector(unlist(dimnames(Y)[3]))))
  COLS <- COLS[order(apply(Y[,CT,], 2, sum), decreasing = F)]
  Y_CT_df <- as.data.frame.table(Y[,CT,])
  # Y_CT_df <- as.data.frame.table(Y[,CT,]/rowSums(Y[,CT,]))
  colnames(Y_CT_df) <- c("time", "clone", "count")

  sc <- (ggplot(Y_CT_df,
                aes(x = time,
                    y = count,
                    group = clone,
                    fill = factor(clone)))
         +geom_area(size=.3, colour="white")
         + theme_bw() + scale_fill_manual(values=COLS)
         + theme_classic() + labs(y= "barcode count", x = "time (months)")
         + theme(legend.position="none",axis.text=element_text(size=20),
                 axis.title=element_text(size=20,face="bold")));
  return(sc)
})

pdf(file = paste(figPath, "/vectorID_", vectorID, "_stackedArea.pdf", sep = ''), height=16, width=8)
ggarrange(plots = sc_CTs, nrow = length(sc_CTs), ncol = 1)
dev.off()


