##Don't have those libraries? Try using 
## install.packages(reshape2)
## and so on, for each of the libraries here. It will be quicker than using the menu.

library(reshape2)
library(data.table)
library(psych)
library(Hmisc)
library(corrplot)
library(apaTables)
library(rtf)
library(plyr)
library(IAT)
library(dplyr)
library(doBy)
library(ggplot2)
library(graphics)
library(gplots)
library(gridExtra)
library(grid)
library(nloptr)
library(car)
library(schoRsch)
library(ez)
library(corrgram)
library(ggrepel)
library(lsmeans)
library(effsize)
library('scales')
library('ggthemes')
library(jmv)
library(lavaan)

options(scipen=999)

#########################################
#########HELPER FUNCTION#################

factorToNumeric <- function (inF)
{
  return (as.numeric(levels(inF))[inF])
}

factorToNumeric.many <- function(df, inVars = NULL, inFormula = NULL)
{
  if (is.null(inVars))
  {
    sdvs <- as.character(inFormula)[2]
    vars <- strsplit(sdvs," \\+ ")[[1]]
  }
  else
  {
    vars = inVars
  }
  df[,vars] = apply(df[,vars], 2, factorToNumeric)
  return(df)
}

charToNumeric.many <- function(df, inVars = NULL, inFormula = NULL)
{
  if (is.null(inVars))
  {
    sdvs <- as.character(inFormula)[2]
    vars <- strsplit(sdvs," \\+ ")[[1]]
  }
  else
  {
    vars = inVars
  }
  df[,vars] = apply(df[,vars], 2, function(x) as.numeric(as.character(x)))
  return(df)
}


#Show table in the plot output.
showTable <- function(theT, inFile=NULL, inWidth=12, inHeight=12)
{
  if (!is.null(inFile))
  {
    pdf(inFile,width=inWidth,height=inHeight)
  }
  
  g <- tableGrob(theT)
  grid.newpage()
  grid.draw(g)

  if (!is.null(inFile))
  {
    dev.off()
  }
}

showTablesList <- function(tables)
{
  gs <- tableGrob(tables[[1]])
  tables[[1]] <- NULL
  for (tbl in tables) 
  {
    gt <- tableGrob(tbl)
    gs <- combine(gs,gt,along=2)
  }
  grid.newpage()
  grid.draw(gs)
  return(gs)
}

showTablesMany <- function(...)
{
  tables <- list(...)
  return(showTablesList(tables))
}

html.showTablesList <- function(tables)
{
  #padding <- "padding-left: .7em; padding-right: .7em;"
  hhh <- as.character(htmlTable::htmlTable(tables[[1]]))
  tables[[1]] <- NULL
  for (tbl in tables) 
  {
    hhh1 <- htmlTable::htmlTable(tbl)
    #hhh <- paste0(hhh,'<br/><br/>', as.character(hhh1))
    hhh <- paste0(hhh,as.character(hhh1))
  }
  vvv <- getOption('viewer')
  htmlfile<-tempfile(fileext='.html')
  cat(hhh,file=htmlfile)
  vvv(htmlfile)
  return(hhh)
}

##Creates a correlation matrix table
cornp <- function(x, addVarNames=T, type='pearson'){ 
  require(Hmisc)
  nums <- sapply(x, is.numeric)
  x <- as.matrix(x[ , nums])
  rc <- rcorr(x, type = type)
  R <- rc$r 
  p <- rc$P 
  n <- rc$n
  
  ## define notions for significance levels; spacing is important.
  ps <- ifelse(p < .001, "< .001", as.character(round(p,3)))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- round(R, 3)
  
  ## build a new matrix that includes the correlations with their apropriate stars 
  Rnew <- matrix(paste(R, ps, n, sep="\n"), ncol=ncol(x)) 
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="") 
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew) 
  
  ## remove last column
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  ##remove the first row
  Rnew <- Rnew[2:nrow(Rnew),]
  if (addVarNames)#Add the rownames as the first column, if needed.
  {
    Rnew$varName <- rownames(Rnew)
    Rnew <- Rnew[,c(ncol(Rnew),1:(ncol(Rnew)-1))]
  }
  #showTable(Rnew)
  return(Rnew) 
}

#Create a correlation table from x and y, with a SAS-like format (r, p, n).
cornpxy <- function(x, y, addVarNames=T, method='pearson'){ 
  require(psych)
  nums <- sapply(x, is.numeric)
  x <- x[ , nums, drop=F]
  nums <- sapply(y, is.numeric)
  y <- y[ , nums, drop=F]

  ct <- corr.test(as.matrix(x), as.matrix(y), method=method, use = 'pairwise', adjust='none')
  
  #rc <- rcorr(x)
  R <- ct$r 
  p <- ct$p 
  n <- ct$n
  
  ## define notions for significance levels; spacing is important.
  ps <- ifelse(p < .001, "< .001", as.character(round(p,3)))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- round(R, 3)
  
  ## build a new matrix that includes the correlations with their apropriate stars 
  Rnew <- matrix(paste(R, ps, n, sep="\n"), ncol=ncol(R)) 
  #diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(y), "", sep="") 
  
  ## remove upper triangle
  Rnew <- as.data.frame(Rnew) 
  if (addVarNames)#Add the rownames as the first column, if needed.
  {
    Rnew$varName <- rownames(Rnew)
    Rnew <- Rnew[,c(ncol(Rnew),1:(ncol(Rnew)-1))]
  }
  
  ## remove last column and return the matrix (which is now a data frame)
  #showTable(Rnew)
  return(Rnew) 
}

#Run cornpxy by variables. 
corby <- function(inData, xs, ys, byvars)
{

  func <- function(xx)
  {
    return(cornpxy(xx[,xs, drop=F], 
                   xx[,ys, drop=F], 
                   addVarNames = T))
  }
  corTTT <- ddply(inData, as.quoted(byvars), func)
  #showTable(corTTT)
  return(corTTT)
}

#Run cornpxy by variables, but show only the correlations. 
corby.ronly <- function(inData, xs, ys, byvars, round.mat=T, addSig=F, sigs=c(.05, .01, .001))
{
  func <- function(xx)
  {
    ttt <- corr.test(as.matrix(xx[,xs, drop=F]),  as.matrix(xx[,ys, drop=F]), 
                     method='pearson', use = 'pairwise', adjust='none')
    if (round.mat)
    {
      ttt <- structure(rapply(ttt,function(x) if(is.numeric(x)) round(x,3) else x,how="replace"),
                class="htest")
    }
    
    if (addSig)
    {
      stars <- ifelse(is.na(ttt$p), '', 
                      ifelse(ttt$p < sigs[3], "***", 
                             ifelse(ttt$p < sigs[2], '**', 
                                    ifelse(ttt$p < sigs[1], '*', 
                                           ''))))
      ttt$r <- matrix(paste0(ttt$r, stars), 
                      nrow=nrow(ttt$r), ncol=ncol(ttt$r), 
                      dimnames = list(row.names(ttt$r), colnames(ttt$r)))
    }
    
    ttt2 <- as.data.frame(ttt$r)
    ttt2$vs <- xs
    return(ttt2)
  }
  corTTT <- ddply(inData, as.quoted(byvars), func)
  #showTable(corTTT)
  return(corTTT)
}

#Create a table of partialled out correlations 
# pout is the variable to partial out
mypcorl <- function(inData, xs, ys, pout, nround = 3)
{
  library(ppcor)
  
  dt <- data.table(inData)
  pVar <- dt[,get(pout)]
  pcs <- list()
  for (x in xs)
  {
    xVar <- dt[,get(x)]
    for (y in ys)
    {
      yVar <- dt[,get(y)]
      pcs[[length(pcs)+1]] <- pcor.test(xVar, yVar, pVar)
      pcs[[length(pcs)]]$x <- x
      pcs[[length(pcs)]]$y <- y
      pcs[[length(pcs)]]$pout <- pout
    }
  }
  outt <- data.frame(matrix(unlist(pcs), nrow=length(pcs), byrow=T),stringsAsFactors=FALSE)
  outt <- outt[, !(colnames(outt) %in% c("X3","X5","X6"))]
  setnames(outt, c('X1', 'X2', 'X4', 'X7', 'X8', 'X9'), c('pcor', 'pval', 'n', 'x', 'y', 'pout'))
  outt <- outt[,c('x', 'y', 'pout', 'pcor', 'pval', 'n')]
  
  if (nround>0)
  {
    outt$pcor <- round(as.numeric(outt$pcor),nround)
    outt$pval <- round(as.numeric(outt$pval),nround)
    outt$pval <- ifelse(outt$pval==0, '<.001', as.character(outt$pval))
  }
  
  return(outt)
}

#Create a SAS-like table of partialled out correlations 
# pout is the variable to partial out
mypcor <- function(inData, xs, ys, pout, nround = 3)
{
  pt <- mypcorl(inData, xs, ys, pout, nround)
  outt <- as.data.frame(xs)
  for (x in xs)
  {
    for (y in ys)
    {
      line <- pt[pt$x==x & pt$y==y,]
      addWhat <- paste(line$pcor, line$pval, line$n, sep='\n')
      colName <- paste(line$y, line$pout, sep='|')
      outt[[colName]][outt$xs==x] <- addWhat
    }
  }
  setnames(outt, c('xs'), c('varName'))
  return(outt)
}

##For creating a correlaiton matrix plot for the corrgram function
panel.shadeNtext <- function (x, y, corr = NULL, col.regions, ...) 
{
  #corr <- cor(x, y, use = "pair")
  rc <- rcorr(as.matrix(x), as.matrix(y))
  corr <- rc$r[1,2]
  n <- rc$n[1,2]
  #results <- cor.test(x, y, alternative = "two.sided")
  #est <- round(results$p.value, 3)
  est <- round(rc$P[1,2], 3)
  stars <- ifelse(est < .001, "< .001", as.character(est))
  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1, 
                                               length = ncol + 1), include.lowest = TRUE))
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind], 
       border = NA)
  box(col = "lightgray")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- formatC(corr, digits = 3, format = "f")
  cex.cor <- 1/strwidth("-X.xxx")
  #cex.cor <- .6/strwidth("-X.xxx")
  fonts <- ifelse(stars != "", 2,1)
  # option 1: stars:
  #text(0.5, 0.4, paste0(r,"\n", stars), cex = cex.cor)
  text(0.5, 0.7, r, cex = cex.cor*0.6)
  text(0.5, 0.4, paste0( '(', stars, ')' ), cex = cex.cor*0.5)
  text(0.5, 0.2, paste0( '(', n, ')' ), cex = cex.cor*0.5)
  # option 2: bolding:
  #text(0.5, 0.5, r, cex = cex.cor, font=fonts)
}

#Create corrgram for all the variables provided in dt
my.corrgram <- function(dt)
{
  return (corrgram(dt, 
           order=TRUE, lower.panel=panel.shadeNtext,
           upper.panel=panel.pts, 
           diag.panel=panel.minmax, 
           main="My matrix"))
}

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="lavender", ...)
}
## put correlations & 95% CIs on the upper panels,
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y,use="complete.obs")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  prefix <- "r = "
  rc <- cor.test(x,y)
  rci <- rc$conf.int
  txt2 <- format(c(rci, 0.123456789), digits=digits)[1]
  txt3 <- format(c(rci, 0.123456789), digits=digits)[2]
  prefix2 <- "\nCI = "
  txt <- paste(prefix, txt, prefix2, txt2, ", ", txt3, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1)
}


my.matrixConf <- function(dt)
{
  return (pairs(dt, lower.panel=panel.smooth, cex = .8, pch = 21, bg="steelblue",
        diag.panel=panel.hist, cex.labels = 1.2, font.labels=2, upper.panel=panel.cor))
}

#Function to compute mean without outliers
meanNoOutliers <- function(xx, nSDs = 3)
{
  zxx <- scale(xx)
  return (mean(xx[abs(zxx)<nSDs], na.rm=TRUE))
}


##For providing basic descriptives
mysum <- function(x, 
                  shouldSE = TRUE, shouldSD = TRUE, shouldT = FALSE, 
                  shouldMed = TRUE, shouldPvalue = FALSE, probs = NULL) { 
  ret <- c(n = sum(!is.na(x)), M = round(mean(x, na.rm=TRUE),3))
  if (shouldSD)
  {
    ret <- c(ret, SD = round(sd(x, na.rm=TRUE),3))
  }
  if (shouldSE)
  {
    ret <- c(ret, SE = round(sqrt(var(x, na.rm=TRUE)/length(!is.na(x))),3))
  }
  if (shouldMed)
  {
    ret <- c(ret, med = round(median(x, na.rm=TRUE),3))
  }
  if (shouldPvalue | shouldT)
  {
    ttt <-  t.test(x, mu=0)
    pVal <- ifelse(ttt$p.value <= .0001, .0001, ttt$p.value)
    tVal <- ttt$statistic
    dfVal <- ttt$parameter
    if (shouldPvalue)
    {
      ret <- c(ret, p = pVal)  
    }
    if (shouldT)
    {
      ret <- c(ret, t = tVal)  
      ret <- c(ret, df = dfVal)  
    }
  }
  if (!is.null(probs))
  {
    ret <- c(ret, round(quantile(x, probs=probs, na.rm=TRUE),3))
  }
  return(ret)
}

#Run mysym by a factor
mysumBy <- function(formula, dt, shortCols=F, bindTables=T, ...)
{
  sdvs = as.character(formula)[2]
  ivs <- as.character(formula)[3]
  dvs <- strsplit(sdvs," \\+ ")[[1]]
  tbls <- list()
  index <- 1
  for (dv in dvs)
  {
    frm <- paste(dv,'~', ivs)
    tbl <- summaryBy(as.formula(frm), data = dt, FUN = mysum, ...)
    #browser()
    if (shortCols | bindTables)
    {
      #nc <- ncol(tbl)
      #st <- nc-4
      #colnames(tbl)[st:nc] <- c('N', 'M', 'SD', 'SE', 'Med')
      colnames(tbl) <- gsub(paste0(dv,"."), "", colnames(tbl), fixed = TRUE)
      tbl <- cbind(var = dv, tbl)
    }
    tbls[[index]] <- tbl
    index <- index+1
  }
  if (bindTables)
  {
    tbls <- ldply(tbls, rbind)
  }
  tbls <- tbls[tbls$n > 0,]
  return (tbls)
}

###Run mysum to a few variables.
##Usage: my.sumMany(formula = DV1 + DV2 + DV3 ~ factor1 + factor2, dt=mydata)
my.sumMany <- function(formula, dt, shortCols=F, bindTables=T, ...)
{
  sdvs = as.character(formula)[2]
  dt$dummy <- as.factor('-')
  frm <- paste0(sdvs, ' ~ dummy')
  ttt <- mysumBy(as.formula(frm), dt, shortCols, bindTables, ...)
  if (bindTables)
  {
    ttt$dummy <- NULL
  }
  else
  {
    ttt <- lapply(ttt, function(x) { x["dummy"] <- NULL; x })
  }
  return(ttt)
}
mysum.many <- my.sumMany

###Run mysum to a few variables and present in a table.
mysumMany.table <- function(formula=NULL, theDVs=NULL, dt, ...)
{
  if (!is.null(formula))
  {
    sdvs = as.character(formula)[2]
    dvs <- strsplit(sdvs," \\+ ")[[1]]
  }
  else
  {
    dvs = theDVs
  }
  tbls <- list()
  index <- 1
  for (dv in dvs)
  {
    tbl <- as.data.frame(t(mysum(dt[,dv], ...)))
    tbl$name <- dv
    tbls[[index]] <- tbl
    index <- index+1
  }
  tbls <- ldply(tbls, rbind)
  return (tbls)
}

###Run a regression analysis from data and formula 
##(standardize all the variables)
myReg <- function (indata, form)
{
  #Get data
  forreg <- indata[,all.vars(form)]
  #Standardize
  forreg.s <- as.data.frame(scale(forreg))

  #Run multiple regression analysis
  fit <- lm(formula = form, data=forreg.s)
  #Get stuff for the table
  sumreg <- summary(fit)
  sumreg$call <- form

  return(sumreg)
}

###Run a regression analysis from data and formula 
##(standardize all the variables)
myReg.lm <- function (indata, form)
{
  #Get data
  forreg <- indata[,all.vars(form)]
  #Standardize
  forreg.s <- as.data.frame(scale(forreg))
  
  #Run multiple regression analysis
  fit <- lm(formula = form, data=forreg.s)
  #Get stuff for the table
  return(fit)
}

###Melt using a formula
mymelt <- function(dt, formula)
{
  allvars = all.vars(formula)
  idvars = all.vars(terms(formula)[[3]])
  return (reshape2::melt(dt[,allvars],id.vars=idvars))  
}

handleOutliersBy <- function(df, targetVars, byVars, cutOffSD=2.5)
{
  #convert to data.table, using only the relevant variables.
  DT <- setDT(df[,c(targetVars, byVars)])
  #Standardize the target variables by the group variables.
  DT <- DT[, paste0(targetVars, ".z") := lapply(.SD, 
                                                function(x) as.vector(scale(x))), by=byVars]
  #Create max and min values, based on the cut-off argument.
  DT <- DT[, paste0(targetVars, ".min") := lapply(.SD, 
                                                  function(x) as.vector(mean(x,na.rm=T)-(cutOffSD*sd(x,na.rm=T)))), by=byVars]
  DT <- DT[, paste0(targetVars, ".max") := lapply(.SD, 
                                                  function(x) as.vector(mean(x,na.rm=T)+(cutOffSD*sd(x,na.rm=T)))), by=byVars]
  #Convert back to data.frame
  dfn <- setDF(DT)

  #Create Winsorized and trimmed variables for all the target variables.
  for (varName in targetVars) {
    #Create some variable names from the target variable.
    varNameW <- paste0(varName, ".w")
    varNameT <- paste0(varName, ".t")
    varNameMax <- paste0(varName, ".max")
    varNameMin <- paste0(varName, ".min")
    #By default the winsorized value is the same as the original value.
    dfn[,c(varNameW)] <- dfn[ ,c(varName)]
    #But if the original value is larger than the maximum value allowed by the cut-off, we cut it into the maximum value.
    condMatMax <- dfn[,c(varName)] > dfn[,c(varNameMax)]
    ##YBYB: For some reason condMatMax sometimes has NA rows. We'll try to ignore them and do nothing with them
    condMatMax <- ifelse(is.na(condMatMax),FALSE, condMatMax)
    dfn[condMatMax, c(varNameW)] <- dfn[condMatMax,c(varNameMax)]
    #But if the original value is smaller than the minimum value allowed by the cut-off, we cut it into the minimum value.
    condMatMin <- dfn[,c(varName)] < dfn[,c(varNameMin)]
    ##YBYB: For some reason condMatMax sometimes has NA rows. We'll try to ignore them and do nothing with them
    condMatMin <- ifelse(is.na(condMatMin),FALSE, condMatMin)
    dfn[condMatMin,c(varNameW)] <- dfn[condMatMin,c(varNameMin)]
    ##The trimmed variables:
    dfn[,c(varNameT)] <- dfn[ ,c(varName)]
    #If the original value is larger than the maximum value allowed we omit it.
    dfn[condMatMax,c(varNameT)] <- NA
    #If the original value is smaller than the minimum value allowed we omit it.
    dfn[condMatMin,c(varNameT)] <- NA
  }
  dfr <- cbind(df, dfn[,c(paste0(targetVars, ".w"), 
                          paste0(targetVars, ".z"), 
                          paste0(targetVars, ".t"))])
  return (dfr)
}

my.t.test.p.value <- function(...) {
      obj<-try(t.test(...), silent=TRUE)
      if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

myDescribeBy <- function(inData, inBy, ranges=FALSE, D=F, nround=3) { 

  #Get the description
  myD <- describeBy(inData, inBy, skew=FALSE, ranges=ranges, mat=TRUE)  
  if (D) {
    myD$D <- myD$mean/myD$sd
  }
  #if (ttest) {
  #  pvals<- summarise_each(inDataBy, funs(my.t.test.p.value(., mu=0, na.rm=TRUE)))
  #}
  
  if (nround>0)
  {
    is.num <- sapply(myD, is.numeric)
    myD[is.num] <- lapply(myD[is.num], round, nround)
  }
  
    
  theNames <- names(inData)
  myD$item <- theNames[myD$vars]
  myD$vars <- NULL
  rownames(myD) <- c()
  return(myD)
}

#Provide a better frequencies table 
my.freq <- function(x, name='value')
{
  t1 <- table(x, exclude = NULL)
  t11 <- as.data.table(t1)
  ttt <- as.data.table(cbind(Freq=t11, CumFreq=cumsum(t1), prcnt=prop.table(t1), cumPrcnt=cumsum(prop.table(t1))))
  setnames(ttt, old=c('Freq.x', 'Freq.N', 'prcnt.N'), new=c(name,'N','Prcnt'))
  ttt$Prcnt <- round(ttt$Prcnt,3)
  ttt$cumPrcnt <- round(ttt$cumPrcnt,3)
  ttt$prcnt.x <- NULL
  return(ttt)
}

my.freq2d <- function(v1, v2, exclude=NULL)
{
  t1 <- table(v1, v2, exclude=exclude)
  t11 <- as.data.table(t1)
  ttt <- as.data.table(cbind(Freq=t11, prcnt=prop.table(t1)))
  ttt <- ttt[,c('Freq.v1','Freq.v2','Freq.N','prcnt.N')]
  names(ttt) <- c('v1','v2','N','Prcnt')
  ttt$Prcnt <- round(ttt$Prcnt,3)
  return(ttt)
}

#Helper to provide mean and N
meanN <- function(x) { 
  c(
    n = sum(!is.na(x)),
    M = mean(x, na.rm=TRUE)
  )
}

#For ggplot2: creates a few statistics that we can use for graphing
min.mean.sd.max <- function(x) {
  #r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  xsd <- sd(x, na.rm=TRUE)
  xn <- sum(!is.na(x))
  #xerr <- qnorm(0.975)*xsd/sqrt(xn)
  xmin = min(x, na.rm = TRUE)
  xmax = max(x, na.rm = TRUE)
  xmean =  mean(x, na.rm=TRUE)
  r <- c(xmin, xmean-xsd, xmean, xmean + xsd, xmax)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#For ggplot2: creates a few statistics that we can use for graphing (I often use this one for creating box-plots)
min.median.sd.max <- function(x) {
  #r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  xsd <- sd(x, na.rm=TRUE)
  xn <- sum(!is.na(x))
  #xerr <- qnorm(0.975)*xsd/sqrt(xn)
  xmin = min(x, na.rm = TRUE)
  xmax = max(x, na.rm = TRUE)
  xmean =  mean(x, na.rm=TRUE)
  xmedian =  median(x, na.rm=TRUE)
  r <- c(xmin, xmean-xsd, xmedian, xmean + xsd, xmax)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#For ggplot2
mean.sd.ci <- function(x) {
  #r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  xsd <- sd(x, na.rm=TRUE)
  #xerr <- qnorm(0.975)*xsd/sqrt(xn)
  xmean =  mean(x, na.rm=TRUE)
  ci <- getCI(x)
  r <- c(xmean, xsd, ci[1], ci[2])
  names(r) <- c("mean", "sd", "ci.low", "ci.high")
  return(r)
}

#For ggplot2
mean.sd.ci.n <- function(x) {
  #r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  r <- mean.sd.ci(x)
  r[length(r)+1] <- length(x)
  names(r)[length(r)] <- 'n'
  return(r)
}

#Compute rowmeans from the left side of a formula
#Uasge: allds$var.mean <- my.rowMean(dt=allds, var1 + var2 ~ .)
my.rowMean <- function(dt, formula=NULL, dvs=NULL, na.rm=T)
{
  if (!is.null(formula))
  {
    sdvs = as.character(formula)[2]
    dvs <- strsplit(sdvs," \\+ ")[[1]]
  }
  return (rowMeans(dt[,dvs], na.rm=na.rm))
}

my.rowMean.z <- function(dt, formula=NULL, dvs=NULL, na.rm=T)
{
  if (!is.null(formula))
  {
    sdvs = as.character(formula)[2]
    dvs <- strsplit(sdvs," \\+ ")[[1]]
  }
  dt.t <- scale(dt[,dvs])
  return (rowMeans(dt.t, na.rm=na.rm))
}


my.rowMean.fisher <- function(dt, formula=NULL, dvs=NULL, na.rm=T)
{
  if (!is.null(formula))
  {
    sdvs = as.character(formula)[2]
    dvs <- strsplit(sdvs," \\+ ")[[1]]
  }
  fdt <- psych::fisherz(dt[,dvs])
  return (psych::fisherz2r(rowMeans(fdt, na.rm=na.rm)))
}


#Compute mean of a row, after reversing reversed items. 
##Returns the data with the reversed variables and with the mean, as the
my.rowMeanRev <- function(dt, rev, no.rev, mVar, deductFrom=NULL, multiplyBy=NULL)
{
  #Slice only the data we want to reverse.
  dt.rev <- dt[,rev]
  
  if (!is.null(deductFrom))
  {#Reverse by deducting (usually, from max-response+1)
    dt.rev <- dt.rev %>% mutate_all(funs(deductFrom-.))
  }
  if (!is.null(multiplyBy))
  {#Reverse by multiplying (usually, by -1)
    dt.rev <- dt.rev %>% mutate_all(funs(multiplyBy*.))
  }
  #Rename the reversed variables
  names(dt.rev) <- paste0('r.', names(dt.rev))
  #Add the reversed variables to the data.frame
  dt <- cbind(dt,dt.rev)
  #Compute the mean after the reverse
  dt[,mVar] <- my.rowMean(dt=dt, dvs = c(names(dt.rev),no.rev))
  return(dt)
}

#Count how many of the variables in the left side of the formula are non-missing
#Uasge: allds$nVars <- my.countValues(dt=allds, var1 + var2 ~ .)
my.countValues <- function(dt, formula=NULL, dvs=NULL)
{
  if (!is.null(formula))
  {
    sdvs = as.character(formula)[2]
    dvs <- strsplit(sdvs," \\+ ")[[1]]
  }
  return (rowSums(!is.na(dt[,dvs])))
}


#Get the bottom and top of a confidence interval of a mean (not in the context of ANOVA). Used for creating box-plot.
##YBYB: Cannot use for publications, unless your test was a t-test because you had only two groups.
getCI <- function(x) {
  #r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  xsd <- sd(x, na.rm=TRUE)
  xn <- sum(!is.na(x))
  xerr <- qnorm(0.975)*xsd/sqrt(xn)
  xmean <-  mean(x, na.rm=TRUE)
  xmin <- xmean-xerr
  xmax <- xmean+xerr
  r <- c(xmin, xmax)
  names(r) <- c("ymin", "ymax")
  r
}

get.poolSD <- function(sds, ns) {
  aboves <- (ns-1)*sds*sds
  above <- sum(aboves)
  belows <- ns-1
  below <- sum(belows)
  varPool <- above/below
  sdPool <- sqrt(varPool)
  return(sdPool)
}

##ggplot2 helper
n_fun <- function(x){
  return(data.frame(y = max(x)*1.1, label = paste0("n = ",sum(!is.na(x)))))
}

#Show a scatter plots by variables (use the wrapper function defined below)
scatterBy <- function(inData, title='A scatter plot', 
                      xTitle = 'IV', yTitle = 'DV', points='points', freeY=T)
{
  if (is.null(inData$group))
  {
    inData$group <- ''
  }
  if (is.null(inData$facetCond1))
  {
    formula <- NULL
    inData$facetCond1 = '1'
    inData$facetCond2 = '1'
  }
  else if (is.null(inData$facetCond2))
  {
    formula <- ~facetCond1
    inData$facetCond2 = '1'
  }
  else
  {
    formula <- facetCond1 ~ facetCond2
  }
  
  facetg <- NULL
  if (!is.null(formula))
  {
    if (freeY)
    {
      facetg <- facet_wrap(formula, scales="free_y") 
    }
    else 
    {
      facetg <- facet_grid(formula) 
    }
  }
  
  outBox <- ggplot(inData, aes(x=xCond, y=DV, shape = group, color = group)) +
    geom_point() +    # Use hollow circles
    geom_smooth(method=lm) + 
    facetg + 
    ggtitle(title) + 
    labs(x=xTitle, y=yTitle)
    
  
  return(outBox)
}

###A useful wrapper for scatterBY
#DV: the y variable, xCond = the x variable, 
#facetCond1 and facetCond2 are two variables for data subsetting.
my.scatterBy <- function(DV, xCond, group=NULL, facetCond1 = NULL, facetCond2 = NULL, 
                         title='A scatter plot', 
                         xTitle = 'IV', yTitle = 'DV', points='points', freeY=T)
{
  dt <- data.frame(DV = DV, xCond = xCond)
  if (!is.null(group))
  {
    dt$group <- group
  }
  if (!is.null(facetCond1))
  {
    dt$facetCond1 <- facetCond1
  }
  if (!is.null(facetCond2))
  {
    dt$facetCond2 <- facetCond2
  }
  
  return (scatterBy(dt, title, xTitle, yTitle, points, freeY))
}

###Create a simple bar plot to show frequencies of a variable
countBar <- function(var, include.na=F, labelsSize=0, yTitle = 'Count', xTitle = 'Var'){
  if (include.na)
  {
     var <- factor(var, exclude = NULL)
  }
  ret <- ggplot(data.frame(var), aes(x=var)) +  geom_bar(fill="#A6ACAF", colour="black") + 
    geom_text(aes(label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust=-1)+ 
    theme_bw(base_size = 23) + 
    labs(x=xTitle, y=yTitle)
  
  if (labelsSize > 0)
  {
    ret <- ret + geom_text(aes(y = 0, label = var), 
                           vjust=-0.9, colour="black", check_overlap = TRUE, size = labelsSize)
  }
  
  return(ret)
}

###countBar by groups.
countBarBy <- function(DV, group, include.na=F, labelsSize=0, yTitle = 'Count', xTitle = 'Var'){
  df <- data.frame(DV=DV, group=group)
  if (include.na)
  {
    df$DV <- factor(df$DV, exclude = NULL)
  }
  ret <- ggplot(df, aes(x=DV, fill=group)) +  geom_bar() + 
    geom_text(aes(label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust=-1)+ 
    theme_bw(base_size = 23) + 
    labs(x=xTitle, y=yTitle)
  
  if (labelsSize > 0)
  {
    ret <- ret + geom_text(aes(y = 0, label = DV), 
                           vjust=-0.9, colour="black", check_overlap = TRUE, size = labelsSize)
  }
  
  return(ret)
}

#facet.na, DV.na = whether to remove missing values for each of these
my.countBarBy2 <- function(DV, facet, yTitle = 'Count', xTitle = 'Var', facet.na=T, DV.na=T){
  df <- data.frame(DV=as.factor(DV), facet=facet)
  if (!facet.na)
  {
    df <- df[!is.na(df$facet),]
  }
  if (!DV.na)
  {
    df <- df[!is.na(df$DV),]
  }
  else
  {
    
  }
  dfr_prop <- df %>% 
    dplyr::count(DV, facet) %>%            # facet_by() & summarise(n = n()) are implicit
    group_by(facet) %>%
    mutate(prop = prop.table(n))
  
  gg_prop <- ggplot(data = data.frame()
                     , aes(x = DV, y = prop)) + 
    #facet_grid(.~facet) +
    #facet_wrap(~ facet, scales = "free")+
    facet_wrap(~ facet)+
    geom_text(aes(y = prop,label = n), vjust = 0, nudge_y = 0.01, nudge_x = 0.05)+
    geom_bar(stat = 'identity', position = 'dodge', alpha = 2/3) +  
    scale_y_continuous(labels = percent, limits = c(0, 1)) +    
    scale_fill_few('Medium', drop = FALSE) +              # keep levels, if data is filtered
    labs(x = xTitle, y = yTitle)  

  ret <- gg_prop %+%   # use %+% to add...
    dfr_prop    # a dataframe
  
  return(ret)
}

my.countBarBy2.many <- function (formula, dt, yTitle = 'Count', xTitle = 'Var', facet.na=T, DV.na=T)
{
  mlt <- mymelt(dt=dt, formula=formula)
  return(my.countBarBy2(DV=mlt$value, facet=mlt$variable, yTitle, xTitle, facet.na, DV.na))
}

forPlot <- function(xCond, group, DV)
{
  if (is.null(group))
  {
    dt <- data.frame(DV = DV, xCond = xCond, group = rep('', length(xCond)))
  }
  else
  {
    dt <- data.frame(DV = DV, xCond = xCond, group = group)
  }

  dt$xCond <- factor(dt$xCond)
  dt$group <- factor(dt$group)
  
  forplot <- ddply(.data=dt,
                   .variables = .(xCond, group),
                   summarise,
                   xsd = sd(DV, na.rm=TRUE),
                   xn = sum(!is.na(DV)),
                   xmin = min(DV, na.rm = TRUE),
                   xmax = max(DV, na.rm = TRUE),
                   xtop = max(DV, na.rm = TRUE)+0.3,
                   xmean =  mean(DV, na.rm=TRUE),
                   xmedian =  median(DV, na.rm=TRUE), 
                   xlower = mean(DV, na.rm=TRUE) - sd(DV, na.rm=TRUE), 
                   xupper = mean(DV, na.rm=TRUE) + sd(DV, na.rm=TRUE), 
                   xerrlow = mean(DV, na.rm=TRUE) - (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))), 
                   xerrup = mean(DV, na.rm=TRUE) + (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))) 
  )
  
  return(forplot)
}

##A function for creating a matrix of bar graphs for a 2 by 2 design.
#Use the wrapper function defined further below.
#The inData argument is the row-per-participant data. 
#That data must have these variable names: 
#DV - the variable for the y-axis
#xCond - one the factor for the x-axis
#group - second factor for the x-axis.
twoByTwoBar <- function(inData, title='A 2 by 2 graph', 
                            xTitle = 'IV', yTitle = 'DV', points='jitter', na.rm=T) {

  inData$xCond <- factor(inData$xCond)
  inData$group <- factor(inData$group)
  
  forplot <- ddply(.data=inData,
                   .variables = .(xCond, group),
                   summarise,
                   xsd = sd(DV, na.rm=TRUE),
                   xn = sum(!is.na(DV)),
                   xmin = min(DV, na.rm = TRUE),
                   xmax = max(DV, na.rm = TRUE),
                   xtop = max(DV, na.rm = TRUE)+0.3,
                   xmean =  mean(DV, na.rm=TRUE),
                   xmedian =  median(DV, na.rm=TRUE), 
                   xlower = mean(DV, na.rm=TRUE) - sd(DV, na.rm=TRUE), 
                   xupper = mean(DV, na.rm=TRUE) + sd(DV, na.rm=TRUE), 
                   xerrlow = mean(DV, na.rm=TRUE) - (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))), 
                   xerrup = mean(DV, na.rm=TRUE) + (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))) 
  )
  
  dodge <- position_dodge(width=0.9)
  outBox <- ggplot(data=forplot, aes(x=factor(xCond), fill=factor(group), y= xmean))+
    geom_bar(position=dodge, stat='identity')+
    geom_errorbar(aes(ymin = xerrlow, ymax = xerrup), width = 0.25, position=dodge) + 
    ggtitle(title) + 
    labs(x=xTitle, y=yTitle)+
    theme_bw(base_size = 19)+
    scale_size_continuous(range = c(1,3))+
    #geom_text(aes(y = xmax,label = paste('n =', xn, sep=' ')), position=dodge)
    geom_text(aes(y = xerrup*1.1,label = paste('n =', xn, sep=' ')), position=dodge)
  
  return(outBox)
}

twoByTwoBarH <- function(DV, xCond, group=NULL, title='A 2 by 2 graph', 
                        xTitle = 'IV', yTitle = 'DV', points='jitter', na.rm=T) {
  if (is.null(group))
  {
    dt <- data.frame(DV = DV, xCond = xCond, group = rep('', length(xCond)))
  }
  else
  {
    dt <- data.frame(DV = DV, xCond = xCond, group = group)
  }
    
  return (twoByTwoBar(dt, title, xTitle, yTitle, points, na.rm))
}
myg.twoByTwoBar <- twoByTwoBarH

  
##A function for creating a matrix of bar graphs for a 2 by 2 design.
#The inData argument is the row-per-participant data. 
#That data must have these variable names: 
#DV - the variable for the y-axis
#xCond - one the factor for the x-axis
#group - second factor for the x-axis.
#facetCond1 - separate the graphs by this variable.
#facetCond2 - separate the graphs also by this variable.
twoByTwoBars <- function(inData, title='A 2 by 2 graph', 
                        xTitle = 'IV', yTitle = 'DV', points='jitter', na.rm=T) {
  if (is.null(inData$facetCond1))
  {
    facetGrid <- NULL
    inData$facetCond1 = '1'
    inData$facetCond2 = '1'
  }
  else if (is.null(inData$facetCond2))
  {
    facetGrid <- facet_grid(.~facetCond1)
    inData$facetCond2 = '1'
  }
  else
  {
    facetGrid <- facet_grid(facetCond1 ~ facetCond2)
  }

   forplot <- ddply(.data=inData,
                   .variables = .(facetCond1, facetCond2, xCond, group),
                   summarise,
                   xsd = sd(DV, na.rm=TRUE),
                   xn = sum(!is.na(DV)),
                   xmin = min(DV, na.rm = TRUE),
                   xmax = max(DV, na.rm = TRUE),
                   xtop = max(DV, na.rm = TRUE)+0.3,
                   xmean =  mean(DV, na.rm=TRUE),
                   xmedian =  median(DV, na.rm=TRUE), 
                   xlower = mean(DV, na.rm=TRUE) - sd(DV, na.rm=TRUE), 
                   xupper = mean(DV, na.rm=TRUE) + sd(DV, na.rm=TRUE), 
                   xerrlow = mean(DV, na.rm=TRUE) - (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))), 
                   xerrup = mean(DV, na.rm=TRUE) + (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))) 
  )
    
  if (na.rm == TRUE)
  {
    forplot <- forplot[!is.na(forplot$xCond) & !is.na(forplot$facetCond1) & !is.na(forplot$facetCond2),]
  }
    
  
  dodge <- position_dodge(width=0.9)
  outBox <- ggplot(data=forplot, aes(x=factor(xCond), fill=factor(group), y= xmean))+
    geom_bar(position=dodge, stat='identity')+
    geom_errorbar(aes(ymin = xerrlow, ymax = xerrup), width = 0.25, position=dodge) + 
    ggtitle(title) + 
    labs(x=xTitle, y=yTitle)+
    theme_bw(base_size = 19)+
    facetGrid+
    scale_size_continuous(range = c(1,3))+
    geom_text(aes(y = xmax,label = paste('n =', xn, sep=' ')), position=dodge)
  
  return(outBox)
}

#Data should have:
#DV: The y-axis variable
#xCond = x-axis variable
#group = how to separate lines.
#facetCond = a variable for separating the graphs 
#facetCond2 = another variable for separating the graphs
lineGraph <- function(inData, title='A line graph', 
                        xTitle = 'IV', yTitle = 'DV', na.rm=T, showN=T, 
                      usePoints=T, pointsSize=4) {
  
  
  if (is.null(inData$facetCond))
  {
    facetGrid <- NULL
    inData$facetCond = '1'
    inData$facetCond2 = '1'
  }
  else if (is.null(inData$facetCond2))
  {
    facetGrid <- facet_grid(.~facetCond)
    inData$facetCond2 = '1'
  }
  else
  {
    facetGrid <- facet_grid(facetCond ~ facetCond2)
  }
  
  geomPoints <- NULL
  if (usePoints)
  {
    geomPoints <- geom_point(aes(y=xmean),color='grey', size=pointsSize)
  }
  
  ngeom <- NULL
  if(showN)
  {
    ngeom <- geom_text_repel(aes(y = xupper,label = paste('n =', xn, sep=' ')))
  }
  
  forplot <- ddply(.data=inData,
                   .variables = .(facetCond, facetCond2, group, xCond),
                   summarise,
                   xsd = sd(DV, na.rm=TRUE),
                   xn = sum(!is.na(DV)),
                   xmin = min(DV, na.rm = TRUE),
                   xmax = max(DV, na.rm = TRUE),
                   xtop = max(DV, na.rm = TRUE)+0.3,
                   xmean =  mean(DV, na.rm=TRUE),
                   xmedian =  median(DV, na.rm=TRUE), 
                   xlower = mean(DV, na.rm=TRUE) - sd(DV, na.rm=TRUE), 
                   xupper = mean(DV, na.rm=TRUE) + sd(DV, na.rm=TRUE), 
                   xerrlow = mean(DV, na.rm=TRUE) - (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))), 
                   xerrup = mean(DV, na.rm=TRUE) + (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV))))) 
    

    if (na.rm == TRUE)
    {
      forplot <- forplot[!is.na(forplot$xCond) & !is.na(forplot$facetCond) & !is.na(forplot$facetCond2),]
    }
    
  outBox <- ggplot(data=forplot, aes(x=factor(xCond), y=xmean, group=group, colour=group, shape=group, linetype=group))+
    geom_line(size=1.5) + 
    geomPoints +
    geom_errorbar(aes(ymin = xerrlow, ymax = xerrup), width = 0.5, linetype = 'longdash', colour = 'black') + 
    facetGrid + 
    theme_bw(base_size = 19)+
    ggtitle(title) + 
    labs(x=xTitle, y=yTitle)+
    scale_size_continuous(range = c(1,3))+
    ngeom
  
  
  return(outBox)
}
my.lineGraph <- function(DV, xCond, group=NULL, facet1=NULL, facet2=NULL, 
                         title='', xTitle = '', yTitle = '', 
                         na.rm=T, showN=T, usePoints=T, pointsSize=4)
{
  dt <- data.frame(DV, xCond)

  if (is.null(group))
  {
    dt$group <- ''
  }
  else
  {
    dt$group <- group
  }
  
  dt$facetCond <- facet1
  dt$facetCond2 <- facet2
  
  return(lineGraph(dt, title, xTitle, yTitle, na.rm, showN, usePoints, pointsSize))
}

#Create a few simple line graphs, with no error bars.
lineGraphs <- function(DV, xVar, groupVar, facet1=NULL, facet2=NULL, 
                       xTitle='', yTitle='', title='')
{
  dt <- data.frame(DV, xVar, groupVar)
  if (is.null(facet1))
  {
    facetGrid <- NULL
  }
  else if (is.null(facet2))
  {
    facetGrid <- facet_grid(.~facet1)
    dt$facet1 <- facet1
  }
  else
  {
    facetGrid <- facet_grid(facet1 ~ facet2)
    dt$facet1 <- facet1
    dt$facet2 <- facet2
  }
  
  gg <- ggplot(data=dt, aes(x=factor(xVar), group=groupVar, y=DV, color=groupVar))+
    geom_point() + facetGrid  + geom_line() + 
    ggtitle(title) + 
    labs(x=xTitle, y=yTitle)
    return(gg)
}
##A function for creating a matrix of boxplot graphs for a 2 by 2 design.
#The inData argument is the row-per-participant data. 
#That data must have these variable names: 
#DV - the variable for the y-axis
#xCond - the factor for the x-axis
#facetCond - the factor for separating the graphs.
twoByTwoBoxPlot <- function(inData, title='', 
                            xTitle = '', yTitle = '', points='jitter', na.rm=T, refactor=T) {
  if (refactor == T)
  {
    inData$xCond <- factor(inData$xCond)
    inData$facetCond <- factor(inData$facetCond)
  }
  
  forplot <- ddply(.data=inData,
                    .variables = .(xCond, facetCond),
                    summarise,
                    xsd = sd(DV, na.rm=TRUE),
                    xn = sum(!is.na(DV)),
                    xmin = min(DV, na.rm = TRUE),
                    xmax = max(DV, na.rm = TRUE),
                    xtop = max(DV, na.rm = TRUE)+0.3,
                    xmean =  mean(DV, na.rm=TRUE),
                    xmedian =  median(DV, na.rm=TRUE), 
                    #xlower = mean(DV, na.rm=TRUE) - sd(DV, na.rm=TRUE), 
                    #xupper = mean(DV, na.rm=TRUE) + sd(DV, na.rm=TRUE), 
                    xlower = quantile(DV, 0.25, na.rm=TRUE), 
                    xupper = quantile(DV, 0.75, na.rm=TRUE), 
                    xerrlow = mean(DV, na.rm=TRUE) - (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))), 
                    xerrup = mean(DV, na.rm=TRUE) + (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))) 
  )
  
  if (na.rm == TRUE)
  {
    forplot <- forplot[!is.na(forplot$xCond) & !is.na(forplot$facetCond),]
  }
  
  geomPts <- NULL
  if (points == 'jitter')
  {
    geomPts <- geom_jitter(data = inData, aes(y=DV, x=factor(xCond), facet=factor(facetCond)), position=position_jitter(width=.8, height=0), shape=4)
  }
  else if (points == 'points')
  {
    geomPts <- geom_count(data = inData, aes(y=DV, x=factor(xCond), facet=factor(facetCond)), shape=4)
  }
  
  outBox <- ggplot(data=forplot, aes(x=factor(xCond), facet=factor(facetCond)))+
    geom_boxplot(aes(lower = xlower, upper = xupper, middle = xmedian, ymin = xmin, ymax = xmax), stat = "identity") + 
    geom_errorbar(aes(ymin = xerrlow, ymax = xerrup), width = 0.9, linetype = 'longdash', colour = 'black', size=1.1) + 
    facet_wrap(~ facetCond)+
    #facet_grid(.~facetCond)+
    #facet_grid(facetCond ~ .)+
    geom_point(aes(y=xmean),color='grey', size=4, shape=17) + 
    ggtitle(title) + 
    labs(x=xTitle, y=yTitle)+
    #xlab("Sorting condition") + 
    #ylab("IAT Score") + 
    
    #scale_fill_brewer(palette="Greys")+
    theme_bw(base_size = 19)+
    geomPts+
    #ifelse(jitter,
    #  geom_jitter(data = inData, aes(y=DV, x=factor(xCond), facet=factor(facetCond)), position=position_jitter(width=.4, height=0))+
    #  geom_count(data = inData, aes(y=DV, x=factor(xCond), facet=factor(facetCond)))+
      #)+
    scale_size_continuous(range = c(1,3))+
    geom_text(aes(y = xmax,label = paste('n =', xn, sep=' ')), vjust = 0, nudge_y = 0.05, nudge_x = 0.05)
  
  return(outBox)
}
myg.twoByTwoBoxPlot <- function(DV, xCond, facetCond, ...)
{
  dt <- data.frame(DV = DV, xCond = xCond, facetCond = facetCond)
  return(twoByTwoBoxPlot(dt, ...))
}

##A function for creating a matrix of boxplot graphs for a 2 by 2 by 2 design.
#***Use the wrapper function defined further below***
#The inData argument is the row-per-participant data. 
#That data must have these variable names: 
#DV - the variable for the y-axis
#xCond - the factor for the x-axis
#facetCond1 - the first factor for separating the graphs.
#facetCond2 - the second factor for separating the graphs.
twoByTwoBoxPlots <- function(inData, title='', 
                            xTitle = '', yTitle = '', points='none') {
  inData$xCond <- factor(inData$xCond)
  inData$facetCond1 <- factor(inData$facetCond1)
  inData$facetCond2 <- factor(inData$facetCond2)
  
  forplot <- ddply(.data=inData,
                   .variables = .(xCond, facetCond1, facetCond2),
                   summarise,
                   xsd = sd(DV, na.rm=TRUE),
                   xn = sum(!is.na(DV)),
                   xmin = min(DV, na.rm = TRUE),
                   xmax = max(DV, na.rm = TRUE),
                   xtop = max(DV, na.rm = TRUE)+0.3,
                   xmean =  mean(DV, na.rm=TRUE),
                   xmedian =  median(DV, na.rm=TRUE), 
                   #xlower = mean(DV, na.rm=TRUE) - sd(DV, na.rm=TRUE), 
                   #xupper = mean(DV, na.rm=TRUE) + sd(DV, na.rm=TRUE), 
                   xlower = quantile(DV, 0.25, na.rm=TRUE), 
                   xupper = quantile(DV, 0.75, na.rm=TRUE), 
                   xerrlow = mean(DV, na.rm=TRUE) - (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))), 
                   xerrup = mean(DV, na.rm=TRUE) + (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))) 
  )

  geomPts <- NULL
  if (points == 'jitter')
  {
    geomPts <- geom_jitter(data = inData, aes(y=DV, x=factor(xCond)), position=position_jitter(width=.8, height=0), shape=4)
  }
  else if (points == 'points')
  {
    geomPts <- geom_count(data = inData, aes(y=DV, x=factor(xCond)), shape=4)
  }
  
  outBox <- ggplot(data=forplot, aes(x=factor(xCond)))+
    geom_boxplot(aes(lower = xlower, upper = xupper, middle = xmedian, ymin = xmin, ymax = xmax), stat = "identity") + 
    geom_errorbar(aes(ymin = xerrlow, ymax = xerrup), width = 0.8, linetype = 'longdash', colour = 'black', size=1.1) + 
    facet_grid(facetCond1 ~ facetCond2)+
    geom_point(aes(y=xmean),color='grey', size=4, shape=17) + 
    ggtitle(title) + 
    geomPts +
    labs(x=xTitle, y=yTitle)+
    #scale_fill_brewer(palette="Greys")+
    theme_bw(base_size = 19)+
    scale_size_continuous(range = c(1,3))+
    geom_text(aes(y = xmax,label = paste('n =', xn, sep=' ')), vjust = 0, nudge_y = 0.05, nudge_x = 0.05)
  
  return(outBox)
}
#Wrapper for twoByTwoBoxPlots
myg.twoByTwoBoxPlots <- function(DV, xCond, facetCond1, facetCond2, 
                                 title='', 
                                 xTitle = '', yTitle = '', points='none')
{
  dt <- data.frame(DV = DV, xCond = xCond, facetCond1 = facetCond1, facetCond2 = facetCond2)
  return(twoByTwoBoxPlots(dt, title, xTitle, yTitle, points))
}

##A function for creating one boxplot (one factor design).
#***Use the wrapper function defined below***
#The inData argument is the row-per-participant data. 
#That data must have these variable names: 
#DV - the variable for the y-axis
#xCond - the factor for the x-axis
oneFactorBoxPlot <- function(inData, title='', 
                            xTitle = '', yTitle = '', jitter='jitter', jitterWidth=.8, 
                            pointShape=4, pointSize=3, meanColor='grey', meanSize=4, meanShape=17) {
  inData$xCond <- factor(inData$xCond)
  
  forplot <- ddply(.data=inData,
                   .variables = .(xCond),
                   summarise,
                   xsd = sd(DV, na.rm=TRUE),
                   xn = sum(!is.na(DV)),
                   xmin = min(DV, na.rm = TRUE),
                   xmax = max(DV, na.rm = TRUE),
                   xtop = max(DV, na.rm = TRUE)+0.3,
                   xmean =  mean(DV, na.rm=TRUE),
                   xmedian =  median(DV, na.rm=TRUE), 
                   #xlower = mean(DV, na.rm=TRUE) - sd(DV, na.rm=TRUE), 
                   #xupper = mean(DV, na.rm=TRUE) + sd(DV, na.rm=TRUE), 
                   xlower = quantile(DV, 0.25, na.rm=TRUE), 
                   xupper = quantile(DV, 0.75, na.rm=TRUE), 
                   xerrlow = mean(DV, na.rm=TRUE) - (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))), 
                   xerrup = mean(DV, na.rm=TRUE) + (qnorm(0.975)*sd(DV, na.rm=TRUE)/sqrt(length(!is.na(DV)))) 
  )

  geomPts <- NULL
  if (jitter == 'points')
  {
    geomPts <- geom_count(data = inData, aes(y=DV, x=factor(xCond), color=factor(xCond)))
  }
  else if (jitter == 'jitter')
  {
    geomPts <- geom_jitter(data = inData, aes(y=DV, x=factor(xCond), color=factor(xCond)), position=position_jitter(width=jitterWidth, height=0), 
                           shape=pointShape, size=pointSize)
  }
  
  outBox <- ggplot(data=forplot, aes(x=factor(xCond)))+
    geom_boxplot(aes(lower = xlower, upper = xupper, middle = xmedian, ymin = xmin, ymax = xmax), stat = "identity") +
    geom_errorbar(aes(ymin = xerrlow, ymax = xerrup), width = 0.5, linetype = 'longdash', colour = 'black', size=1.1) + 
    #facet_grid(.~facetCond)+
    geom_point(aes(y=xmean), color=meanColor, size=meanSize, shape=meanShape) + 
    ggtitle(title) + 
    labs(x=xTitle, y=yTitle, color=xTitle)+
    #xlab("Sorting condition") + 
    #ylab("IAT Score") + 
    #scale_fill_brewer(palette="Greys")+
    theme_bw(base_size = 19)+
    geomPts+
    scale_size_continuous(range = c(1,3))+
    geom_text(aes(y = xmax,label = paste('n =', xn, sep=' ')), vjust = 0, nudge_y = 0.05, nudge_x = 0.05)
  
  return(outBox)
}
#wrapper for oneFactorBoxPlot
myg.oneFactorBoxPlot <- function(DV, xCond, 
                                 title='', 
                                 xTitle = '', yTitle = '', jitter='jitter', jitterWidth=.8, 
                                 pointShape=4, pointSize=3, meanColor='grey', meanSize=4, meanShape=17)
{
  dt <- data.frame(DV = DV, xCond = xCond)
  return(oneFactorBoxPlot(dt, title, xTitle, yTitle, jitter, jitterWidth, pointShape, pointSize, meanColor, meanSize, meanShape))
}

my.bubble <- function(DV, xCond, sizeVar=NULL, fillVar=NULL, xTitle='', yTitle='', title='')
{
  dt <- data.frame(xCond = xCond, DV = DV)
  if(is.null(sizeVar))
  {
    dt <- setDF(setDT(dt)[, sizeVar := .N, by = .(xCond, DV)])
  }
  else
  {
    dt$sizeVar = sizeVar
  }
  
  if(is.null(fillVar))
  {
    dt$fillVar <- dt$sizeVar
  }
  else
  {
    dt$fillVar <- fillVar
  }
  
  outg <- ggplot(data=dt, aes(x = xCond, y = DV, size = sizeVar, fill = fillVar)) +
    geom_point(shape = 21) +
    ggtitle(title) +
    labs(x = xTitle, y = yTitle)
  
  return(outg)
}

##A function for creating many boxplot (one factor design) - one for each variable in the data.
#The inData argument is a "melted" dataset. It must have the variables: 
#value - the variable for the y-axis
#facetCond - the name of the DV variable (usually outputs as the column 'variable' from the melt function)
#xCond - the variable for the x-axis.
manyBoxes <- function(inData, title='', xTitle = '', yTitle = '', points='jitter', freeY=TRUE, meanSize=4, nSize=2)
{
  inData$xCond <- factor(inData$xCond)
  inData$facetCond <- factor(inData$facetCond)
  
  geomPts <- NULL
  if (points=='points')
  {
    geomPts <- geom_count(shape=1)
  }
  else if (points=='jitter')
  {
    geomPts <- geom_jitter(position=position_jitter(width=.8, height=0), shape=4)
  }

  facetW <- facet_wrap( ~ facetCond) 
  if (freeY)
  {
    facetW <- facet_wrap( ~ facetCond, scales="free_y") 
  }
  
  outG <- ggplot(data = inData, aes(x=xCond, y=value, facet=facetCond)) + 
    #geom_boxplot(aes(fill=sortCond))+
    scale_fill_brewer(palette="Greys")+
    stat_summary(fun.data = min.median.sd.max, geom = "boxplot") + 
    stat_summary(fun.data = getCI, geom = "errorbar", width=0.5, linetype = 'longdash', colour = 'black') + 
    stat_summary(fun.y="mean", geom="point", shape=17, size=meanSize, color="grey") +
    
    #So far, I couldn't figure out how to add the n of each bar above the bar.
    stat_summary(fun.data = n_fun, geom = "text", size=nSize, color="black") + 
    #geom_jitter(position=position_jitter(width=.8, height=0), shape=1)+
    #geom_text(aes(y = "max",label = getN), vjust = 0, nudge_y = 0.05, nudge_x = 0.05)+
    geomPts+
    ggtitle(title) + 
    labs(x=xTitle, y=yTitle)+
    facetW 
  
  return(outG)
}

my.manySingleBoxes.Facets <- function(formula, dt, title='', xTitle = '', yTitle = '', points='jitter', freeY=TRUE, meanSize=4, nSize=2)
{
  mlt <- mymelt(dt=dt, formula=formula)
  mlt$xCond = 'all'
  mlt$facetCond = mlt$variable
  return(manyBoxes(inData = mlt, title=title, xTitle = xTitle, yTitle = yTitle, points=points, freeY=freeY, meanSize=meanSize, nSize=nSize))
}


my.manySingleBoxes <- function(formula, dt, title='', 
                               xTitle = '', yTitle = '', jitter='jitter', jitterWidth=.8, 
                               pointShape=4, pointSize=3, meanColor='grey', meanSize=4, meanShape=17)
{
  mlt <- mymelt(dt=dt, formula=formula)
  return(myg.oneFactorBoxPlot(DV=mlt$value, xCond = mlt$variable, title=title, 
                              xTitle = xTitle, yTitle = yTitle, jitter=jitter, jitterWidth=jitterWidth, 
                              pointShape=pointShape, pointSize=pointSize, meanColor=meanColor, meanSize=meanSize, meanShape=meanShape))
}

###Create a violin graph
#DV = the y-variable, xFactor = the x-variable, 
#fillFactor = a factor for creating violins with a different fill
#facet1 and facet2: factors for dividing the violins, if need
my.violin <- function(DV, xFactor, fillFactor=NULL, shapeFactor=NULL, colorFactor=NULL,
                      scales='fixed', yTitle='', xTitle='', title='', fillTitle='',
                      fillColors=c("#EEEEEE", "#CCCCCC", "#999999", '#676767', '#343434'), 
                      facet1=NULL, facet2=NULL, points=NULL, x.na=T, 
                      jitterWidth=.4, jitterHeight=0, no.violin=F, dotSize=0.3)
{
    dt <- data.frame(DV=DV, xFactor=as.factor(xFactor))
    if (!x.na)
    {
      #browser()
      dt <- dt[!is.na(dt$xFactor),]
      dt$xFactor <- factor(dt$xFactor)
    }
    dt$facet1 <- facet1
    dt$facet2 <- facet2
    dt$fillFactor <- fillFactor
    my.aes <- aes(x=xFactor, fill=fillFactor, y= DV)
    theViolin <- geom_violin(aes(fill=fillFactor), position=position_dodge(.9))
    theMean <- stat_summary(fun.y="mean", colour="black", size=4, shape=1, geom="point", aes(group=fillFactor), position=position_dodge(.9))
    theMedian <- stat_summary(fun.y="median", colour="wheat4", size=3, shape=0, geom="point", aes(group=fillFactor), position=position_dodge(.9))
    theError <- stat_summary(fun.data="mean_cl_normal", geom="errorbar", aes(group=fillFactor), position=position_dodge(.9))
    thePoint <- geom_count(aes(colour=fillFactor))
    if (is.null(fillFactor))
    {
      my.aes <- aes(x=xFactor, y= DV)
      theViolin <- geom_violin()
      theMean <- stat_summary(fun.y="mean", colour="black", geom="point")
      theError <- stat_summary(fun.data="mean_cl_normal", geom="errorbar")
    }

    if (is.null(facet1))
    {
      my.facet.wrap <- NULL
    }
    else if (is.null(facet2))
    {
      my.facet.wrap <- facet_wrap(~ facet1,scales=scales)
    }
    else 
    {
      my.facet.wrap <- facet_wrap(facet1 ~ facet2,scales=scales)
    }
    
    colFactor = dt$xFactor
    if (!is.null(colorFactor))
    {
      colFactor <- colorFactor
    }
    else if (!is.null(fillFactor))
    {
      colFactor <- fillFactor
    }
    
    theJitter <- geom_jitter(position=position_jitter(width=jitterWidth, height=jitterHeight), 
                             shape=4, aes(color = colFactor))
    if (!is.null(shapeFactor))
    {
      theJitter <- geom_jitter(position=position_jitter(width=jitterWidth, height=jitterHeight), 
                               aes(color = colFactor, shape=shapeFactor))
    }
    
    theDots <- geom_dotplot(binaxis = 'y', dotsize = dotSize, stackdir = 'center')
    
    if (is.null(points))
    {
      geomPts <- NULL
    }
    else if (points == 'jitter')
    {
      geomPts <- theJitter
    }
    else if (points == 'points')
    {
      geomPts <- thePoint
    }
    else if (points == 'dotplot')
    {
      geomPts <- theDots
    }
    
    if (no.violin)
    {
      theViolin = NULL
    }
    

    #browser()
    outG <- ggplot(data=dt, my.aes)+
      theViolin+
      theMean+
      theMedian+
      theError+
      my.facet.wrap+
      scale_fill_manual(values=fillColors)+
      ggtitle(title)+
      geomPts+
      labs(x=xTitle, y=yTitle, fill=fillTitle, colour=fillTitle)+
      theme_bw(base_size = 19) + 
      scale_color_hue(l=50, c=55)
    
    return(outG)
}

###Show density plots for DV by xFactor
##Best to use with melted data in which DV is the value and xFactor is the name of variable.
my.densities <- function(DV, xFactor)
{
  df <- data.frame(DV=DV, xFactor=xFactor)
  ggg <- ggplot(df, aes(DV)) + geom_density() + 
    facet_wrap(~ xFactor, scales = "free") 
  return(ggg)
}

plot.density <- function(dv, probs=c(.1, .2, .3, .4, .5, .6, .7, .8, .9, 1))
{
  dv <- dv[!is.na(dv)]
  dens <- density(dv)
  df <- data.frame(xx=dens$x, yy=dens$y)
  quantiles <- quantile(dv, prob=probs)
  df$quant <- factor(findInterval(df$xx,quantiles))
  gg <- ggplot(df, aes(xx,yy)) + 
    geom_line() + 
    geom_ribbon(aes(ymin=0, ymax=yy, fill=quant)) + 
    scale_x_continuous(breaks=quantiles) + 
    scale_fill_brewer(guide="none")
  return(gg)
}

plot.percentile <- function(dv, probs=seq(0.1, 0.9, by= 0.1), labelsSize = 3, yTitle = '', xTitle = '')
{
  dv <- dv[!is.na(dv)]
  pcs <- quantile(dv, probs = seq(0, 1, by= 0.01))
  df <- data.frame(x=c(0:100), p = round(pcs,2))
  quantiles <- quantile(df$p, prob=probs)
  df$quants <- factor(findInterval(df$p,quantiles)/10)
  
  xn <- length(dv)
  xmax <- max(dv)
  xmin <- min(dv)
  
  gg <- ggplot(df, aes(x,p, label=p)) + 
    geom_line() + 
    geom_ribbon(aes(ymin=0, ymax=p, fill=quants)) + 
    scale_x_continuous(breaks=seq(0,100,by=5)) + 
    scale_fill_manual(values=c("#A6ACAF", "#283747", "#B03A2E", "#6C3483", "#1F618D", "#117A65", "#B9770E", "#A04000", "#641E16", "#4A235A", "#1B4F72", "#0B5345", "#145A32", "#7E5109", "#626567", "#17202A")) + 
    labs(x=xTitle, y=yTitle) + 
    annotate("text", label = paste0("n = ", xn), x = 10, y = (xmax-xmin)*0.85, size = 5, colour = "black")
  
    
  if (labelsSize>0)
  {
    gg <- gg + geom_text(vjust=-0.9, colour="black", check_overlap = TRUE, size = labelsSize)
  }
  return(gg)
}

###Show a table in the HTML viewer
#Replaces \n with <br/> == perfect for showing my correlation matrix created by my correlation functions
my.htmlTable <- function(ttt)
{
  htt <- as.data.frame(lapply(ttt, FUN= function(x)gsub('\n','<br/>',x)))
  #YBYB: Added after a new problem with presentation of the table. Might need to remove in the future
  hhtt <- htmlTable::htmlTable(htt)
  class(hhtt) <- c("htmlTable", "html", "character")
  return (hhtt)
}

####Examines correlation with and without outliers using a few different methods.
#****After running: click back on the plot area to see a few interesting plots
###inSID = subject id
##x and y are the variable to correlate
my.bivariate <- function(inSID, x, y)
{
  retList <- list()
  forc <- data.frame(session_id = inSID, x = x, y = y)
  retList$corMat.orig <- cornp(forc[,c('x', 'y')])
  retList$scatter.orig <- my.scatterBy(x, y)
  
  library(mvoutlier)
  retList$corrPlot <- corr.plot(x, y)
  
  library(aplpack)
  retList$bagPlot <- bagplot(cbind(x,y),pch=16,cex=2)
  
  mod <- lm(y ~ x, data=forc)
  forc$cooksd <- cooks.distance(mod)
  mcooksd <- mean(forc$cooksd, na.rm=T)
  hcooksd <- 4*mcooksd
  
  plot(forc$cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
  abline(h = hcooksd, col="red")  # add cutoff line
  text(x=1:length(forc$cooksd)+1, y=forc$cooksd, labels=ifelse(forc$cooksd>hcooksd,forc$session_id,""), col="red") 
  forc$influential <- forc$cooksd >= hcooksd
  
  retList$corMat.postCook <- cornp(forc[!forc$influential,c('x', 'y')])
  retList$dt <- forc
  return(retList)
}

my.densMany <- function(dt, sid='session_id', inVars=NULL, formula=NULL)
{
  if (is.null(inVars))
  {
    sdvs <- as.character(formula)[2]
    dvs <- strsplit(sdvs," \\+ ")[[1]]
    vars <- c(sid, dvs)
  }
  else
  {
    vars = c(sid, inVars)
  }
  
  df <- melt(allds[,vars], id.vars = sid, variable.name = 'vars')
  ggg <- ggplot(df, aes(value)) + geom_density() + facet_wrap(~ vars, scales = "free")
  return(ggg)
}

my.alphaBy <- function(inData, formula=NULL, alphaVars=NULL, byVars=NULL)
{
  if (is.null(alphaVars))
  {
    alphaVars <- as.character(formula)[2]
  }
  if (is.null(byVars))
  {
    byVars <- as.character(formula)[3]
  }
  
  
  get.alpha <- function(xx)
  {
    aaa <- psych::alpha(as.data.frame(xx[,alphaVars]), check.keys=FALSE)
    return(aaa$total$std.alpha)
  }
  aby <- ddply(inData, as.quoted(byVars), get.alpha)
  return(aby)
}

missingNames <- function(names.formula, the.names=NULL, dt)
{
  if (is.null(the.names))
  {
    sdvs <- as.character(names.formula)[2]
    the.names <- strsplit(sdvs," \\+ ")[[1]]
  }
  return (the.names[!(the.names %in% names(dt))])
}

round.df <- function(x, digits=3) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, class) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  return(x)
}

mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}


iat.scr.dbl <- function(dt, v_error=2, v_extreme=2, v_std=1) 
  #v_error=2 means recode error latency to m+600
  {
  dt1 <- dt
  dt1$blockName <- ifelse(dt1$pairing=='cong', 'B3', 'B6')
  dt2 <- dt 
  dt2$blockName <- ifelse(dt2$pairing=='cong', 'B4', 'B7')
  iatscr <- cleanIAT(rbind(dt1,dt2), block_name="blockName", 
                     trial_blocks = c("B3", "B4", "B6", "B7"), 
                     session_id="session_id", 
                            trial_latency="trial_latency",
                            trial_error = "trial_error", 
                            v_error=v_error, v_extreme=v_extreme, v_std=v_std) 
  return(iatscr)
}

#Will return a scored data with session_id, subexcl, iat1, iat2... by parts.
iat.parcel <- function(dt, parts=6, prefix='iat.part.', subVar='session_id', trialNumVar='trial_number',
                       rtVar='trial_latency', pairingVar='pairing', errVar='trial_error', parseBy='seq') {
  #create the data.frame
  the.d <- dt[,c(subVar, trialNumVar, rtVar, pairingVar, errVar)]
  setnames(the.d, old = c(subVar, trialNumVar, rtVar, pairingVar, errVar), 
           new = c('session_id', 'trial_number', 'trial_latency', 'pairing', 'trial_error'))
  #Sort by trial_number  
  the.d <- the.d[order(the.d$session_id, the.d$pairing, the.d$trial_number),]
  
  setDT(the.d)
  #Count number of trials for each session within each pairing condition, 
  #and also create a running number within each session-pairing subset.
  the.d <- the.d[ , `:=`( trialCount = .N , trialNum = 1:.N ), by = c("session_id", 'pairing') ]
  #For each trial, set the part it belongs to.
  if (parseBy == 'seq')
  {#Parse sequentially
    the.d$part <- floor((the.d$trialNum-1) / (the.d$trialCount / parts))+1
  }
  else
  {#Parse by Modulo 
    the.d$part <- (the.d$trialNum %% parts) + 1
  }
  #compute an IAT score for each part [the IAT scoring expects two parts, so the iat.scr.dbl function duplicated the data]
  dtm <- the.d[,iat.scr.dbl(.SD), by='part']
  setDF(the.d)

  #Within each part in the dtm (the dataset with scores), IAT1 and IAT2 were computed from the same data, so we need only IAT1
  outdt <- dcast(dtm, session_id ~ part, value.var = c('IAT1'))
  names(outdt)[2:(parts+1)] <- paste0(prefix, names(outdt)[2:(parts+1)])
  return(setDF(outdt))
}

my.mosaic <- function(var1, var2)
{
  library("grid"); library("vcd")
  dt <- data.frame(varA=var1, varB=var2)
  # prepare data for plot as a "Structured Contingency Table"
  form <- structable(varA ~ varB, dt)
  
  #Note that the arguments to the prop.table are the opposite order of the arguments to structable
  labeling_values1 <- table(dt$varB, dt$varA)
  labeling_values2 <- round(prop.table(table(dt$varB, dt$varA)),2)*100
  labeling_values <- matrix(paste0(labeling_values1,' (',labeling_values2,'%)'), ncol = ncol(labeling_values1))
  
  # basic plot
  dimnames(labeling_values)<-dimnames(form)
  mosaic(form, pop=FALSE)
  labeling_cells(text = as.table(labeling_values), margin = FALSE)(as.table(form))
}

my.barplot <- function(y, x, inLegend=c('1st','2nd','3rd','4th'), 
                       inColors=c("grey","skyblue","royalblue","darkblue"), 
                       inTitle='', axisSize=1.5, namesSize=1.5, pctSize=1.5)
{
  bartable = table(y, x)  ## get the cross tab
  ppp = round(prop.table(bartable,2)*100,1) # percent
  rotulo=paste(ppp,'%')
  barras <- barplot(bartable, beside = TRUE, legend = inLegend,
                    col=inColors, main=inTitle, cex.axis=axisSize, cex.names=namesSize)## plot 
  text(barras, 0, rotulo,cex=pctSize,pos=3, col ="#ff0000")
}

#create nice anova output (from Tzipi)
my.niceAnova <- function(aovOut)
{
  output <- nice(aovOut, es = "pes") %>%
    kable() %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F, position = "left")
  
  return(output)
}
