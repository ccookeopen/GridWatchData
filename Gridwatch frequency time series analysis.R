remove(list = ls())
library("lubridate")
library("plyr")
library("moments")
library("fitdistrplus")
library("sn")
library("rmutil")
library("gnorm")
library("zoo")
library("extraDistr")
library("fGarch")

pdf("grid-frequency-analysis.pdf")  #Save as a pdf file

# Read in Gridwatch Data

GridwatchData <- read.csv(file="gridwatch-fulldataset.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
attach(GridwatchData)

# Remove the zero-frequency data - should just replace with NA but this what I did first
GridwatchDataFreqNonZero <- GridwatchData[which(GridwatchData$frequency > 0.0),]
length(GridwatchDataFreqNonZero$frequency)  # Number left
GridwatchDataFreqZero <- GridwatchData[which(GridwatchData$frequency == 0.0),]
length(GridwatchDataFreqZero$frequency) # Number removed



# Summary statistics for the frequency
descdist(GridwatchDataFreqNonZero$frequency)

# Iterate over the years 2012 - 2019

GridwatchDataFreqNonZeroYear <- list()
GridwatchDataFreqNonZeroYearFreqMean <- list()
GridwatchDataFreqNonZeroYearFreqSD <- list()
GridwatchDataFreqNonZeroYearFreqSkewness <- list()
GridwatchDataFreqNonZeroYearFreqKurtosis <- list()
Years <- list()
index <- 0
for (j in 2012:2019) {
  index <- index+1
  print("Data for year")
  print(j)
  Years[[index]] <- j
  GridwatchDataFreqNonZeroYear[[j]]  <- GridwatchDataFreqNonZero[which(year(ymd_hms(GridwatchDataFreqNonZero$timestamp)) == j),]
  print(length(GridwatchDataFreqNonZeroYear[[j]]$frequency))
  print(descdist(GridwatchDataFreqNonZeroYear[[j]]$frequency))
  print("-------------------------")
  GridwatchDataFreqNonZeroYearFreqMean[[index]] <- mean(GridwatchDataFreqNonZeroYear[[j]]$frequency)
  GridwatchDataFreqNonZeroYearFreqSD[[index]] <- sd(GridwatchDataFreqNonZeroYear[[j]]$frequency)
  GridwatchDataFreqNonZeroYearFreqSkewness[[index]] <- skewness(GridwatchDataFreqNonZeroYear[[j]]$frequency)
  GridwatchDataFreqNonZeroYearFreqKurtosis[[index]] <- kurtosis(GridwatchDataFreqNonZeroYear[[j]]$frequency)
  hist(GridwatchDataFreqNonZeroYear[[j]]$frequency, freq=FALSE, breaks = "FD", xlim = c(49.7, 50.3),main=paste(Years[index],"GB Grid frequency distribution", sep = " "), xlab="Grid frequency (Hz)")
}

#Plot of summary stats by year

plot(Years, GridwatchDataFreqNonZeroYearFreqMean, xlab=" ", ylab="Mean", main="GB Grid frequency mean, 2012 - 2019")
plot(Years, GridwatchDataFreqNonZeroYearFreqSD, xlab=" ", ylab="Standard Deviation", main="GB Grid frequency standard deviation, 2012 - 2019")
plot(Years, GridwatchDataFreqNonZeroYearFreqSkewness, xlab=" ", ylab="Skewness", main="GB Grid frequency skewness, 2012 - 2019")
plot(Years, GridwatchDataFreqNonZeroYearFreqKurtosis, xlab=" ", ylab="Kurtosis", main="GB Grid frequency kurtosis, 2012 - 2019")


# Counts of 'extreme events' by year

GridwatchDataFreqNonZeroSigma <- list()
GridwatchDataFreqNonZeroSigmaYear <- list()
GridwatchDataFreqNonZeroSigmaYears <- list()
GridwatchDataFreqNonZeroFreqMean <- mean(GridwatchDataFreqNonZero$frequency)
GridwatchDataFreqNonZeroFreqSD <- sd(GridwatchDataFreqNonZero$frequency)


for(i in 1:20) {
  GridwatchDataFreqNonZeroSigma[[i]] <- GridwatchDataFreqNonZero[which(abs(GridwatchDataFreqNonZero$frequency - GridwatchDataFreqNonZeroFreqMean)  > i*GridwatchDataFreqNonZeroFreqSD),]
  GridwatchDataFreqNonZeroSigmaYear[[i]] <- list()
  GridwatchDataFreqNonZeroSigmaYears[[i]] <- list()
  for (j in 2012:2019){
    GridwatchDataFreqNonZeroSigmaYear[[i]][[j]] <- GridwatchDataFreqNonZero[which((abs(GridwatchDataFreqNonZero$frequency - GridwatchDataFreqNonZeroFreqMean)  > i*GridwatchDataFreqNonZeroFreqSD) & (year(ymd_hms(GridwatchDataFreqNonZero$timestamp)) == j)),]
    GridwatchDataFreqNonZeroSigmaYears[[i]] <- c(GridwatchDataFreqNonZeroSigmaYears[[i]], nrow(GridwatchDataFreqNonZeroSigmaYear[[i]][[j]]))
  }
  print(GridwatchDataFreqNonZeroSigmaYears[[i]])
  print(plot(Years,GridwatchDataFreqNonZeroSigmaYears[[i]], xlab="Calendar year", ylab="Frequency of events", main=paste(toString(i),"sigma GB grid frequency events, \n 2012 - 2019",sep =" ")))
}

# Fits

hist(GridwatchDataFreqNonZero$frequency, breaks=500, freq=FALSE, xlim = c(49.75,50.25), xlab = "Grid frequency (Hz)", main = "Distribution of GB Grid Frequency\n May 2011 - Jan 2020\n fitted with a normal distribution")
fitnormal <- fitdist(GridwatchDataFreqNonZero$frequency,
                 distr = 'norm',
                 method = 'mge',
                 start = list(mean=50,sd=1))
fitnormal
x <- seq(from=49,to=51,length.out = 500)
freqnrm <- dnorm(x,mean=fitnormal$estimate[1],sd=fitnormal$estimate[2])
lines(x,freqnrm,lwd=2, col="blue")


hist(GridwatchDataFreqNonZero$frequency, breaks=500, freq=FALSE, xlim = c(49.75,50.25), xlab = "Grid frequency (Hz)", main = "Distribution of GB Grid Frequency\n May 2011 - Jan 2020\n fitted with a skew Cauchy distribution")
fitsc <- fitdist(GridwatchDataFreqNonZero$frequency,
                 distr = 'sc',
                 method = 'mge',
                 start = list(xi=50,omega=1,alpha=0))

fitsc
x <- seq(from=49,to=51,length.out = 500)
freqnrm <- dsc(x,xi=fitsc$estimate[1],omega=fitsc$estimate[2], alpha=fitsc$estimate[3])
lines(x,freqnrm,lwd=2, col="blue")


# hist(GridwatchDataFreqNonZero$frequency, breaks=500, freq=FALSE, xlim = c(49.75,50.25), xlab = "Grid frequency (Hz)", main = "Distribution of GB Grid Frequency\n May 2011 - Jan 2020")
# fitlaplace <- fitdist(GridwatchDataFreqNonZero$frequency,
#                 distr = 'laplace',
#                 method = 'mge',
#                 start = list(m=50,s=0.05))
#fitlaplace  # fit fails - not sure why
# Instead we fit by eye
hist(GridwatchDataFreqNonZero$frequency, breaks=500, freq=FALSE, xlim = c(49.75,50.25), xlab = "Grid frequency (Hz)", main = "Distribution of GB Grid Frequency\n May 2011 - Jan 2020\n with an approximate Laplace distribution")
x<- seq(from=49,to=51,length.out = 500)
freqnrm <- dlaplace(x,m=50,s=0.035)
lines(x,freqnrm,lwd=2, col="yellow")

hist(GridwatchDataFreqNonZero$frequency, breaks=500, freq=FALSE, xlim = c(49.75,50.25), xlab = "Grid frequency (Hz)", main = "Distribution of GB Grid Frequency\n May 2011 - Jan 2020\n fitted with a logistic distribution")
fitlogis <- fitdist(GridwatchDataFreqNonZero$frequency,
                 distr = 'logis',
                 method = 'mge',
                 start = list(location=50,scale=1))
fitlogis

x <- seq(from=49,to=51,length.out = 500)
freqnrm <- dlogis(x,location=fitlogis$estimate[1],scale=fitlogis$estimate[2])
lines(x,freqnrm,lwd=2, col="blue")




hist(GridwatchDataFreqNonZero$frequency, breaks=500, freq=FALSE, xlim = c(49.75,50.25), xlab = "Grid frequency (Hz)", main = "Distribution of GB Grid Frequency\n May 2011 - Jan 2020\n fitted with a generalized normal distribution")
fitgnorm <- fitdist(GridwatchDataFreqNonZero$frequency,
                 distr = 'gnorm',
                 method = 'mge',
                 start = list(mu=50,alpha=sqrt(2), beta=2))
fitgnorm
x <- seq(from=49,to=51,length.out = 500)
freqnrm <- dgnorm(x,mu=fitgnorm$estimate[1],alpha=fitgnorm$estimate[2], beta=fitgnorm$estimate[3])
lines(x,freqnrm,lwd=2, col="blue")


# time series analysis

# Convert frequency into a timeseries

freq <- GridwatchData$frequency
freqNA <- freq 
freqNA[freqNA == 0] <- NA # Replace the zero values with NA
freq.ts <- ts(freq, start =c(2011,5,27,15,50), frequency = 365*24*12)  #Only approximate timing - need to get this right
freqNA.ts <- na.approx(ts(freqNA, start =c(2011,5,27,15,50), frequency = 365*24*12)) 
#Only approximate timing - need to get this right, also replacing NsA with approx values

plot(freqNA.ts, ylab="Grid frequency (Hz)", xlab="", main="GB Grid frequency to 2019")
plot(freq.ts,ylab="Grid frequency (Hz)", xlab="", main= paste(c("GB Grid frequency to 2019", "showing missing data as zero")))
plot(freq.ts, ylim=c(48, 50.5), ylab="Grid frequency (Hz)", xlab="", main= paste(c("GB Grid frequency to 2019", "showing missing data as zero")))

# Decompose time series 
freqNA.ts.decomp <- decompose(freqNA.ts)
plot(freqNA.ts.decomp) #graph not brilliant - not sure how to fix

#Separate out parts
freqNA.ts.trend <- freqNA.ts.decomp$trend
plot(freqNA.ts.trend, ylab="Grid frequency (Hz)", main="Trend for GB Grid frequency to 2019")
freqNA.ts.seasonal <- freqNA.ts.decomp$seasonal
plot(freqNA.ts.seasonal, ylab="Grid frequency (Hz)", main="Seasonal component for GB Grid frequency to 2019")
freqNA.ts.random <- freqNA.ts.decomp$random
plot(freqNA.ts.random, ylab="Grid frequency (Hz)", main="Random component for GB Grid frequency to 2019")


# Analysis of whole time series 
freqNA.ts.vector = as.vector(freqNA.ts)
hist(freqNA.ts.vector, breaks="FD", freq=FALSE, xlab = "Frequency (Hz)", main = "Histogram of GB Grid frequency to 2019")
hist(freqNA.ts.vector, breaks=500, freq=FALSE, xlim= c(49.8, 50.2), xlab = "Frequency (Hz)", main = "Histogram of GB Grid frequency to 2019")

hist(freqNA.ts.vector, breaks=500, freq=FALSE,xlim= c(49.8, 50.2), xlab = "Frequency (Hz)", main = paste(c("Histogram of GB Grid frequency to 2019", "with fitted normal distribution")))
fitnormal <- fitdist(freqNA.ts.vector,
                     distr = 'norm',
                     method = 'mge',
                     start = list(mean=50,sd=0.05))

fitnormal
x <- seq(from=50-0.2,to=50+0.2,length.out = 500)
freqnrm <- dnorm(x,mean=fitnormal$estimate[1],sd=fitnormal$estimate[2])
lines(x,freqnrm,lwd=2, col="blue")

hist(freqNA.ts.vector, breaks=500, freq=FALSE, xlim= c(49.8, 50.2), xlab = "Frequency (Hz)", main = paste(c("Histogram of GB Grid frequency to 2019", "with fitted skew normal distribution")))
fitsnorm <- fitdist(freqNA.ts.vector,
                    distr = 'snorm',
                    method = 'mge',
                    start = list(mean=50, sd=0.055, xi=0.4))

fitsnorm
x <- seq(from=50-0.2,to=50+0.2,length.out = 500)
freqnrm <- dsnorm(x,mean=fitsnorm$estimate[1],sd=fitsnorm$estimate[2], xi=fitsnorm$estimate[3])
lines(x,freqnrm,lwd=2, col="blue")


# Analysis of random part

length(freqNA.ts.random) - length(na.omit(freqNA.ts.random))   # Lots of NAs - not all time period calculated
freqNA.ts.random.vector <- as.vector(na.omit(freqNA.ts.random))  # Write as vector for fitting
hist(freqNA.ts.random.vector, breaks="FD", freq=FALSE, ylim=c(0,8), xlab = "Frequency (Hz)", main = "Histogram of random component of GB Grid frequency")

hist(freqNA.ts.random.vector, breaks="FD", freq=FALSE, ylim=c(0,8), xlab = "Frequency (Hz)", main = paste(c("Histogram of random component of GB Grid frequency to 2019", "with fitted normal distribution")))
fitnormal <- fitdist(freqNA.ts.random.vector,
                     distr = 'norm',
                     method = 'mge',
                     start = list(mean=0,sd=0.05))

fitnormal
x <- seq(from=-0.2,to=0.2,length.out = 500)
freqnrm <- dnorm(x,mean=fitnormal$estimate[1],sd=fitnormal$estimate[2])
lines(x,freqnrm,lwd=2, col="blue")


hist(freqNA.ts.random.vector, breaks="FD", freq=FALSE, ylim=c(0,8), xlab = "Frequency (Hz)", main = paste(c("Histogram of random component of GB Grid frequency to 2019", "with fitted generalized normal distribution")))
fitgnorm <- fitdist(freqNA.ts.random.vector,
                    distr = 'gnorm',
                    method = 'mge',
                    start = list(mu=0,alpha=1, beta=1))

fitgnorm
x <- seq(from=-0.2,to=0.2,length.out = 500)
freqnrm <- dgnorm(x,mu=fitgnorm$estimate[1],alpha=fitgnorm$estimate[2], beta=fitgnorm$estimate[3])
lines(x,freqnrm,lwd=2, col="blue")


hist(freqNA.ts.random.vector, breaks="FD", freq=FALSE, ylim=c(0,8), xlab = "Frequency (Hz)", main = paste(c("Histogram of random component of GB Grid frequency to 2019", "with fitted skew normal distribution")))
fitsnorm <- fitdist(freqNA.ts.random.vector,
                    distr = 'snorm',
                    method = 'mge',
                    start = list(mean=-0.01, sd=0.055, xi=0.4))

fitsnorm
x <- seq(from=-0.2,to=0.2,length.out = 500)
freqnrm <- dsnorm(x,mean=fitsnorm$estimate[1],sd=fitsnorm$estimate[2], xi=fitsnorm$estimate[3])
lines(x,freqnrm,lwd=2, col="blue")


# Overlay normal and skew normal distributions
hist(freqNA.ts.random.vector, breaks="FD", freq=FALSE, ylim=c(0,8), xlab = "Frequency (Hz)", main = paste(c("Histogram of random component of GB Grid frequency to 2019", "with fitted normal (blue) and skew normal (red) distributions")))
x <- seq(from=-0.2,to=0.2,length.out = 500)
freqnrm <- dgnorm(x,mu=fitgnorm$estimate[1],alpha=fitgnorm$estimate[2], beta=fitgnorm$estimate[3])
lines(x,freqnrm,lwd=2, col="blue")
freqnrm <- dsnorm(x,mean=fitsnorm$estimate[1],sd=fitsnorm$estimate[2], xi=fitsnorm$estimate[3])
lines(x,freqnrm,lwd=2, col="red")


# Analysis of seasonal component of the time series 




freqNA.ts.seasonal.vector <- as.vector(na.omit(freqNA.ts.seasonal))  # Write as vector for fitting
hist(freqNA.ts.seasonal.vector, breaks="FD", freq=FALSE, xlab = "Frequency (Hz)", main = "Histogram of seasonal component of GB Grid frequency")


hist(freqNA.ts.seasonal.vector, breaks="FD", freq=FALSE, xlab = "Frequency (Hz)", main = paste(c("Histogram of seasonal component of GB Grid frequency to 2019", "with fitted normal distribution")))
fitnormal <- fitdist(freqNA.ts.seasonal.vector,
                     distr = 'norm',
                     method = 'mge',
                     start = list(mean=0,sd=0.05))

fitnormal
x <- seq(from=-0.2,to=0.2,length.out = 500)
freqnrm <- dnorm(x,mean=fitnormal$estimate[1],sd=fitnormal$estimate[2])
lines(x,freqnrm,lwd=2, col="blue")

dev.off() # Stop saving as a pdf file
