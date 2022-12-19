# Any problems let me know and I can sort it! - Tyson
# Requires specific layout of dataset, see qGene Analysis Script (Dec 2022) Sample Data.xslx for details
# I have confirmed all analysis is identical to that of old Excel spreadsheet

# Standard Package to install
install.packages("dplyr")
library(dplyr)

# First Step: Import Dataset -> Your Data using template, then proceed
# Before analysis confirm Amplification efficiencies and standard deviation of samples is usable

# Change this to your Dataset name
Data <- read.csv("homozygotesData.csv")

#Sets % SEM cutoff for samples (default 25%)
SEMCutoff = 25

# Do you wish to automatically censor data over the above SEM% threshold? 0 = no, 1 = yes (default yes)
UseCutoff = 0

# Change std Amp. efficiencies if necessary
TarAmp = 2.07 # SINV
RefAmp = 1.92 # RpL32

# If you need to delete the dataframe \/ (do at start of each analysis)
rm(Normalised)
rm(list)

# Run me (housekeeping stuff to sort the data/prepare for analysis)
Normalised <- data.frame("Sample","MNE","SE of MNE", "SE of MNE as % (x<=25%)")
i=0
val = 1:nrow(Data)
val <- val[seq(1,(length(val)), 2)]
w = unlist(c(Data[4]))
x = unlist(c(Data[3]))
y = unlist(c(Data[2]))
s = matrix(x,nrow=1,ncol=nrow(Data))
v = matrix(w,nrow=1,ncol=nrow(Data))
t = matrix(y,nrow=nrow(Data),ncol=1)
list <- data.frame("Censored Samples","MNE","SE of MNE", "SE of MNE as % (x>25%)")

# For loop processes all the data and adds it to the dataframe
for (i in val) {
  
  # Calculates MNE of Duplicates
  NE1 = (RefAmp^v[1,i])/(TarAmp^s[1,i])
  NE2 = (RefAmp^v[1,i+1])/(TarAmp^s[1,i+1])
  MNE2 = ((NE1+NE2)/2)
  
  # Calculates SE of MNE
  b = (RefAmp^v[1,i])/(TarAmp^s[1,i])
  d = (RefAmp^v[1,i+1])/(TarAmp^s[1,i+1])
  SEofMNE = (sd(c(b,d)))/sqrt(2)
  
  # Calculates SE of MNE as %
  SEofMNEasP = (SEofMNE/MNE2)*100
  
  # Adds each set of data to the data frame (you may wish to change the digits= option to be smaller if you wish)
  if (UseCutoff == 1) {
    if (SEofMNEasP > SEMCutoff) {
      cat("Sample:",t[i,1],"was censored due to exceeding SEM cutoff of",SEMCutoff, "with a value of",(round(SEofMNEasP,digits=2)),"\n")
      list[nrow(list) + 1,] = c(t[i,1], round(MNE2,digits=10), round(SEofMNE,digits=5), round(SEofMNEasP,digits=2)) }
    else {
      Normalised[nrow(Normalised) + 1,] = c(t[i,1], round(MNE2,digits=10), round(SEofMNE,digits=5), round(SEofMNEasP,digits=2))}
    }
  else {
    Normalised[nrow(Normalised) + 1,] = c(t[i,1], round(MNE2,digits=10), round(SEofMNE,digits=5), round(SEofMNEasP,digits=2))}
  }

# Visualize the processed Dataset
Normalised

# Visualize censored data (if applicable)
list

# Write to your drive as an excel file, you will have to edit below file pathway
write.csv(Normalised,"homozygotesAnalysis.csv", row.names = FALSE)


# Write censored data to your drive (if you wish)
write.csv(list,"homozygotesCensored.csv", row.names = FALSE)
write.csv(list,"sampleCensored.csv", row.names = FALSE)



