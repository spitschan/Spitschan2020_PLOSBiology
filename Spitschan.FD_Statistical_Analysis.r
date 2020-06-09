# Load libraries
library(openxlsx)
library(BayesFactor)
library(pracma)
library(ggplot2)
library(gridExtra)

# Load Excel spreadsheet
fd_data <- read.xlsx('Spitschan.FD_Sample_Protocol.xlsx')

# Extract only those rows that actually have measurements in them
fd_data <- fd_data[!is.na(fd_data$Measurement), ]

# Iterate parameters for one participant
NSpacing = 12
ratio = logspace(-1.5, 1, NSpacing)
#ratio = 0.5
NSims = 1000

# Number of participants
NParticipantsMax = 8

# Set up some empty arrays to capture data in the loop
BFs_h0 = array(dim = c(length(ratio), NParticipantsMax, NSims))
BFs_h1 = array(dim = c(length(ratio), NParticipantsMax, NSims))

p1 <- array()
sdev = 1

for (a in 1:length(ratio)) {
  # SNR
  for (p in 1:NParticipantsMax) {
    for (s in 1:NSims) {
      # Resample data N times
      
      Meas_h0 = vector()
      Meas_h1 = vector()
      participantID = vector()
      
      # Create some empty variables
      for (r in 1:p) {
        # Generate the data for p participants
        # Simulate data under H1
        Meas_h1_baseline <- rnorm(nrow(fd_data), 0, sdev)
        Meas_h1_sinewave_A <-
          rnorm(nrow(fd_data), sdev * ratio[a], 0.1) # Generate a sine wave with noisy amplitude
        Meas_h1_sinewave <-
          Meas_h1_sinewave_A * fd_data$Sine.component # Multiply with sine wave
        Meas_h1 <-
          append(Meas_h1, Meas_h1_baseline + Meas_h1_sinewave)
        
        # Simlate data under H0
        Meas_h0 <-
          append(Meas_h0, rnorm(nrow(fd_data), 0, sdev)) # This is simply the noisy distribution
        
        # Add participant ID here
        participantID <-
          append(participantID, rep(r, nrow(fd_data)))
      }
      
      fd_data1 <-
        data.frame(
          rep(fd_data$`Time.since.awake.[hour]`, p),
          rep(fd_data$`Time.since.entering.FD.[hour]`, p),
          rep(fd_data$`Sine.component`, p),
          rep(fd_data$`Cosine.component`, p),
          Meas_h0,
          Meas_h1,
          participantID
        )
      colnames(fd_data1) <-
        c(
          'Time.since.awake.[hour]',
          'Time.since.entering.FD.[hour]',
          'Sine.component',
          'Cosine.component',
          'Meas_h0',
          'Meas_h1',
          'ParticipantID'
        )
      
      # Set up a linear model for the full model and the null model
      fd_model_full_h0 <-
        lmBF(
          `Meas_h0` ~ `Time.since.awake.[hour]` + `Time.since.entering.FD.[hour]` + `Sine.component` + `Cosine.component` + `ParticipantID`,
          data = fd_data1
        )
      fd_model_null_h0 <-
        lmBF(
          `Meas_h0` ~ `Time.since.awake.[hour]` + `Time.since.entering.FD.[hour]` + `ParticipantID`,
          data = fd_data1
        )
      
      # Calculate the Bayes Factor
      BF_h0 <- fd_model_full_h0 / fd_model_null_h0
      BFs_h0[a, p, s] <- extractBF(BF_h0)$bf
      
      # Set up a linear model for the full model and the null model
      fd_model_full_h1 <-
        lmBF(
          `Meas_h1` ~ `Time.since.awake.[hour]` + `Time.since.entering.FD.[hour]` + `Sine.component` + `Cosine.component`,
          data = fd_data1
        )
      fd_model_null_h1 <-
        lmBF(
          `Meas_h1` ~ `Time.since.awake.[hour]` + `Time.since.entering.FD.[hour]`,
          data = fd_data1
        )
      
      # Calculate the Bayes Factor
      BF_h1 <- fd_model_full_h1 / fd_model_null_h1
      BFs_h1[a, p, s] <- extractBF(BF_h1)$bf
      #BFs_error[a, s] <- extractBF(BF)$error
      
      # Calculate posterior. Adding code here for demonstration
      #posteriorcalcs = posterior(BF_h0, iterations = 1000)
      #summary(posteriorcalcs)
    }
  }
}

falsePositive <- array(dim = c(length(ratio), NParticipantsMax))
truePositive <- array(dim = c(length(ratio), NParticipantsMax))
ratios <- array(dim = c(length(ratio), NParticipantsMax))
NParticipants <- array(dim = c(length(ratio), NParticipantsMax))
thresh <- 10
for (a in 1:length(ratio)) {
  for (p in 1:NParticipantsMax) {
    ratios[a, p] <- ratio[a]
    NParticipants[a, p] <- p
    falsePositive[a, p] <-
      length(BFs_h0[a, p, BFs_h0[a, p, ] > thresh]) / length(BFs_h0[a, p, ])
    truePositive[a, p] <-
      length(BFs_h1[a, p, BFs_h1[a, p, ] > thresh]) / length(BFs_h1[a, p, ])
  }
}


# Turn them all into long vectors
falsePositive <- c(falsePositive)
truePositive <- c(truePositive)
ppv <- (truePositive) / (falsePositive + truePositive)
ratios <- c(ratios)
NParticipants <- c(NParticipants)

# Assemble as a data frame
roc <-
  data.frame(ratios, falsePositive, truePositive, NParticipants, ppv)

# Plot the data
p1 <-
  ggplot(roc, aes(
    x = log10(ratios),
    y = truePositive,
    color = NParticipants
  )) + geom_point() + geom_line(aes(group = NParticipants)) + ylim(0, 1) + coord_fixed(ratio = 2) + ggtitle('True positives')
p2 <-
  ggplot(roc, aes(
    x = log10(ratios),
    y = falsePositive,
    color = NParticipants
  )) + geom_point() + geom_line(aes(group = NParticipants)) + ylim(0, 1) + coord_fixed(ratio = 2) + ggtitle('False positives')
p3 <-
  ggplot(roc, aes(
    x = log10(ratios),
    y = ppv,
    color = NParticipants
  )) + geom_point() + geom_line(aes(group = NParticipants)) + ylim(0, 1) + coord_fixed(ratio = 2) + ggtitle('Sensitivity')

# Grab the legend
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x)
    x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Assemble the data into a grid
aleg <- g_legend(p1)
p1 <-
  p1 + theme_classic() + theme(legend.position = "none") + labs(x = "Log10 a/sigma", y = "True positive rate")
p2 <-
  p2 + theme_classic() + theme(legend.position = "none") + labs(x = "Log10 a/sigma", y = "False positive rate")
p3 <-
  p3 + theme_classic() + theme(legend.position = "none") + labs(x = "Log10 a/sigma", y = "Sensitivity")

# Save analysis
pdf("figures/Spitschan.FD_BF_10.pdf",
    width = 8,
    height = 8)

# Arrange the data
grid.arrange(p1, p2, p3, aleg, nrow = 2, ncol = 2)

# Close the pdf file
dev.off() 
