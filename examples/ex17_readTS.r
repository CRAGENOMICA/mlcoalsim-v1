tableR <- read.table(file="ex17_recpost.out", header = T, sep = "", quote = "", dec = ".",row.names=NULL, as.is = FALSE, na.strings = "NA",colClasses = NA, nrows = 100,skip = 0, check.names = TRUE,strip.white = FALSE, blank.lines.skip = TRUE,comment.char = "#", allowEscapes = FALSE);
tableT <- read.table(file="ex17_thetapost.out", header = T, sep = "", quote = "", dec = ".",row.names=NULL, as.is = FALSE, na.strings = "NA",colClasses = NA, nrows = 100,skip = 0, check.names = TRUE,strip.white = FALSE, blank.lines.skip = TRUE,comment.char = "#", allowEscapes = FALSE);
nloci <- 1;
for(x in 1:nloci) {
	par(mfrow=c(2,2))
	hist(tableR[,x]);
	hist(tableT[,x]);
	plot(tableR[,x],tableT[,x]);
	reg <- lm(tableR[,x]~tableT[,x]);
	reg;
	RT <- tableR[,x]/tableT[,x];
	hist(RT);
}