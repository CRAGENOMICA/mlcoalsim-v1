#theta and rec ex01 (real theta=10, rec=10)
tableT <- read.table(file="./ex01par_mcmc.out", header = T, sep = "", quote = "", dec = ".", as.is = FALSE, na.strings = "na",colClasses = NA, nrows = 5000,skip = 0, check.names = TRUE,strip.white = FALSE, blank.lines.skip = TRUE,comment.char = "#", allowEscapes = FALSE);

file = "./plot_parameters_sum.pdf";
pdf(file);
par(mfrow=c(1,3));
for(x in 1:2) {
		hist(tableT[,x],main = paste("Histogram parameters ",colnames(tableT)[x]),br=c(seq(from=0,to=25,by=0.5)),xlim=c(0,25));
}
plot(tableT[,1],tableT[,2],xlab=colnames(tableT)[1],ylab=colnames(tableT)[2])
dev.off();

