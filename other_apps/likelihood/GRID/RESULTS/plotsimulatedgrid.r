#theta and rec ex01 (real theta=10, rec=10)
tableP <- read.table(file="./ex01par_grid.out", header = T, sep = "", quote = "", dec = ".", as.is = FALSE, na.strings = "na",colClasses = NA, nrows = 5000,skip = 0, check.names = TRUE,strip.white = FALSE, blank.lines.skip = TRUE,comment.char = "#", allowEscapes = FALSE);

tableT <- read.table(file="./ex01grid.out", header = T, sep = "", quote = "", dec = ".", as.is = FALSE, na.strings = "na",colClasses = NA, nrows = 5000,skip = 0, check.names = TRUE,strip.white = FALSE, blank.lines.skip = TRUE,comment.char = "#", allowEscapes = FALSE);

file = "./plot_surface_grid.pdf";
pdf(file);
par(mfrow=c(1,2));
#levelplot(tableT$Total ~ tableP$thetaw * tableP$recombination,main="Grid: Likelihood",xlab="thetaw",ylab="recombination",region=T,col.regions=gray((16:0)/16));
contourplot(tableT$Total ~ tableP$thetaw * tableP$recombination,main="Grid: Likelihood",xlab="thetaw",ylab="recombination",region=T,col.regions=gray((16:0)/16),labels=FALSE,cuts=7);
dev.off();

