library(MASS)

vecOfFiles=character()
args= commandArgs(trailingOnly = T)
myargs =args[length(args)]
print(typeof(args))
if(length(args)==0){
	vecOfFiles=append(vecOfFiles,"peaksForR.txt")
	print("no args passed to R")
	quit()
}else{
	for(i in 1:length(args)){
		vecOfFiles=append(vecOfFiles,args[i])
		print(args[i])
	}
}
print(vecOfFiles[[1]])

sink("r.out")

vOfPvalues = numeric()
testFunc <- function(x) {
	v = strsplit(x, '\t')
	y =as.integer(strsplit(as.character(v[3]), ',', fixed=TRUE)[[1]])
	results = tryCatch({
		ff = fitdistr(y, "Negative Binomial")
		est_size = ff$estimate[1]
		est_mu = ff$estimate[2]
		#v[2] is height
		cat(
		sprintf("%s\t%s\t%s\t%e\n", v[1], v[2], v[3], 1-pnbinom(q=(as.integer(v[2])-1), size=est_size, mu=est_mu, log.p=F) )
		)
	}, error = function(e) {
		cat(
		sprintf("%s\t%s\t%s\t%e\n", v[1], v[2], v[3], 0 )
		)
	})
}

for (i in 1:length(vecOfFiles)) { 
	t = read.table(file=vecOfFiles[i], sep="\t")
 	tmp = apply(t, 1, testFunc )
}

sink()
