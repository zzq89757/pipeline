#!/usr/bin/env Rscript
### VERSION: 3.0
### AUTHOR: Fang Peng 
### Updated date: 2019-05

### Librarys required: 'getopt','seqTools' <installed automatically>
### export R_LIBS=$HOME/.R_lib:$R_LIBS


version <- '     VERSION: 3.0
      AUTHOR: Fang Peng 
Updated date: 2019-05
'
usage <- 'VERSION: 1.0 by Fang Peng
options:
  -v | --verbose
  -h | --help
  -i | --input <dir>
     Path to FASTQ file.
  -c | --cpu <num>
     CPUs, default: 20.
  -m | --mode <character>
     10x FASTQ as "10x" or "none", default: "none".
  -o | --output <dir>
      default: "result.txt".
'
spec = matrix(c(
  'version', 'v', 0, "logical",
  'help'   , 'h', 0, "logical",
  'input'  , 'i', 2, "character",
  'output'  , 'o', 2, "character",
  'mode'  , 'm', 2, "character",
  'cpu'  , 'c', 2, "character"
), byrow=TRUE, ncol=4);

lib <- .libPaths()
lib.acc <- lib[which(file.access(lib,2) != -1)]
if(length(lib.acc) == 0){
	if(file.access('~/.R_lib', 0) == -1){
		system('mkdir -p ~/.R_lib')
	}
	.libPaths(c('~/.R_lib', .libPaths()))
}else{
	.libPaths(c(lib.acc, lib[!lib %in% lib.acc]))
}
Sys.setenv(R_LIBS = paste(paste(Sys.getenv('HOME'),'/.R_lib',sep=""),Sys.getenv('R_LIBS'),sep=':'))
if(!suppressMessages(suppressWarnings(require(getopt)))){
	repos <- 'http://mirrors.tuna.tsinghua.edu.cn/CRAN/'
	install.packages("getopt", repos = repos, quiet = T)
}
suppressMessages(library(getopt))
if(!suppressMessages(suppressWarnings(require(seqTools)))){
	source("http://bioconductor.org/biocLite.R")
	biocLite('seqTools', suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
}
suppressMessages(library(seqTools))
suppressMessages(library(parallel))

opt = getopt(spec);
if (length(names(opt)) == 1){cat(usage); q(status=1)}
if (!is.null(opt$help)){cat(usage); q(status=1)}
if (!is.null(opt$version)){cat(version); q(status=1)}
#if(is.null(opt$version )){opt$version = FALSE }
if(is.null(opt$input)){cat(usage); q(status=1)}
if(is.null(opt$output)){opt$output = "result.txt"}
if(is.null(opt$cpu)){opt$cpu = 20}
if(is.null(opt$mode)){opt$mode = 'none'}

## function 
quality_ratio <- function(fq_2, cut=20){
	require(seqTools)
	quality_ratio <- character()
	for(i in 1:nFiles(fq_2)){
		quality_dis_2 <- phred(fq_2,i)
		colnames(quality_dis_2) <- 1:ncol(quality_dis_2)
		quality_ratio_2 <- format(sum(as.numeric(quality_dis_2[which(as.numeric(rownames(quality_dis_2)) > cut),])) / sum(as.numeric(quality_dis_2)) * 100, nsmall=2, digits=2)
		quality_ratio[i] <- quality_ratio_2
	}
	return(quality_ratio)
}

merge_quality_ratio <- function(fq_1, fq_2, cut=20){
	
}

nucl_freq <- function(fq_2){
	require(seqTools)
	bases <- c('A','T','G','C','N')
	nucl_ratio <- data.frame(Total_bases=numeric(), a=numeric(), c=numeric(), g=numeric(), t=numeric(), n=numeric(),stringsAsFactors=F)
	nucl_freq <- list()
	for(i in 1:nFiles(fq_2)){
		nucl <- nucFreq(fq_2,i)
		nucl_ratio[i,] <- c(sum(as.numeric(nucl)),rowSums(nucl[which(toupper(rownames(nucl)) %in% bases),]) / sum(as.numeric(nucl))* 100)
		nucl_freq[[i]] <- nucl[which(toupper(rownames(nucl)) %in% bases),]
	}
	nucl_ratio$gc <- nucl_ratio$g + nucl_ratio$c
	nucl_ratio[,2:ncol(nucl_ratio)] <- format(nucl_ratio[,2:ncol(nucl_ratio)], nsmall=2, digits=0)
	return(nucl_ratio)
}

## main
if(toupper(opt$mode) == "10X"){
	files <- list.files(path=opt$input, pattern="*(_)?(R)?2(_001)?.f(ast)?q(.gz)?$")
}else{
	files <- list.files(path=opt$input, pattern="*_(R)?[12](_001)?.f(ast)?q(.gz)?$")
}
if(length(files) == 0){
	warning(paste("  Can not find FASTQ files in dir: ", opt$input,"", sep="'"))
	q(save = 'no')
}

files <- files[file.info(file.path(opt$input, files))$size > 20]

cores <- detectCores(logical = FALSE)
if(cores < opt$cpu){opt$cpu <- cores}
if(length(files) < opt$cpu){opt$cpu <- length(files)}
fq <- mclapply(file.path(opt$input, files), fastqq, mc.cores = opt$cpu)

quality20 <- character()
for(i in 1:length(fq)){
	quality20 <- c(quality20, quality_ratio(fq[[i]],cut=20))
}
quality30 <- character()
for(i in 1:length(fq)){
	quality30 <- c(quality30, quality_ratio(fq[[i]],cut=30))
}

nucl_ratio <- data.frame()
total_reads <- numeric()
for(i in 1:length(fq)){
	nucl_ratio <- rbind(nucl_ratio, nucl_freq(fq[[i]]))
	total_reads <- c(total_reads, nReads(fq[[i]]))
}
#nucl <- nucl_freq(fq2)

q20 <- quality20
q30 <- quality30
#nucl_ratio <- nucl$nucl_ratio
#total_reads <- nReads(fq1) + nReads(fq2)
sample <- gsub('(_S\\d+)?(_L\\d+)?(_(R)?[12])?(_001)?.f(ast)?q(.gz)?','',files, perl=T)
data <- data.frame(Sample=sample, Total_Reads=total_reads, nucl_ratio, Q20=q20, Q30=q30, stringsAsFactors=F)

data.final <- data.frame(stringsAsFactors=FALSE)
if(length(which(duplicated(data$Sample))) > 0){
	for(i in data$Sample[which(duplicated(data$Sample))]){
		data.final[i,1] <- i
		for(j in 2:3){
			data.final[i,j] <- sum(data[data$Sample %in% i,j])
		}
		for(j in 4:11){
			data.final[i,j] <- paste(data[data$Sample %in% i,j], collapse=';') 
		}
	}
	unpaired_sample <- data$Sample[which(!data$Sample %in% data$Sample[which(duplicated(data$Sample))])]
	if(length(unpaired_sample) > 0){
		colnames(data.final) <- colnames(data)
		data.final <- rbind(data.final, data[data$Sample %in% unpaired_sample,])
	}
}else{
	data.final <- data
}
colnames(data.final) <- c('Sample','Total_Reads','Total_bases','A (%)','C (%)','G (%)','T (%)','N (%)','GC (%)','Q20 (%)','Q30 (%)')

if(toupper(opt$mode) == "10X"){
	samplename <- do.call(rbind,strsplit(data.final$Sample,"-"))
	sample <- unique(samplename[,1])
	data.final.10x <- data.frame(stringsAsFactors=FALSE)
	for(i in sample){
		nucl <- matrix(as.numeric(as.matrix(data.final[,-c(1:3)])),ncol=ncol(data.final[,-c(1:3)]))
		data.final.10x <- rbind(data.final.10x, cbind(Sample=i, t(colSums(data.final[samplename[,1] == i,c(2,3)])) , t(format(colSums(nucl[samplename[, 1] == i,] * data.final[samplename[,1] == i,3]) / sum(data.final[samplename[,1] == i,3]), nsmall=2, digits=0))))
	}
	colnames(data.final.10x) <- c('Sample','Total_Reads','Total_bases','A (%)','C (%)','G (%)','T (%)','N (%)','GC (%)','Q20 (%)','Q30 (%)')
	data.final.10x$Total_Reads <-  prettyNum(data.final.10x$Total_Reads, big.mark=',')
	data.final.10x$Total_bases <-  prettyNum(data.final.10x$Total_bases, big.mark=',')
	write.table(data.final.10x, file=paste(sub('.txt$',"",opt$output), ".10X.txt",sep=""), quote=F,sep="\t",row.names =F)
}

# output
data.final$Total_Reads <-  prettyNum(data.final$Total_Reads, big.mark=',')
data.final$Total_bases <-  prettyNum(data.final$Total_bases, big.mark=',')
write.table(data.final, file=opt$output, quote=F,sep="\t",row.names =F)


### quality score distribution along reads
# files <- c("yangpinS-AACCGTAA_HL3HCCCXY_L8_1.fq.gz","yangpinS-AACCGTAA_HL3HCCCXY_L8_2.fq.gz")
# r1 <- fastqq(files[1])
# r2 <- fastqq(files[2])
# r1_q <- phred(r1,1)
# r2_q <- phred(r2,1)
# r1_q_d <- apply(r1_q, 2, function(x, quality=0:93){
	# names(x) <- quality
	# x <- x[x != 0]
	# sum(x * as.numeric(names(x))) / sum(x)
# })
# r2_q_d <- apply(r2_q, 2, function(x, quality=0:93){
	# names(x) <- quality
	# x <- x[x != 0]
	# sum(x * as.numeric(names(x))) / sum(x)
# })
# quality <- data.frame(x=c(1:r1@maxSeqLen, (r1@maxSeqLen + 1): (r1@maxSeqLen +r2@maxSeqLen)), y=c(r1_q_d, r2_q_d))
# ggplot(data = quality, aes(x=x, y=y)) + 
	# geom_point(colour ="springgreen4",size=0.5) + ylim(0,45) + 
	# scale_x_continuous(breaks=seq(0,300,20)) + 
	# xlab("Position along reads") + ylab("Quality score") + 
	# ggtitle(paste0("Quality score distribution along reads (",sub("_1.fq.gz","",r1@filenames),")")) + 
	# geom_vline(xintercept = 150,linetype=2, colour = 'skyblue3')



