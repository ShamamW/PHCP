#!/usr/bin/Rscript
time = Sys.time()
library(dplyr)
####read in the parameters
args = commandArgs(trailingOnly=TRUE)
hh <- paste(unlist(args),collapse=' ')
listoptions <- unlist(strsplit(hh,'--'))[-1]
options.args <- sapply(listoptions,function(x){
  unlist(strsplit(x, ' '))[-1]
})
options.names <- sapply(listoptions,function(x){
  option <-  unlist(strsplit(x, ' '))[1]
})

names(options.args) <- unlist(options.names)

if(!all(c("out","target","chr","ancient","modern")%in% options.names)) {stop("arguments are missing with no default")}

donors = F # prefix of phased donors data
theta <- 0.00138488 # Mutation rate
N <- 64.5698 # Ne

OUTPUT = options.args["out"]
ANCIENT_GENOME = options.args["target"] #id of the ancient sample
num = as.numeric(options.args["chr"]) # chromosome number
prefix_ancient <- options.args["ancient"] # prefix of ancient genomes data
prefix_modern = options.args["modern"]
for (i in c("donors","theta","N")){
  if(i %in% options.names) assign(i,options.args[i])
}

##### organize the input data
modern_data <- read.table(gzfile(paste(prefix_modern,num,"haps","gz",sep = ".")), stringsAsFactors = F) #modern donors data, haps format
ancient_data <- read.table(gzfile(paste(paste(prefix_ancient,num,sep = "_"),"gz", sep = ".")), header=F, stringsAsFactors = F) # ancient genomes data
ID_modern <- read.table(paste(prefix_modern,"ind",sep = "."), stringsAsFactors = F, col.names = c("id","pop")) # ID's of modern data. First coulumn with samples' ids, second column with populations name.
ID_ancient <- read.table(paste(prefix_ancient,"ind", sep = "."), header = F, stringsAsFactors = F, col.names = c("id","pop"))
colnames(ancient_data) <- c("chr","bp","a1","a2","snp_id","location",ID_ancient$id)
filterd_ancient_data <- ancient_data %>% select(1:6, one_of(ANCIENT_GENOME))

# Keep only the chosen donors data.
# First create character vector with the samples' ids and assign it as the names of the haps file coulumns
# afterwards, merge the haps file with the genetic map

if (donors!=F) {
  donors_pop <- read.table(donors, header=F, stringsAsFactors = F, col.names = "pop")
  ID_modern$index <- 1:nrow(ID_modern)
  donors_ord <- semi_join(ID_modern, donors_pop, by = "pop") %>% pull(index)
  donors_index <- sort(c(5+donors_ord*2-1, 5+donors_ord*2))
  modern_data <- modern_data[,c(1:5,donors_index)]
  ID_modern  = ID_modern[donors_ord, ]
} 

ID.haps.names <- ID_modern$id[sort(rep(1:nrow(ID_modern),2))]
ID.haps.names <- paste(ID.haps.names, c("A","B"), sep = "_")
colnames(modern_data) <- c("chr","snp_id","bp","a1","a2",ID.haps.names)

all_genomes <- inner_join(filterd_ancient_data, modern_data, by = c("chr","bp","a1","a2")) %>% arrange(chr, bp)

ancient_genomes <- as.matrix(all_genomes[,ANCIENT_GENOME, drop = F])
mode(ancient_genomes) <- "logical"
reference_genomes <- as.matrix(all_genomes[,ID.haps.names, drop = F])
mode(reference_genomes) <- "numeric"
mode(reference_genomes) <- "logical"
genetic_location <- all_genomes$location
J <- ncol(reference_genomes) #number of donors haplotypes
prior <- rep(1/J,J) #vector of a-priori copying probabilities from one haplotype
prior.prior <- outer(prior,prior)



ancient <- ancient_genomes[ ,ANCIENT_GENOME]
#remove SNPs that are missing in the ancient genome from both donors data and ancient genome

informative_snps <- !is.na(ancient)
df <- reference_genomes[informative_snps,] 
ancient <- ancient[informative_snps]

L <- nrow(df) #number of SNPs

## calculate the genetic distance between each pair of following SNPs
rho <- N * c(diff(genetic_location[informative_snps]),Inf)

#forward algorithm (with scaling)

Theta <- matrix(theta,J,J)
Prior.A.B <- outer(prior,prior)
df.snp <- df[1,1:J]
ancient.equal.df <- (ancient[1] == df.snp)*(1-2*theta)
A.ne.B <- outer(df.snp, df.snp, function(x,y) x != y)
emission_l <- sweep(Theta, 1, ancient.equal.df, "+")
emission_l[A.ne.B] <- 0.5
alpha_hat <- emission_l*prior.prior
ct <- numeric(L) #scaling vector
ct[1] <- 1/sum(alpha_hat)

alpha_hat <- alpha_hat*ct[1]

prob_array <- array(dim = c(J,J,L))
prob_array[ , ,1] <- alpha_hat

for(l in 2:L) {
  expr <- exp(-rho[l-1])
  recomb2 <- (1-expr)^2*(1/(J^2))
  recomb1 <- expr*(1-expr)*(1/J)*rowSums(alpha_hat)
  recomb02 <- exp(-2*rho[l-1])*alpha_hat + recomb2
  recomb021 <- sweep(recomb02, 1, recomb1, "+") 
  tot_rec <- sweep(recomb021, 2, recomb1, "+") 
  
  df.snp <- df[l,1:J]
  ancient.equal.df <- (ancient[l] == df.snp)*(1-2*theta)
  A.ne.B <- outer(df.snp, df.snp, function(x,y) x != y)
  emission_l <- sweep(Theta, 1, ancient.equal.df, "+")
  emission_l[A.ne.B] <- 0.5
  
  alpha_hat <- emission_l * tot_rec
  ct[l] <- 1/sum(alpha_hat)
  alpha_hat <- ct[l]*alpha_hat
  prob_array[ , ,l] <- alpha_hat
}

print("finish forward")
print(Sys.time() - time)

#### backward
beta_hat <- ct[L]*matrix(1, nrow = J, ncol = J)
xk <- matrix(nrow = J, ncol = J)
Sigma <- 0
len <- 0

for(l in (L-1):1) {
  expr <- exp(-rho[l])
  
  df.snp <- df[l+1,1:J]
  ancient.equal.df <- (ancient[l+1] == df.snp)*(1-2*theta)
  A.ne.B <- outer(df.snp, df.snp, function(x,y) x != y)
  emission_l <- sweep(Theta, 1, ancient.equal.df, "+")
  emission_l[A.ne.B] <- 0.5
  v <- beta_hat*prob_array[, , (l+1)]/ct[l+1]
  u <- (prob_array[ , ,l]*beta_hat*(exp(-2*rho[l]) + 2*expr*(1-expr)*(1/J) + (1-expr)^2*(1/J^2))) * emission_l
  Sigma<-Sigma + v - u
  
  b <- beta_hat*emission_l
  recomb2 <- (1-expr)^2*(1/(J^2))*sum(b)
  recomb1 <- expr*(1-expr)*(1/J)*rowSums(b)
  recomb02 <- exp(-2*rho[l])*b + recomb2
  recomb021 <- sweep(recomb02, 1, recomb1, "+") 
  beta_hat <- sweep(recomb021, 2, recomb1, "+") 
  
  beta_hat <- beta_hat*ct[l]
  len <- len + 0.5* (rho[l]/N) * (prob_array[ , ,l]*beta_hat/ct[l] + v) #length been copied
}
xk1 <- prob_array[ , ,1]*beta_hat/ct[l]
xk <- xk1 + Sigma # number of chunks
colnames(xk) <- ID.haps.names
rownames(xk) <- ID.haps.names
colnames(len) <- ID.haps.names
rownames(len) <- ID.haps.names
print("finish HMM")
print(Sys.time() - time)

#### export the output
chrom_num <- paste("chrom",num,sep="")
output_count<-paste(chrom_num,ANCIENT_GENOME,"chunkscount",OUTPUT,sep="_")
write.table(xk,output_count)
output_len<-paste(chrom_num,ANCIENT_GENOME,"length",OUTPUT,sep="_")
write.table(len,output_len)


print(Sys.time() - time)