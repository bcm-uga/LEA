struct2geno <- function(input.file, ploidy, FORMAT = 1, extra.row = 0, extra.column = 0){
  
  # check input.file 
  if(missing(input.file)) {
    stop("'input.file' argument is missing.")}
  else if (!is.character(input.file)){
    stop("'input.file' argument must be of type character.")}
  
  if (!file.exists(input.file)) stop("Input file not found.")
  
  # check ploidy
  ploidy = test_integer("ploidy", ploidy, NULL)
  if (ploidy != 1 & ploidy !=2) 
    {stop("Ploidy must be equal to 1 (haploids) or 2 (diploids).")}
  diploid <- (ploidy == 2)
  
  # check Format
  FORMAT = test_integer("FORMAT", FORMAT, NULL)  
  if (FORMAT != 1 & FORMAT !=2) 
    {stop("FORMAT must be equal to 1 (1 row per individual) or 2 (2 rows per individual).")}
  
  if (!diploid & FORMAT == 2) stop("FORMAT = 2 can be used with diploids only.")    

  # check extra.row/extra.col
  extra.row = test_integer("extra.row", extra.row, NULL) 
  extra.column = test_integer("extra.column", extra.column, NULL)   

  dat = read.table(input.file)

    if (extra.row > 0) dat <- dat[-(1:extra.row),]
    if (extra.column > 0) dat <-  dat[,-(1:extra.column)]
    n = nrow(dat)
    L = ncol(dat)
    
    if (FORMAT == 1 & diploid == FALSE) {n.ind = n; n.loc = L}
    if (FORMAT == 1 & diploid == TRUE) {n.ind = n; n.loc = L/2}
    if (FORMAT == 2 & diploid == TRUE) {n.ind = n/2; n.loc = L}
    cat("Input file in the STRUCTURE format. The genotypic matrix has", n.ind, "individuals and", n.loc,
        "markers.","\n")
    cat("The number of extra rows is", extra.row,"and the number of extra columns is",extra.column,".\n")
  

  dat = as.matrix(dat)

  unique.dat = unique(as.numeric(dat))
  missing.dat = unique.dat[unique.dat < 0]

  if (length(missing.dat) == 0)  cat("The input file contains no missing genotypes.","\n")
  if (length(missing.dat) == 1)  cat("Missing alleles are encoded as",missing.dat,", converted as 9.\n")
  if (length(missing.dat) > 1) stop("Multiple values for missing data.","\n")



  # Convert allelic data into absence/presence data at each locus
  # Results are stored in the "dat.binary" object

  if (FORMAT == 1 & diploid == FALSE) {
    dat.binary = NULL
    for (j in 1:L){
      allele = sort(unique(dat[,j]))
      for (i in allele[allele >= 0]) dat.binary=cbind(dat.binary, dat[,j]== i)
      LL = ncol(dat.binary)
      ind = which(dat[,j] < 0)
      if (length(ind) != 0){dat.binary[ind, (LL - length(allele) + 2):LL] = -9}
     }
    }

  if (FORMAT == 1 & diploid == TRUE) {
    dat.2 = matrix(NA, ncol = L/2, nrow = 2*n)
    for (ii in 1:n){
      dat.2[2*ii-1,] = dat[ii, seq(1,L,by = 2)]
      dat.2[2*ii,] = dat[ii,seq(2,L,by = 2)]
      }
    L = ncol(dat.2)

    dat.binary = NULL

    for (j in 1:L){
      allele = sort(unique(dat.2[,j]))
      for (i in allele[allele >= 0]) dat.binary=cbind(dat.binary, dat.2[,j]==i)
      LL = dim(dat.binary)[2]
      ind = which(dat.2[,j] < 0)
      if (length(ind) != 0){dat.binary[ind, (LL - length(allele) + 2):LL] = -9}
      }
    }



  if (FORMAT == 2 & diploid == TRUE) {
    dat.binary = NULL
    for (j in 1:L){
      allele = sort(unique(dat[,j]))
      for (i in allele[allele >= 0]) dat.binary=cbind(dat.binary, dat[,j]==i)
      LL = dim(dat.binary)[2]
      ind = which(dat[,j] < 0)
      if (length(ind) != 0){dat.binary[ind, (LL - length(allele) + 2):LL] = -9}
    }}

  # Compute a genotype count for each allele (0,1,2 or 9 for a missing value)
  # The results are stored in 'genotype'

  n = nrow(dat.binary)
  L <- ncol(dat) 
  LL <- ncol(dat.binary)
  
  if (diploid == TRUE){
    n = n/2
    genotype <- matrix(NA, nrow = n, ncol = LL)
    for(i in 1:n){
      genotype[i,] <- dat.binary[2*i-1, ] + dat.binary[2*i, ]
      genotype[i, (genotype[i,] < 0)] <-  NA
    }
    if (LL == L & FORMAT == 1) {genotype <- genotype[ , seq(2, LL, by = 2)]}
    if (LL == 2*L & FORMAT == 2) {genotype <- genotype[ , seq(2, LL, by = 2)]}    
   }
  
 
  
  if (FORMAT == 1 & diploid == FALSE){
    genotype = dat.binary
    for(i in 1:n){
      genotype[i, (genotype[i,] < 0)] = NA
    }
    if (LL == 2*L) {genotype <- genotype[ , seq(2, LL, by = 2)]}
    }
  
  genotype[is.na(genotype)] <- 9
  lst.monomorphic <- apply(genotype, 2, FUN = function(x) {length(unique(x[x != 9]))})
  
  if (sum(lst.monomorphic == 1) > 0) {
    warning("Monomorphic alleles generated during conversion were removed. \n")
    genotype <- genotype[,lst.monomorphic > 1] 
    }


  
  write.lfmm(R = genotype, output.file = paste(input.file,".lfmm", sep=""))
  write.geno(R = genotype, output.file = paste(input.file,".geno", sep=""))
  cat("Output files:", paste(input.file,".geno  .lfmm.", sep=""),"\n")
}
