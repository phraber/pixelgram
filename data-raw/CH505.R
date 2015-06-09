aas_file <- system.file("extdata", "CH505-gp160.fasta", package="pixgramr")
tre_file <- system.file("extdata", "CH505-gp160.tre", package="pixgramr")

CH505 <- list()
CH505$aas <- seqinr::as.matrix.alignment(seqinr::read.alignment(aas_file, 
	         "fasta"))
CH505$tre <- ape::read.tree(tre_file)
CH505$wpi <- sapply(1:length(CH505$tre$tip.label), function(i)
	         unlist(strsplit(CH505$tre$tip.label[i], '[.]'))[[1]])

devtools::use_data(CH505, overwrite=T)
