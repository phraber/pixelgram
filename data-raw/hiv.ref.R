nts_file <- system.file("extdata", "hiv-ref.fasta", package="pixgramr")
tre_file <- system.file("extdata", "hiv-ref.tre", package="pixgramr")

hiv.ref <- list()
hiv.ref$nts <- seqinr::as.matrix.alignment(seqinr::read.alignment(nts_file, "fasta"))
hiv.ref$tre <- ape::read.tree(tre_file)
hiv.ref$clade <- 
    gsub("^Ref-", "", gsub("[12]$", "",
            sapply(1:length(hiv.ref$tre$tip.label), function(i)
		unlist(strsplit(hiv.ref$tre$tip.label[i], '[.]'))[[1]])))

devtools::use_data(hiv.ref, overwrite=T)
