# CHECK THE FILE PATHS BEFORE RUNNING THIS SCRIPT
# THE SET PATHS ARE FOR greedy based alignment algorithm ALIGNMENTS --> basic variant (SAlign)

library(GOSemSim)
terms <- c('MF', 'BP')

sp1_set = c('mouse','mouse','yeast','mouse','mouse')
sp2_set = c('human','yeast','human','fly','worm')

args <- commandArgs(trailingOnly = TRUE)
cat(args, sep = "\n")

sp1 <- args[1]
sp2 <- args[2]
stage <- args[3] # t or b

	for(term_i in 1:2){ # MF and BP

	  term <- terms[term_i]
	  goterm <- tolower(term)

	  align_file <- paste(c("alignments/",sp1, '-', sp2, '-', stage, ".alignment"),collapse="")
	  alignment_file <- scan(align_file, what="character", sep='\n',quiet=TRUE)
	  len <- length(alignment_file)
	  d <- godata('org.Hs.eg.db', ont=term, computeIC=FALSE)

          fileConn <- paste(c("evaluation/semantic_files/",sp1,'-',sp2,'-',stage,'-',term,".txt"),collapse="")

	  path1 <- paste(c('evaluation/go-ids/', sp1) ,collapse="")
	  path2 <- paste(c('evaluation/go-ids/', sp2) ,collapse="")

	  match1 <- paste(c(path1, '/', goterm, '-go.txt'), collapse='')
	  match2 <- paste(c(path2, '/', goterm, '-go.txt'), collapse='')
	  main_stream <- ''
          
	  for(i in (1:len)){
	    local_stream <- ''
	    p<- strsplit(alignment_file[i], '\t')
	    p1 <- p[[1]][1]
	    p2 <- p[[1]][2]
	    line1 <- grep(p1, readLines(match1), value = TRUE)
	    line2 <- grep(p2, readLines(match2), value = TRUE)
	    val1 <- unlist(strsplit(line1, '\t')) # //'\t'
	    val2 <- unlist(strsplit(line2, '\t')) # //'\t'
	    if(is.null(val1[2])) {
	      print(c(paste(p1, "\t", sp1, "\t", "skipped")))
	      next
	    }
	    if(is.null(val2[2])) {
	      print(c(paste(p2, "\t", sp2, "\t", "skipped")))
	      next
	    }
	    go1 <- strsplit(val1[2], ',')
	    go2 <- strsplit(val2[2], ',')
	    sm <- 0
	    sm <- mgoSim(go1[[1]], go2[[1]], semData=d, measure="Wang", combine="BMA")
	    local_stream <- c(paste(p1, "\t", p2, "\t", sm, sep = ""))
	    main_stream <- c(trimws(paste(main_stream,'\n', local_stream, sep = "")))
	  }
	  writeLines(main_stream,fileConn) 
	  #close(fileConn)
}
print(c(paste("Semantic Similarities are computed and stored in evaluation directory.")))
print(c(paste("Run average_semantic_similarity file to get the average results.")))

