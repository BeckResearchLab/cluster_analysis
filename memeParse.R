# to be sourced

# if training.set is F, the training set will be loaded, but not validated
# if training.set is a data frame with one column length and row names the
# sequence names, the data frame will be compared to the one generated from
# the file, the value of training.set in the result list will be NULL
memeParse <- function(meme.filename, training.set = F) {
	# load the meme output file
	if (file.exists(meme.filename)) {
		meme.output <- readLines(meme.filename)
	} else {
		stop(sprintf("unable to open meme output file: %s", meme.filename))
	}

	# return data
	meme.data <- list()

	# parse out the list of sequences in the training set so that
	# it can be validated against the expected
	# this will help catch syncronization failures between the 
	# file system meme output and the clustering results
	ts.lines <- grep("^TRAINING SET", meme.output)
	if (length(ts.lines) <= 0) {
		stop(sprintf("unable to find TRAINING SET text in meme output file: %s", meme.filename))
	}
	if (length(ts.lines) != 1) {
		stop(sprintf("multiple TRAINING SET headers in meme output file: %s", meme.filename))
	}
	ts.start.line <- ts.lines[1] + 6
	# note the fixed = T because we are grepping on a regex special character '*'
	ts.end.line <- grep("********************************************************************************",
		fixed = T,
		meme.output[ts.start.line:length(meme.output)]
	)[1] - 2 + ts.start.line
	file.training.set <- do.call(rbind,
		lapply(strsplit(meme.output[ts.start.line:ts.end.line], "[\\t\\s]+", perl = T),
			function (x) {
			# there will either be two triplets on a line or a single triplet
				if (length(x) == 6) {
					return(rbind(x[1:3], x[4:6]))
				} else if (length(x) == 3) {
					return(x[1:3])
				} else {
					stop(
					sprintf("unexpected number of fields in file '%s' on training set line like:\n",
							meme.filename, paste(x, collapse=" ")
						)
					)
				}
			}
		)
	)
	file.training.set <- data.frame(length = as.integer(file.training.set[,3]),
		row.names = file.training.set[,1])
	# compare the two data sets, if the call supplied a reference or expected training set
	# data frame with a single column length
	if (!identical(training.set, F) && !identical(training.set, file.training.set)) {
		stop (
			sprintf("the expected training set and the training set in meme output are different for file: %s",
				meme.filename
			)
		)
	}
	
	# now find the location of each of the motif header lines
	lines <- grep("^MOTIF\\s+\\d", meme.output, perl = T)
	motifs <- length(lines)
	if (motifs <= 0) {
		stop(sprintf("unable to find MOTIF header lines in meme output file: %s", meme.filename))
	}
	# split the header lines
	header.split <- strsplit(meme.output[lines], "[\\t\\s]+", perl = T)

	# parse the position specific score matrices out of the file for all motifs at once
	pssms <- list()
	for (i in 1:motifs) {
		pssm.lines <- grep(sprintf("Motif %d position-specific probability matrix", i), meme.output)
		if (length(pssm.lines) > 0) {
			pssm.header <- strsplit(meme.output[pssm.lines + 2], " ")[[1]]
			pssm.w <- as.numeric(pssm.header[6])
			# extract the nucleotides on columns 2,3,4,5
			pssm <- do.call(rbind, 
				strsplit(meme.output[pssm.lines + 2 + 1:pssm.w], "\\s+", perl = T))[, 2:5]
			# convert to a matrix
			pssm <- matrix(as.numeric(pssm), nrow = pssm.w, ncol = 4, byrow = F)
			pssms[[i]] <- list(pssm = pssm)
		} else {
			stop(sprintf("unable to parse PSSM for motif %d in file: %s", i, meme.filename))
		}
	}

	# iterate over the motifs and extract the header, the consensus and the sites
	for (i in 1:motifs) {
		splt <- header.split[[i]]
		motif <- as.integer(splt[2])
		width <- as.integer(splt[6])
		sites <- as.integer(splt[9])
		llr <- as.integer(splt[12])
		# fix up for exponent in meme output via the sub
		e.value <- as.numeric(sub("\\+", "", splt[15]))

		# pull the first line of the multilevel consensus
		description.start.line <- grep(sprintf("Motif %d Description", motif), meme.output)
		consensus.start.line <- grep("Multilevel", 
			meme.output[(description.start.line+1):length(meme.output)])[1] + description.start.line
		consensus <- strsplit(meme.output[consensus.start.line], "[\\t\\s]+", perl=T)[[1]][2]

		# pull the sites matrix
		positions.start.line <- grep(sprintf("Motif %d sites sorted by position p-value", motif), 
			meme.output) + 4
		positions.end.line <- grep("--------------------------------------------------------------------------------", 
			meme.output[(positions.start.line + 1):length(meme.output)])[1] + positions.start.line - 1
		# note, we only want the sequence name, start, p-value and center site column
		# i.e. ignore the flanking sequences reported by meme
		# first regularize the split rows
		regular.rows <- lapply(strsplit(meme.output[positions.start.line:positions.end.line], "[\\t\\s]+", perl = T),
			function (x) {
				if (length(x) == 4) {
					# the 5' flanking portion of the site sequence is not present,
					# this can happen if the site appears at the very begining of the sequence
					return(x[c(1:4)])
				} else if (length(x) == 5 || length(x) == 6) {
					# either the line is complete or the 3' flanking region is absent
					return(x[c(1:3, 5)])
				} else {
					stop(
						sprintf("the sites postion for motif %d in file %s has an unexpected number of fields: %d", 
							i,
							meme.filename,
							length(x)
						)
					)
				}
			}
		)
		positions <- do.call(rbind, regular.rows)
		colnames(positions) <- c("gene", "start", "p.value", "site")
		# note that we don't use gene as a rowname because it may appear more than once
		# if a motif occurs more than once in a single sequence
		positions <- data.frame(gene = positions[, "gene"], 
			start = as.integer(positions[, "start"]), 
			p.value = as.numeric(positions[, "p.value"]), 
			site = positions[, "site"])

		# assemble the full record
		meme.data[[motif]] <- list(training.set = training.set, width = width, sites = sites, llr = llr, e.value = e.value, consensus = consensus, pssm = pssms[[motif]]$pssm, positions = positions)
	}

	return(meme.data)
}
