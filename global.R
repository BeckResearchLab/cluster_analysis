source("env.R")
source("utilities.R")
source("plotting.R")

options(error = recover)

env.assemble()

getConnection <- function(dbname) {
	if (!exists('.connection', where=.GlobalEnv)) {
		.connection <<- dbConnect(RMySQL::MySQL(), dbname = dbname)
	} else if (class(try(dbGetQuery(.connection, "SELECT 1"))) == "try-error") {
		dbDisconnect(.connection)
		.connection <<- dbConnect(RMySQL::MySQL(), dbname = dbname)
	}

  return(.connection)
}


