source("env.R")
source("utilities.R")
source("plotting.R")

options(error = recover)

env.assemble()

get.connection <- function(dbname) {
	if (!exists('.connection', where=.GlobalEnv)) {
		.connection <<- dbConnect(RMySQL::MySQL(), dbname = dbname)
	} else if (class(try(dbGetQuery(.connection, "SELECT 1"))) == "try-error") {
		dbDisconnect(.connection)
		.connection <<- dbConnect(RMySQL::MySQL(), dbname = dbname)
	}

  return(.connection)
}

dir.my.cluster <- function(dir.output, dir.my.cluster, instance.pid, instance.time, session.id) {
	paste(dir.output, dir.my.cluster, sprintf("%d.%d.%d", instance.pid, instance.time, session.id), sep="/")
}
