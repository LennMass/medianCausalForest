#' Save HTML Index Embedded in Shiny App
#'
#' @inheritParams saveFiles
#' @returns No return value. It is used to save necessary files temporarily
#' to run Shiny App.

saveInd <- function(filePath){

  indstr <- "<!doctype html>
<html>
<head>
	<meta charset='utf8'>
	<meta name='viewport' content='width=device-width'>

	<title>Svelte app</title>

	<link rel='icon' type='image/png' href='/favicon.png'>
	<link rel='stylesheet' href='/global.css'>
  <link rel='stylesheet' href='/bundle.css'>
  <link href='https://fonts.googleapis.com/css?family=Montserrat:400,700&display=swap' rel='stylesheet'>
</head>

<body>
	<script src='/bundle.js'></script>
</body>
</html>"
  indsplt <- strsplit(indstr, "\n")
  indvec <- unlist(indsplt)

  fileConn <- file(paste(filePath, "/shinyapp/www/index.html", sep=""))
  writeLines(indvec, fileConn)
  close(fileConn)

}
