options(
  HTTPUserAgent =
    sprintf(
      "R/%s R (%s)",
      getRversion(),
      paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])
    )
)

install.packages("arrow", repos = "https://packagemanager.rstudio.com/all/__linux__/focal/latest")