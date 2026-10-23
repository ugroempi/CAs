# 1. Create a hidden internal environment to store the data safely
.pkg_internal_env <- new.env(parent = emptyenv())
.pkg_internal_env$smolenBigFrame <- NULL # Default fallback

.onLoad <- function(libname, pkgname) {
  # 2. Wrap the loading process in a tryCatch block
  tryCatch({
    # Replace with your actual path or URL
    pfadGithubSmolen <-
      "https://raw.githubusercontent.com/chrispbsmolen/ca-tools/main/table/ca_table.csv"
    # 2. Build the request and configure a strict timeout (e.g., 3 seconds)
    req <- httr2::request(pfadGithubSmolen)
    req <- httr2::req_timeout(req, 3)

    # 3. Perform the request (will throw an error if offline, 404, or 500)
    resp <- httr2::req_perform(req)
    resp_text <- httr2::resp_body_string(resp)
    .pkg_internal_env$smolenBigFrame <- utils::read.table(text = resp_text,
                        header=TRUE, sep=",", quote = "\"", dec = ".",
                        fill = TRUE, comment.char = "")
  })
}

function (file, header = TRUE, sep = ",", quote = "\"", dec = ".",
          fill = TRUE, comment.char = "", ...)
  read.table(file = file, header = header, sep = sep, quote = quote,
             dec = dec, fill = fill, comment.char = comment.char, ...)
.onAttach <- function(libname, pkgname){
  if (is.null(smolenBigFrame))
    packageStartupMessage(
      "Note: smolenBigFrame could not be loaded.",
      "eCAN and eCAK will use frozen Colbourn tables.")
}
makeActiveBinding("smolenBigFrame",
    function() .pkg_internal_env$smolenBigFrame,
    env = asNamespace("CAs"))
