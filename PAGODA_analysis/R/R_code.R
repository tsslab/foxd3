# load library
if (!require("scde")) {
require(devtools)
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)
}
library(scde)

app = readRDS("pagoda_app_86_cells")

# show app in the browser (port 1468)
show.app(app, “foxd3”, browse = TRUE, port = 1468)
