read_LIS_data <-  function(indices = NA, variable = "dhi", raw = FALSE){

## Store files names only
files <- list.files( "/DATA/", pattern = "\\.dta$", full.names = TRUE)

## survey number
if(is.na(indices)) {indices = 1:length(files)}
#c(3, 7, 12, 14, 22, 24, 26, 27, 29, 33) ##

## store results
results = list()

######################### Implement on all Countries #################3
for (i in indices) { ##

  # Read and process each file (your existing code)
  file <- files[i]
  data <- haven::read_dta(file)

  microdata = data[,c("hwgt", variable)]
  microdata= microdata[complete.cases(microdata),]

  ## Expand data by HH weight
  if(raw){
  microdata$hwgt <- round(microdata$hwgt*10, 0)  ## multiply HHweight by 10 to minimize expansion error
  microdata=replicate_rows(microdata)
  }

}
}
