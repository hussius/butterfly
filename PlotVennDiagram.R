install.packages('venneuler')
require(venneuler)


plotVenn <- function(A,B,COMMON,Names){
  info = c((A-COMMON), (B-COMMON), COMMON)
  Names = c(Names,paste(Names[1],Names[2],sep="&"))
  names(info) <- Names
  v <- venneuler(info)
  plot(v)
  
}



Names = c("blast2GO","ARGOT2")
A=40745
B=94512
COMMON=19860
