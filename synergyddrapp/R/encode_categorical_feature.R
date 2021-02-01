#' Encode drug names
#' @param x name of drug
#' @return numeric encode
#' @export

drug2idx <- function(x){
	load(file = "../data/drug_summary.rda")   ### remove "MSC" fron drug names!
	drug2idx <- list()
	for(i in seq_along(unique(drug_summary[,1]))){
		drug2idx[[drug_summary[i,1]]] <- i
	}
	drug2idx[[x]]	
}

#' Encode mode of actions
#' @param x mode of action
#' @return numeric encode
#' @export
moa2idx <- function(x){
        load(file = "../data/drug_summary.rda")
        moa2idx <- list()
        for(i in seq_along(unique(drug_summary[,2]))){
                moa2idx[[drug_summary[i,2]]] <- i
        }
	moa2idx[[x]]
}
