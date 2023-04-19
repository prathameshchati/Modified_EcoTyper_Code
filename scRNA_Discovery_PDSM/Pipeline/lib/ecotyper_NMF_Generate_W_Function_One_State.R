library(NMF)

NMFGenerateW_One_State <- function(cellTypeBinaryH, scaledExpMat) {
    
    #trainig_gene_set = gsub("__pos", "", rownames(W)[grepl("__pos", rownames(W))])
    #new_data = new_data[match(trainig_gene_set, rownames(new_data)),]
    #rownames(new_data) = trainig_gene_set
        
    cellTypeBinaryH = cellTypeBinaryH[colnames(scaledExpMat)]
    scaledExpMat[is.na(scaledExpMat)] = 0
    
    cellTypeBinaryH = as.matrix(cellTypeBinaryH)
    scaledExpMat = as.matrix(scaledExpMat)
    
    # best method
    
    # CHANGE FILTERING IF NECESSARY
#     scaledExpMatF = scaledExpMat[,colSums(scaledExpMat) > 0]
#     scaledExpMatF = scaledExpMatF[rowSums(scaledExpMatF) > 0,]

    # scaledExpMatF = scaledExpMat[rowSums(scaledExpMat) > 0,]
    # scaledExpMatF = scaledExpMatF[,colSums(scaledExpMatF) > 0]

    to_predict = posneg(scaledExpMat)

    to_predict[is.na(to_predict)] = 0

#     to_predict = to_predict[,colSums(to_predict) > 0]
#     to_predict = to_predict[rowSums(to_predict) > 0,]

    cellTypeBinaryHF = cellTypeBinaryH[,colnames(to_predict)]
    cellTypeBinaryHF = as.matrix(cellTypeBinaryHF)
    cellTypeBinaryHF = t(cellTypeBinaryHF)
    rownames(cellTypeBinaryHF) = c(cell_type)

    my_method <- function (i, v, x, copy = FALSE, eps = .Machine$double.eps, ...) {
        w <- .basis(x)
        h <- .coef(x)
        nb <- nbterms(x)
        nc <- ncterms(x)
    #   h <- NMF:::std.divergence.update.h(v, w, h, nbterms = nb, ncterms = nc, copy = copy)
        w <- NMF:::std.divergence.update.w(v, w, h, nbterms = nb, ncterms = nc, copy = copy)
        if (i%%10 == 0) {
    #        h <- pmax.inplace(h, eps, icterms(x))
             w <- pmax.inplace(w, eps, ibterms(x))
        }
        if (copy) {
            .basis(x) <- w
    #       .coef(x) <- h
        }
        return(x)
    }

    dummy = rnmf(nrow(to_predict), H = cellTypeBinaryHF)

#     dummyW = dummy@W
#     dummyH = dummy@H
    
     dummyWF = dummy@W
     dummyHF = dummy@H


#     dummyW[is.na(dummyW)] = 0
#     dummyH[is.na(dummyH)] = 0

#     dummyWF = dummyW[,colSums(dummyW) > 0]
#     dummyHF = dummyH[,colSums(dummyH) > 0]

#     dummyWF = dummyWF[rowSums(dummyWF) > 0,]
#     dummyHF = dummyHF[rowSums(dummyHF) > 0,]

#     dummyHF = dummyHF[,colnames(to_predict)]
    
    
    # NEW CODE ADDED - TOGGLE THIS LINE AS NEEDED
    # dummyHF = dummyHF + 0.0001
    dummyHF = dummyHF + .Machine$double.eps

    my.seeding.method <- function(model, target){
        basis(model) <- dummyWF #estim.r@fit@W
        # initialize H randomly
        coef(model) <- dummyHF
        # return updated object
        return(model)
    }

    nmf_method <- NMFStrategy('my-method', 'brunet', Update = my_method, objective = 'KL', Stop='connectivity')

    new_nmf = nmf(to_predict, nrow(cellTypeBinaryHF), nrun = 1, method = nmf_method, seed = my.seeding.method, .opt='P1')

    W = new_nmf@fit@W

    W[is.na(W)] = 0

#     rownames(W) = c(paste0(rownames(scaledExpMatF), "__pos"),paste0(rownames(scaledExpMatF), "__neg"))
    rownames(W) = rownames(scaledExpMat)

    colnames(W) = rownames(cellTypeBinaryHF)

    W
    
#     WF = W[rowSums(W) > 0,]

#     WF


}



