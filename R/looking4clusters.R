### run clusters ###

## Input Vars
# selectedk: Number of clusters selected by the user.
# If 'NULL', default values are used (kmin=2,kmax=10), else, a rank of 5 units
# around 'selectedk' is calculated
##

run_kmeans <- function(object,selectedk=NULL){

data <- object$data

if(is.null(selectedk)){
    kmin<-2
    kmax<-10
}else{
    kmin<-(selectedk-5)
    kmax<-(selectedk+5)
}

if(kmin <= 1)
    kmin <- 2
if(kmax >= nrow(data))
    kmax <- nrow(data)-1
if(kmin >= kmax)
    kmin <- kmax-1

iter<-kmax-kmin+1

# kmeans (applied for all selected 'k' values)
all.kmeans<-lapply(kmin:kmax, function(k){
    kmeans(data, k, nstart=1 ,iter.max = 15)
})

if(is.null(selectedk)){
    numberClustersKmeans <- calinski_harabasz_index(iter,
        all.kmeans,nrow(data),kmin)
}else{
    numberClustersKmeans <- selectedk
}

# Clustering table

for(i in seq_len(iter)){
    clusters <- as.factor(paste0("kmeans_",all.kmeans[[i]]$cluster))
    optim_cluster <- FALSE
    if(length(levels(clusters))==numberClustersKmeans){
        optim_cluster <- TRUE
    }
    object <- addcluster(object,clusters,"kmeans",optim_cluster=optim_cluster)
}

return(object)

# 'numberClustersKmeans' is the default number of clusters values to represent
# in the graphs if automatic mode is select, else, 'selectedk' is used. 
}

calinski_harabasz_index <- function(iter,all.kmeans,N,kmin){
    # Calinski-Harabasz Index
    SSw <- vapply(seq_len(iter), function(x){
        all.kmeans[[x]]$tot.withinss
    },numeric(1))
    SSb <- vapply(seq_len(iter), function(x){
        all.kmeans[[x]]$betweens
    },numeric(1))
    chIndex <- vapply(seq_len(iter), function(x){
        (SSb[x]*(N-(x+kmin-1)))/(SSw[x]*((x+kmin-1)-1))
    },numeric(1))

    # Number of optimum clusters
    return(which(chIndex==max(chIndex))+kmin-1)
}

## Input Vars
# distance: the distance measure to be used. This must be one of "euclidean",
# "maximum", "manhattan", "canberra", "binary" or "minkowski".
# agglomeration: the agglomeration method to be used.  This should be one of
# "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or
# "centroid".
# selectedk: Number of clusters selected by the user. If 'NULL', default values
# are used (kmin=2,kmax=10), else, a rank of 5 units around 'selectedk' is
# calculated .
##

run_pam_hclust <- function(object, distance="euclidean",
        agglomeration="complete", selectedk=NULL, threads=NULL){
    if(requireNamespace("fpc",quietly=TRUE)){

data <- object$data

if(is.null(selectedk)){
    kmin<-2
    kmax<-10
}else{
    kmin<-(selectedk-5)
    kmax<-(selectedk+5)
}

if(kmin <= 1){
    kmin <- 2
}
if(kmax >= nrow(data)){
    kmax <- nrow(data)-1
}
if(kmin >= kmax){
    kmin <- kmax-1
}

iter <- kmax-kmin+1

dissimilarity <- get_dissimilarity(distance,data,threads)

pamresults <- run_pam(selectedk,dissimilarity,kmin,kmax,iter,data)
allClassfPam <- pamresults[[1]]
numberClustersPAM <- pamresults[[2]]

hclustresults <- run_hclust(agglomeration,dissimilarity,
    kmin,kmax,selectedk,iter)
Hclusters <- hclustresults[[1]]
numberClustersHclust <- hclustresults[[2]]

object <- clustering_tables(object,iter,allClassfPam,numberClustersPAM,
    Hclusters,numberClustersHclust)

    }else{
        warning("Install 'fpc' to get pam and hclust clusters.")
    }
    return(object)
}

run_pam <- function(selectedk,dissimilarity,kmin,kmax,iter,data){
## PAM:
allClassfPam <- NULL
numberClustersPAM <- NULL
if(requireNamespace("cluster",quietly=TRUE)){
    if(is.null(selectedk)){
        # Calinski Harabasz
        pamk<-fpc::pamk(dissimilarity, krange=kmin:kmax, criterion="ch",
            usepam=TRUE, scaling=FALSE, diss=TRUE)
        numberClustersPAM<-pamk$nc
    }else{
        numberClustersPAM <- selectedk
    }

    # All clusters 
    all.pam <- lapply(kmin:kmax, function(k){
        cluster::pam(dissimilarity, k, diss=TRUE, cluster.only=FALSE,
            keep.diss=FALSE, keep.data=FALSE)
    })

    allClassfPam <- vapply(seq_len(iter),function(k){
        all.pam[[k]]$clustering
    }, integer(nrow(data)))

# Explanation of any default parameters: 
# criterion: "ch" = calinski Harabasz index
# usepam: If TRUE, PAM clustering method is applied, else, CLARA clustering
# method is computed. 
# diss: Are input data a dissimilarity matrix?
# scaling: After calculate the dissimilarity matrix, Want to scale it?
# cluster.only: If TRUE, only the clustering will be computed and returned
}else{
    warning("Install 'cluster' to get pam clusters.")
}
return(list(allClassfPam,numberClustersPAM))
}

get_dissimilarity <- function(distance,data,threads){
    # Calculate dissimilarity matrix (parallel)
    if(!(distance %in% c("euclidean", "maximum", "manhattan",
            "canberra", "binary", "minkowski"))){
        distance <- "euclidean"
    }
    if(requireNamespace("parallelDist",quietly=TRUE)){
        dissimilarity <- parallelDist::parDist(data,
            method=distance, threads=threads)
    }else{
        dissimilarity <- dist(data,method=distance)
    }
    return(dissimilarity)
}

run_hclust <- function(agglomeration,dissimilarity,kmin,kmax,selectedk,iter){
## HCLUST:
Hclusters <- NULL
numberClustersHclust <- NULL
if(requireNamespace("dendextend",quietly=TRUE)){

if(!(agglomeration %in% c("ward.D", "ward.D2", "single", "complete",
        "average", "mcquitty", "median", "centroid")))
    agglomeration <- "complete"
hierarchical<-hclust(dissimilarity,method = agglomeration)
Hclusters<-dendextend::cutree(hierarchical,k=seq(kmin,kmax,1))

if(is.null(selectedk)){
    # Calinski Harabasz
    chIndex<-numeric()
    for(i in seq_len(iter))
        chIndex[i] <- fpc::cluster.stats(d = dissimilarity,
            clustering = Hclusters[,i])$ch
    numberClustersHclust<-(which(chIndex==max(chIndex))+kmin-1)
}else{
    numberClustersHclust <- selectedk
}

}else{
    warning("Install 'dendextend' to get hclust clusters.")
}
return(list(Hclusters,numberClustersHclust))
}

clustering_tables <- function(object,iter,allClassfPam,numberClustersPAM,
    Hclusters,numberClustersHclust){
# Clustering tables

for(i in seq_len(iter)){
    if(length(allClassfPam)){
        clusters <- as.factor(paste0("pam_",allClassfPam[,i]))
        optim_cluster <- FALSE
        if(length(levels(clusters))==numberClustersPAM){
            optim_cluster <- TRUE
        }
        object <- addcluster(object, clusters,
            "pam", optim_cluster = optim_cluster)
    }
    if(length(Hclusters)){
        clusters <- as.factor(paste0("hclust_",Hclusters[,i]))
        optim_cluster <- FALSE
        if(length(levels(clusters))==numberClustersHclust){
            optim_cluster <- TRUE
        }
        object <- addcluster(object, clusters, "hclust",
            optim_cluster = optim_cluster)
    }
}

# 'numberClustersPAM' and 'numberClustersHclust' are the default number of
# clusters values to represent in the graphs if automatic mode is select, else,
# 'selectedk' is used.
return(object)
}

### run dimensional reductions ###

run_pca <- function(object){

    PCAcomponents <- prcomp(object$data, scale=FALSE)
    pca<-PCAcomponents$x[,seq_len(2)]

    return(addreduction(object,pca,"pca"))
}

## Input Vars
# perplex: numeric. Perplexity parameter (should not be bigger than
# 3 * perplexity < nrow(X)-1). This value effectively controls how many nearest
# neighbours are taken into account when constructing the embedding in the
# low-dimensional space (default: 30)
# maxIter: integer. Number of iterations (default: 1000)
##

run_tsne <- function(object,perplex=30,maxIter=1000){
    if(requireNamespace("Rtsne",quietly=TRUE)){
        data <- object$data
        if(3 * perplex > nrow(data)-1){
            perplex <- (nrow(data)-1)/3
        }
        tSNEcomponents <- Rtsne::Rtsne(data, dims=2, perplexity=perplex,
            verbose=FALSE, check_duplicates=FALSE, max_iter=maxIter)
        tsne <- tSNEcomponents$Y
        colnames(tsne) <- c("tSNE1","tSNE2")
        object <- addreduction(object,tsne,"tsne")
    }else{
        warning("Install 'Rtsne' to get tsne dimensionality reduction.")
    }
    return(object)
}

run_mds <- function(object,threads=NULL){

    data <- object$data

    # Distance matrix
    if(requireNamespace("parallelDist",quietly=TRUE)){
        distMatrix <- parallelDist::parDist(data,
            method="canberra", threads=threads)
    }else{
        message("installing 'parallelDist' can improve performance")
        distMatrix <- dist(data,method="canberra")
    }

    # Classical Multi Dimensional Scaling (cMDS)
    cmds <- cmdscale(distMatrix, k=2)
    colnames(cmds) <- c("cMDS1","cMDS2")

    return(addreduction(object,cmds,"cmds"))
}

run_nmf <- function(object){
    if(requireNamespace("NMF",quietly=TRUE)){
        data <- object$data
        if(!sum(colSums(data)==0)){
            if(!sum(rowSums(data)==0)){
                if(!sum(data<0)){
                    NMFcomponents <- NMF::nmf(data, rank=2, method="brunet")
                    nmf <- NMF::basis(NMFcomponents)
                    colnames(nmf) <- c("NMF1","NMF2")
                    object <- addreduction(object,nmf,"nmf")
                }else{
                    warning(
"Your data contains some negative entries, this is not supported for nmf."
                    )
                }
            }else{
                warning(
"Your data has rows that are all zeros, this is not supported for nmf."
                )
            }
        }else{
            warning(
"Your data has columns that are all zeros, this is not supported for nmf."
            )
        }
    }else{
        warning("Install 'NMF' to get nmf dimensionality reduction.")
    }
    return(object)
}

run_umap <- function(object){
    if(requireNamespace("uwot", quietly=TRUE)){
        data <- object$data
        n_neighbors <- 15
        if(n_neighbors>(nrow(data)/3)){
            n_neighbors <- floor(nrow(data)/3)
        }
        if(n_neighbors < 2){
            n_neighbors <- 2
        }
        UMAPcomponents <- uwot::umap(data, n_neighbors = n_neighbors)
        colnames(UMAPcomponents) <- c("UMAP1","UMAP2")
        object <- addreduction(object,UMAPcomponents,"umap")
    }else{
        warning("Install 'uwot' to get umap dimensionality reduction.")
    }
    return(object)
}

### looking4clusters object management ###
l4c <- function(data, groups = NULL, components = FALSE, running_all = TRUE,
    distance = "euclidean", agglomeration = "complete", selectedk = NULL,
    perplex = 30, maxIter = 1000, threads = NULL, force_execution = FALSE){

    if(is.null(dim(data))){
        stop("data: incorrect dimensions. A matrix like object is required.")
    }

    if(running_all && (ncol(data) > 5000 || nrow(data) > 5000)){
        if(force_execution){
            message(
"Too large matrix, could cause performance problems with some methods"
            )
        }else{
            message(
"Too large matrix, could cause performance problems with some methods,
they will be omitted"
            )
        }
    }

    object <- create_l4c(data,components)

    if(!is.null(groups)){
        object <- addcluster(object,groups,myGroups=TRUE)
    }

    if(running_all){
        object <- running_clusters(object, selectedk, data, force_execution,
            distance, agglomeration, threads)
        object <- running_reductions(object, data, force_execution,
            threads, perplex, maxIter)
    }

    return(object)
}

running_clusters <- function(object, selectedk, data, force_execution,
        distance, agglomeration, threads){

    message("Running kmeans...")
    object <- run_kmeans(object,selectedk)
    if(force_execution || !(nrow(data) > 5000)){
        message("Running pam and hclust...")
        object <- run_pam_hclust(object, distance, agglomeration,
            selectedk, threads)
    }

    return(object)
}

running_reductions <- function(object, data, force_execution,
        threads, perplex, maxIter){

    message("Running pca...")
    object <- run_pca(object)
    message("Running tsne...")
    object <- run_tsne(object,perplex,maxIter)
    if(force_execution || !(nrow(data) > 5000)){
        message("Running mds...")
        object <- run_mds(object,threads)
    }
    if(force_execution || !(ncol(data) > 5000)){
        message("Running nmf...")
        object <- run_nmf(object)
    }
    message("Running umap...")
    object <- run_umap(object)

    return(object)
}

Ncomponents <- function(variability){
    nComp <- 4
    divisor <- 10
    if(length(variability)*0.9>1){
        variability <-
            variability[seq_len(as.integer(length(variability)*0.9))]
    }

    change <- variability[-length(variability)] - variability[-1]
    propChange <-
        (mean(change)/mean(sort(change,decreasing=TRUE)[seq_len(3)]))*100
    if(propChange>10 && propChange<=15){
        divisor <- 5
    }
    if(propChange>15){
        return(nComp)
    }
    minVar <- max(change)/divisor
    CutoffCriterion <- (change<minVar)
    position <- numeric()
    for(i in seq_len(length(variability)-1)){
        if(CutoffCriterion[i]==TRUE){
            k <- k+1
        }else{
            k <- 0
        }
        position <- c(position,k)
        if(k>7 || (i>10 && k>3)){
            n <- (max(which(position==1))-1)
            if(n>1){
                nComp <- n
            }
        }
    }
    return(nComp)
}

create_l4c <- function(data,components=FALSE){

    samples <- rownames(data)
    if(is.null(samples)){
        samples <- paste0("sample_",seq_len(nrow(data)))
    }

    if(components){
        if(nrow(data)<5){
            stop(
"You cannot apply to components for less than five samples."
            )
        }
        if(ncol(data)<5){
            stop(
"You cannot apply to components for less than five variables."
            )
        }

        ldata <- log2(data+1)
        lpca <- prcomp(ldata)
        var <- summary(lpca)$importance[2,]
        nc <- Ncomponents(var)
        pca <- prcomp(data,scale=FALSE)
        data <- pca$x[,seq_len(nc)]
    }

    return(structure(list(data=data.matrix(data), samples=samples,
        variables=colnames(data), options=list()), class="looking4clusters"))
}

addreduction <- function(object, data, name=NULL){
    if(!inherits(object,"looking4clusters")){
        stop("object: must be a 'looking4clusters' object")
    }
    if(nrow(data)!=length(object$samples)){
        stop("data: there must be one row per sample")
    }
    if(!length(object$reductions)){
        object$reductions <- list()
    }
    if(is.null(name)){
        name <- paste0("reduction_",length(object$reductions)+1)
    }else{
        name <- clean_names(name)
    }
    object$reductions[[name]] <- data[,seq_len(2)]
    return(object)
}

addcluster <- function(object, data, name=NULL, groupStatsBy=FALSE,
        myGroups=FALSE, optim_cluster=FALSE){
    if(!inherits(object,"looking4clusters")){
        stop("object: must be a 'looking4clusters' object")
    }
    if(length(data)!=length(object$samples)){
        stop("data: there must be one per sample")
    }
    if(!length(object$clusters)){
        object$clusters <- list()
    }
    if(is.null(name)){
        name <- paste0("cluster_",length(object$clusters)+1)
    }else{
        name <- clean_names(name)
    }
    data <- as.factor(data)
    if(length(object$clusters[[name]])){
        if(is.factor(object$clusters[[name]])){
            object$clusters[[name]] <-
                data.frame(V1=object$clusters[[name]],V2=data)
            colnames(object$clusters[[name]]) <-
                vapply(object$clusters[[name]],function(x){
                    return(paste0("levels_",length(levels(x))))
                }, character(1))
            attr(object$clusters[[name]],"optim_cluster") <-
                length(levels(object$clusters[[name]][[1]]))
        }else{
            object$clusters[[name]][[paste0("levels_",
                length(levels(data)))]] <- data
        }
        if(optim_cluster){
            attr(object$clusters[[name]], "optim_cluster") <-
                length(levels(data))
        }
    }else{
        object$clusters[[name]] <- data
    }
    if(groupStatsBy){
        object$options$groupStatsBy <- c(object$options$groupStatsBy,name)
    }
    if(myGroups){
        object$options$myGroups <- name
    }
    return(object)
}

clean_names <- function(name){
    return(gsub("[^A-Za-z0-9]","_",name))
}

print.looking4clusters <- function(x, ...){
    cat("An object of class looking4clusters\n")
    cat(paste0(length(x$variables)," variables across ",
        length(x$samples)," samples\n"))
    if(length(x$clusters)){
        cat(paste0(length(x$clusters)," clusters added: ",
            paste0(names(x$clusters),collapse=", "),"\n"))
    }
    if(length(x$reductions)){
        cat(paste0(length(x$reductions)," dimensional reductions added: ",
            paste0(names(x$reductions),collapse=", "),"\n"))
    }
}

looking4clusters <- function(data, groups = NULL,
    components = FALSE, running_all = TRUE, distance = "euclidean",
    agglomeration = "complete", selectedk = NULL, perplex = 30,
    maxIter = 1000, threads = NULL, force_execution = FALSE){
    l4c(data, groups, components, running_all,
        distance, agglomeration, selectedk,
        perplex, maxIter, threads, force_execution)
}
