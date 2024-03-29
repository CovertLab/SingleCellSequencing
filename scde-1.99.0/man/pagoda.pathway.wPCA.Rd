% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/functions.R
\name{pagoda.pathway.wPCA}
\alias{pagoda.pathway.wPCA}
\title{Run weighted PCA analysis on pre-annotated gene sets}
\usage{
pagoda.pathway.wPCA(varinfo, setenv, n.components = 2,
  n.cores = detectCores(), min.pathway.size = 10, max.pathway.size = 1000,
  n.randomizations = 10, n.internal.shuffles = 0, n.starts = 10,
  center = TRUE, batch.center = TRUE, proper.gene.names = NULL,
  verbose = 0)
}
\arguments{
\item{varinfo}{adjusted variance info from pagoda.varinfo() (or pagoda.subtract.aspect())}

\item{setenv}{environment listing gene sets (contains variables with names corresponding to gene set name, and values being vectors of gene names within each gene set)}

\item{n.components}{number of principal components to determine for each gene set}

\item{n.cores}{number of cores to use}

\item{min.pathway.size}{minimum number of observed genes that should be contained in a valid gene set}

\item{max.pathway.size}{maximum number of observed genes in a valid gene set}

\item{n.randomizations}{number of random gene sets (of the same size) to be evaluated in parallel with each gene set (can be kept at 5 or 10, but should be increased to 50-100 if the significance of pathway overdispersion will be determined relative to random gene set models)}

\item{n.internal.shuffles}{number of internal (independent row shuffles) randomizations of expression data that should be evaluated for each gene set (needed only if one is interested in gene set coherence P values, disabled by default; set to 10-30 to estimate)}

\item{n.starts}{number of random starts for the EM method in each evaluation}

\item{center}{whether the expression matrix should be recentered}

\item{batch.center}{whether batch-specific centering should be used}

\item{proper.gene.names}{alternative vector of gene names (replacing rownames(varinfo$mat)) to be used in cases when the provided setenv uses different gene names}

\item{verbose}{verbosity level}
}
\value{
a list of weighted PCA info for each valid gene set
}
\description{
For each valid gene set (having appropriate number of genes) in the provided environment (setenv),
the method will run weighted PCA analysis, along with analogous analyses of random gene sets of the
same size, or shuffled expression magnitudes for the same gene set.
}
\examples{
\donttest{
data(pollen)
cd <- pollen
cd <- cd[,colSums(cd>0)>1.8e3]
cd <- cd[rowSums(cd)>10,]
cd <- cd[rowSums(cd>0)>5,]
knn <- knn.error.models(cd, k=ncol(cd)/4, n.cores=10, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = FALSE)
# create go environment
library(org.Hs.eg.db)
ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids)
names(rids) <- ids
go.env <- eapply(org.Hs.egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))
go.env <- go.env[unlist(lapply(go.env, length))>5]
library(GO.db)
desc <- unlist(lapply(mget(names(go.env), GOTERM, ifnotfound = NA), function(x) if(is.logical(x)) { return("") } else { slot(x, "Term")}))
names(go.env) <- paste(names(go.env), desc)  # append description to the names
go.env <- list2env(go.env)  # convert to an environment
pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=10, n.internal.shuffles=50)
}
}

