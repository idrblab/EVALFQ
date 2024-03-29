% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evalfqspikedpart.R
\name{lfqspikepart}
\alias{lfqspikepart}
\title{Conduct LFQ and assess performance for part of LFQ workflows by collectively considering the spiked proteins.}
\usage{
lfqspikepart(data_s, spiked, selectFile,
Ca="1", Cb="1", Cc="1", Cd="1", Ce="1")
}
\arguments{
\item{data_s}{This input file should be numeric type except the first and second column containing the names and label (control or case) of the studied samples, respectively. The intensity data should be provided in this input file with the following order: samples in row and proteins/peptides in column. Missing value (NA) of protein intensity are allowed.}

\item{spiked}{The file should provide the concentrations of known proteins (such as spiked proteins). This file is required, if the user want to conduct assessment using criteria (e) This file should contain the class of samples and the Sample ID. The Sample ID should be unique and defined by the preference of EVALFQ users, and the class of samples refers to the group of Sample ID. The ID of the spiked proteins should be consistent in both “data_s" and "spiked”. Detail information are described in the online "Example".}

\item{selectFile}{Input the name of your prefered strategies. Sample data of this data type is in the working directory (in github) “idrblab/EVALFQ/data/selectworkflows.rda”.}

\item{Ca}{Criterion (a): precision of LFQ based on the proteomes among replicates (Proteomics. 15:3140-51, 2015). If set 1, the user chooses to assess LFQ workflows using Criterion (a). If set 0, the user excludes Criterion (a) from performance assessment. The default setting of this value is “1”.}

\item{Cb}{Criterion (b): classification ability of LFQ between distinct sample groups (Nat Biotechnol. 28:83-9, 2010). If set 1, the user chooses to assess LFQ workflows using Criterion (b). If set 0, the user excludes Criterion (b) from performance assessment. The default setting of this value is “1”.}

\item{Cc}{Criterion (c): differential expression analysis by reproducibility-optimization (Nat Biotechnol. 32:896-902, 2014). If set 1, the user chooses to assess LFQ workflows using Criterion (c). If set 0, the user excludes Criterion (c) from performance assessment. The default setting of this value is “1”.}

\item{Cd}{Criterion (d): reproducibility of the identified protein markers among different datasets (Mol Biosyst. 11:1235-40, 2015). If set 1, the user chooses to assess LFQ workflows using Criterion (d). If set 0, the user excludes Criterion (d) from performance assessment. The default setting of this value is “1”.}

\item{Ce}{Criterion (e): accuracy of LFQ based on spiked and background proteins (Nat Biotechnol. 34:1130-6, 2016). If set 1, the user chooses to assess LFQ workflows using Criterion (e). If set 0, the user excludes Criterion (e) from performance assessment. The default setting of this value is “1”.}
}
\value{
preprocessed spiked matrix
}
\description{
The EVALFQ enables the label-free quantification of proteomic
data and the performance assessment of each LFQ workflow from multiple perspectives. This
tool function gives the results based on each LFQ workflow and the criteria preferred and selected by the users.
For function definitions and descriptions please use "??EVALFQ" command in R.
}
\examples{
library(EVALFQ)
\donttest{data_s <- read.csv(file = "ANPELA-Spiked.INPUT.txt", sep="\t", header = TRUE, row.names = 1)}
\donttest{lfqspikepart(data_s = data_s,
spiked = "ANPELA-Spiked.Protein.txt", selectFile = "selectedmethods.csv",
assum_a="Y", assum_b="Y", assum_c="Y", Ca="1", Cb="1", Cc="1", Cd="1", Ce="1")}
}
