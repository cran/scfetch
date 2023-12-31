% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bam.R
\name{Bam2Fastq}
\alias{Bam2Fastq}
\title{Convert bam files to fastq files.}
\usage{
Bam2Fastq(
  bam.folder = NULL,
  bam.path = NULL,
  bam.type = c("10x", "other"),
  pair.end = NULL,
  bamtofastq.path = NULL,
  bamtofastq.paras = "--nthreads 4",
  sort.name = FALSE,
  sort.thread = 4
)
}
\arguments{
\item{bam.folder}{Folder contains bam files, obtained from \code{DownloadSRA}. Default: NULL.}

\item{bam.path}{Paths of bams. \code{bam.folder} and \code{bam.path} cannot be both NULL. Default: NULL.}

\item{bam.type}{The source of bam files, choose from 10x (e.g. CellRanger) or other. Default: 10x.}

\item{pair.end}{The bam files are pair-end or single-end, used when \code{bam.type} is other. Default: NULL (identify with flag).}

\item{bamtofastq.path}{Path to 10x bamtofastq (\code{bam.type} is 10x) or samtools (\code{bam.type} is other).
Default: NULL (conduct automatic detection with bamtofastq_linux/samtools).}

\item{bamtofastq.paras}{Parameters for \code{bamtofastq.path}. Default: "--nthreads 4".}

\item{sort.name}{Logical value, whether the bam files are sorted by name, required when \code{bam.type} is other. Default: FALSE.}

\item{sort.thread}{The number of threads for bam sorting, used when \code{bam.type} is other. Default: 4.}
}
\value{
NULL or paths of failed bams.
}
\description{
Convert bam files to fastq files.
}
\examples{
\dontrun{
# need users to provide prefetch.path and bamtofastq.path
GSE138266.runs <- ExtractRun(acce = "GSE138266", platform = "GPL18573")
GSE138266.down <- DownloadBam(
  gsm.df = GSE138266.runs, prefetch.path = "/path/to/prefetch",
  out.folder = "/path/to/output"
)
GSE138266.convert <- Bam2Fastq(
  bam.folder = "/path/to/output",
  bamtofastq.path = "/path/to/bamtofastq_linux or samtools",
  bamtofastq.paras = "--nthreads 4"
)
}
}
