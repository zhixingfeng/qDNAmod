## Copyright (c) 2010, Pacific Biosciences of California, Inc.

## All rights reserved.
 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted (subject to the limitations in the
## disclaimer below) provided that the following conditions are met:
 
##     * Redistributions of source code must retain the above copyright
##        notice, this list of conditions and the following disclaimer.
 
##     * Redistributions in binary form must reproduce the above
##        copyright notice, this list of conditions and the following
##        disclaimer in the documentation and/or other materials provided
##        with the distribution.
  
##     * Neither the name of Pacific Biosciences nor the names of its
##        contributors may be used to endorse or promote products derived
##        from this software without specific prior written permission.
 
## NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
## GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
## BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
## WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
## MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS CONTRIBUTORS
## BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
## BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
## WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
## OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
## IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## ###################################################################
##
## This file manages data access to the files and their object
## model. More complicated functionality, i.e., manipulation of
## alignments, joining data together is deal with in the *Util.R
## files.
##
## R is 1-based. h5 is zero-based. The h5r package expects 1-based
## indexing. In the various files, we store index values. These are
## 0-based as one would expect - hence the necessity to add 1 to
## index values.
##
## ###################################################################

##
## Class definitions
##
setClassUnion("stringOrNull", c("character", "NULL"))
setClassUnion("H5ObjOrNull", c("H5Obj", "NULL"))
setClassUnion("matrixOrNull", c("matrix", "NULL"))

setClass("PacBioDataFile", contains = "H5File",
         representation = representation(
           version  = "stringOrNull"),
         prototype = prototype(version = NULL))

setClass("PacBioCmpH5", contains = "PacBioDataFile",
         representation = representation(
           AlnIndex = "data.frame",
           AlnGroup = "data.frame",
           RefGroup = "data.frame",
           MovieInfo = "data.frame",
           RefInfo = "data.frame",
           isSorted = "logical",
           mask = "logical"))

setClass("HasBaseAccessors")

setClass("PacBioBasH5", contains = c("PacBioDataFile", "HasBaseAccessors"), 
         representation = representation(
           baseEvents = "data.frame",
           baseCallsG = "H5ObjOrNull",
           ccsBaseEvents = "data.frame",
           ccsBaseCallsG = "H5ObjOrNull"),
         prototype = prototype(baseEvents = NULL, baseCallsG = NULL,
           ccsBaseEvents = NULL, ccsBaseCallsG = NULL))

setClass("PacBioPlsH5", contains = c("PacBioBasH5"),
         representation = representation(
           pulseEvents = "data.frame",
           pulseCallsG = "H5ObjOrNull"),
         prototype = prototype(pulseEvents = NULL, pulseCallsG = NULL))

setClass("MultiPart", contains = c("PacBioDataFile", "HasBaseAccessors"),
         representation = representation(parts = "list",
           backingClass = "character",
           holeMap = "integer"),
         prototype = prototype(backingClass = NULL))

## infrequently used classes. 
setClass("PacBioTrcH5", contains = "PacBioDataFile",
         representation = representation())

setClass("PacBioAlnH5", contains = "PacBioDataFile",
         representation = representation(
           ZMWs = "matrix",
           alignment = "H5Obj"))

##
## Generics - These generics represent functions which make sense for
## various types of h5 files, e.g., aln.h5 files have alignments as do
## cmp.h5 files.
##
setGeneric("getAlignments", function(h5Obj, ...) {
  standardGeneric("getAlignments")
})

setGeneric("getStartFrame", function(h5Obj, ...) {
  standardGeneric("getStartFrame")
})

setGeneric("getWidthInFrames", function(h5Obj, ...) {
  standardGeneric("getWidthInFrames")
})

setGeneric("getMovieName", function(h5Obj, ...) {
  standardGeneric("getMovieName")
})

setGeneric("getSNR", function(h5Obj, ...) {
  standardGeneric("getSNR")
})

setGeneric("getBpzvar", function(h5Obj, ...) {
  standardGeneric("getBpzvar")
})

setGeneric("getBpzvarw", function(h5Obj, ...) {
  standardGeneric("getBpzvarw")
})

setGeneric("getPkzvar", function(h5Obj, ...) {
  standardGeneric("getPkzvar")
})

setGeneric("getPkzvarw", function(h5Obj, ...) {
  standardGeneric("getPkzvarw")
})

setGeneric("getReadScore", function(h5Obj, ...) {
  standardGeneric("getReadScore")
})

setGeneric("getHoleNumbers", function(h5Obj, ...) {
  standardGeneric("getHoleNumbers")
})

## The QualityValues are in both basH5 and cmpH5 files.
setGeneric("getQualityValue", function(h5Obj, ...) {
  standardGeneric("getQualityValue")
})
setGeneric("getDeletionQV", function(h5Obj, ...) {
  standardGeneric("getDeletionQV")
})
setGeneric("getDeletionTag", function(h5Obj, ...) {
  standardGeneric("getDeletionTag")
})
setGeneric("getInsertionQV", function(h5Obj, ...) {
  standardGeneric("getInsertionQV")
})
setGeneric("getMergeQV", function(h5Obj, ...) {
  standardGeneric("getMergeQV")
})
setGeneric("getSubstitutionQV", function(h5Obj, ...) {
  standardGeneric("getSubstitutionQV")
})
setGeneric("getSubstitutionTag", function(h5Obj, ...) {
  standardGeneric("getSubstitutionTag")
})
setGeneric("getBaseEvents", function(h5Obj, ...) {
  standardGeneric("getBaseEvents")
})
setGeneric("getPulseEvents", function(h5Obj, ...) {
  standardGeneric("getPulseEvents")
})
setGeneric("getBaseDataset", function(h5Obj, ...) {
  standardGeneric("getBaseDataset")
})
setGeneric("getPulseDataset", function(h5Obj, ...) {
  standardGeneric("getPulseDataset")
})
setGeneric("getCCSDataset", function(h5Obj, ...) {
  standardGeneric("getCCSDataset")
})

## XXX : possibly need a check here to see if you have been called
## without file.
setMethod("initialize", "PacBioDataFile", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)

  if (h5AttributeExists(.Object, "Version")) {
    version <- getH5Attribute(.Object, "Version")[]
  } else {
    version <- ""
  }
  .Object@version <- version
  return(.Object)
})

setMethod("show", "PacBioDataFile", function(object) {
  callNextMethod()
  cat("Version:", object@version, "\n")
})

getVersion <- function(pacBioDataFile) {
  return(pacBioDataFile@version)
}

## #############################################################################
##
## Compare H5 Files
##
## #############################################################################
PacBioCmpH5 <- function(fileName) {
  obj <- new("PacBioCmpH5", fileName = fileName)

  .initMap <- function(grpName, names = c("ID", "Path")) {
    group <- getH5Group(obj, grpName)
    names(names) <- names
    dta <- data.frame(lapply(names, function(nm) {
      if(h5DatasetExists(group, nm)) {
        dset <- getH5Dataset(group, nm)[]
      } else {
        NA
      }
    }), stringsAsFactors = FALSE)
    rownames(dta) <- dta$ID
    return(dta)
  }
  AlnGroup <- .initMap("AlnGroup")
  RefGroup <- .initMap("RefGroup", c("ID", "Path", "RefInfoID"))
  RefInfo <- .initMap("RefInfo", c("ID", "FullName", "Length", "MD5"))
  MovieInfo <- .initMap("MovieInfo", c("ID", "Name", "Exp", "Run", "FrameRate"))
  
  ## add some things to these data structures.
  AlnGroup$RefGroupName <- sapply(strsplit(AlnGroup$Path, "/"), "[[", 2)
  
  if (h5DatasetExists(obj, "RefGroup/OffsetTable")) {
    refBlocks <- getH5Dataset(obj, "RefGroup/OffsetTable")[]
    if (is.null(dim(refBlocks))) {
      dim(refBlocks) <- c(1, 3)
    }
    refBlocks <- refBlocks[,-1, drop = FALSE]
    colnames(refBlocks) <- c("offsetBegin", "offsetEnd")
    refBlocks[,1] <- refBlocks[,1] + 1
    RefGroup <- cbind(RefGroup, refBlocks)
    isSorted <- TRUE
  } else {
    isSorted <- FALSE
  }
  RefGroup <- merge(RefGroup, RefInfo, by.x = "RefInfoID", by.y = "ID")
  RefGroup$Name <- gsub("^/", "", RefGroup$Path)

  obj@isSorted <- isSorted
  obj@AlnGroup <- AlnGroup
  obj@RefGroup <- RefGroup
  obj@RefInfo <- RefInfo
  obj@MovieInfo <- MovieInfo
  
  aI <- getH5Dataset(obj, "AlnInfo/AlnIndex")
  ## necessary when the dim is <= 1.
  oDims <- dim(aI)
  aI <- aI[]
  dim(aI) <- oDims
  
  obj@AlnIndex <-
    data.frame("ID"            = aI[,1],
               "alnGroupPath"  = AlnGroup[match(aI[,2], AlnGroup$ID),  "Path"],
               "movieName"     = MovieInfo[match(aI[,3], MovieInfo$ID), "Name"],
               "refName"       = RefGroup[match(aI[,4], RefGroup$ID),  "Name"],
               "fullRefName"   = RefGroup[match(aI[,4], RefGroup$ID),  "FullName"],
               "tStart"        = aI[,5],
               "tEnd"          = aI[,6],
               "alignedStrand" = aI[,7],
               "holeNumber"    = aI[,8],
               "setNumber"     = aI[,9],
               "strobeNumber"  = aI[,10],
               "moleculeID"    = aI[,11],
               "rStart"        = aI[,12],
               "rEnd"          = aI[,13],
               "mapQV"         = aI[,14],
               "nMatches"      = aI[,15],
               "nMisMatches"   = aI[,16],
               "nInsertions"   = aI[,17],
               "nDeletions"    = aI[,18],
               "offsetBegin"   = aI[,19],
               "offsetEnd"     = aI[,20],
               "nBackRead"     = aI[,21],
               "nOverlap"      = aI[,22], stringsAsFactors = FALSE)
  
  ## adjusting for 0-v-1 based.
  obj@AlnIndex[, "tStart"] <- as.integer(obj@AlnIndex[, "tStart"] + 1)
  obj@AlnIndex[, "rStart"] <- as.integer(obj@AlnIndex[, "rStart"] + 1)
  obj@AlnIndex[, "offsetBegin"] <- as.integer(obj@AlnIndex[, "offsetBegin"] + 1)

  ## Deterministically clean up all the garbage we created in building up this
  ## PacBioCmpH5 object, before returning.
  unneeded <- setdiff(ls(), c("obj"))
  rm(list=unneeded)
  gc()
  
  return(obj)
}

alnIndex  <- function(cmpH5) cmpH5@AlnIndex
alnGroup  <- function(cmpH5) cmpH5@AlnGroup
refGroup  <- function(cmpH5) cmpH5@RefGroup
refInfo   <- function(cmpH5) cmpH5@RefInfo
movieInfo <- function(cmpH5) cmpH5@MovieInfo
isSorted  <- function(cmpH5) cmpH5@isSorted

##
## This function gets the ref elt by Name, FullName, or ID.
##
.getRefElt <- function(cmpH5, refSeq, which) {
  index <- if (is.character(refSeq)) {
    tmp <- which(refSeq == refGroup(cmpH5)$Name)
    if (length(tmp) == 0) {
      which(refSeq == refGroup(cmpH5)$FullName)
    } else {
      tmp
    }
  } else if (is.numeric(refSeq)) {
    which(refSeq == refGroup(cmpH5)$ID)
  } else {
    stop("refSeq needs to be either integer (ID) or character (Name or FullName).")
  }
  if (length(index) != 1) {
    stop(paste("refSeq:", refSeq, "not found or not unique:", index))
  }
  refGroup(cmpH5)[index, which]
}
getRefPath <- function(cmpH5, refSeq) {
  .getRefElt(cmpH5, refSeq, "Path")
}
getRefName <- function(cmpH5, refSeq) {
  .getRefElt(cmpH5, refSeq, "Name")
}
getRefFullName <- function(cmpH5, refSeq) {
  .getRefElt(cmpH5, refSeq, "FullName")
}
getRefLength <- function(cmpH5, refSeq) {
  .getRefElt(cmpH5, refSeq, "Length")
}
getRefMD5 <- function(cmpH5, refSeq) {
  .getRefElt(cmpH5, refSeq, "MD5")
}
.getRefBlock <- function(cmpH5, refSeq) {
  if (! isSorted(cmpH5)) {
    stop("cmpH5 file not sorted.")
  }
  as.integer(.getRefElt(cmpH5, refSeq, c("offsetBegin", "offsetEnd")))
}

setMethod("show", "PacBioCmpH5", function(object) {
  callNextMethod(object)
  cat("N Alignments:", nrow(alnIndex(object)), "\n")
  cat("N Bases:", format(sum(getReadLength(object)), digits = 4), "\n")
  cat("N ReadGroups:", nrow(alnGroup(object)), "\n")
  cat("N RefSeqs:",    nrow(refGroup(object)), "\n")
})

setMethod("show", "MultiPart", function(object) {
  callNextMethod(object)
  cat("Backing Class:", object@backingClass, "\n")
  cat("Parts:\n")
  cat(paste("\t", names(object@parts), sep = "", collapse = "\n"), "\n")
})

setMethod("summary", "PacBioCmpH5", function(object) {
  show(object)
  print(refGroup(object))
  print(movieInfo(object))
})

setMethod("nrow", "PacBioCmpH5", function(x) {
  nrow(alnIndex(x))
})

setMethod("head", "PacBioCmpH5", function(x, ...) {
  head(alnIndex(x), ...)
})

setMethod("colnames", "PacBioCmpH5", function(x) {
  colnames(alnIndex(x))
})

setMethod("$", "PacBioCmpH5", function(x, name) {
  alnIndex(x)[,name]
})

setMethod("[", c("PacBioCmpH5", "ANY", "ANY", "ANY"), function(x, i, j) {
  alnIndex(x)[i, j, drop = FALSE]
})

setMethod("getMovieName", "PacBioCmpH5", function(h5Obj, idx = seq.int(1, nrow(h5Obj))) {
  h5Obj$movieName[idx]
})

##
## This stores the mapping from base pairs to characters.
##
.NUMBERS <- c(17,18,20,24,16,31,33,34,36,40,32,47,65,66,68,72,64,79,
              129,130,132,136,128,143,1,2,4,8,0,15,241,242,244,248,240,255)
.PAIRS   <- c("AA", "AC", "AG", "AT", "A-", "AN", "CA", "CC", "CG", "CT", "C-",
              "CN", "GA", "GC", "GG", "GT", "G-", "GN", "TA", "TC", "TG", "TT",
              "T-", "TN", "-A", "-C", "-G", "-T", "--", "-N", "NA", "NC", "NG",
              "NT", "N-", "NN")
.bMap                 <- matrix(character(256*2), ncol = 2)
colnames(.bMap)       <- c("read", "reference")
.bMap[]               <- NA
.bMap[.NUMBERS + 1, ] <- do.call(rbind, strsplit(.PAIRS, ""))

.bMapC <- function(idx) .bMap[idx + 1, ]

.geq.1.3.1 <- function(cmpH5) {
  bits <- as.integer(strsplit(getVersion(cmpH5), "\\.")[[1]][1:3])
  return(sum(bits * 100^(3:1)) >= sum(c(1, 3, 1) * 100^(3:1)))
}

.toNA <- function(cmpH5, v, naValue) {
  if (.geq.1.3.1(cmpH5)) {
    if (is.na(naValue)) {
      v
    } else {
      if (is.nan(naValue)) {
        lapply(v, function(e) ifelse(is.nan(e), NA, e))
      } else {
        lapply(v, function(e) ifelse(e == naValue, NA, e))
      }
    }
  } else {
    v
  }
}
  
.hasAlignmentDataset <- function(cmpH5, dsName) {
  ## The assumption here is that all groups either have or don't have
  ## a particular dataset.
  h5DatasetExists(cmpH5, base:::paste(cmpH5$alnGroupPath[1], dsName, sep = '/'))
}

## 
## This reads a dataset which is structured like the alignments. It
## looks overly complicated because we read the file in an organized
## fashion to avoid reading the file too much
##
.getDatasetByIdxFast <- function(cmpH5, idx = seq.int(1, nrow(cmpH5)), dsName,
                                 convert = NULL, naValue = NA) {
  if (missing(idx))
    idx <- seq.int(1, nrow(cmpH5))
  if (length(idx) == 0)
    stop("Index length must be > 0.")
  
  aln <- alnIndex(cmpH5)
  dsetNames <- base:::paste(aln$alnGroupPath[idx], dsName, sep = '/')
  
  udsets <- dsetNames[!duplicated(dsetNames)]
  dsets <- lapply(udsets, function(a) getH5Dataset(cmpH5, a))
  names(dsets) <- udsets
  dfactor <- factor(dsetNames, udsets)
  sidxs <- split(idx, dfactor)
  iidxs <- split(seq.int(length(idx)), dfactor)
  res <- vector("list", length(idx))
  
  for ( j in 1:length(dsets) ) {
    x <- read1DSlabs(dsets[[j]], aln$offsetBegin[sidxs[[j]]],
                     aln$offsetEnd[sidxs[[j]]] - aln$offsetBegin[sidxs[[j]]] + 1)
    if (! is.null(convert)) {
      x <- lapply(x, convert)
    }
    res[iidxs[[j]]] <- x
  }
  return(.toNA(cmpH5, res, naValue))
}

##
## This was the old function which still works for datasets with
## dimensions, hence it is preserved.
## 
.getDatasetByIdx <- function(cmpH5, idx = seq.int(1, nrow(cmpH5)), dsName,
                             convert = NULL, naValue = NA) {
  if (missing(idx))
    idx <- seq.int(1, nrow(cmpH5))
  if (length(idx) == 0)
    stop("Index length must be > 0.")

  
  aln <- alnIndex(cmpH5)[idx, c("alnGroupPath", "offsetBegin", "offsetEnd"),
                         drop = FALSE]

  rgp <- aln[,1]
  sqa <- seq.int(1, nrow(aln))
  rgp <- factor(rgp, rgp[!duplicated(rgp)])
  lst <- split(sqa, rgp)

  a <- lapply(lst, function(j) {
    gg <- getH5Group(cmpH5, aln[j[1], "alnGroupPath"])

    if (! h5DatasetExists(gg, dsName)) {
      stop(sprintf("dataset: %s doesn't exist.", dsName))
    }
    aa <- getH5Dataset(gg, dsName)

    if (is.null(dim(aa))) {
      v <- read1DSlabs(aa, aln$offsetBegin[j],
                       aln$offsetEnd[j] - aln$offsetBegin[j] + 1)
      if (! is.null(convert)) {
        lapply(v, convert)
      } else {
        return ( v )
      }
    } else {
      lapply(j, function(i) {
        hs <- hSlab(c(aln[i,"offsetBegin"], 1), end = c(aln[i,"offsetEnd"],
                                                  ncol(aa)))
        if (is.null(convert))
          aa[hs]
        else 
          convert(aa[hs])
      })
    }
  })
  ## make sure that the order of the results is same as the idx.
  a <- do.call(c, a)
  b <- character(length(a))
  b[do.call(c, lst)] <- names(a)
  a <- a[b]
  names(a) <- NULL

  return(.toNA(cmpH5, a, naValue))
}

.inFrames <- function(cmpH5) {
  ("FrameRate" %in% colnames(movieInfo(cmpH5))) && all(! is.na(movieInfo(cmpH5)$FrameRate))
}

.convertToSeconds <-  function(cmpH5, idx, elts) {
  if (.inFrames(cmpH5)) {
    mapply(movieInfo(cmpH5)$FrameRate[match(getMovieName(cmpH5, idx), movieInfo(cmpH5)$Name)],
           elts, FUN = function(fr, elt) {
             elt/fr
           }, SIMPLIFY = FALSE)
  } else {
    elts
  }
}

.getFrameBasedDataset <- function(cmpH5, idx, name, unit = c("seconds", "frames"),
                                  naValue = 0x7fffffff) {
  d <- .getDatasetByIdxFast(cmpH5, idx, name, naValue = naValue)
  if (match.arg(unit) == "seconds") {
    .convertToSeconds(cmpH5, idx, d)
  } else {
    d
  }
}

## XXX : IRanges has a generic called 'restrict' -- my preferred name
## -- and I don't want to deal with overwriting it and load order,
## etc.
restrictCmpH5 <- function(cmpH5, mask = rep(TRUE, nrow(cmpH5))) {
  ## A mask is used instead of the more common permutation (idx) because I
  ## want to preserve order.
  if (any(! is.finite(mask)) || length(mask) != nrow(cmpH5)) {
    stop("mask must not contain missing values and, nrow(cmpH5) == length(mask), must hold.")
  }
  nalnIdx <- cmpH5@AlnIndex[mask, ]
  nrefGroup <- refGroup(cmpH5)
  ntbl <- table(factor(nalnIdx[,"refName"],
                       levels = nrefGroup$Name[order(nrefGroup$offsetBegin)]))

  ## The key thing to be careful about is that the offsets aren't in
  ## any particular order.
  csum <- cumsum(ntbl)
  noffsets <- cbind(c(1, csum[-length(csum)] + 1), csum)
  rownames(noffsets) <- names(ntbl)
  nrefGroup[match(rownames(noffsets),
                  nrefGroup$Name), c("offsetBegin", "offsetEnd")] <- noffsets[,c(1,2)]
  ## And you can lose references.
  nrefGroup <- nrefGroup[nrefGroup[,"offsetEnd"]-nrefGroup[,"offsetBegin"] > 0,]

  cmpH5@AlnIndex <- nalnIdx
  cmpH5@RefGroup <- nrefGroup
  cmpH5@mask <- mask
  
  ## the other thing that I could consider doing is consider fixing up
  ## the nBack* columns, but it just reverts to linear search (it
  ## should unless there are bugs.)
  return(cmpH5)
}

getIPD <- function(cmpH5, idx = seq.int(1, nrow(cmpH5)), unit = c("seconds", "frames")) {
  .getFrameBasedDataset(cmpH5, idx, "IPD", unit, naValue = 0xffff)
}

getPulseWidth <- function(cmpH5, idx = seq.int(1, nrow(cmpH5)), unit = c("seconds", "frames")) {
  .getFrameBasedDataset(cmpH5, idx, "PulseWidth", unit, naValue = 0xffff)
}

getStartTime <- function(cmpH5, idx = seq.int(1, nrow(cmpH5)), unit = c("seconds", "frames")) {
  if (.hasAlignmentDataset(cmpH5, "StartTime")) 
    .getDatasetByIdxFast(cmpH5, idx, "StartTime", naValue = 0x7fffffff)
  else if (.hasAlignmentDataset(cmpH5, "StartFrame")) {
    .getFrameBasedDataset(cmpH5, idx, "StartFrame", unit)
  } else {
    stop("No appropriate dataset to compute StartTime")
  }
}

setMethod("getStartFrame", "PacBioCmpH5",
          function(h5Obj, idx = seq.int(1, nrow(h5Obj))) {
            .getDatasetByIdxFast(h5Obj, idx, "StartFrame", naValue = 0x7fffffff)
          })

getCumulativeAdvanceTime <- function(cmpH5, idx = seq.int(1, nrow(cmpH5)),
                                     unit = c("seconds", "frames")) {
  mapply(getIPD(cmpH5, idx, unit), getPulseWidth(cmpH5, idx, unit),
         FUN = function(ipd, pw) {
           v <- ipd + pw
           v <- cumsum(ifelse(is.na(v), 0, v))
           ## it actually makes a difference if you use pw versus ipd;
           ## unfortunately, the first IPD is sometimes NA - this is a bug in
           ## loadPulses.
           v[is.na(pw)] <- NA
           v
         }, SIMPLIFY = FALSE)
}

getPkmid <- function(cmpH5, idx) .getDatasetByIdxFast(cmpH5, idx, "pkmid")
getAlignmentsRaw <- function(cmpH5, idx) .getDatasetByIdxFast(cmpH5, idx, "AlnArray")

setMethod("getAlignments", "PacBioCmpH5", function(h5Obj, idx) {
  .getDatasetByIdxFast(h5Obj, idx, "AlnArray", convert = .bMapC)
})
setMethod("getQualityValue", "PacBioCmpH5", function(h5Obj, idx) {
  .getDatasetByIdxFast(h5Obj, idx, "QualityValue", naValue = 255)
})
setMethod("getDeletionQV", "PacBioCmpH5", function(h5Obj, idx) {
  .getDatasetByIdxFast(h5Obj, idx, "DeletionQV", naValue = 255)
})
setMethod("getDeletionTag", "PacBioCmpH5", function(h5Obj, idx) {
  .getDatasetByIdxFast(h5Obj, idx, "DeletionTag")
})
setMethod("getInsertionQV", "PacBioCmpH5", function(h5Obj, idx) {
  .getDatasetByIdxFast(h5Obj, idx, "InsertionQV", naValue = 255)
})
setMethod("getMergeQV", "PacBioCmpH5", function(h5Obj, idx) {
  .getDatasetByIdxFast(h5Obj, idx, "MergeQV", naValue = 255)
})
setMethod("getSubstitutionQV", "PacBioCmpH5", function(h5Obj, idx) {
  .getDatasetByIdxFast(h5Obj, idx, "SubstitutionQV", naValue = 255)
})
setMethod("getSubstitutionTag", "PacBioCmpH5", function(h5Obj, idx) {
  .getDatasetByIdxFast(h5Obj, idx, "SubstitutionTag")
})

## #############################################################################
##
## Base H5 Files
##
## #############################################################################
.ofBasecalls <- c("Basecall", "DeletionQV", "DeletionTag",
                  "InsertionQV", "PreBaseFrames", "PulseIndex",
                  "QualityValue", "SubstitutionQV", "SubstitutionTag",
                  "WidthInFrames")

setMethod("getMovieName", "PacBioBasH5", function(h5Obj) {
  getH5Attribute(getH5Group(h5Obj, "ScanData/RunInfo"), "MovieName")[]
})

setMethod("getMovieName", "MultiPart", function(h5Obj) {
  getMovieName(part(h5Obj))
})

.initBasEvents <- function(.Object, gname) {
  baseCalls <- getH5Group(.Object, gname)
  bcZMW <- do.call(cbind, lapply(c("HoleNumber", "NumEvent", "HoleXY", "HoleStatus"), function(nm) {
    getH5Dataset(getH5Group(baseCalls, "ZMW"), nm)[]
  }))
  bcZMW <- cbind(bcZMW, "offset" = cumsum(c(1, bcZMW[-nrow(bcZMW), 2])))
  colnames(bcZMW) <- c("holeNumber", "numEvent", "x", "y", "holeStatus", "offset")
  status <- getH5Attribute(getH5Dataset(getH5Group(baseCalls, "ZMW"), "HoleStatus"),
                           "LookupTable")[]
  bcZMW <- as.data.frame(bcZMW)
  bcZMW$holeStatus <- status[match(bcZMW$holeStatus + 1, seq_along(status))]
  list(events = bcZMW, group = baseCalls)
}

setMethod("initialize", "PacBioBasH5", function(.Object, fileName = NULL) {
  .Object <- callNextMethod(.Object, fileName = fileName)

  l <- .initBasEvents(.Object, "PulseData/BaseCalls")
  .Object@baseEvents  <- l$events
  .Object@baseCallsG  <- l$group
  
  if (h5GroupExists(.Object, "PulseData/ConsensusBaseCalls")) {
    ## XXX : These should always be there, but old files might not
    ## have them - really I'm just delaying the breaking to the call
    ## which needs them.
    l <- .initBasEvents(.Object, "PulseData/ConsensusBaseCalls")
    .Object@ccsBaseEvents  <- l$events
    .Object@ccsBaseCallsG  <- l$group
  }
  return(.Object)
})

isMultiPart <- function(fileName) {
  h5GroupExists(H5File(fileName), 'MultiPart')
}

part <- function(multiPart, i = 1) {
  multiPart@parts[[i]]
}
parts <- function(multiPart) {
  multiPart@parts
}

.MultiPart <- function(fileName, klassName = c("PacBioBasH5", "PacBioPlsH5")) {
  .Object <- new("MultiPart", fileName = fileName)
  .Object@backingClass = match.arg(klassName)

  parts <- getH5Dataset(.Object, "MultiPart/Parts")[]
  names(parts) <- parts
  .Object@parts <- lapply(parts, function(n) {
    k <- call(klassName)
    k[[2]] <- gsub(basename(fileName), n, fileName, fixed = TRUE)
    eval(k)
  })
  
  ## holeLookup is going to be defined as a vector mapping v[hn] =>
  ## part. 
  holeLookup <- getH5Dataset(.Object, "MultiPart/HoleLookup")[]
  m <- as.integer(rep(NA, max(holeLookup[,1])))
  m[holeLookup[,1] + 1] = holeLookup[,2] 
  .Object@holeMap <- m

  return(.Object)
}

.getParts <- function(multiPart, holeNumbers) {
  multiPart@holeMap[holeNumbers + 1]
}

PacBioBasH5 <- function(fileName) {
  if (isMultiPart(fileName)) {
    .MultiPart(fileName, "PacBioBasH5")
  }
  else {
    new("PacBioBasH5", fileName = fileName)
  }
}

## #############################################################################
##
## Pulse H5 Files
##
## #############################################################################
.ofPulses <- c("Channel", "Chi2", "ClassifierQV", "IsPulse", "MaxSignal",
               "MeanSignal", "MidSignal", "MidStdDev",
               "StartFrame", "WidthInFrames")

PacBioPlsH5 <- function(fileName) {
  if (isMultiPart(fileName)) {
    .MultiPart(fileName, "PacBioPlsH5")
  } else {
    .Object <- new("PacBioPlsH5", fileName = fileName)
    l <- .initBasEvents(.Object, "PulseData/PulseCalls")
    .Object@pulseEvents <- l$events
    .Object@pulseCallsG <- l$group
    return(.Object)
  }
}

setMethod("getPulseEvents", "PacBioPlsH5", function(h5Obj) {
  h5Obj@pulseEvents
})
setMethod("getPulseEvents", "MultiPart", function(h5Obj) {
  do.call(rbind, lapply(parts(h5Obj), getPulseEvents))
})

setMethod("getBaseEvents", "PacBioBasH5", function(h5Obj) {
  h5Obj@baseEvents
})
setMethod("getBaseEvents", "MultiPart", function(h5Obj) {
  do.call(rbind, lapply(parts(h5Obj), getBaseEvents))
})

setGeneric("getCCSEvents", function(h5Obj, ...) {
  standardGeneric("getCCSEvents")
})
setMethod("getCCSEvents", "PacBioBasH5", function(h5Obj) {
  h5Obj@ccsBaseEvents
})
setMethod("getCCSEvents", "MultiPart", function(h5Obj) {
  do.call(rbind, lapply(parts(h5Obj), getCCSEvents))
})

setMethod("getHoleNumbers", "PacBioBasH5", function(h5Obj) {
  h5Obj@baseEvents[, "holeNumber"]
})
setMethod("getHoleNumbers", "MultiPart", function(h5Obj) {
  which(!is.na(h5Obj@holeMap)) - 1
})

setMethod("getBaseDataset", "PacBioBasH5",
          function(h5Obj, dsName, holeNumbers = getHoleNumbers(h5Obj),
                   convert = NULL, ...) {
            .getFromBasLikeFile(h5Obj@baseCallsG, getBaseEvents(h5Obj),
                                dsName, holeNumbers, convert, ...)
          })

setMethod("getBaseDataset", "MultiPart",
          function(h5Obj, dsName, holeNumbers = getHoleNumbers(h5Obj),
                   convert = NULL, ...) {
            .getFromMultiPartBasLikeFile(h5Obj, dsName, holeNumbers, convert,
                                         getBaseDataset, ...)
          })

setMethod("getPulseDataset", "PacBioPlsH5",
          function(h5Obj, dsName, holeNumbers = getHoleNumbers(h5Obj),
                   convert = NULL, ...) {
            .getFromBasLikeFile(h5Obj@pulseCallsG, getPulseEvents(h5Obj),
                                dsName, holeNumbers, convert, ...)
          })

setMethod("getPulseDataset", "MultiPart",
          function(h5Obj, dsName, holeNumbers = getHoleNumbers(h5Obj),
                   convert = NULL, ...) {
            .getFromMultiPartBasLikeFile(h5Obj, dsName, holeNumbers, convert,
                                         getPulseDataset, ...)
          })

setMethod("getCCSDataset", "PacBioBasH5",
          function(h5Obj, dsName, holeNumbers = getHoleNumbers(h5Obj),
                   convert = NULL, ...) {
            .getFromBasLikeFile(h5Obj@ccsBaseCallsG, getCCSEvents(h5Obj),
                                dsName, holeNumbers, convert, ...)
          })

setMethod("getCCSDataset", "MultiPart",
          function(h5Obj, dsName, holeNumbers = getHoleNumbers(h5Obj),
                   convert = NULL, ...) {
            .getFromMultiPartBasLikeFile(h5Obj, dsName, holeNumbers,
                                         convert, getCCSDataset, ...)
          })

.getFromMultiPartBasLikeFile <- function(h5Obj, dsName, holeNumbers, convert, method,
                                         leftPads = rep(0, length(holeNumbers)),
                                         rightPads = rep(0, length(holeNumbers))) {
  partIdx <- split(seq.int(1, length(holeNumbers)),
                   factor(.getParts(h5Obj, holeNumbers),
                          1:length(h5Obj@parts)))
  res <- vector("list", length(holeNumbers))
  names(res) <- holeNumbers
  for (i in seq_along(partIdx)) {
    h <- holeNumbers[partIdx[[i]]]
    if (length(h) > 0) {
      res[match(h, holeNumbers)] <- method(part(h5Obj, i), dsName, holeNumbers = h,
                                           convert = convert,
                                           leftPads = leftPads[partIdx[[i]]],
                                           rightPads = rightPads[partIdx[[i]]])
    }
  }
  return(res)
}

.getFromBasLikeFile <- function(group, events, dsName, holeNumbers,
                                convert = NULL,
                                leftPads = rep(0, length(holeNumbers)),
                                rightPads = rep(0, length(holeNumbers))) {
  dset <- getH5Dataset(group, dsName)
  rows <- match(holeNumbers, events$holeNumber)
  evts <- as.matrix(events[, c("offset", "numEvent"), drop = FALSE])
  if (any(is.na(rows)) || length(rows) <= 0)
    stop("Invalid ZMW: ", paste(holeNumbers, collapse = ", "))

  if (any((leftPads + rightPads) > evts[rows, 2])) 
    stop("rightPads + leftPads must not be larger than the number of ZMW events")
    
  if (is.null(dim(dset))) {
    offset <- evts[rows,1,drop=FALSE] + leftPads
    nevts  <- evts[rows,2,drop=FALSE] - rightPads - leftPads 
    dta <- read1DSlabs(dset, offset, nevts)
    if (! is.null(convert)) {
      dta <- lapply(dta, convert)
    } 
  } else {
    evts <- evts[rows,,drop = FALSE]
    dta <- lapply(1:nrow(evts), function(i) {
      if (evts[i,2] <= 0) {
        return(NULL)
      }
      rng <- (evts[i,1] + leftPads[[i]]):(evts[i,1] + evts[i,2] - rightPads[[i]] - 1)
      d <- dset[rng,]
      if (! is.null(convert)) {
        d <- convert(d)
      }
      d
    })
  }
  names(dta) <- holeNumbers
  return(dta)
}

.convertToLettersFromAscii <- function(a) {
  m <- c(65, 67, 71, 84, 45, 97, 99, 103, 116)
  b <- c("A", "C", "G", "T", "-", "A", "C", "G", "T")
  b[match(a, m)]
}

getBasecalls <- function(basH5, holeNumbers = getHoleNumbers(basH5), ...)
  getBaseDataset(basH5, "Basecall", holeNumbers, .convertToLettersFromAscii, ...)

getPreBaseFrames <- function(basH5, holeNumbers = getHoleNumbers(basH5), ...)
  getBaseDataset(basH5, "PreBaseFrames", holeNumbers, ...)

getPulseIndex <- function(basH5, holeNumbers = getHoleNumbers(basH5), ...)
  getBaseDataset(basH5, "PulseIndex", holeNumbers, ...)

setMethod("getQualityValue", "HasBaseAccessors",
          function(h5Obj, holeNumbers = getHoleNumbers(h5Obj), ...) {
            getBaseDataset(h5Obj, "QualityValue", holeNumbers, ...)
          })
setMethod("getDeletionQV", "HasBaseAccessors",
          function(h5Obj, holeNumbers = getHoleNumbers(h5Obj), ...) {
            getBaseDataset(h5Obj, "DeletionQV", holeNumbers, ...)
          })
setMethod("getDeletionTag", "HasBaseAccessors",
          function(h5Obj, holeNumbers = getHoleNumbers(h5Obj), ...) {
            getBaseDataset(h5Obj, "DeletionTag", holeNumbers, convert = .convertToLettersFromAscii, ...)
          })
setMethod("getInsertionQV", "HasBaseAccessors",
          function(h5Obj, holeNumbers = getHoleNumbers(h5Obj), ...) {
            getBaseDataset(h5Obj, "InsertionQV", holeNumbers, ...)
          })
setMethod("getMergeQV", "HasBaseAccessors",
          function(h5Obj, holeNumbers = getHoleNumbers(h5Obj), ...) {
            getBaseDataset(h5Obj, "MergeQV", holeNumbers, ...)
          })
setMethod("getSubstitutionQV", "HasBaseAccessors",
          function(h5Obj, holeNumbers = getHoleNumbers(h5Obj), ...) {
            getBaseDataset(h5Obj, "SubstitutionQV", holeNumbers, ...)
          })
setMethod("getSubstitutionTag", "HasBaseAccessors",
          function(h5Obj, holeNumbers = getHoleNumbers(h5Obj), ...) {
            getBaseDataset(h5Obj, "SubstitutionTag", holeNumbers,
                           convert = .convertToLettersFromAscii, ...)
          })

getChannel <- function(plsH5, holeNumbers = getHoleNumbers(plsH5), ...)
  getPulseDataset(plsH5, "Channel", holeNumbers, ...)
getChi2 <- function(plsH5, holeNumbers = getHoleNumbers(plsH5), ...)
  getPulseDataset(plsH5, "Chi2", holeNumbers, ...)
getIsPulse <- function(plsH5, holeNumbers = getHoleNumbers(plsH5), ...)
  getPulseDataset(plsH5, "IsPulse", holeNumbers, ...)
getLabelQV <- function(plsH5, holeNumbers = getHoleNumbers(plsH5), ...)
  getPulseDataset(plsH5, "LabelQV", holeNumbers, ...)
getMaxSignal <- function(plsH5, holeNumbers = getHoleNumbers(plsH5), ...)
  getPulseDataset(plsH5, "MaxSignal", holeNumbers, ...)
getMeanSignal <- function(plsH5, holeNumbers = getHoleNumbers(plsH5), ...)
  getPulseDataset(plsH5, "MeanSignal", holeNumbers, ...)
getMidSignal <- function(plsH5, holeNumbers = getHoleNumbers(plsH5), ...)
  getPulseDataset(plsH5, "MidSignal", holeNumbers, ...)
getMidStdDev <- function(plsH5, holeNumbers = getHoleNumbers(plsH5), ...)
  getPulseDataset(plsH5, "MidStdDev", holeNumbers, ...)


## XXX : this dataset exists in the PulseCalls as well. If a user
## wants it, they'll have to use: getPulseDataset( .. ).
setMethod("getWidthInFrames", "HasBaseAccessors",
          function(h5Obj, holeNumbers = getHoleNumbers(basH5), ...) {
            getBaseDataset(h5Obj, "WidthInFrames", holeNumbers, ...)
          })

setMethod("getStartFrame", "HasBaseAccessors",
          function(h5Obj, holeNumbers = getHoleNumbers(h5Obj), ...) {
            getPulseDataset(h5Obj, "StartFrame", holeNumbers, ...)
          })

getCCSBasecalls <- function(basH5, holeNumbers = getHoleNumbers(basH5), ...)
  getCCSDataset(basH5, "Basecall", holeNumbers = holeNumbers,
                convert = .convertToLettersFromAscii, ...)


## #############################################################################
##
## Align H5 Files
##
## #############################################################################
PacBioAlnH5 <- function(fileName, whAlignment = c("Alignments", "CRFAlignments",
                                    "GlobalAlignments")) {
  obj <- new("PacBioAlnH5", fileName = fileName)
  whAlignment <- match.arg(whAlignment)
    
  ## initialize groups.
  alignments <- getH5Group(obj, whAlignment)
  
  ZMW <- sapply(c("HoleNumber", "ReadLength", "ReadStartBase",
                  "ReadString.Index", "TemplateLength", "TemplateStartBase",
                  "TemplateString.Index"), function(nm) {
                    getH5Dataset(alignments, nm)[]
                  })
  cn <- colnames(ZMW)
  ZMW <- cbind(ZMW, getH5Dataset(alignments, "Accuracy")[])
  colnames(ZMW) <- c(cn, "Accuracy")

  ZMW <- ZMW[order(ZMW[,"HoleNumber"]), ]
    
  obj@ZMWs <- ZMW 
  obj@alignment <- alignments
  return(obj)
}

getZMWs <- function(alnH5) {
  alnH5@ZMWs
}

setMethod("nrow", "PacBioAlnH5", function(x) {
  nrow(x@ZMWs)
})

.getAlnH5String <- function(alnH5, zmws, wh = c("ReadString", "TemplateString"),
                            convert = .convertToLettersFromAscii) {
  wh   <- match.arg(wh)
  whS  <- paste(wh, "Index", sep = ".")
  cs   <- cbind(start  = cumsum(c(1, getZMWs(alnH5)[-nrow(alnH5), whS])),
                offset = getZMWs(alnH5)[, whS])

  if (missing(zmws))
    zmws <- getZMWs(alnH5)[, "HoleNumber"]
  
  mIdx <- match(zmws, getZMWs(alnH5)[,"HoleNumber"])
  cs <- cs[mIdx,]
  rownames(cs) <- zmws
  cs <- cs[cs[,2] > 0, ]
  cs[,2] <- cs[,2]
  
  ds <- getH5Dataset(alnH5@alignment, wh)
  
  apply(cs, 1, function(slc) {
    convert(ds[hSlab(slc[1], width = slc[2])])
  })
}

setMethod("getAlignments", "PacBioAlnH5", function(h5Obj, zmws, ...) {
  mapply(function(r, t) {
    cbind("read" = r, "reference" = t)
  }, .getAlnH5String(h5Obj, zmws = zmws, wh = "ReadString", ...),
         .getAlnH5String(h5Obj, zmws = zmws, wh = "TemplateString", ...),
         SIMPLIFY = FALSE)
})

## #############################################################################
##
## Trace H5 Files
##
## #############################################################################
PacBioTrcH5 <- function(fileName) {
  obj <- new("PacBioTrcH5", fileName = fileName)  
  return(obj)
}

getTraces <- function(trcH5) {
  getH5Dataset(trcH5, "TraceData/Traces", inMemory = FALSE)
}

