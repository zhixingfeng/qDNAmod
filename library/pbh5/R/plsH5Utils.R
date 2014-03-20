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

setGeneric("getRegionsTable", function(h5Obj, ...) {
  standardGeneric("getRegionsTable")
})
setGeneric("getAcqParams", function(h5Obj, ...) {
  standardGeneric("getAcqParams")
})

setGeneric("getBaselineSigma", function(h5Obj, ...) {
  standardGeneric("getBaselineSigma")
})
setGeneric("getBaselineLevel", function(h5Obj, ...) {
  standardGeneric("getBaselineLevel")
})

setMethod("getRegionsTable", "PacBioBasH5", function(h5Obj) {
  tbl           <- getH5Dataset(h5Obj, "PulseData/Regions")
  regionTypes   <- getH5Attribute(tbl, "RegionTypes")[]
  regionsDescs  <- getH5Attribute(tbl, "RegionDescriptions")[]
  regionSources <- getH5Attribute(tbl, "RegionSources")[]
  columnNames   <- getH5Attribute(tbl, "ColumnNames")[]
  
  mat <- as.data.frame(tbl[])
  colnames(mat) <- c("holeNumber", "type", "start", "end", "score")

  mat$type  <- regionTypes[mat[,"type"] + 1]
  mat$start <- mat$start + 1
  mat[mat$type != "GlobalAccuracy", ]
})
setMethod("getRegionsTable", "MultiPart", function(h5Obj) {
  do.call(rbind, lapply(parts(h5Obj), getRegionsTable))
})

setMethod("getAcqParams", "PacBioBasH5", function(h5Obj) {
  names <- listH5Contents(g <- getH5Group(h5Obj, "ScanData/AcqParams"))$"."$attributes
  s <- sapply(names, function(a) {
    getH5Attribute(g, a)[]
  })
  s$MovieTime <- s$NumFrames / s$FrameRate
  return(s)
})
setMethod("getAcqParams", "MultiPart", function(h5Obj) {
  getAcqParams(part(h5Obj))
})

getNumEvent <- function(h5Obj,  holeNumbers = getHoleNumbers(h5Obj)) {
  getBaseEvents(h5Obj)[match(holeNumbers, getHoleNumbers(h5Obj))]
}

setGeneric(".getHoleDataset", function(h5Obj, ...) {
  standardGeneric(".getHoleDataset")
})
setGeneric(".getDatasetFromBasLike", function(h5Obj, ...) {
  standardGeneric(".getDatasetFromBasLike")
})
setMethod(".getDatasetFromBasLike", "PacBioBasH5", function(h5Obj, name) {
  getH5Dataset(h5Obj, name)
})
setMethod(".getDatasetFromBasLike", "MultiPart", function(h5Obj, name) {
  getH5Dataset(h5Obj@parts[[1]], name)
})

setGeneric(".getGroupFromBasLike", function(h5Obj, ...) {
  standardGeneric(".getGroupFromBasLike")
})
setMethod(".getGroupFromBasLike", "PacBioBasH5", function(h5Obj, name) {
  getH5Group(h5Obj, name)
})
setMethod(".getGroupFromBasLike", "MultiPart", function(h5Obj, name) {
  getH5Group(h5Obj@parts[[1]], name)
})


setMethod(".getHoleDataset", "PacBioBasH5", function(h5Obj, name, holeNumbers) {
  stopifnot(h5DatasetExists(h5Obj, name))
  getH5Dataset(h5Obj, name)[match(holeNumbers, getHoleNumbers(h5Obj))]
})
setMethod(".getHoleDataset", "MultiPart", function(h5Obj, name, holeNumbers) {
  byPart <- split(holeNumbers, factor(.getParts(h5Obj, holeNumbers), 1:length(h5Obj@parts)))

  ## determine the dimensions of the data.
  isV <- is.null(bigd <- dim(getH5Dataset(part(h5Obj), name)))
  if (isV) {
    res <- vector("numeric", length(holeNumbers))
  } else {
    ## For 3-dimensional data, return it as 3-d array.
    if (length(bigd) == 2) {
        res <- matrix(NA, nrow = length(holeNumbers), ncol = bigd[2])
    } else {
        res <- array(NA, dim = c(length(holeNumbers), bigd[3], bigd[2]))
    }
  }
  for (i in seq_along(byPart)) {
    h <- byPart[[i]]
    if (length(h) > 0) {
      rm <- match(h, getHoleNumbers(part(h5Obj, i)))
      lm <- match(h, holeNumbers)
      if (isV) {
        res[lm] <- getH5Dataset(part(h5Obj, i), name)[rm]
      } else {
        if (length(bigd) == 2) {
          res[lm,] <- getH5Dataset(part(h5Obj, i), name)[rm,]
        } else {
          for (d in 1:bigd[2]) {
            res[lm,,d] <- getH5Dataset(part(h5Obj, i), name)[rm,d,]    
          }
        }
      }
    }
  }
  return(res)
})
          
setMethod("getReadScore", "HasBaseAccessors", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getHoleDataset(h5Obj, "PulseData/BaseCalls/ZMWMetrics/ReadScore", holeNumbers)
})

setGeneric("getChannelToBaseMap", function(h5Obj) {
  standardGeneric("getChannelToBaseMap")
})
setMethod("getChannelToBaseMap", "PacBioBasH5", function(h5Obj) {
  sapply(0:3, function(d) {
    getH5Attribute(getH5Group(h5Obj, sprintf("ScanData/DyeSet/Analog[%d]", d)), "Base")[]
  })
})
setMethod("getChannelToBaseMap", "MultiPart", function(h5Obj) {
  sapply(0:3, function(d) {
    getH5Attribute(getH5Group(part(h5Obj), sprintf("ScanData/DyeSet/Analog[%d]", d)), "Base")[]
  })
})

.getDatasetWithChannelCols <- function(h5Obj, name, holeNumbers) {
  d <- .getHoleDataset(h5Obj, name, holeNumbers)
  if (is.null(dim(d))) {
    stopifnot(length(d) == 4)
    dim(d) <- c(1, 4)
  }
  colnames(d) <- getChannelToBaseMap(h5Obj)
  d
}

.getBpzvarw <- function(basH5, holeNumbers = getHoleNumbers(basH5)) {
  d <- .getHoleDataset(basH5, "PulseData/BaseCalls/ZMWMetrics/HQRegionBpzvarw", holeNumbers)
  g <- .getGroupFromBasLike(basH5, "PulseData/BaseCalls/ZMWMetrics")
  zvarwidthcuts <- getH5Attribute(g, "ZvarWidthCuts")[]
  bpzvarw <- lapply(1:length(zvarwidthcuts), function(z) {
    df <- as.data.frame(d[1:length(holeNumbers),,z])
    colnames(df) <- getChannelToBaseMap(basH5)
    df
  }) 
  names(bpzvarw) <- zvarwidthcuts
  return(bpzvarw)
}

.getPkzvarw <- function(basH5, holeNumbers = getHoleNumbers(basH5)) {
  d <- .getHoleDataset(basH5, "PulseData/BaseCalls/ZMWMetrics/HQRegionPkzvarw", holeNumbers)
  g <- .getGroupFromBasLike(basH5, "PulseData/BaseCalls/ZMWMetrics")
  zvarwidthcuts <- getH5Attribute(g, "ZvarWidthCuts")[]
  pkzvarw <- lapply(1:length(zvarwidthcuts), function(z) {
    df <- as.data.frame(d[1:length(holeNumbers),,z])
    colnames(df) <- getChannelToBaseMap(basH5)
  }) 
  names(pkzvarw) <- zvarwidthcuts
  return(pkzvarw)
}

setMethod("getBpzvarw", "PacBioBasH5", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getHoleDataset(h5Obj, "PulseData/BaseCalls/ZMWMetrics/HQRegionBpzvarw", holeNumbers)
})
setMethod("getBpzvarw", "MultiPart", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getBpzvarw(h5Obj, holeNumbers)
})

setMethod("getPkzvarw", "PacBioBasH5", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getHoleDataset(h5Obj, "PulseData/BaseCalls/ZMWMetrics/HQRegionPkzvarw", holeNumbers)
})
setMethod("getPkzvarw", "MultiPart", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getPkzvarw(h5Obj, holeNumbers)
})

setMethod("getBpzvar", "PacBioBasH5", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getDatasetWithChannelCols(h5Obj, "PulseData/BaseCalls/ZMWMetrics/HQRegionBpzvar", holeNumbers)
})
setMethod("getBpzvar", "MultiPart", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getDatasetWithChannelCols(h5Obj, "PulseData/BaseCalls/ZMWMetrics/HQRegionBpzvar", holeNumbers)
})

setMethod("getPkzvar", "PacBioBasH5", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getDatasetWithChannelCols(h5Obj, "PulseData/BaseCalls/ZMWMetrics/HQRegionPkzvar", holeNumbers)
})
setMethod("getPkzvar", "MultiPart", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getDatasetWithChannelCols(h5Obj, "PulseData/BaseCalls/ZMWMetrics/HQRegionPkzvar", holeNumbers)
})

setMethod("getSNR", "PacBioBasH5", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getDatasetWithChannelCols(h5Obj, "PulseData/BaseCalls/ZMWMetrics/HQRegionSNR", holeNumbers)
})
setMethod("getSNR", "MultiPart", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getDatasetWithChannelCols(h5Obj, "PulseData/BaseCalls/ZMWMetrics/HQRegionSNR", holeNumbers)
})

setMethod("getBaselineSigma", "PacBioPlsH5", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getDatasetWithChannelCols(h5Obj, "PulseData/PulseCalls/ZMW/BaselineSigma", holeNumbers)
})
setMethod("getBaselineSigma", "MultiPart", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getDatasetWithChannelCols(h5Obj, "PulseData/PulseCalls/ZMW/BaselineSigma", holeNumbers)
})

setMethod("getBaselineLevel", "PacBioPlsH5", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getDatasetWithChannelCols(h5Obj, "PulseData/PulseCalls/ZMW/BaselineLevel", holeNumbers)
})
setMethod("getBaselineLevel", "MultiPart", function(h5Obj, holeNumbers = getHoleNumbers(h5Obj)) {
  .getDatasetWithChannelCols(h5Obj, "PulseData/PulseCalls/ZMW/BaselineLevel", holeNumbers)
})

.makeZMWMetricFx <- function(metric) {
  force(metric)

  function(basH5, holeNumbers = getHoleNumbers(basH5)) {
    .getHoleDataset(basH5, sprintf("PulseData/BaseCalls/ZMWMetrics/%s", metric), 
                    holeNumbers)
  }
}

getHoleStatus <- function(basH5, holeNumbers = getHoleNumbers(basH5)) {
  dsName <- sprintf("PulseData/BaseCalls/%s/%s", "ZMW", "HoleStatus")
  ds <- .getDatasetFromBasLike(basH5, dsName)
  attrs <- getH5Attribute(ds, "LookupTable")[]
  vals <- .getHoleDataset(basH5, dsName, holeNumbers)
  attrs[match(vals+1, 1:length(attrs))]
}

getHoleXY <- function(basH5, holeNumbers = getHoleNumbers(basH5)) {
  .getHoleDataset(basH5, sprintf("PulseData/BaseCalls/%s/%s", "ZMW", "HoleXY"),
                  holeNumbers)
}

## ##
## ## these are auto-generated doing something like:
## ## nms <- ls(getH5Group(fetchPlsH5(175171)[[1]], "PulseData/BaseCalls/ZMWMetrics"))[-1]
## ## cat(sprintf("%s <- .makeZMWMetricFx(\"%s\")\n", paste("get", nms, sep = ""), nms))
## ##
getProductivity       <- .makeZMWMetricFx("Productivity")
getBaseFraction       <- .makeZMWMetricFx("BaseFraction")
getBaseRate           <- .makeZMWMetricFx("BaseRate")
getBaseRateVsT        <- .makeZMWMetricFx("BaseRateVsT")
getBaseWidth          <- .makeZMWMetricFx("BaseWidth")
getCmBasQv            <- .makeZMWMetricFx("CmBasQv")
getCmDelQv            <- .makeZMWMetricFx("CmDelQv")
getCmInsQv            <- .makeZMWMetricFx("CmInsQv")
getCmSubQv            <- .makeZMWMetricFx("CmSubQv")
getDarkBaseRate       <- .makeZMWMetricFx("DarkBaseRate")
getHQRegionEndTime    <- .makeZMWMetricFx("HQRegionEndTime")
getHQRegionStartTime  <- .makeZMWMetricFx("HQRegionStartTime")
getHQRegionDuration   <- function(basH5, holeNumbers = getHoleNumbers(basH5)) {
  getHQRegionEndTime(basH5, holeNumbers) - getHQRegionStartTime(basH5, holeNumbers)
}
getLocalBaseRate      <- .makeZMWMetricFx("LocalBaseRate")
getNumBaseVsT         <- .makeZMWMetricFx("NumBaseVsT")
getNumPauseVsT        <- .makeZMWMetricFx("NumPauseVsT")
getPausiness          <- .makeZMWMetricFx("Pausiness")
getRmBasQv            <- .makeZMWMetricFx("RmBasQv")
getRmDelQv            <- .makeZMWMetricFx("RmDelQv")
getRmInsQv            <- .makeZMWMetricFx("RmInsQv")
getRmSubQv            <- .makeZMWMetricFx("RmSubQv")

getHQLength           <- function(basH5, holeNumbers = getHoleNumbers(basH5)) {
  rt <- getRegionsTable(basH5)
  rt <- rt[rt$type == "HQRegion",,drop = FALSE]
  rt <- rt[match(holeNumbers, rt$holeNumber), c("start", "end")]
  rt[,2] - rt[,1]
}

getReadType <- function(basH5, holeNumbers = getHoleNumbers(basH5)) {
  readType <- .makeZMWMetricFx("ReadType")(basH5, holeNumbers)
  ds <- .getDatasetFromBasLike(basH5, "PulseData/BaseCalls/ZMWMetrics/ReadType")
  attr <- getH5Attribute(ds, "UnitsOrEncoding")[]

  nAndv <- do.call(rbind, lapply(strsplit(attr, ","), function(z) {
    sapply(strsplit(z, ":"), '[', 1:2)
  }))
  nAndv[2, match(readType, as.integer(nAndv[1,]))]
}

associateZMWMetric <- function(cmpH5, plsH5s, zmwMetric = getProductivity) {
  doWithPlsAndCmp(cmpH5, plsH5s, function(cmpH5, plsH5, idx) {
    z <- zmwMetric(plsH5)
    m <- match(cmpH5$holeNumber[idx], getHoleNumbers(plsH5))
    if (is.null(dim(z)))
      z[m]
    else if (length(dim(z)) == 3) {
      z[m,,]
    } else {
      z[m,,drop = FALSE]
    }
  })
}


