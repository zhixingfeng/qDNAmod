require(pbh5)

if (! require(pbls)) {
  ## Essentially an internal test because of the lack of available
  ## "small" Multi-stream files.
  quit(save = "no")
}

source("tinyTestHarness.R")

plsH5 <- fetchPlsH5('2311641-0001')[[1]]
basH5 <- fetchBasH5('2311641-0001')[[1]]

plsH5
basH5

TH <- TestHarness("MultiStream")
tf <- function(basH5, f) {
  holes <- getHoleNumbers(basH5)
  parts <- pbh5:::.getParts(basH5, holes)
  byPart <- split(holes, parts)
  all(sapply(seq_along(byPart), function(i) {
    part <- pbh5:::part(basH5, i)
    holes <- sample(byPart[[i]], size = 500)
    all.equal(f(part, holes),
              f(basH5, holes))
  }))
}

TH("Basecalls Match", {
  tf(basH5, getBasecalls)
})

basFunctions <- c("getPreBaseFrames", "getPulseIndex", "getQualityValue",
                  "getDeletionQV", "getDeletionTag", "getInsertionQV",
                  "getMergeQV", "getSubstitutionQV", "getSubstitutionTag")
invisible(lapply(basFunctions, function(f) {
  TH(f, {
    tf(basH5, get(f))
  })
}))

plsFunctions <- c("getChannel", "getChi2", "getIsPulse", "getMaxSignal",
                  "getMeanSignal", "getMidSignal", "getMidStdDev",
                  "getStartFrame")
invisible(lapply(plsFunctions, function(f) {
  TH(f, {
    tf(plsH5, get(f))
  })
}))

zmwMetrics <- c('getProductivity',
                'getBaseFraction',     
                'getBaseRate',          
                'getBaseRateVsT',       
                'getBaseWidth',         
                'getCmBasQv',           
                'getCmDelQv',           
                'getCmInsQv',           
                'getCmSubQv',           
                'getDarkBaseRate',      
                'getHQRegionEndTime',   
                'getHQRegionStartTime', 
                'getHQRegionDuration',
                'getLocalBaseRate',   
                'getNumBaseVsT',      
                'getNumPauseVsT',     
                'getPausiness',       
                'getRmBasQv',          
                'getRmDelQv',          
                'getRmInsQv',          
                'getRmSubQv',
                'getHQLength',         
                'getSNR',
                'getReadScore')

invisible(lapply(zmwMetrics, function(f) {
  TH(f, {
    tf(basH5, get(f))
  })
}))

plsZMWMetrics <- c('getBaselineSigma', 'getBaselineLevel')

invisible(lapply(plsZMWMetrics, function(f) {
  TH(f, {
    tf(plsH5, get(f))
  })
}))


TH("Multi-Column Dataset Match", {
  ms <- getMeanSignal(plsH5, hn <- sample(getHoleNumbers(plsH5), size = 1000))
  ch <- getChannel(plsH5, hn)
  w <- mapply(ms, ch, FUN = function(a,b) ((is.null(a) && is.null(b)) || (nrow(a) == length(b))))
  all(w)
})

TH("Multi-Column Dataset Match, with pads", {
  pevts <- getPulseEvents(plsH5)
  hn <- sample(pevts$holeNumber[pevts$numEvent > 200], size = 1000)
  lp <- floor(runif(1000, 0, 100))
  rp <- floor(runif(1000, 0, 100))
  ms <- getMeanSignal(plsH5, holeNumbers = hn, leftPads = lp, rightPads = rp)
  ch <- getChannel(plsH5, hn, leftPads = lp, rightPads = rp)
  w <- mapply(ms, ch, FUN = function(a,b) (nrow(a) - length(b)))
  all(w == 0)
})

TH(action = 'print')
