# t1 <- antsImageRead(getANTsRData('r64'), 3)
# antsImageWrite(t1, 'r64.nii.gz')
# t1 <- antsImageRead(getANTsRData('r16'), 3)
# antsImageWrite(t1, 'r16.nii.gz')
library(ANTsR)
set.seed(1213) # for reproducibility

### Set up parameters ####
nloc <- 1000 # number of locations for added activity
noise <- c(150, 1500)
area <- 37 # brodmann area to take
randsrc.amp <- 15
blob.amp <- 15
speckle.amp <- 15
specklegm.amp <- 20
smooth.sigma <- 2


recompute <- F
aslres <- F
for(specklegm.amp in c(40)){
if(recompute){
  t1 <- antsImageRead(getANTsRData('mni'), 3)
  lab <- antsImageRead(getANTsRData('mnib'), 3)
  blob <- new('antsImage', 'float', 3)
  ThresholdImage(3, lab, blob, area, area)
  cerebrum <- antsImageRead('data/simulation/cerebrum.nii.gz', 3)
  cortex <- antsImageClone(lab)
  ThresholdImage(3, cortex, cortex, 1, 999) 
  gm_prior <- antsImageRead('data/simulation/gm_prior.nii.gz', 3)
  wm_prior <- antsImageRead('data/simulation/wm_prior.nii.gz', 3)
  csf_prior <- antsImageRead('data/simulation/csf_prior.nii.gz', 3)  
  if(aslres){
    for (img in c(cortex, blob, cerebrum)){
      ResampleImageBySpacing(list(3, img, img, 3.25, 3.25, 6))
      ThresholdImage(3, img, img, 0.5, 999)
    }
    for ( img in c(gm_prior, wm_prior, csf_prior, t1)){
      ResampleImageBySpacing(list(3, img, img, 3.25, 3.25, 6))
    }
  }
  myseg <- Atropos(d=3, a=t1,  m = "[0.1,1x1x1]", c = "[25,0]", 
                   x = cerebrum, i=list(csf_prior, gm_prior, wm_prior), priorweight=0.8)
  antsImageWrite(t1, 'data/simulation/mni.nii.gz')
  antsImageWrite(blob, 'data/simulation/blob.nii.gz')
  antsImageWrite(myseg$segmentation, 'data/simulation/seg.nii.gz')
  antsImageWrite(myseg$probabilityimages[[1]], 'data/simulation/csf.nii.gz')
  antsImageWrite(myseg$probabilityimages[[2]], 'data/simulation/gm.nii.gz')
  antsImageWrite(myseg$probabilityimages[[3]], 'data/simulation/wm.nii.gz')
} 
t1 <- antsImageRead('data/simulation/mni.nii.gz', 3)
mymask <- getMask(t1)
csf <- antsImageRead('data/simulation/csf.nii.gz', 3)
gm <- antsImageRead('data/simulation/gm.nii.gz', 3)
wm <- antsImageRead('data/simulation/wm.nii.gz', 3)
seg <- antsImageRead('data/simulation/seg.nii.gz', 3)
gm.mask <- maskImage(seg, seg, 2)
gm.mask[gm.mask>0] <- 1
antsImageWrite(gm.mask, 'data/simulation/gm_mask.nii.gz')

mynoise <- rep(0, length(gm.mask[gm.mask>0]))
locs <- round(runif(nloc, min=0, max=length(mynoise)))
mynoise[locs] <- runif(length(locs), min=noise[1], max=noise[2])
mynoise.img <- antsImageClone(gm.mask)
mynoise.img[mymask>0] <- 0
mynoise.img[gm.mask>0] <- mynoise
SmoothImage(3, mynoise.img, 1, mynoise.img)
ImageMath(3, mynoise.img, 'm', mynoise.img, 5)
antsImageWrite(mynoise.img, 'data/simulation/noise.nii.gz')
system(paste("SurfaceBasedSmoothing data/simulation/noise.nii.gz 1", 
             "data/simulation/gm_mask.nii.gz data/simulation/noise_smooth.nii.gz 40"))
mynoise.img <- antsImageRead('data/simulation/noise_smooth.nii.gz', 3)
myblob <- antsImageRead('data/simulation/blob.nii.gz', 3)
system(paste("SurfaceBasedSmoothing data/simulation/blob.nii.gz 1", 
             "data/simulation/gm_mask.nii.gz data/simulation/blob_smooth.nii.gz 10"))
myblob <- antsImageRead('data/simulation/blob_smooth.nii.gz', 3)

func <- antsImageClone(t1)

func[mymask>0] <- csf[mymask>0] * 5 +
  gm[mymask>0] * 100 +  wm[mymask>0] * 0.4 * 100 + 
  + myblob[mymask>0] * blob.amp + mynoise.img[mymask>0] * randsrc.amp + 
  rnorm(length(mymask[mymask>0]), sd=10) * speckle.amp
func[gm.mask>0] <- func[gm.mask>0] + rnorm(length(gm.mask[gm.mask>0]), sd=10) * specklegm.amp
SmoothImage(3, func, smooth.sigma, func)
antsImageWrite(func, paste('data/simulation/func', specklegm.amp, '.nii.gz', sep=''))

if(aslres){
  ResampleImageBySpacing(list(3, func, func, 1, 1, 1))
  antsImageWrite(func, 'data/simulation/func_upsample.nii.gz')
}


}