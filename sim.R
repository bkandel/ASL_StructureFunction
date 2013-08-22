# t1 <- antsImageRead(getANTsRData('r64'), 3)
# antsImageWrite(t1, 'r64.nii.gz')
# t1 <- antsImageRead(getANTsRData('r16'), 3)
# antsImageWrite(t1, 'r16.nii.gz')
library(ANTsR)
set.seed(1213) # for reproducibility

nloc <- 400 # number of locations for added activity
noise <- c(150, 1500)

recompute <- F
if(recompute){
t1 <- antsImageRead(getANTsRData('mni'), 3)
antsImageWrite(t1, 'data/simulation/mni.nii.gz')
lab <- antsImageRead(getANTsRData('mnib'), 3)
antsImageWrite(lab, 'data/simulation/mnib.nii.gz')
cerebrum <- antsImageRead('data/simulation/cerebrum.nii.gz', 3)
mymask <- getMask(lab) # only cortical areas are labeled
antsImageWrite(mymask, 'data/simulation/mask.nii.gz')
cortex <- antsImageClone(lab)
cortex[cortex != 0] <- 1
antsImageWrite(cortex, 'data/simulation/cortex.nii.gz')
t1.gm <- maskImage(t1, cortex)
antsImageWrite(t1.gm, 'data/simulation/t1_cortex.nii.gz')
myseg <- Atropos(d=3, a=t1,  m = "[0.1,1x1x1]", c = "[25,0]", 
                  x = cerebrum, i="kmeans[3]")
antsImageWrite(myseg$segmentation, 'data/simulation/seg.nii.gz')
antsImageWrite(myseg$probabilityimages[[1]], 'data/simulation/csf.nii.gz')
antsImageWrite(myseg$probabilityimages[[2]], 'data/simulation/gm.nii.gz')
antsImageWrite(myseg$probabilityimages[[3]], 'data/simulation/wm.nii.gz')
} 
lab <- antsImageRead('data/simulation/mnib.nii.gz', 3)
t1 <- antsImageRead('data/simulation/mni.nii.gz', 3)
mymask <- getMask(t1)
csf <- antsImageRead('data/simulation/csf.nii.gz', 3)
gm <- antsImageRead('data/simulation/gm.nii.gz', 3)
wm <- antsImageRead('data/simulation/wm.nii.gz', 3)
seg <- antsImageRead('data/simulation/seg.nii.gz', 3)
cortex <- antsImageRead('data/simulation/cortex.nii.gz', 3)
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
tmp <- rep(0, length(mymask[mymask>0]))
tmp1 <- lab[mymask >0]
ant.left <- which(tmp1==37)
tmp[ant.left] <- 30
myblob <- antsImageClone(t1)
myblob[mymask > 0] <- 0
myblob[mymask>0][ant.left] <- 20
antsImageWrite(myblob, 'data/simulation/blob.nii.gz')
system(paste("SurfaceBasedSmoothing data/simulation/blob.nii.gz 1", 
             "data/simulation/cortex.nii.gz data/simulation/blob_smooth.nii.gz 10"))
myblob <- antsImageRead('data/simulation/blob_smooth.nii.gz', 3)

func <- antsImageClone(t1)

func[mymask>0] <- csf[mymask>0] * 5 +
  gm[mymask>0] * 100 +  wm[mymask>0] * 0.2 * 100 + 
  + myblob[mymask>0] * 1 + mynoise.img[mymask>0] * 5 + 
  rnorm(length(mymask[mymask>0]), sd=10) * 5
func[gm.mask>0] <- func[gm.mask>0] + rnorm(length(gm.mask[gm.mask>0]), sd=10) * 30
SmoothImage(3, func, 2, func)
antsImageWrite(func, 'data/simulation/func.nii.gz')
