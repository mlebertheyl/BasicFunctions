# save plots for gif in gif_r folder
#####
# THE TRUE SIMULATED SPATIAL FIELD
#####

nfrac <- 20
###############
for (r in 1:nfrac) {
  png(filename = )
  f <- r/20
  png.name <- paste0(paste0("gif_r/r",f),".png")
  png(filename = png.name)
  mydiff.r(mesh, fem = mat.plus, barrier.triangles = list(bar1, bar2),
           range.fraction = c(0.05, f))
  dev.off()
}

################ dir_out
# something like this but it needs to be fixed
dir_out <- file.path(tempdir(), "gif_r")
dir.create(dir_out, recursive = TRUE)

for (r in 1:nfrac) {
  png(filename = )
  f <- r/20
  png.name <- paste0("r",f)
  fpp <- file.path(dir_out, paste0(png.name, ".png"))
  png(filename = fpp)
  mydiff.r(mesh, fem = mat.plus, barrier.triangles = list(bar1, bar2),
           range.fraction = c(0.05, f))
  dev.off()
}
################END dir_out

library(magick)
## list file names and read in
imgs <- list.files("gif_r", full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "diff_r20.gif")
##################

#-------------
#different frac

nfrac <- 10

for (r in 1:nfrac) {
  f <- r/10
  png.name <- paste0(paste0("gif_r/r",f),".png")
  png(filename = png.name)
  mydiff.r(mesh, fem = mat.plus, barrier.triangles = list(bar1, bar2),
           range.fraction = c(0.05, f))
  dev.off()
}


library(magick)
## list file names and read in
imgs <- list.files("gif_r", full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "diff_r.gif")



# ----------------------------------
#####
# SPATIAL POSTERIOR FOR BARRIER MODEL
#####

nfrac <- 20
################ dir_out
dir_out <- file.path(tempdir(), "posmean_gif")
dir.create(dir_out, recursive = TRUE)

for(r in 1:nfrac) {
  f <- r/20
  mats <- mydiff.r(mesh, fem = mat.plus, barrier.triangles = list(bar1, bar2),
                   range.fraction = c(0.05, f), 
                   return.list = TRUE, return.plot = FALSE)
  
  fpp <- file.path(dir_out, paste0(f, ".png"))
  png(fpp)
  my_diff.pos(mesh, loc.data, 
              mats[[2]],
              mats$sample, 
              range.fraction = c(0.05, f),
              return.list = FALSE, return.plot = TRUE)
  dev.off()
}

## list file names and read in
imgs <- list.files(dir_out, full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "posmean_r20.gif")

