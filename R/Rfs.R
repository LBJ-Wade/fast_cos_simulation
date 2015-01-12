dyn.load("Rfs.so")

#
# Function definitions
#
fs.init <- function(omega.m0) {
  .C("rfs_init", as.double(omega.m0))
}

read.powerspectrum <- function(filename)
  .Call("rfs_read_powerspectrum", filename)


lpt.init <- function(nc, boxsize) {
  .C("rfs_lpt_init", as.integer(nc), as.double(boxsize))
  NULL
}

set.lpt <- function(seed, a) {
  .C("rfs_set_LPT", as.integer(seed), as.double(a))
  NULL
}

particles <- function(arg) {
  # particles(0:9)              => paticles with indices 0..9
  # particles(z_max)            => paticles with 0 <= z < z_max
  # particles(c(z_min, z_max))  => particles with z_min <= z < z_max
  # particles(c(x_min, x_max, y_min, y_max, z_min, z_max))
  .Call("rfs_particles", arg)
}

#
# Initial setup
#
omega.m0 <- 0.273
fs.init(omega.m0)

ps <- read.powerspectrum("camb_matterpower.dat")

nc <- 16
boxsize <- 100

lpt.init(nc, boxsize)

d <- particles(1:(nc^2-1))
