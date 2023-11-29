
mydiff.r <- function(mesh, fem, barrier.triangles, 
                     prior.range = c(1.44, 0.5), prior.sigma = c(0.7, 0.5),
                     range.fraction,
                     range = 20,
                     set.inla.seed = 8,
                     return.list = FALSE,
                     return.plot = TRUE) {
  
  barrier.model <- inla.barrier.pcmatern.plus(mesh, fem, barrier.triangles, 
                                              prior.range, prior.sigma, 
                                              range.fraction)
  # range fraction has to have the same length has barrier triangles
  Q <- inla.rgeneric.q(barrier.model, "Q", theta = c(0, log(range)))
  u <- suppressWarnings(inla.qsample(n=1, Q=Q, seed = set.inla.seed))
  u <- u[ ,1]
  
  if (return.plot == TRUE) {
  local.plot.field(u, main="True spatial field",
                   sub = paste("Barrier Model, range.fraction = c(", 
                               range.fraction[1], ", ", range.fraction[2], ")")) 
  }
  if (return.list == TRUE) {
    return(list(barrier.model = barrier.model, Q = Q, sample = u))
  }
}

##################

my_diff.pos <- function(mesh, loc.data, Q, u,
         sigma.u = 1, sigma.epsilon = 0.2,
         range.fraction = c(0.1, 0.1),
         pos = res.barrier$summary.random$s$mean,
         return.list = FALSE,
         return.plot = TRUE) {
  
  A.data <- inla.spde.make.A(mesh, loc.data)
  #Q from mydiff.r function
  u.data <- A.data %*% u
  df <- data.frame(loc.data)
  names(df) <- c('locx', 'locy')

  df$y <- drop(sigma.u * u.data + sigma.epsilon * rnorm(nrow(df)))
  
  stk <- inla.stack(data=list(y=df$y), 
                    A=list(A.data, 1),
                    effects=list(s=1:mesh$n, 
                                 intercept=rep(1, nrow(df))), 
                    remove.unused = FALSE, 
                    tag='est')
  formula <- y ~ 0 + intercept + f(s, model = barrier.model)
  # - The spatial model component is different from stationary
  # - The rest of the model setup is the same as in the stationary case!
  # - - e.g. the inla(...) call below is the same, 
  #     only this formula is different
  
  res.barrier <- inla(formula, data = inla.stack.data(stk),
                        control.predictor = list(A = inla.stack.A(stk)),
                        family = 'gaussian',
                        control.family = list(hyper = list(
                          prec = list(prior = "pc.prec", fixed = FALSE, param = c(0.2,0.5)))),
                        control.mode=list(restart=T, theta=c(3.2, 0.4, 1.6)))
  if (return.plot == TRUE) {
  local.plot.field(pos, 
                   main="Spatial posterior for Barrier model",
                   sub = paste("Barrier Model, range.fraction = c(", 
                               range.fraction[1], ", ", range.fraction[2], ")")) 
  plot(poly.bar.orginal, add=T)
  }
  
  if (return.list == TRUE) {
    return(list(pos = quote(pos), df = df, summ.df = summary(df), 
                res = res.barrier, summ.res = summary(res.barrier)))
  }
}


