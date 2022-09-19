initParms <- function(newParms = NULL) {
  parms <- c(
    half_life = 70.0,
    F_m = 0.094,
    F_milk = 0.154,
    r_f_m = 0.35,
    F_abs = 0.9,
    n_i = 10,
    food_dose = 0,
    k = 0.0
  )

  if (!is.null(newParms)) {
    if (!all(names(newParms) %in% c(names(parms)))) {
      stop("illegal parameter name")
    }
    parms[names(newParms)] <- newParms
  }

  parms <- within(as.list(parms), {
  })
  out <- .C("getParms",  as.double(parms),
            out=double(length(parms)),
            as.integer(length(parms)))$out
  names(out) <- names(parms)
  out
}

Outputs <- c(
    "C_m",
    "C_i",
    "M_mf",
    "M_i",
    "D_m",
    "D_i",
    "R_milk",
    "r_m_mf",
    "A_bal"
)

initStates <- function(parms, newStates = NULL) {
  Y <- c(
    A_mf = 0.0,
    A_i = 0.0,
    AUC_m = 0.0,
    AUC_i = 0.0,
    d_m = 0.0,
    d_i = 0.0,
    T_in = 0.0,
    T_out = 0.0
  )

  Y <- within(c(as.list(parms),as.list(Y)), {
    Y["A_mf"] <- 0.0 
    Y["A_i"] <- 0.0 
    Y["AUC_m"] <- 0.0 
    Y["AUC_i"] <- 0.0 
    Y["d_m"] <- 0.0 
    Y["d_i"] <- 0.0 
    Y["T_in"] <- 0.0 
    Y["T_out"] <- 0.0 

  })$Y

  if (!is.null(newStates)) {
    if (!all(names(newStates) %in% c(names(Y)))) {
      stop("illegal state variable name in newStates")
    }
    Y[names(newStates)] <- newStates
  }

.C("initState", as.double(Y));
Y
}
