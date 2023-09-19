is.Raster <- function(x) {
  return((
    class(x)[1] == "RasterLayer" ||
      class(x)[1] == "RasterBrick" ||
      class(x)[1] == "RasterStack" || class(x)[1] == "SpatRaster"
  ))
}

get_wr2_ss <- function(y, y_pred, w) {
  ss_residual <- sum(w * (y - y_pred)^2)
  ss_total <- sum(w * (y - weighted.mean(y, w))^2)
  return(1 - ss_residual / ss_total)
}

add_terms <-
  function(av) {
    reduce(av, function(x, y) {
      str_c(x, " + ", y)
    })
  }

object_size <- function(aobject) {
  aobject %>%
    object.size() %>%
    print(unit = "auto")
}

if (Sys.info()[1] == "Windows") {
  qn <- 5
} else {
  qn <- 3
}

trim_xy <- function(adata_xy) {
  adata_xy %>%
    mutate(across(c(x, y), ~ round(.x, 4)))
}

check_join <- function(adata) {
  targeted <- rast("data/outputs/masks/mask_extraction.tif") %>%
    sum(na.rm = T) %>%
    as.data.frame(xy = T) %>%
    distinct(x, y) %>%
    trim_xy()

  de <- adata %>%
    distinct(x, y)

  joined <- inner_join(targeted, de)
  antied <- anti_join(targeted, de)

  tibble(
    mask_grid = nrow(targeted),
    test_grid = nrow(de),
    inner_joined = nrow(joined),
    anti_joined = nrow(antied)
  )
}

arrow <-
  arrow(
    length = unit(0.015, "npc"),
    ends = "last",
    type = "open"
  )

get_data <- function(varss =
                       c(
                         "GOSIF_sum",
                         "cloud",
                         "AOD",
                         "PM25",
                         "W126_1",
                         str_c("bin", 0:39),
                         str_c("step", 1:10),
                         str_c("root_", 1:9)
                       )) {
  qread("data/tidied.qs", nthreads = qn) %>%
    lazy_dt() %>%
    mutate(
      crop_parent = fifelse(str_detect(crop, "Rice"), "Rice", crop),
      x_y = str_c(crop, x_y),
      across(
        c(
          starts_with("GOSIF"),
          starts_with("Wenetal"),
          starts_with("RTSIF"),
          starts_with("CSIF")
        ),
        log
      ),
      DF = DR / GR,
      DirR = GR - DR,
      dif = diff / down,
      AOT40 = fcase(
        crop_parent == "Maize",
        AOT40_6,
        crop_parent == "Rice",
        AOT40_7,
        crop_parent == "Wheat",
        AOT40_7 - AOT40_1
      ),
      W126 = fcase(
        crop_parent == "Maize",
        W126_6,
        crop_parent == "Rice",
        W126_7,
        crop_parent == "Wheat",
        W126_7 - W126_1
      ),
      O3 = fcase(
        crop_parent == "Maize",
        O3_6,
        crop_parent == "Rice",
        O3_7,
        crop_parent == "Wheat",
        (O3_7 * (MA - `GR&EM` + 1) - O3_1) / (MA - `GR&EM` + 0)
      )
    ) %>%
    drop_na(all_of(varss)) %>%
    group_by(crop, x, y) %>%
    filter(n() >= 10) %>%
    ungroup() %>%
    as_tibble() %>%
    group_by(crop) %>%
    mutate(across(
      c(AOD, NO2, starts_with("W126"), starts_with("AOT40")),
      ~ DescTools::Winsorize(.x, probs = c(0, 0.995), na.rm = T)
    )) %>%
    ungroup() %>%
    arrange(crop) %>%
    inner_join(read_xlsx("regions.xlsx") %>% mutate(region = str_c("R", region)),
      by = join_by(province, crop_parent)
    )
}

order <- c(
  "North China",
  "Northeast China",
  "East China",
  "China",
  "South China",
  "Central China",
  "Southwest China",
  "Northwest China"
)

relpred <-
  function(object,
           newdata,
           baseline = NULL,
           level = 0.90) {
    # object: lm or fit_feols fitted model
    # newdata: data.frame (or matrix) of X variables with names that correspond to fitted model (names that don't fit are ignored)
    # baseline: list with same length as the number of columns in newdata (or, eventually, matching names). without this, baseline is assumed to be zero.
    # object <- lm(y ~ x + I(x^2) + group, data = data)
    # object <- feols(y ~ x + I(x^2) | group, data = data)
    # newdata <- newdata; level = 0.95
    # baseline <- NULL
    # baseline <- list(x = 1)
    # END DEBUG

    # TODO: Add support for factor variables in newdata later. Requirements
    # - Must either check or convert variables to match levels in object.
    # - Decide how baseline should be chosen for factor variables if baseline is explicitly supplied.
    if (any(sapply(newdata, function(x) {
      !is.numeric(x)
    }))) {
      stop("For now, only numeric variables are supported in newdata.")
    }

    if (!is.null(baseline) & ncol(newdata) != length(baseline)) {
      stop("baseline vector length must be equal to number of columns in newdata.")
    }

    # Populate a data.frame with all variables used in object estimation and return a design matrix that matches that object.
    newdata_filled <- fill_missing_vars(object, newdata)

    # If baseline isn't passed, assume zero
    if (is.null(baseline)) {
      baseline_df <- newdata
      baseline_df[] <- 0
    } else {
      baseline_df <- as.data.frame(baseline)
    }

    baseline_filled <- fill_missing_vars(object, baseline_df)

    # Subtract baseline from newdata to get X
    X <- as.matrix(newdata_filled - baseline_filled)

    # Compute coefficients and standard errors ----
    B <- as.numeric(coef(object))

    # Compute degrees of freedom
    if (inherits(object, "fixest")) {
      df <- attributes(vcov(object, attr = T))$G
    } else {
      df <- object$df.residual
    }

    # Drop any variables in X that are not the fitted model object
    # This will handle a variety of issues, including a missing intercept term
    X <- X[, names(coef(object))]

    fit <- data.frame(fit = X %*% B)
    sig <- vcov(object)
    se <- apply(X,
      MARGIN = 1,
      FUN = get_se,
      sig = sig
    )

    t_val <- qt((1 - level) / 2 + level, df = df)
    fit$lwr <- fit$fit - t_val * se
    fit$upr <- fit$fit + t_val * se

    fit %>% as_tibble()
  }

fill_missing_vars <- function(object, X) {
  # Populate a data.frame with all variables used in object estimation and return a design matrix with the same coefficients.
  # Factor variables filled in using the first factor level, numeric filled with zeroes

  orig_vars <- all.vars(formula(object))[-1]

  # Add any factor variables that are in the RHS but not newdata with the baseline value
  # TODO: How to add back a missing factor variable for feols? We can do it with lm using object$xlevels
  # TODO: If the user passes, e.g., factor(group) this will not work. Not sure how to resolve.
  fact_vars <- names(object$xlevels)
  for (f in fact_vars) {
    if (!(f %in% names(X))) {
      X[[f]] <- factor(object$xlevels[[f]][1],
        levels = object$xlevels[[f]]
      )
    } else {
      # This is risky, but will ensure levels in object match those in X...
      X[[f]] <- factor(X[[f]],
        levels = object$xlevels[[f]]
      )
    }
  }

  # If there are any other variables in the formula but not in X, add them with value = 0
  for (v in orig_vars) {
    if (!(v %in% names(X))) {
      X[[v]] <- 0
    }
  }

  # Extract terms object, removing response variable
  tt <- terms(object)
  Terms <- delete.response(tt)

  as.data.frame(model.matrix(Terms, data = X))
}

get_se <- function(r, sig) {
  # Compute linear combination, helper function for predict_partial
  # Given numeric vector r (the constants) and vcov sig, compute SE
  r <- matrix(r, nrow = 1)
  sqrt(r %*% sig %*% t(r))
}
