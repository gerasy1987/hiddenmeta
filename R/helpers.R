## Helper functions

#' Generate sizes of groups in population
#'
#' @param N sample size
#' @param prev_K numeric vector of prevalence for each group with last group being hidden population
#' @param rho_K numeric vector that gives lower and upper triangles of correlation matrix between group memberships. If single number all group membership correlations are assumed to be equal
#' @param .ord character vector of ordered names of all possible combinations of group memberships. Default is \code{NULL} to generate the order automatically
#'
#' @return Vector of group sizes
#' @export
#'
#' @examples
#' \dontrun{
#' gen_group_sizes(N = 1000,
#'                 prev_K = c(known = .3, hidden = .1),
#'                 rho_K = .05)
#' }
#'
#' @import dplyr
#' @importFrom magrittr `%>%`
#' @importFrom purrr quietly
#' @importFrom MultiRNG draw.correlated.binary
gen_group_sizes <- function(N, prev_K, rho_K, .ord = NULL) {

  if (!(length(rho_K) %in% c(1, sum((length(prev_K) - 1):1))))
    stop("Length of correlations vector should be either 1 or sum((K - 1):1)")

  rho <- diag(rep(1, length(prev_K)))

  rho[lower.tri(rho)] <- rho_K
  rho[upper.tri(rho)] <- t(rho)[upper.tri(rho)]

  if (is.null(.ord)) {
    .ord <-
      seq_along(prev_K) %>%
      lapply(X = ., function(x) 0:1) %>%
      expand.grid() %>%
      apply(1, paste0, collapse = "")
  }

  MultiRNG::draw.correlated.binary(no.row = N,
                                   d = length(prev_K),
                                   prop.vec = prev_K,
                                   corr.mat = rho) %>%
    dplyr::as_tibble(.name_repair = function(names) paste0("V", 1:length(names))) %>%
    dplyr::group_by_all() %>%
    purrr::quietly(dplyr::summarise)(n = dplyr::n()) %>%
    .$result %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(group = paste0(dplyr::c_across(dplyr::starts_with("V")), collapse = "")) %>%
    dplyr::select(.data$group, .data$n) %>%
    {`names<-`(.$n, .$group)} %>%
    .[.ord]

  # # can only handle 2 groups at the moment !!!
  # {rho_K * sqrt(prod(prev_K * (1 - prev_K))) + prod(prev_K)} %>%
  #   round(.,floor(log10(N))) %>%
  #   {N * c(1 - sum(prev_K) + ., prev_K - ., .)} %>%
  #   `names<-`(
  #     apply(
  #       expand.grid(
  #         lapply(
  #           seq_along(prev_K), function(x) 0:1)), 1, paste0, collapse = ""))
}


#' Generate block matrix of link probabilities
#'
#' Returns full block matrix of link probabilities given the group memberships based on probabilities of links within each of the groups and between in- and out-group members
#'
#' @param p_edge_within named list of numeric vectors giving probability of link between in-group members and out-group members for each of groups
#' @param p_edge_between named list of numeric values giving probability of link between in- and out-group member for each of groups
#' @param .ord character vector of ordered names of all possible combinations of group memberships. Default is \code{NULL} to generate the order automatically
#'
#' @return \code{2^K} by \code{2^K} matrix of probabilities of establishing a link
#' @export
#'
#' @examples
#' \dontrun{
#' gen_block_matrix(p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.9)),
#'                  p_edge_between = list(known = 0.05, hidden = 0.01))
#' }
#'
#' @importFrom magrittr `%>%`
gen_block_matrix <- function(p_edge_within, p_edge_between, .ord = NULL) {

  if (length(p_edge_within) != length(p_edge_between))
    stop("p_edge_within needs to be a list of length K")

  if (!all(sapply(p_edge_within, function(x) length(x) == 2)))
    stop("Each element of p_edge_within should be of length 2 (in and out group probability)")

  within_mat <-
    lapply(rev(p_edge_within), function(x) 1 - diag(1 - x)) %>%
    Reduce(kronecker, x = .)

  between_mat <-
    lapply(rev(p_edge_between), function(x) get_perms(c(1,x))) %>%
    Reduce(kronecker, x = .)

  out_mat <- within_mat * between_mat

  if (is.null(.ord)) {
    colnames(out_mat) <- rownames(out_mat) <-
      lapply(1:length(p_edge_within), function(x) 0:1) %>%
      do.call(expand.grid, .) %>%
      apply(1, paste0, collapse = "")
  } else {
    colnames(out_mat) <- rownames(out_mat) <- .ord
  }

  return(out_mat)
}

#' Random generation from Beta distribution
#'
#' Reparametrization of base R rbeta to take location and scale rather than 2 shape parameters
#'
#' @param n number of observations
#' @param mu location parameter
#' @param sd scale parameter
#'
#' @return vector of draws from Beta distribution
#' @export
#'
#' @examples
#' \dontrun{
#' rbeta_mod(4, .5, .1)
#' }
#'
#' @importFrom stats rbeta
rbeta_mod <- function(n, mu, sd) {
  alpha <- ((1 - mu) / sd^2 - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(rbeta(n = n, shape1 = alpha, shape2 = beta))
}

# HIDDEN HELPERS ------------------------------------------------------------------------------



merge_vertices <- function(graph_list) {
  if (!all(c("edges", "vertices") %in% names(graph_list)))
    stop("Graph list has to have edges and vertices elements.")

  graph_list %>%
    {
      left_join(
        left_join(.$edges, cbind(.$vertices, from = 1:N),  by = "from"),
        cbind(.$vertices, to = 1:N), by = "to", suffix = c("_from", "_to")
      )
    }

}

# generates vector of marginal type names
gen_marginal_types <- function(K) {
  sapply(0:K,
         function(x) paste0(paste0(rep("*", times = x), collapse = ""),
                            "1",
                            paste0(rep("*", times = K - x), collapse = "")))
}

# helper for population simulations
mutate_visibility_tidy <- function(data, vars, p_visible, beta_sd = 0.05) {
  for (var in vars) {
    if (p_visible[[var]] == 1) {
      data[,paste0("p_visible_", var)] <- 1
    } else {
      data[,paste0("p_visible_", var)] <-
       rbeta_mod(nrow(data), mu = p_visible[[var]], sd = beta_sd)
    }

    data[,paste0(var, "_visible")] <-
      rbinom(nrow(data), 1, unlist(data[,paste0("p_visible_", var)])) * data[,var]
  }
  return(data)
}

# helper for population simulations using data.table
mutate_visibility <- function(dat, vars, p_visible, beta_sd = 0.05) {
  # dat <- copy(dat)
  for (var in vars) {
    if (p_visible[[var]] == 1) {
      dat[,paste0("p_visible_", var) := 1]
    } else {
      dat[,paste0("p_visible_", var) := rbeta_mod(.N, mu = p_visible[[var]], sd = beta_sd)]
    }
    dat[, paste0(var, "_visible") := rbinom(.N, 1, dat[[paste0("p_visible_", var)]]) * dat[[var]]]
  }
  return(dat)
}

# helper for population simulations
mutate_add_groups_tidy <- function(data, add_groups) {
  for (var in seq_along(add_groups)) {
    if (class(add_groups[[var]]) == "numeric") {
      data %<>%
        dplyr::mutate("{ names(add_groups)[var] }" := rbinom(n(), 1, add_groups[[var]]))
    } else if (class(add_groups[[var]]) == "character" & names(add_groups)[var] != "") {
      data %<>%
        dplyr::mutate("{ names(add_groups)[var] }" := eval(parse(text = add_groups[[var]])))
    } else if (class(add_groups[[var]]) == "character" & names(add_groups)[var] == "") {
      data %<>%
        dplyr::mutate(eval(parse(text = add_groups[[var]])))
    }
  }
  return(data)
}

# helper for population simulations using data.table package
mutate_add_groups <- function (dat, add_groups) {
  # dat <- data.table::copy(dat)
  for (var in seq_along(add_groups)) {
    if (class(add_groups[[var]]) == "numeric") {
      dat[, names(add_groups)[var] := rbinom(.N, 1, add_groups[[var]])]
    }
    else if (class(add_groups[[var]]) == "character" & names(add_groups)[var] != "") {
      dat <- dat[, `:=`(names(add_groups)[var], eval(parse(text = add_groups[[var]])))]
    }
    else if (class(add_groups[[var]]) == "character" & names(add_groups)[var] == "") {
      dat <- dat[, eval(parse(text = add_groups[[var]]))]
    }
  }
  return(dat)
}

# permutations helper from https://stackoverflow.com/a/34287541
get_perms <- function(x) {
  if (length(x) == 1) {
    return(x)
  }
  else {
    res <- matrix(nrow = 0, ncol = length(x))
    for (i in seq_along(x)) {
      res <- rbind(res, cbind(x[i], Recall(x[-i])))
    }
    return(res)
  }
}

# generate all possible group membership types based on number of groups k
gen_group_types <- function(K) {
  lapply(1:K, function(x) 0:1) %>%
    expand.grid() %>%
    apply(X = ., MARGIN = 1, FUN = paste0, collapse = "")
}

# helper to create data.table with characteristics of new eligible vertices
get_new_eligible <- function(sampled_df, to, parent, coup) {
  if (length(to) == 1) {
    data.table::data.table(
      from = parent,
      to = as.integer(to),
      own_coupon =
        sampled_df[name == parent][[paste0("coupon_", sample(1:coup, size = 1))]]
    )
  } else if (length(to) > 1) {
    data.table::data.table(
      from = parent,
      to = sample(x = as.integer(to),
                  size = min(length(to), coup),
                  replace = FALSE),
      own_coupon =
        as.character(
          sampled_df[
            name == parent,
            paste0("coupon_", sample(1:coup, size = min(length(to), coup))),
            with = FALSE
          ]
        )
    )
  } else {
    data.table::data.table(from = integer(),
                           to = integer(),
                           own_coupon = character())
  }
}

# get graph from adjacency list transformed to character vector
retrieve_graph <- function(adj_vec, split_vec = ";", list_mode = "all") {
  lapply(adj_vec, function(x) as.numeric(base::strsplit(x, split = split_vec)[[1]])) %>%
    igraph::graph.adjlist(mode = list_mode)
}

# get adjacency list from the one transformed to character vector
retrieve_adjlist <- function(adj_vec, split_vec = ";", list_mode = "all") {
   lapply(adj_vec, function(x) as.numeric(base::strsplit(x, split = split_vec)[[1]]))
}

# extract psu and strata variables from survey design formula
extract_terms <- function(fmla) {
  term_labels <- attr(terms(fmla), "term.labels")
  strata_ids <- grepl("strata\\(", term_labels)
  psu <- term_labels[!grepl("strata\\(", term_labels)]
  if (sum(strata_ids) == 1) {
    strata_text <- gsub(replacement = "", x = term_labels[strata_ids], pattern = "strata\\(|\\)")
    strata <- attr(terms(as.formula(paste("~ ", strata_text))), "term.labels")
  } else if (sum(strata_ids) == 0) {
    strata <- NULL
  } else if (sum(strata_ids) > 1) {
    stop("Design formula cannot have more than one strata() specification")
  }

  return(
    list(psu = psu, strata = strata)
  )
}

## -----------------------------------------------------------------------------

# helper to avoid naming conflicts between variables and column names for data.table

check_naming_conflict <- function(data,var,name) {
  names <- colnames(data)
  check_name <- var == name

  if(!is.null(var)) {
    if(any(check_name)){
      warning(paste(name, " column cannot be named '", name, "', renaming column", sep = ""))
      name_dup <- TRUE
      i <- 1
      name_new <- ""
      while(name_dup) {
        name_new <- paste(name,i, sep = "")
        if(!name_new %in% names) {
          name_dup <- FALSE
        }
        i <- i + 1
      }
      names[names == name] <- name_new
      colnames(data) <- names
      var[check_name] <- name_new
    }
  }
  return(list(colnames = names, var = var))
}















