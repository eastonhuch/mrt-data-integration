final_sandwich_step <- function(hessian, meat, n_users, d) {
  half_sandwich <- solve(hessian, t(chol(meat)), tol=1e-50)
  sandwich <- tcrossprod(half_sandwich) * n_users / (n_users-d)
  sandwich
}

construct_sandwich_balanced <- function(scores, hessian, n_users, t_max, d) {
  scores_agg <- apply(
    aperm(
      array(scores, dim = c(t_max, n_users, d)),
      c(2,1,3)),
    MARGIN=c(1,3), FUN=sum)
  meat <- crossprod(scores_agg)
  final_sandwich_step(hessian, meat, n_users, d)
}

construct_sandwich_unbalanced <- function(scores, hessian, user_ids, n_users, d){
  scores_agg <- as.matrix(aggregate(scores ~ user_ids, FUN=sum)[, -1])
  colnames(scores_agg) <- NULL
  meat <- crossprod(scores_agg)
  final_sandwich_step(hessian, meat, n_users, d)
}