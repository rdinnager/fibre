set_env <- function(env_from, env_to = rlang::global_env()) {
  rlang::env_unbind(env_to, rlang::env_names(env_from))
  rlang::env_coalesce(env_to, env_from)
}