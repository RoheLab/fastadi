p_omega_f_norm_ut_impl <- function(U, d, V, row, col, num_threads) {
    RcppParallel::setThreadOptions(numThreads = num_threads)
    res <- p_omega_f_norm_ut_impl_cpp(U, d, V, row, col)
    RcppParallel::defaultNumThreads()
    return(res)
}

p_u_ztx_impl <- function(U, d, V, x, num_threads) {
    RcppParallel::setThreadOptions(numThreads = num_threads)
    res <- p_u_ztx_impl_cpp(U, d, V, x)
    RcppParallel::defaultNumThreads()
    return(res)
}

p_u_zx_impl <- function(U, d, V, x, num_threads) {
    RcppParallel::setThreadOptions(numThreads = num_threads)
    res <- p_u_zx_impl_cpp(U, d, V, x)
    RcppParallel::defaultNumThreads()
    return(res)
}

relative_f_norm_change_impl <- function(new_U, new_d, new_V, U, d, V, num_threads) {
    RcppParallel::setThreadOptions(numThreads = num_threads)
    res <- relative_f_norm_change_impl_cpp(new_U, new_d, new_V, U, d, V)
    RcppParallel::defaultNumThreads()
    return(res)
}