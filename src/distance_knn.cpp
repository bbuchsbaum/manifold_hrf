#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// compute full pairwise Euclidean distance matrix for columns of M
// M is p x N (features x samples)
// [[Rcpp::export]]
arma::mat pairwise_distances_cpp(const arma::mat& M) {
    arma::rowvec ss = arma::sum(arma::square(M), 0);
    arma::mat dist2 = arma::repmat(ss.t(), 1, M.n_cols) + arma::repmat(ss, M.n_cols, 1) - 2 * M.t() * M;
    dist2.elem(arma::find(dist2 < 0)).zeros();
    arma::mat D = arma::sqrt(dist2);
    D.diag().zeros();
    return D;
}

// compute k nearest neighbors of query columns against data columns
// returns list with idx (k x n_query) and dist (k x n_query)
// [[Rcpp::export]]
Rcpp::List knn_search_cpp(const arma::mat& data, const arma::mat& query, int k) {
    arma::rowvec data_ss = arma::sum(arma::square(data), 0);
    arma::rowvec query_ss = arma::sum(arma::square(query), 0);
    arma::mat gram = data.t() * query;
    arma::mat dist2 = arma::repmat(data_ss.t(), 1, query.n_cols) + arma::repmat(query_ss, data.n_cols, 1) - 2 * gram;
    dist2.elem(arma::find(dist2 < 0)).zeros();
    int n_query = query.n_cols;
    arma::umat idx(k, n_query);
    arma::mat dists(k, n_query);
    for(int i = 0; i < n_query; ++i) {
        arma::vec col = dist2.col(i);
        arma::uvec ord = arma::sort_index(col);
        arma::uvec top = ord.head(k);
        idx.col(i) = top;
        dists.col(i) = arma::sqrt(col.elem(top));
    }
    idx += 1; // convert to 1-based indexing
    return Rcpp::List::create(Rcpp::Named("idx") = idx, Rcpp::Named("dist") = dists);
}

