#ifndef algebra_hpp
#define algebra_hpp

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

std::vector<double> operator+(const std::vector<double>& vec1, const double& value) {
    int size = vec1.size();
    std::vector<double> vec(size, 0.);
    for(int i = 0; i < size; i++) vec[i] = vec1[i] + value;
    return vec;
}

std::vector<double> operator+(const double& value, const std::vector<double>& vec1);

std::vector<double> operator+(const std::vector<double>& vec1, const std::vector<double>& vec2);

std::vector<std::vector<double>> operator+(const std::vector<std::vector<double> >& mat1, const std::vector<std::vector<double> >& mat2);

std::vector<std::vector<double>> operator+(const double& value, const std::vector<std::vector<double> >& mat1);

void output(const std::vector<std::vector<double> >& A);

std::vector<double> operator+(const std::vector<double>& vec1, const double& value){
    int size = vec1.size();
    std::vector<double> vec(size, 0.);
    for(int i = 0; i < size; i++) vec[i] = vec1[i] + value;
    return vec;
}

std::vector<double> operator+(const double& value, const std::vector<double>& vec1){
    int size = vec1.size();
    std::vector<double> vec(size, 0.);
    for(int i = 0; i < size; i++) vec[i] = value + vec1[i];
    return vec;
}

std::vector<double> operator+(const std::vector<double>& vec1, const std::vector<double>& vec2){
    int size = vec1.size();
    std::vector<double> vec(size, 0.);
    for(int i = 0; i < size; i++) vec[i] = vec1[i] + vec2[i];
    return vec;
}

std::vector<std::vector<double> > operator+(const std::vector<std::vector<double> >& mat1, const std::vector<std::vector<double> >& mat2){
    int row = mat1.size();
    int col = mat1.front().size();
    std::vector<std::vector<double> > mat(row, std::vector<double>(col, 0));
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            mat[i][j] = mat1[i][j] + mat2[i][j];
        }
    }
    return mat;
}

std::vector<std::vector<double> > operator+(const double& value, const std::vector<std::vector<double> >& mat1){
    int row = mat1.size();
    int col = mat1.front().size();
    std::vector<std::vector<double> > mat(row, std::vector<double>(col, 0));
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            mat[i][j] = value + mat1[i][j];
        }
    }
    return mat;
}

std::vector<std::vector<double> > operator+(const std::vector<std::vector<double> >& mat1, const double& value){
    int row = mat1.size();
    int col = mat1.front().size();
    std::vector<std::vector<double> > mat(row, std::vector<double>(col, 0));
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            mat[i][j] = mat1[i][j] + value;
        }
    }
    return mat;
}


std::vector<double> operator-(const std::vector<double>& vec1, const double& value){
    int size = vec1.size();
    std::vector<double> vec(size, 0.);
    for(int i = 0; i < size; i++) vec[i] = vec1[i] - value;
    return vec;
}

std::vector<double> operator-(const double& value, const std::vector<double>& vec1){
    int size = vec1.size();
    std::vector<double> vec(size, 0.);
    for(int i = 0; i < size; i++) vec[i] = value - vec1[i];
    return vec;
}

std::vector<double> operator-(const std::vector<double>& vec1, const std::vector<double>& vec2){
    int size = vec1.size();
    std::vector<double> vec(size, 0.);
    for(int i = 0; i < size; i++) vec[i] = vec1[i] - vec2[i];
    return vec;
}

std::vector<std::vector<double> > operator-(const std::vector<std::vector<double> >& mat1, const std::vector<std::vector<double> >& mat2){
    int row = mat1.size();
    int col = mat1.front().size();
    std::vector<std::vector<double> > mat(row, std::vector<double>(col, 0));
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            mat[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
    return mat;
}

std::vector<std::vector<double> > operator-(const double& value, const std::vector<std::vector<double> >& mat1){
    int row = mat1.size();
    int col = mat1.front().size();
    std::vector<std::vector<double> > mat(row, std::vector<double>(col, 0));
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            mat[i][j] = value - mat1[i][j];
        }
    }
    return mat;
}

std::vector<std::vector<double> > operator-(const std::vector<std::vector<double> >& mat1, const double& value){
    int row = mat1.size();
    int col = mat1.front().size();
    std::vector<std::vector<double> > mat(row, std::vector<double>(col, 0));
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            mat[i][j] = mat1[i][j] - value;
        }
    }
    return mat;
}

std::vector<double> operator*(const double& value, const std::vector<double>& vec1){
    int size = vec1.size();
    std::vector<double> vec(size, 0.);
    for(int i = 0; i < size; i++) vec[i] = value * vec1[i];
    return vec;
}

std::vector<double> operator*(const std::vector<double>& vec1, const double& value){
    int size = vec1.size();
    std::vector<double> vec(size, 0.);
    for(int i = 0; i < size; i++) vec[i] = vec1[i] * value;
    return vec;
}

std::vector<std::vector<double> > operator*(const std::vector<double>& vec1, const std::vector<double>& vec2){
    int row = vec1.size();
    int col = vec2.size();
    std::vector<std::vector<double> > mat(row, std::vector<double>(col, 0));
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            mat[i][j] = vec1[i] * vec2[j];
        }
    }
    return mat;
}

std::vector<std::vector<double> > operator*(const std::vector<std::vector<double> >& mat1, const std::vector<std::vector<double> >& mat2){
    int row1 = mat1.size();
    int col1 = mat1.front().size();
    int col2 = mat2.front().size();
    std::vector<std::vector<double> > mat(row1, std::vector<double>(col2, 0));
    for(int i = 0; i < row1; i++){
        for(int j = 0; j < col1; j++){
            for(int k = 0; k < col2; k++){
                mat[i][k] += mat1[i][j] * mat2[j][k];
            }
        }
    }
    return mat;
}

std::vector<std::vector<double> > operator*(const double& value, const std::vector<std::vector<double> >& mat1){
    int row = mat1.size();
    int col = mat1.front().size();
    std::vector<std::vector<double> > mat(row, std::vector<double>(col, 0));
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            mat[i][j] = value * mat1[i][j];
        }
    }
    return mat;
}

std::vector<std::vector<double> > operator*(const std::vector<std::vector<double> >& mat1, const double& value){
    int row = mat1.size();
    int col = mat1.front().size();
    std::vector<std::vector<double> > mat(row, std::vector<double>(col, 0));
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            mat[i][j] = mat1[i][j] * value;
        }
    }
    return mat;
}

std::vector<double> operator*(const std::vector<std::vector<double> >& mat1, const std::vector<double>& vec1){
    int row = mat1.size();
    int col = vec1.size();
    std::vector<double> vec(row, 0.);
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            vec[i] += mat1[i][j] * vec1[j];
        }
    }
    return vec;
}

std::vector<double> operator*(const std::vector<double>& vec1, const std::vector<std::vector<double> >& mat1){
    int row = mat1.size();
    int col = vec1.size();
    std::vector<double> vec(col, 0.);
    for(int i = 0; i < col; i++){
        for(int j = 0; j < row; j++){
            vec[i] += vec1[j] * mat1[j][i];
        }
    }
    return vec;
}

std::vector<double> operator/(const std::vector<double>& vec1, const double& value){
    int size = vec1.size();
    std::vector<double> vec(size, 0.);
    for(int i = 0; i < size; i++) vec[i] = vec1[i] / value;
    return vec;
}

std::vector<double> operator/(const double& value, const std::vector<double>& vec1){
    int size = vec1.size();
    std::vector<double> vec(size, 0.);
    for(int i = 0; i < size; i++) vec[i] = value / vec1[i];
    return vec;
}

void operator+=(std::vector<double>& vec1, const std::vector<double>& vec2){
    for(int i = 0; i < vec1.size(); i++) vec1[i] += vec2[i];
}

void operator+=(std::vector<double>& vec1, const double& value){
    for(int i = 0; i < vec1.size(); i++) vec1[i] += value;
}

void operator+=(std::vector<std::vector<double> >& mat1, const std::vector<std::vector<double> >& mat2){
    for(int i = 0; i < mat1.size(); i++){
        for(int j = 0; j < mat1[i].size(); j++){
            mat1[i][j] += mat2[i][j];
        }
    }
}

void operator-=(std::vector<double>& vec1, const std::vector<double>& vec2){
    for(int i = 0; i < vec1.size(); i++) vec1[i] -= vec2[i];
}

void operator-=(std::vector<double>& vec1, const double& value){
    for(int i = 0; i < vec1.size(); i++) vec1[i] -= value;
}

void operator-=(std::vector<std::vector<double> >& mat1, const std::vector<std::vector<double> >& mat2){
    for(int i = 0; i < mat1.size(); i++){
        for(int j = 0; j < mat1[i].size(); j++){
            mat1[i][j] -= mat2[i][j];
        }
    }
}

void operator*=(std::vector<double>& vec1, const double& value){
    for(int i = 0; i < vec1.size(); i++) vec1[i] *= value;
}

void operator/=(std::vector<double>& vec1, const double& value){
    for(int i = 0; i < vec1.size(); i++) vec1[i] /= value;
}

std::vector<double> zeros(const int& n){
    return std::vector<double>(n, 0.);
}

std::vector<double> ones(const int& n){
    return std::vector<double>(n, 1.);
}

std::vector<double> seq(int begin, int end){
    std::vector<double> vec;
    if(begin <= end){
        vec.resize(end - begin + 1);
        for(int i = begin; i <= end; i++) vec[i - begin] = i;
    }else{
        vec.resize(begin - end + 1);
        for(int i = begin; i >= end; i--) vec[begin - i] = i;
    }
    return vec;
}

std::vector<double> log(const std::vector<double>& vec1){
    std::vector<double> vec(vec1);
    for(int i = 0; i < vec.size(); i++) vec[i] = std::log(vec[i]);
    return vec;
}

std::vector<double> sqrt(const std::vector<double>& vec1){
    std::vector<double> vec(vec1);
    for(int i = 0; i < vec.size(); i++) vec[i] = std::sqrt(vec[i]);
    return vec;
}

double sum(const std::vector<double>& vec){
    double sum = 0.;
    for(int i = 0; i < vec.size(); i++) sum += vec[i];
    return sum;
}

double sum(const std::vector<double>& vec, const int& p){
    double sum = 0.;
    for(int i = 0; i < vec.size(); i++) sum += std::pow(vec[i], p);
    return sum;
}

double inner_product(const std::vector<double>& vec1, const std::vector<double>& vec2){
    double value = 0.;
    for(int i = 0; i < vec1.size(); i++) value += vec1[i] * vec2[i];
    return value; 
}

double norm(const std::vector<double>& vec){
    double norm = 0.;
    for(int i = 0; i < vec.size(); i++) norm += vec[i] * vec[i];
    return std::sqrt(norm);
}

std::vector<double> reciprocal(const std::vector<double>& vec){
    std::vector<double> v(vec.size());
    for(int i = 0; i < v.size(); i++) v[i] = 1. / vec[i];
    return v;
}

void normalization(std::vector<double>& vec){
    vec /= norm(vec);
}

std::vector<std::vector<double> > zeros(int row, int col){
    return std::vector<std::vector<double> >(row, std::vector<double>(col, 0));
}

std::vector<std::vector<double> > ones(int row, int col){
    return std::vector<std::vector<double> >(row, std::vector<double>(col, 1));
}

std::vector<double> col(const std::vector<std::vector<double> >& mat, const int& c){
    int col_vector_size = mat.size();
    std::vector<double> v(col_vector_size);
    for(int i = 0; i < col_vector_size; i++) v[i] = mat[i][c];
    return v;
}

std::vector<std::vector<double> > identity(int n){
    std::vector<std::vector<double> > A(n, std::vector<double>(n, 0));
    for(int i = 0; i < n; i++) A[i][i] = 1.;
    return A;
}

std::vector<std::vector<double> > trans(const std::vector<std::vector<double> >& A){
    int row = A.size();
    int col = A.front().size();
    std::vector<std::vector<double> > At(col, std::vector<double>(row, 0));
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            At[j][i] = A[i][j];
        }
    }
    return At;
}

void abs(std::vector<std::vector<double> >& A){
    int row = A.size();
    int col = A.front().size();
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            A[i][j] = std::fabs(A[i][j]);
        }
    }
}

void diag_zero(std::vector<std::vector<double> >& A){
    int n = A.size();
    for(int i = 0; i < n; i++){
        A[i][i] = 0;
    }
}

double non_diag_abs_max_value(const std::vector<std::vector<double> >& A){
    int n = A.size();
    double max = -1e8;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(max < std::fabs(A[i][j])& & i != j){
                max = std::fabs(A[i][j]);
            }
        }
    }
    return max;
}

std::pair<int, int> non_diag_abs_max_index(const std::vector<std::vector<double> >& A){
    int n = A.size();
    double max = -1e8;
    int p, q;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(max < std::fabs(A[i][j])& & i != j){
                max = std::fabs(A[i][j]);
                p = i, q = j;
            }
        }
    }
    return std::make_pair(p, q);
}

std::vector<std::vector<double> > inverse(const std::vector<std::vector<double> >& mat){

    const int n = mat.size();
    std::vector<std::vector<double> > A(mat);
    std::vector<std::vector<double> > Ai(identity(n));

    for(int r = 0; r < n; r++){

        int max_row;
        double max = 0.;
        for(int i = r; i < n; i++){
            if(max < std::fabs(A[i][r])){
                max_row = i;
                max = std::fabs(A[i][r]);
            }
        }

        if(r != max_row){
            for(int c = 0; c < n; c++){
                std::swap(A[r][c], A[max_row][c]);
                std::swap(Ai[r][c], Ai[max_row][c]);
            }
        }

        double a = A[r][r];
        for(int c = 0; c < n; c++){
            A[r][c] /= a;
            Ai[r][c] /= a;
        }

        for(int i = 0; i < n; i++){
            if(i == r) continue;
            a = A[i][r];
            for(int c = 0; c < n; c++){
                A[i][c] -= a * A[r][c];
                Ai[i][c] -= a * Ai[r][c];
            }
        }
    }
    return Ai;
}

std::vector<std::vector<double> > upper_triangle(const std::vector<std::vector<double> >& A, const int& p = 0){
    int n = A.size();
    std::vector<std::vector<double> > mat(zeros(n, n));
    for(int i = 0; i < n; i++){
        for(int j = i + p; j < n; j++){
            mat[i][j] = A[i][j];
        }
    }
    return mat;
}

std::vector<std::vector<double> > lower_triangle(const std::vector<std::vector<double> >& A, const int& p = 0){
    int n = A.size();
    std::vector<std::vector<double> > mat(zeros(n, n));
    for(int i = p; i < n; i++){
        for(int j = 0; j <= i - p; j++){
            mat[i][j] = A[i][j];
        }
    }
    return mat;
}

std::vector<double> diagonalization_component(const std::vector<std::vector<double> >& A){
    std::vector<double> v(zeros(A.size()));
    for(int i = 0; i < v.size(); i++) v[i] = A[i][i];
    return v;
}

std::vector<std::vector<double> > diagonalization_matrix(std::vector<double>& v){
    int n = v.size();
    std::vector<std::vector<double> > M(zeros(n, n));
    for(int i = 0; i < n; i++) M[i][i] = v[i];
    return M;
}

void LU_decomposition(const std::vector<std::vector<double> >& M, std::vector<std::vector<double> >& L, std::vector<std::vector<double> >& U){

    const int n = M.size();
    std::vector<std::vector<double> > A(M);
    U = identity(n);

    for(int i = 0; i < n; i++){

        int N = n - i - 1;

        double l0 = L[i][i] = A[0][0];

        std::vector<double> l1(zeros(n));
        for(int j = 0; j != N; ++j){
            L[j + i + 1][i] = l1[j] = A[j + 1][0];
        }

        std::vector<double> u1(zeros(n));
        for(int j = 0; j != N; ++j){
            U[i][j + i + 1] = u1[j] = A[0][j + 1] / l0;
        }

        std::vector<std::vector<double> > lu(zeros(n, n));
        for(int j = 0; j != N; ++j){
            for(int k = 0; k != N; ++k){
                lu[j][k] = l1[j] * u1[k];
            }
        }

        std::vector<std::vector<double> > A1(zeros(n, n));
        for(int j = 0; j != N; ++j){
            for(int k = 0; k != N; ++k){
                A1[j][k] = A[j + 1][k + 1] - lu[j][k];
            }
        }

        A = A1;
    }
}

std::vector<std::vector<double> > hessenberg(const std::vector<std::vector<double> >& mat){

    int n = mat.size();
    std::vector<std::vector<double> > H(mat);

    for(int k = 1; k < n - 1; k++){
        std::vector<double> u(zeros(n));
        for(int i = k; i != n; ++i) u[i] = H[i][k - 1];

        double ss = 0;
        for(int i = k + 1; i < n; i++) ss += u[i] * u[i];
        if(std::fabs(ss) <= 0.) continue;
        double s = sqrt(ss + u[k] * u[k]);
        if(u[k] > 0.) s = -s;

        u[k] -= s;
        double uu = std::sqrt(ss + u[k] * u[k]);
        for(int i = k; i != n; ++i) u[i] /= uu;

        std::vector<double> f(zeros(n)), g(zeros(n));
        for(int i = 0; i < n; i++){
            for(int j = k; j < n; j++){
                f[i] += H[i][j] * u[j];
                g[i] += H[j][i] * u[j];
            }
        }
        double gamma = inner_product(u, g);

        for(int i = 0; i < n; i++){
            f[i] -= gamma * u[i];
            g[i] -= gamma * u[i];
        }

        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                H[i][j] = H[i][j] - 2. * u[i] * g[j] - 2. * f[i] * u[j];
            }
        }
    }
    return H;
}

std::vector<double> eigen_values(const std::vector<std::vector<double> >& mat){

    int n = mat.size();
    double eps = 1e-8;
    std::vector<std::vector<double> > M(mat);
    std::vector<std::vector<double> > H = hessenberg(M);
    std::vector<double> s(n), c(n);
    for(int m = n; m >= 2;){
        if(std::fabs(H[m - 1][m - 2]) < eps){
            m--;
            continue;
        }
        double shift = H[m - 1][m - 1];
        for(int i = 0; i < m; i++) H[i][i] -= shift;
        for(int k = 0; k < m - 1; k++){
            double a = H[k][k], b = H[k + 1][k], r = std::sqrt(a * a + b * b);
            s[k] = r == 0. ? 0. : b / r;
            c[k] = r == 0. ? 0. : a / r;
            for(int j = k; j < m; j++){
                double x = H[k][j], y = H[k + 1][j];
                H[k][j] = c[k] * x + s[k] * y;
                H[k + 1][j] = -s[k] * x + c[k] * y;
            }
        }
        for(int k = 0; k < m - 1; k++){
            for(int i = 0; i < k + 2; i++){
                double x = H[i][k], y = H[i][k + 1];
                H[i][k] = c[k] * x + s[k] * y;
                H[i][k + 1] = -s[k] * x + c[k] * y;
            }
        }
        for(int i = 0; i < m; i++) H[i][i] += shift;
    }

    std::vector<double> lambda(zeros(n));
    for(int i = 0; i != n; ++i) lambda[i] = H[i][i];

    std::sort(lambda.begin(), lambda.end());

    return lambda;
}

std::vector<double> eigen_vector(const std::vector<std::vector<double> >& mat, const double& lambda){

    double eps = 1e-6;
    const int n = mat.size();
    std::vector<std::vector<double> > A(mat);
    std::vector<double> eigenvec(zeros(n)); eigenvec[0] = 1.;
    for(int i = 0; i < n; i++) A[i][i] -= lambda;

    double muo = 0., mu = 0.;
    do{
        muo = mu;
        std::vector<double> next = inverse(A) * eigenvec;
        mu = inner_product(next, eigenvec);
        normalization(next);
        eigenvec = next;    
    }while(std::fabs((mu - muo) / mu) > eps);

    return eigenvec;
}

void diagonalization(const std::vector<std::vector<double> >& mat, std::vector<std::vector<double> >& V, std::vector<std::vector<double> >& D){
    
    const int n = mat.size();
    std::vector<double> lambda = eigen_values(mat);
    D = diagonalization_matrix(lambda); 
    if(lambda == diagonalization_component(mat)){
        V = identity(n);
    }else{
        for(int i = 0; i < n; i++){
            std::vector<double> ev = eigen_vector(mat, lambda[i]);
            for(int j = 0; j < n; j++) V[j][i] = ev[j];
        }
    }
}

double similar_trans(std::vector<std::vector<double> >& A, const std::pair<int, int>& index){

    double eps = 1e-8;
    double phi = 0.;
    int p = index.first, q = index.second;
    double App = A[p][p], Apq = A[p][q], Aqq = A[q][q];

    if(std::fabs(App - Aqq) < eps){
        phi = M_PI / 4.;
    }else{
        phi = 0.5 * std::atan((2. * Apq) / (App - Aqq));
    }

    double cos = std::cos(phi), sin = std::sin(phi), coscos = cos * cos, sinsin = sin * sin, sincos = sin * cos;

    A[p][p] = (coscos * App) + (2. * sincos * Apq) + (sinsin * Aqq);
    A[q][q] = (sinsin * App) - (2. * sincos * Apq) + (coscos * Aqq);
    A[p][q] = 0.;
    A[q][p] = 0.;

    for(int i = 0; i < A.size(); i++){
        if(i != p& & i != q){
            double Api = A[p][i], Aqi = A[q][i];
            A[i][p] = (cos * Api) + (sin * Aqi);
            A[p][i] = A[i][p];
            A[i][q] = (cos * Aqi) - (sin * Api);
            A[q][i] = A[i][q];
        }
    }

    return phi;
}

void givens_rot(std::vector<std::vector<double> >& G, const std::pair<int, int>& index, const double& phi){

    int p = index.first, q = index.second;

    for(int i = 0; i < G.size(); i++){
        G[i][p] = std::cos(phi) * G[p][i] + std::sin(phi) * G[q][i];
        G[i][q] = std::cos(phi) * G[q][i] - std::sin(phi) * G[p][i];
    }
}

void jacobi(const std::vector<std::vector<double> >& A, std::vector<std::vector<double> >& R, std::vector<std::vector<double> >& D){
    
    double eps = 1e-12;
    D = A;
    R = identity(R.size());
    while(non_diag_abs_max_value(D) >= eps){
        std::pair<int, int> max_index = non_diag_abs_max_index(D);
        double phi = similar_trans(D, max_index);
        givens_rot(R, max_index, phi);
    }
}

void output(const std::vector<double>& v){
    std::cout << std::endl;
    for(int i = 0; i < v.size(); i++){
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

void output(const std::vector<std::vector<double> >& A){
    std::cout << std::endl;
    for(int i = 0; i < A.size(); i++){
        for(int j = 0; j < A.front().size(); j++){
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

#endif