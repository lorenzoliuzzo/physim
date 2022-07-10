
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Library defining possible operations among std::vector<double> and std::vector<std::vector<double>>.
// last updated:    08/07/2022


#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>


// =============================================                                                                                         
// Utilities
// =============================================  

std::vector<double> zeros(const int& n) { return std::vector<double>(n, 0.); }

std::vector<std::vector<double>> zeros(const int& n_rows, const int& n_cols) { return std::vector<std::vector<double>>(n_rows, zeros(n_cols)); }

void print(const std::vector<double>& vec) {
    for (auto i : vec) std::cout << "[" << i << "]\t"; 
}

void print(const std::vector<std::vector<double>>& mat) {
    for (int i{}; i < mat.size(); i++) {
        for (int j{}; j < mat.front().size(); j++) {
            std::cout << "[" << mat[i][j] << "]\t";
        }
        std::cout << std::endl; 
    }   
}

// =============================================                                                                                         
// Sum
// =============================================  

// sum of a vector and a scalar
std::vector<double> operator+(const std::vector<double>& vec1, const double& value) {
    std::vector<double> vec = zeros(vec1.size());
    for (unsigned int i{}; i < vec1.size(); i++) vec[i] = vec1[i] + value;
    return vec;
}

// sum of a scalar and a vector  
std::vector<double> operator+(const double& value, const std::vector<double>& vec1) {
    std::vector<double> vec = zeros(vec1.size());
    for (unsigned int i{}; i < vec1.size(); i++) vec[i] = value + vec1[i];
    return vec;
}

// sum of two vectors
std::vector<double> operator+(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<double> vec = zeros(vec1.size());
    for (int i{}; i < vec1.size(); i++) vec[i] = vec1[i] + vec2[i];
    return vec;
}

// sum of a multidimentional vector and a scalar
std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& mat1, const double& value) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat1.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = mat1[i][j] + value;
        }
    }
    return mat;
}

// sum of a scalar and a multidimentional vector 
std::vector<std::vector<double>> operator+(const double& value, const std::vector<std::vector<double>>& mat1) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat1.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = value + mat1[i][j];
        }
    }
    return mat;
}

// sum of two multidimentional vectors
std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat1.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = mat1[i][j] + mat2[i][j];
        }
    }
    return mat;
}

// increase vector with a scalar
void operator+=(std::vector<double> vec1, const double& value) {
    for (int i{}; i < vec1.size(); i++) vec1[i] += value;
}

// increase vector with a vector
void operator+=(std::vector<double> vec1, const std::vector<double>& vec2) {
    for (int i{}; i < vec1.size(); i++) vec1[i] += vec2[i];
}

// increase multidimentional vector with a scalar
void operator+=(std::vector<std::vector<double>> mat1, const double& value) {
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1[i].size(); j++) {
            mat1[i][j] += value;
        }
    }
}

// increase multidimentional vector with a multidimentional vector
void operator+=(std::vector<std::vector<double>> mat1, const std::vector<std::vector<double>>& mat2) {
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1[i].size(); j++) {
            mat1[i][j] += mat2[i][j];
        }
    }
}


// =============================================                                                                                         
// Subtraction
// =============================================  

// subtraction of a vector and a scalar 
std::vector<double> operator-(const std::vector<double>& vec1, const double& value) {
    std::vector<double> vec = zeros(vec1.size());
    for (int i{}; i < vec1.size(); i++) vec[i] = vec1[i] - value;
    return vec;
}

// subtraction of a scalar and a vector 
std::vector<double> operator-(const double& value, const std::vector<double>& vec1) {
    std::vector<double> vec = zeros(vec1.size());
    for (int i{}; i < vec1.size(); i++) vec[i] = value - vec1[i];
    return vec;
}

// subtraction of two vectors
std::vector<double> operator-(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<double> vec = zeros(vec1.size());
    for (int i{}; i < vec1.size(); i++) vec[i] = vec1[i] - vec2[i];
    return vec;
}

// subtraction of a multidimentional vector and a scalar  
std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& mat1, const double& value) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat1.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = mat1[i][j] - value;
        }
    }
    return mat;
}

// subtraction of a scalar and a multidimentional vector 
std::vector<std::vector<double>> operator-(const double& value, const std::vector<std::vector<double>>& mat1) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat1.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = value - mat1[i][j];
        }
    }
    return mat;
}

// subtraction of two multidimentional vectors
std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat1.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
    return mat;
}

// decrease vector with a scalar
void operator-=(std::vector<double> vec1, const double& value) {
    for (int i{}; i < vec1.size(); i++) vec1[i] -= value;
}

// decrease vector with a vector
void operator-=(std::vector<double> vec1, const std::vector<double>& vec2) {
    for (int i{}; i < vec1.size(); i++) vec1[i] -= vec2[i];
}

// decrease multidimentional vector with a scalar
void operator-=(std::vector<std::vector<double>> mat1, const double& value) {
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1[i].size(); j++) {
            mat1[i][j] -= value;
        }
    }
}

// decrease multidimentional vector with a multidimentional vector
void operator-=(std::vector<std::vector<double>> mat1, const std::vector<std::vector<double>>& mat2) {
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1[i].size(); j++) {
            mat1[i][j] -= mat2[i][j];
        }
    }
}


// =============================================                                                                                         
// Moltiplication
// =============================================  

// moltiplication of a vector and a scalar 
std::vector<double> operator*(const std::vector<double>& vec1, const double& value) {
    std::vector<double> vec = zeros(vec1.size());
    for (int i{}; i < vec1.size(); i++) vec[i] = vec1[i] * value;
    return vec;
}

// moltiplication of a scalar and a vector 
std::vector<double> operator*(const double& value, const std::vector<double>& vec1) {
    std::vector<double> vec = zeros(vec1.size());
    for (int i{}; i < vec1.size(); i++) vec[i] = value * vec1[i];
    return vec;
}

// moltiplication of a vector and a scalar
void operator*=(std::vector<double> vec1, const double& value) {
    for (int i{}; i < vec1.size(); i++) vec1[i] *= value;
}

// moltiplication of two vectors
std::vector<std::vector<double>> operator*(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<std::vector<double>> mat = zeros(vec1.size(), vec2.size()); 
    for (int i{}; i < vec1.size(); i++) {
        for (int j{}; j < vec2.size(); j++) {
            mat[i][j] = vec1[i] * vec2[j];
        }
    }
    return mat;
}

// moltiplication of a multidimentional vector and a scalar
std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& mat1, const double& value) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat1.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = mat1[i][j] * value;
        }
    }
    return mat;
}

// moltiplication of a scalar and a multidimentional vector 
std::vector<std::vector<double>> operator*(const double& value, const std::vector<std::vector<double>>& mat1) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat1.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = value * mat1[i][j];
        }
    }
    return mat;
}

// moltiplication of a multidimentional vector and a scalar
void operator*=(std::vector<std::vector<double>> mat1, const double& value) {
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat1[i][j] *= value;
        }
    }
}

// moltiplication of two multidimentional vectors
std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat2.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            for (int k{}; k < mat2.front().size(); k++) {
                mat[i][k] += mat1[i][j] * mat2[j][k];
            }
        }
    }
    return mat;
}

// moltiplication of a multidimentional vector and a vector
std::vector<double> operator*(const std::vector<std::vector<double>>& mat1, const std::vector<double>& vec1) {
    std::vector<double> vec = zeros(mat1.size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < vec1.size(); j++) {
            vec[i] += mat1[i][j] * vec1[j];
        }
    }
    return vec;
}

// moltiplication of a vector and a multidimentional vector 
std::vector<double> operator*(const std::vector<double>& vec1, const std::vector<std::vector<double>>& mat1) {
    std::vector<double> vec = zeros(vec1.size());
    for (int i{}; i < vec1.size(); i++) {
        for (int j{}; j < mat1.size(); j++) {
            vec[i] += vec1[j] * mat1[j][i];
        }
    }
    return vec;
}


// =============================================                                                                                         
// Division
// =============================================  

// division of a vector and a scalar 
std::vector<double> operator/(const std::vector<double>& vec1, const double& value) {
    std::vector<double> vec = zeros(vec1.size());
    for (int i{}; i < vec1.size(); i++) vec[i] = vec1[i] / value;
    return vec;
}

// division of a scalar and a vector 
std::vector<double> operator/(const double& value, const std::vector<double>& vec1) {
    std::vector<double> vec = zeros(vec1.size());
    for (int i{}; i < vec1.size(); i++) vec[i] = value / vec1[i];
    return vec;
}

// division of a vector and a scalar
void operator/=(std::vector<double> vec1, const double& value) {
    for (int i{}; i < vec1.size(); i++) vec1[i] /= value;
}

// division of two vectors
std::vector<std::vector<double>> operator/(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<std::vector<double>> mat = zeros(vec1.size(), vec2.size());
    for (int i{}; i < vec1.size(); i++) {
        for (int j{}; j < vec2.size(); j++) {
            mat[i][j] = vec1[i] / vec2[j];
        }
    }
    return mat;
}

// division of a multidimentional vector and a scalar
std::vector<std::vector<double>> operator/(const std::vector<std::vector<double>>& mat1, const double& value) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat1.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = mat1[i][j] / value;
        }
    }
    return mat;
}

// division of a scalar and a multidimentional vector 
std::vector<std::vector<double>> operator/(const double& value, const std::vector<std::vector<double>>& mat1) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat1.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = value / mat1[i][j];
        }
    }
    return mat;
}

// division of a multidimentional vector and a scalar
void operator/=(std::vector<std::vector<double>> mat1, const double& value) {
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat1[i][j] /= value;
        }
    }
}

// division of two multidimentional vectors
std::vector<std::vector<double>> operator/(const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat2.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            for (int k{}; k < mat2.front().size(); k++) {
                mat[i][k] += mat1[i][j] / mat2[j][k];
            }
        }
    }
    return mat;
}

// division of a multidimentional vector and a vector
std::vector<double> operator/(const std::vector<std::vector<double>>& mat1, const std::vector<double>& vec1) {
    std::vector<double> vec = zeros(mat1.size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < vec1.size(); j++) {
            vec[i] += mat1[i][j] / vec1[j];
        }
    }
    return vec;
}

// reciprocal of a vector
std::vector<double> reciprocal(const std::vector<double>& vec) {
    std::vector<double> v = zeros(vec.size());
    for (int i{}; i < v.size(); i++) v[i] = 1. / vec[i];
    return v;
}

// reciprocal of a multidimentional vector
std::vector<std::vector<double>> reciprocal(const std::vector<std::vector<double>>& mat1) {
    std::vector<std::vector<double>> mat = zeros(mat1.size(), mat1.front().size());
    for (int i{}; i < mat1.size(); i++) {
        for (int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = 1. / mat1[i][j];
        }
    }
    return mat;
}



// =============================================                                                                                         
// Other functions
// =============================================  



// std::vector<double> col(const std::vector<std::vector<double>>& mat, const int& c) {
//     std::vector<double> v(mat.size());
//     for (int i{}; i < mat.size(); i++) v[i] = mat[i][c];
//     return v;
// }

// void abs(std::vector<std::vector<double> >& A) {
//     int row = A.size();
//     int col = A.front().size();
//     for (int i{}; i < row; i++) {
//         for (int j{}; j < col; j++) {
//             A[i][j] = std::fabs(A[i][j]);
//         }
//     }
// }

// void normalization(std::vector<double>& vec) {
//     vec /= norm(vec);
// }

// double inner_product(const std::vector<double>& vec1, const std::vector<double>& vec2) {
//     double value = 0.;
//     for (int i{}; i < vec1.size(); i++) value += vec1[i] * vec2[i];
//     return value; 
// }

// std::vector<double> log(const std::vector<double>& vec1) {
//     std::vector<double> vec(vec1);
//     for (int i{}; i < vec.size(); i++) vec[i] = std::log(vec[i]);
//     return vec;
// }

// std::vector<double> sqrt(const std::vector<double>& vec1) {
//     std::vector<double> vec(vec1);
//     for (int i{}; i < vec.size(); i++) vec[i] = std::sqrt(vec[i]);
//     return vec;
// }

// double sum(const std::vector<double>& vec) {
//     double sum = 0.;
//     for (int i{}; i < vec.size(); i++) sum += vec[i];
//     return sum;
// }

// double sum(const std::vector<double>& vec, const int& p) {
//     double sum = 0.;
//     for (int i{}; i < vec.size(); i++) sum += std::pow(vec[i], p);
//     return sum;
// }

