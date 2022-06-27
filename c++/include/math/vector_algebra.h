
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Library defining possible operations among std::vector<>.
// last updated:    19/06/2022
                                                                                         
#pragma once
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>


// =============================================                                                                                         
// Utilities
// =============================================  

std::vector<double> zeros(const int& n) {
    return std::vector<double>(n, 0.);
}

std::vector<std::vector<double>> zeros(const int& n_rows, const int& n_cols) {
    return std::vector<double>(n_rows, std::vector<double>(n_cols, 0.));
}

std::vector<double> ones(const int& n) {
    return std::vector<double>(n, 1.);
}

std::vector<std::vector<double>> ones(const int& n_rows, const int& n_cols) {
    return std::vector<double>(n_rows, std::vector<double>(n_cols, 1.));
}


// =============================================                                                                                         
// Sum
// =============================================  

// sum of a vector and a scalar
std::vector<double> operator+(const std::vector<double>& vec1, const double& value) {
    std::vector<double> vec(vec1.size(), 0.);
    for(unsigned int i{}; i < vec1.size(); i++) vec[i] = vec1[i] + value;
    return vec;
}

// sum of a scalar and a vector  
std::vector<double> operator+(const double& value, const std::vector<double>& vec1) {
    std::vector<double> vec(vec1.size(), 0.);
    for(unsigned int i{}; i < vec1.size(); i++) vec[i] = value + vec1[i];
    return vec;
}

// sum of two vectors
std::vector<double> operator+(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<double> vec(vec1.size(), 0.);
    for(int i{}; i < vec1.size(); i++) vec[i] = vec1[i] + vec2[i];
    return vec;
}

// sum of a multidimentional vector and a scalar
std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& mat1, const double& value) {
    std::vector<std::vector<double>> mat(mat1.size(), std::vector<double>(mat1.front().size(), 0));
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = mat1[i][j] + value;
        }
    }
    return mat;
}

// sum of a scalar and a multidimentional vector 
std::vector<std::vector<double>> operator+(const double& value, const std::vector<std::vector<double>>& mat1) {
    std::vector<std::vector<double>> mat(mat1.size(), std::vector<double>(mat1.front().size(), 0));
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = value + mat1[i][j];
        }
    }
    return mat;
}

// sum of two multidimentional vectors
std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2) {
    std::vector<std::vector<double>> mat(mat1.size(), std::vector<double>(mat1.front().size(), 0));
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = mat1[i][j] + mat2[i][j];
        }
    }
    return mat;
}

// increase vector with a scalar
void operator+=(const std::vector<double>& vec1, const double& value) {
    for(int i = 0; i < vec1.size(); i++) vec1[i] += value;
}

// increase vector with a vector
void operator+=(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    for(int i = 0; i < vec1.size(); i++) vec1[i] += vec2[i];
}

// increase multidimentional vector with a scalar
void operator+=(const std::vector<std::vector<double>>& mat1, const double& value) {
    for(int i = 0; i < mat1.size(); i++) {
        for(int j = 0; j < mat1[i].size(); j++) {
            mat1[i][j] += value;
        }
    }
}

// increase multidimentional vector with a multidimentional vector
void operator+=(const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2) {
    for(int i = 0; i < mat1.size(); i++) {
        for(int j = 0; j < mat1[i].size(); j++) {
            mat1[i][j] += mat2[i][j];
        }
    }
}


// =============================================                                                                                         
// Subtraction
// =============================================  

// subtraction of a vector and a scalar 
std::vector<double> operator-(const std::vector<double>& vec1, const double& value) {
    std::vector<double> vec(vec1.size(), 0.);
    for(int i{}; i < vec1.size(); i++) vec[i] = vec1[i] - value;
    return vec;
}

// subtraction of a scalar and a vector 
std::vector<double> operator-(const double& value, const std::vector<double>& vec1) {
    std::vector<double> vec(vec1.size(), 0.);
    for(int i{}; i < vec1.size(); i++) vec[i] = value - vec1[i];
    return vec;
}

// subtraction of two vectors
std::vector<double> operator-(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<double> vec(vec1.size(), 0.);
    for(int i{}; i < vec1.size(); i++) vec[i] = vec1[i] - vec2[i];
    return vec;
}

// subtraction of a multidimentional vector and a scalar  
std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& mat1, const double& value) {
    std::vector<std::vector<double>> mat(mat1.size(), std::vector<double>(mat1.front().size(), 0));
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = mat1[i][j] - value;
        }
    }
    return mat;
}

// subtraction of a scalar and a multidimentional vector 
std::vector<std::vector<double>> operator-(const double& value, const std::vector<std::vector<double>>& mat1) {
    std::vector<std::vector<double>> mat(mat1.size(), std::vector<double>(mat1.front().size(), 0));
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = value - mat1[i][j];
        }
    }
    return mat;
}

// subtraction of two multidimentional vectors
std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2) {
    std::vector<std::vector<double>> mat(mat1.size(), std::vector<double>(mat1.front().size(), 0));
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
    return mat;
}

// decrease vector with a scalar
void operator-=(const std::vector<double>& vec1, const double& value) {
    for(int i = 0; i < vec1.size(); i++) vec1[i] -= value;
}

// decrease vector with a vector
void operator-=(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    for(int i = 0; i < vec1.size(); i++) vec1[i] -= vec2[i];
}

// decrease multidimentional vector with a scalar
void operator-=(const std::vector<std::vector<double>>& mat1, const double& value) {
    for(int i = 0; i < mat1.size(); i++) {
        for(int j = 0; j < mat1[i].size(); j++) {
            mat1[i][j] -= value;
        }
    }
}

// decrease multidimentional vector with a multidimentional vector
void operator-=(const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2) {
    for(int i = 0; i < mat1.size(); i++) {
        for(int j = 0; j < mat1[i].size(); j++) {
            mat1[i][j] -= mat2[i][j];
        }
    }
}


// =============================================                                                                                         
// Moltiplication
// =============================================  

// moltiplication of a vector and a scalar 
std::vector<double> operator*(const std::vector<double>& vec1, const double& value) {
    std::vector<double> vec(vec1.size(), 0.);
    for(int i{}; i < vec1.size(); i++) vec[i] = vec1[i] * value;
    return vec;
}

// moltiplication of a scalar and a vector 
std::vector<double> operator*(const double& value, const std::vector<double>& vec1) {
    std::vector<double> vec(vec1.size(), 0.);
    for(int i{}; i < vec1.size(); i++) vec[i] = value * vec1[i];
    return vec;
}

// // moltiplication of a vector and a scalar
// void operator*=(const std::vector<double>& vec1, const double& value) {
//     for(int i = 0; i < vec1.size(); i++) vec1[i] *= value;
// }

// moltiplication of two vectors
std::vector<std::vector<double>> operator*(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<std::vector<double>> mat(vec1.size(), std::vector<double>(vec2.size(), 0));
    for(int i{}; i < vec1.size(); i++) {
        for(int j{}; j < vec2.size(); j++) {
            mat[i][j] = vec1[i] * vec2[j];
        }
    }
    return mat;
}

// moltiplication of a multidimentional vector and a scalar
std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& mat1, const double& value) {
    std::vector<std::vector<double>> mat( mat1.size(), std::vector<double>(mat1.front().size(), 0));
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = mat1[i][j] * value;
        }
    }
    return mat;
}

// moltiplication of a scalar and a multidimentional vector 
std::vector<std::vector<double>> operator*(const double& value, const std::vector<std::vector<double>>& mat1) {
    std::vector<std::vector<double>> mat( mat1.size(), std::vector<double>(mat1.front().size(), 0));
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = value * mat1[i][j];
        }
    }
    return mat;
}

// // moltiplication of a multidimentional vector and a scalar
// void operator*=(const std::vector<std::vector<double>>& mat1, const double& value) {
//     for(int i{}; i < mat1.size(); i++) {
//         for(int j{}; j < mat1.front().size(); j++) {
//             mat1[i][j] *= value;
//         }
//     }
// }

// moltiplication of two multidimentional vectors
std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2) {
    std::vector<std::vector<double>> mat(mat1.size(), std::vector<double>(mat2.front().size(), 0));
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < mat1.front().size(); j++) {
            for(int k{}; k < mat2.front().size(); k++) {
                mat[i][k] += mat1[i][j] * mat2[j][k];
            }
        }
    }
    return mat;
}

// moltiplication of a multidimentional vector and a vector
std::vector<double> operator*(const std::vector<std::vector<double>>& mat1, const std::vector<double>& vec1) {
    std::vector<double> vec(mat1.size(), 0.);
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < vec1.size(); j++) {
            vec[i] += mat1[i][j] * vec1[j];
        }
    }
    return vec;
}

// moltiplication of a vector and a multidimentional vector 
std::vector<double> operator*(const std::vector<double>& vec1, const std::vector<std::vector<double>>& mat1) {
    std::vector<double> vec(vec1.size(), 0.);
    for(int i{}; i < vec1.size(); i++) {
        for(int j{}; j < mat1.size(); j++) {
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
    std::vector<double> vec(vec1.size(), 0.);
    for(int i{}; i < vec1.size(); i++) vec[i] = vec1[i] / value;
    return vec;
}

// division of a scalar and a vector 
std::vector<double> operator/(const double& value, const std::vector<double>& vec1) {
    std::vector<double> vec(vec1.size(), 0.);
    for(int i{}; i < vec1.size(); i++) vec[i] = value / vec1[i];
    return vec;
}

// // division of a multidimentional vector and a scalar
// void operator/=(const std::vector<std::vector<double>>& mat1, const double& value) {
//     for(int i{}; i < mat1.size(); i++) {
//         for(int j{}; j < mat1.front().size(); j++) {
//             mat1[i][j] /= value;
//         }
//     }
// }

// division of two vectors
std::vector<std::vector<double>> operator/(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<std::vector<double>> mat(vec1.size(), std::vector<double>(vec2.size(), 0));
    for(int i{}; i < vec1.size(); i++) {
        for(int j{}; j < vec2.size(); j++) {
            mat[i][j] = vec1[i] / vec2[j];
        }
    }
    return mat;
}

// division of a multidimentional vector and a scalar
std::vector<std::vector<double>> operator/(const std::vector<std::vector<double>>& mat1, const double& value) {
    std::vector<std::vector<double>> mat( mat1.size(), std::vector<double>(mat1.front().size(), 0));
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = mat1[i][j] / value;
        }
    }
    return mat;
}

// division of a scalar and a multidimentional vector 
std::vector<std::vector<double>> operator/(const double& value, const std::vector<std::vector<double>>& mat1) {
    std::vector<std::vector<double>> mat( mat1.size(), std::vector<double>(mat1.front().size(), 0));
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < mat1.front().size(); j++) {
            mat[i][j] = value / mat1[i][j];
        }
    }
    return mat;
}

// // division of a multidimentional vector and a scalar
// void operator/=(const std::vector<std::vector<double>>& mat1, const double& value) {
//     for(int i{}; i < mat1.size(); i++) {
//         for(int j{}; j < mat1.front().size(); j++) {
//             mat1[i][j] /= value;
//         }
//     }
// }

// division of two multidimentional vectors
std::vector<std::vector<double>> operator/(const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2) {
    std::vector<std::vector<double>> mat(mat1.size(), std::vector<double>(mat2.front().size(), 0));
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < mat1.front().size(); j++) {
            for(int k{}; k < mat2.front().size(); k++) {
                mat[i][k] += mat1[i][j] / mat2[j][k];
            }
        }
    }
    return mat;
}

// division of a multidimentional vector and a vector
std::vector<double> operator/(const std::vector<std::vector<double>>& mat1, const std::vector<double>& vec1) {
    std::vector<double> vec(mat1.size(), 0.);
    for(int i{}; i < mat1.size(); i++) {
        for(int j{}; j < vec1.size(); j++) {
            vec[i] += mat1[i][j] / vec1[j];
        }
    }
    return vec;
}

// division of a vector and a multidimentional vector 
std::vector<double> operator/(const std::vector<double>& vec1, const std::vector<std::vector<double>>& mat1) {
    std::vector<double> vec(vec1.size(), 0.);
    for(int i{}; i < vec1.size(); i++) {
        for(int j{}; j < mat1.size(); j++) {
            vec[i] += vec1[j] / mat1[j][i];
        }
    }
    return vec;
}

// other functions to be implemented
// std::vector<double> reciprocal(const std::vector<double>& vec) {
//     std::vector<double> v(vec.size());
//     for(int i = 0; i < v.size(); i++) v[i] = 1. / vec[i];
//     return v;
// }

// std::vector<double> col(const std::vector<std::vector<double> >& mat, const int& c) {
//     std::vector<double> v(mat.size());
//     for(int i = 0; i < mat.size(); i++) v[i] = mat[i][c];
//     return v;
// }

// void abs(std::vector<std::vector<double> >& A) {
//     int row = A.size();
//     int col = A.front().size();
//     for(int i = 0; i < row; i++) {
//         for(int j = 0; j < col; j++) {
//             A[i][j] = std::fabs(A[i][j]);
//         }
//     }
// }

// void normalization(std::vector<double>& vec) {
//     vec /= norm(vec);
// }

// double inner_product(const std::vector<double>& vec1, const std::vector<double>& vec2) {
//     double value = 0.;
//     for(int i = 0; i < vec1.size(); i++) value += vec1[i] * vec2[i];
//     return value; 
// }

// std::vector<double> log(const std::vector<double>& vec1) {
//     std::vector<double> vec(vec1);
//     for(int i = 0; i < vec.size(); i++) vec[i] = std::log(vec[i]);
//     return vec;
// }

// std::vector<double> sqrt(const std::vector<double>& vec1) {
//     std::vector<double> vec(vec1);
//     for(int i = 0; i < vec.size(); i++) vec[i] = std::sqrt(vec[i]);
//     return vec;
// }

// double sum(const std::vector<double>& vec) {
//     double sum = 0.;
//     for(int i = 0; i < vec.size(); i++) sum += vec[i];
//     return sum;
// }

// double sum(const std::vector<double>& vec, const int& p) {
//     double sum = 0.;
//     for(int i = 0; i < vec.size(); i++) sum += std::pow(vec[i], p);
//     return sum;
// }

