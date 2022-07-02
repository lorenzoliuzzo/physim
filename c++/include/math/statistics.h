
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Some statistical functions for data analysis. 
// last updated:    18/06/2022


#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>


// =============================================
// mean
// =============================================

template <typename T, typename IT>
T mean(IT begin, IT end) {
    if (end - begin == 0) return 0; 
    return std::accumulate(begin, end, 0.0) / (end - begin); 
}

template <typename T, typename K>
T mean(const std::vector<K>& v) {
    T accu{}; 
    if (v.size() == 0) return accu;
    for (T x : v) accu += x;
    return accu / v.size(); 
}


// =============================================
// median
// =============================================

template <typename T, typename IT>
T median(IT begin, IT end) {
    sort(begin, end); 
    if ((end - begin) % 2 != 0) return *(begin + (end - begin) / 2);
    else return ((*(begin + (end - begin) / 2) + *(begin + (end - begin) / 2 - 1)) / 2);
}

template <typename T, typename K>
T median(std::vector<K> v) {
    sort(v.begin(), v.end()); 
    if (v.size()%2 != 0) return v[v.size() / 2];
    else return (v[v.size() / 2] + v[(v.size() / 2)-1]) / 2; 
}


// =============================================
// variance
// =============================================

template <typename T, typename IT>
T var(IT begin, IT end) {
    if (end - begin == 0) return 0; 
    T accu{}; 
    for (IT i{begin}; i < end; i++) accu += pow(*i - mean<T>(begin, end), 2);
    return accu / (end - begin);
}

template <typename T, typename K>
T var(std::vector<K> v) {
    T accu{}; 
    if (v.size() == 0) return accu;
    for (auto x : v) accu += pow(x - mean<T>(v), 2);
    return accu / v.size();
}


// =============================================
// standard deviation
// =============================================

template <typename T, typename IT>
T sd(IT begin, IT end) {
    return sqrt(var<T>(begin, end));
}

template <typename T, typename K>
T sd(std::vector<K> v) {
    return sqrt(var<T>(v));
}


// =============================================
// standard deviation of mean
// =============================================

template <typename T, typename IT>
T sdom(IT begin, IT end) {
    return sd(begin, end) / sqrt(end - begin);
}

template <typename T, typename K>
T sdom(std::vector<K> v) {
    return sd(v) / sqrt(v.size());
}


// =============================================
// chi squared
// =============================================

template <typename T, typename IT, typename J>
T chi_sq(IT begin, IT end, J expected_value) {
    T accu{}; 
    for (auto x : end) accu += pow(x - expected_value, 2) / sd(begin, end); 
    return accu; 
}

template <typename T, typename K, typename J>
T chi_sq(std::vector<K> v, J expected_value) {
    T accu{}; 
    for (auto x : v) accu += pow(x - expected_value, 2) / sd(v); 
    return accu; 
}

