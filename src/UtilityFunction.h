#ifndef UTILITY_FUNCTION_H
#define UTILITY_FUNCTION_H

#include <cstdlib>
#include <vector>

int PickNumber(const int, const int);
std::size_t PickSizeT(const std::size_t, const std::size_t);
double PickProbabilty();
double PickDouble(const double, const double);
void CalculateMaxAndAvg(const std::vector<double> &, double &, double &);
void CalculateMinMaxAvg(const std::vector<std::size_t> &, std::size_t &, std::size_t &, std::size_t &);
void InsertionSort(std::vector<double> &, std::vector<std::size_t> &);
void RunTournament(std::vector<double> &, std::vector<std::size_t> &);
void InsertionSortIncrease(std::vector<std::size_t> &, std::vector<std::size_t> &);
void CumSum(const double * input, double * cum_sum, const int end_index);
int FindFirstValueGreaterOrEqual(const std::vector<double> &, const std::vector<double> &);

#endif