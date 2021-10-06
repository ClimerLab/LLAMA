#include "UtilityFunction.h"
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <assert.h>

int PickNumber(const int min, const int max)
{
	if (min == max)
	{
		return min;
	}
	else if (max < min)
	{
		printf("ERROR in PickNumber - Max value is smaller than min value.\n");
		return -1;
	}
	else
	{
		return rand() % (max + 1 - min) + min;
	}
}

std::size_t PickSizeT(const std::size_t min, const std::size_t max)
{
	assert(min <= max);

	if(min == max)
	{
		return min;		
	}
	else
	{
		return (std::size_t)rand()% (max+1-min) + min;
	}	
}

double PickProbabilty()
{
	return (double)rand() / (double)RAND_MAX;
}

double PickDouble(const double min, const double max)
{
	return PickProbabilty() * (max - min) + min;
}

void CalculateMaxAndAvg(const std::vector<double> &vec, double &max, double &avg)
{
	// Initialize local variables
	double sum = 0;
	max = vec[0];

	// Loop through all values in array
	for (std::size_t i = 0; i < vec.size(); ++i)
	{
		// Add the value i-th element to the sum
		sum += vec[i];

		// Check if the i-th element is larger than the current max
		if (vec[i] > max)
			max = vec[i];
	}

	avg = sum / vec.size();
}

void CalculateMinMaxAvg(const std::vector<std::size_t> &vec, std::size_t &min, std::size_t &max, std::size_t &avg)
{
	// Initialize local variables
	int sum = 0;
	max = vec[0];

	// Loop through all values in array
	for (std::size_t i = 0; i < vec.size(); ++i)
	{
		// Add the value i-th element to the sum
		sum += vec[i];

		// Check if the i-th element is larger than the current max
		if (vec[i] > max)
			max = vec[i];
		
		// Check if the i-th element is smaller than the current min
		if (vec[i] < min)
			min = vec[i];		
	}

	avg = sum / vec.size();
}

void InsertionSort(std::vector<double> &A, std::vector<std::size_t> &B)
{
	assert(A.size() == B.size());

	std::size_t i, j, tmpB;
	double tmpA;

	i = 0;
	while (i < A.size())
	{
		j = i;
		while ((j > 0) && (A[j - 1] < A[j]))
		{
			tmpA = A[j];
			A[j] = A[j - 1];
			A[j - 1] = tmpA;

			tmpB = B[j];
			B[j] = B[j - 1];
			B[j - 1] = tmpB;
			j--;
		}
		i++;
	}
}

void RunTournament(std::vector<double> &fitness, std::vector<std::size_t> &tour_pool)
{
	// Declare local variables
	int tmp_index, num_picked;
	bool match_found = false;
	std::vector<double> tour_fitness(tour_pool.size());

	// Initialize number of tournament participants picked
	num_picked = 0;

	// Pick values from array A at random until tournament pool is full
	while (num_picked < tour_pool.size())
	{
		// Reset flag
		match_found = false;

		// Pick index at random
		tmp_index = PickNumber(0, fitness.size() - 1);

		// Loop through the first 'num_picked' indicies
		for (int i = 0; i < num_picked; i++)
		{
			// Check if the i-th values matched the tmp_indes
			if (tour_pool[i] == tmp_index) // Match found
			{
				match_found = true; // Set flag
				break; // Break and pick new number
			}
		}

		// Check if no matches were found
		if (!match_found)
		{
			tour_pool[num_picked] = tmp_index; // Add picked value to pool
			tour_fitness[num_picked] = fitness[tmp_index]; // Add the fitness of the corresponding index to the array
			num_picked++; // Increment the number of tournament participants picked
		}
	}

	// Sort tournament pool
	InsertionSort(tour_fitness, tour_pool);
}

void InsertionSortIncrease(int * A, int * B, const int array_length)
{
	int i, j, tmpA, tmpB;

	i = 0;
	while (i < array_length)
	{
		j = i;
		while ((j > 0) && (A[j - 1] > A[j]))
		{
			tmpA = A[j];
			A[j] = A[j - 1];
			A[j - 1] = tmpA;

			tmpB = B[j];
			B[j] = B[j - 1];
			B[j - 1] = tmpB;
			j--;
		}
		i++;
	}
}

void CumSum(const double *input, double *cum_sum, const int end_index)
{
	cum_sum[0] = input[0];

	for (int i = 1; i <= end_index; i++)
		cum_sum[i] = cum_sum[i - 1] + input[i];
}

int FindFirstValueGreaterOrEqual(const double in_value, const double *in_array, const int array_length)
{
	int index = -1;

	for (int i = 0; i < array_length; i++)
	{
		if (in_value <= in_array[i])
		{
			index = i;
			break;
		}
	}

	return index;
}
