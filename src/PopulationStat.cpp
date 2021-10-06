#include "PopulationStat.h"
#include <cstdio>
#include <assert.h>

PopulationStat::PopulationStat() :	max_num_gen(0)
{	
}

PopulationStat::PopulationStat(const std::size_t num_gen) :	max_num_gen(num_gen)
{
	max_fitness.resize(num_gen);
	avg_fitness.resize(num_gen);
	min_num_clust.resize(num_gen);
	max_num_clust.resize(num_gen);
	avg_num_clust.resize(num_gen);	
}

PopulationStat::~PopulationStat()
{	
}

void PopulationStat::setMaxFitness(const std::size_t index, const double fitness)
{
	assert(index < max_num_gen);
	max_fitness[index] = fitness;
}

void PopulationStat::setAvgFitness(const std::size_t index, const double fitness)
{
	assert(index < max_num_gen);
	avg_fitness[index] = fitness;
}

void PopulationStat::setMinNumClust(const std::size_t index, const std::size_t num_clust)
{
	assert(index < max_num_gen);
	min_num_clust[index] = num_clust;
}

void PopulationStat::setMaxNumClust(const std::size_t index, const std::size_t num_clust)
{
	assert(index < max_num_gen);
	max_num_clust[index] = num_clust;
}

void PopulationStat::setAvgNumClust(const std::size_t index, const std::size_t num_clust)
{
	assert(index < max_num_gen);
	avg_num_clust[index] = num_clust;
}

void PopulationStat::setMaxNumGen(const std::size_t num_gen)
{
	max_num_gen = num_gen;
}

void PopulationStat::setFitness(const std::size_t index, const double _max_fitness, const double _avg_fitness)
{
	assert(index < max_num_gen);
	max_fitness[index] = _max_fitness;
	avg_fitness[index] = _avg_fitness;
}

void PopulationStat::setClusters(const std::size_t index, const std::size_t min_clust, const std::size_t max_clust, const std::size_t avg_clust)
{
	assert(index < max_num_gen);
	min_num_clust[index] = min_clust;
	max_num_clust[index] = max_clust;
	avg_num_clust[index] = avg_clust;
}

double PopulationStat::getMaxFitness(const std::size_t index) const
{
	assert(index < max_num_gen);
	return max_fitness[index];
}

double PopulationStat::getAvgFitness(const std::size_t index) const
{
	assert(index < max_num_gen);
	return avg_fitness[index];
}

std::size_t PopulationStat::getMinNumClust(const std::size_t index) const
{
	assert(index < max_num_gen);
	return min_num_clust[index];
}

std::size_t PopulationStat::getMaxNumClust(const std::size_t index) const
{
	assert(index < max_num_gen);
	return max_num_clust[index];
}

std::size_t PopulationStat::getAvgNumClust(const std::size_t index) const
{
	assert(index < max_num_gen);
	return avg_num_clust[index];
}

std::size_t PopulationStat::getMaxNumGen() const
{
	return max_num_gen;
}

void PopulationStat::initialize(const std::size_t num_gen)
{
	max_fitness.resize(num_gen);
	avg_fitness.resize(num_gen);
	min_num_clust.resize(num_gen);
	max_num_clust.resize(num_gen);
	avg_num_clust.resize(num_gen);
	max_num_gen = num_gen;
}

bool PopulationStat::recordStats(const std::string & stat_filename, const std::size_t start_index, const std::size_t end_index, const bool print_header) const
{
	assert(start_index < max_num_gen);
	assert(end_index < max_num_gen);
	// Declare local variables
	FILE *output;

	// Try to open file
	if((output = fopen(stat_filename.c_str(), "a+")) == NULL)
	{
		printf("ERROR in PopulationStat::RecordStats - Could not open output file(%s).\n", stat_filename.c_str());
		return false;
	}

	// Check if header row should be printed
	if (print_header)
		fprintf(output, "Generation,Maximum Fitness,Average Fitness,Minimum Number of Clusters,Maximum Number of Clusters,Average Number of Clusters\n");
	
	// Loop through all entries in the arrays
	for (std::size_t i = start_index; i <= end_index; ++i)
		fprintf(output, "%lu,%0.17f,%f,%lu,%lu,%lu\n", i, max_fitness[i], avg_fitness[i], min_num_clust[i], max_num_clust[i], avg_num_clust[i]);
	
	fclose(output);

	return true;
}