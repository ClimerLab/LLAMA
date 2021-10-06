#include "IslandStat.h"
#include "Chromosome.h"
#include <string>
#include <cstdio>

IslandStat::IslandStat() :	min_num_clust(0),
							max_num_clust(0),
							init_time(0.0),
							ga_time(0.0),
							fitness(0.0),
							num_init_clust(0),
							num_final_clust(0),
							num_split(0),
							num_merge(0),
							num_redis(0),
							num_node_cross(0),
							num_clust_cross(0),
							num_point_mut(0),
							first_gen(0),
							gen_sec_best_chrom(0),
							init_type(0)
{	
}

IslandStat::~IslandStat()
{
}

void IslandStat::recordStats(const std::string &stat_filename, const std::size_t island_num, const bool print_header) const
{
	// Declare local variables
	FILE *output;

	// Try to open file
	if((output = fopen(stat_filename.c_str(), "a+")) == NULL)
	{
		printf("ERROR in IslandStat::RecordStats - Could not open output file(%s).\n", stat_filename.c_str());
		return;
	}

	// Check if header row should be printed
	if (print_header)
	{
		// Print header row
		fprintf(output, "Island,Min Clust,Max Clust,Init Time(s),Generation Time(s),Fitness,Init Type,");
		fprintf(output, "First Generation of Best Chrom,First Generation of 2nd Best Chrom,Num Init Clust,Num Final Clust,Num Splits,Num Merges,Num Redistributions,");
		fprintf(output, "Num Node Crossovers,Num Clust Crossovers,Num Point Mutations\n");
	}
	
	// Print statistics
	fprintf(output, "%lu,%lu,%lu,%f,%f,%f,%lu,", island_num, min_num_clust, max_num_clust, init_time, ga_time, fitness,init_type);
	fprintf(output, "%lu,%lu,%lu,%lu,%lu,%lu,%lu,", first_gen, gen_sec_best_chrom, num_init_clust, num_final_clust, num_split, num_merge, num_redis);
	fprintf(output, "%lu,%lu,%lu\n", num_node_cross, num_clust_cross, num_point_mut);

	// Close file
	fclose(output);
}

void IslandStat::updateMinMax(const std::size_t min, const std::size_t max)
{
	min_num_clust = min;
	max_num_clust = max;
}

void IslandStat::updateTimes(const clock_t ga_start, const clock_t ga_mid, const clock_t ga_end)
{
	setInitTime((float)(ga_mid - ga_start) / CLOCKS_PER_SEC);
	setGA_Time((float)(ga_end - ga_mid) / CLOCKS_PER_SEC);
}

void IslandStat::updateFromChromosome(Chromosome *the_chrom)
{
	fitness = the_chrom->getFitness();
	num_init_clust = the_chrom->getNumInitClusters();
	num_final_clust = the_chrom->getNumClusters();
	num_split = the_chrom->getNumSplits();
	num_merge = the_chrom->getNumMerges();
	num_redis = the_chrom->getNumRedis();
	num_node_cross = the_chrom->getNumNodeCrossovers();
	num_clust_cross = the_chrom->getNumClusterCrossovers();
	num_point_mut = the_chrom->getNumPointMutations();
	first_gen = the_chrom->getCreatedGeneration();
	init_type = the_chrom->getInitType();
}

void IslandStat::setMinNumClust(const std::size_t min_clust)
{
	min_num_clust = min_clust;
}

void IslandStat::setMaxNumClust(const std::size_t max_clust)
{
	max_num_clust = max_clust;
}

void IslandStat::setInitTime(const float time)
{
	init_time = time;
}

void IslandStat::setGA_Time(const float time)
{
	ga_time = time;
}

void IslandStat::setFitness(const float fit)
{
	fitness = fit;
}

void IslandStat::setInitNumClust(const std::size_t int_num_clust)
{
	num_init_clust = int_num_clust;
}

void IslandStat::setFinalNumClust(const std::size_t final_num_clust)
{
	num_final_clust = final_num_clust;
}

void IslandStat::setNumSplits(const std::size_t _num_splits)
{
	num_split = _num_splits;
}

void IslandStat::setNumMerges(const std::size_t _num_merges)
{
	num_merge = _num_merges;
}

void IslandStat::setNumRedis(const std::size_t _num_redis)
{
	num_redis = _num_redis;
}

void IslandStat::setNumNodeCross(const std::size_t _num_node_cross)
{
	num_node_cross = _num_node_cross;
}

void IslandStat::setNumClustCross(const std::size_t _num_clust_cross)
{
	num_clust_cross = _num_clust_cross;
}

void IslandStat::setNumPointMut(const std::size_t _num_point_mut)
{
	num_point_mut = _num_point_mut;
}

void IslandStat::setFirstGen(const std::size_t _first_gen)
{
	first_gen = _first_gen;
}

void IslandStat::setGenSecBestChrom(const std::size_t gen)
{
	gen_sec_best_chrom = gen;
}

std::size_t IslandStat::getMinNumClust() const
{
	return min_num_clust;
}

std::size_t IslandStat::getMaxNumClust() const
{
	return max_num_clust;
}

float IslandStat::getInitTime() const
{
	return init_time;
}

float IslandStat::getGA_Time() const
{
	return ga_time;
}

float IslandStat::getFitness() const
{
	return fitness;
}

std::size_t IslandStat::getInitNumClust() const
{
	return num_init_clust;
}

std::size_t IslandStat::getFinalNumClust() const
{
	return num_final_clust;
}

std::size_t IslandStat::getNumSplits() const
{
	return num_split;
}

std::size_t IslandStat::getNumMerges() const
{
	return num_merge;
}

std::size_t IslandStat::getNumRedis() const
{
	return num_redis;
}

std::size_t IslandStat::getNumNodeCross() const
{
	return num_node_cross;
}

std::size_t IslandStat::getNumClustCross() const
{
	return num_clust_cross;
}

std::size_t IslandStat::getNumPointMut() const
{
	return num_point_mut;
}

std::size_t IslandStat::getFirstGen() const
{
	return first_gen;
}

std::size_t IslandStat::getGenSecBestChrom() const
{
	return gen_sec_best_chrom;
}
