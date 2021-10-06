#ifndef ISLAND_STAT_H
#define ISLAND_STAT_H

#include "Chromosome.h"
#include <ctime>
class Chromosome;

class IslandStat
{
public:
	IslandStat();
	~IslandStat();

	void setMinNumClust(const std::size_t);
	void setMaxNumClust(const std::size_t);
	void setInitTime(const float);
	void setGA_Time(const float);
	void setFitness(const float);
	void setInitNumClust(const std::size_t);
	void setFinalNumClust(const std::size_t);
	void setNumSplits(const std::size_t);
	void setNumMerges(const std::size_t);
	void setNumRedis(const std::size_t);
	void setNumNodeCross(const std::size_t);
	void setNumClustCross(const std::size_t);
	void setNumPointMut(const std::size_t);
	void setFirstGen(const std::size_t);
	void setMaxIslands(const std::size_t);
	void setGenSecBestChrom(const std::size_t);

	std::size_t getMinNumClust() const;
	std::size_t getMaxNumClust() const;
	float getInitTime() const;
	float getGA_Time() const;
	float getFitness() const;
	std::size_t getInitNumClust() const;
	std::size_t getFinalNumClust() const;
	std::size_t getNumSplits() const;
	std::size_t getNumMerges() const;
	std::size_t getNumRedis() const;
	std::size_t getNumNodeCross() const;
	std::size_t getNumClustCross() const;
	std::size_t getNumPointMut() const;
	std::size_t getFirstGen() const;
	std::size_t getMaxIslands() const;
	std::size_t getGenSecBestChrom() const;

	void recordStats(const std::string &, const std::size_t, const bool) const;
	void updateMinMax(const std::size_t, const std::size_t);
	void updateTimes(const clock_t, const clock_t, const clock_t);
	void updateFromChromosome(Chromosome *the_chrom);

private:
	std::size_t min_num_clust;
	std::size_t max_num_clust;
	float init_time;
	float ga_time;
	float fitness;
	std::size_t num_init_clust;
	std::size_t num_final_clust;
	std::size_t num_split;
	std::size_t num_merge;
	std::size_t num_redis;
	std::size_t num_node_cross;
	std::size_t num_clust_cross;
	std::size_t num_point_mut;
	std::size_t first_gen;
	std::size_t gen_sec_best_chrom;
	std::size_t init_type;
};

#endif