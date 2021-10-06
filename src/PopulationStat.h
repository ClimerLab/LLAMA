#ifndef POPILATION_STAT_H
#define POPILATION_STAT_H

#include <string>
#include <vector>

class PopulationStat
{
public:
	PopulationStat();
	PopulationStat(const std::size_t);
	~PopulationStat();

	void setMaxFitness(const std::size_t, const double);
	void setAvgFitness(const std::size_t, const double);
	void setMinNumClust(const std::size_t, const std::size_t);
	void setMaxNumClust(const std::size_t, const std::size_t);
	void setAvgNumClust(const std::size_t, const std::size_t);
	void setMaxNumGen(const std::size_t);
	void setFitness(const std::size_t, const double, const double);
	void setClusters(const std::size_t, const std::size_t, const std::size_t, const std::size_t);

	double getMaxFitness(const std::size_t) const;
	double getAvgFitness(const std::size_t) const;
	std::size_t getMinNumClust(const std::size_t) const;
	std::size_t getMaxNumClust(const std::size_t) const;
	std::size_t getAvgNumClust(const std::size_t) const;
	std::size_t getMaxNumGen() const;

	void initialize(const std::size_t);
	bool recordStats(const std::string &, const std::size_t, const std::size_t, const bool) const;	

private:
	std::vector<double> max_fitness;
	std::vector<double> avg_fitness;
	std::vector<std::size_t> min_num_clust;
	std::vector<std::size_t> max_num_clust;
	std::vector<std::size_t> avg_num_clust;
	std::size_t max_num_gen;
};

#endif