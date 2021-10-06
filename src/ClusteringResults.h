#ifndef CLUSTERING_RESULTS_H
#define CLUSTERING_RESULTS_H

#include "Chromosome.h"
#include "ConfigParser.h"
#include "PopulationStat.h"
#include "SplitMem.h"
#include <vector>


const double FITNESS_TOL = 0.00001;
const std::size_t STEADY_STATE_DURATION = 500;

enum InitType
{
	unweighted_split = 0,
	weighted_split = 1,
	label_propagation = 2,
	random_walk = 3,
	combo = 4
};

class ClusteringResults
{
private:
	const ConfigParser parser;
	std::size_t min_num_clusters;
	std::size_t max_num_clusters;
	bool cross_type;
	bool using_migration;
	bool max_gen_reached;

	const std::size_t POP_SIZE;
	const std::size_t MAX_GEN;
	const int FIT_EQ;
	const int INIT_MODE;
	const double PERC_RND_WALK;
	const std::size_t TOUR_SIZE;
	const std::size_t NUM_FROM_TOUR;
	const std::size_t NUM_ELITE;
	const double PERC_CHILD_CROSS;
	const double PERC_NODE_CROSS;
	const double PERC_CHILD_CROSS_AND_MUT;
	const std::size_t SMALL_GRAGH_SIZE;
	const std::size_t SMALL_GEN;
	const std::size_t MED_GRAPH_SIZE;
	const std::size_t MED_GEN;
	const std::size_t MIGRATION_GEN;
	const bool START_CROSS_TYPE;
	const double CROSS_TYPE_CHANGE_PERC;
	const bool STATUS_UPDATES;

public:
	ClusteringResults(ConfigParser &, const std::size_t);
	~ClusteringResults();

	void initializePopulation(Graph &);
	void runGA_Core(Graph &, const std::string &);
	Chromosome * getSortedChromosome(const bool, const std::size_t);

	bool readPopulation(Graph &, const std::string &);
	bool recordPopulation(const std::string &);
	void readMigrant(Graph &, const std::string &);
	
	void setMinNumClusters(const std::size_t min_num_clust);
	void setMaxNumClusters(const std::size_t max_num_clust);
	void setCrossType(const bool cross_type);
	void setUsingMigration(const bool using_migration);
	void setMaxGenReached(const bool max_gen_reached);

	std::size_t getMinNumClusters() const;
	std::size_t getMaxNumCusters() const;	
	bool getCrossType() const;
	bool getUsingMigration() const;
	bool getMaxGenReached() const;

	std::size_t getSmallGraphSize() const;

	std::vector<std::vector<Chromosome>> pop;// [2] [POP_SIZE]
	std::vector<double> fitness;
	std::vector<std::size_t> num_clusters;
	std::vector<std::size_t> chrom_index;
	bool parent_gen;
	std::size_t generation;
	
	PopulationStat my_pop_stats;
	SplitMem mem;


};
#endif