#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include "Cluster.h"
#include "ConfigParser.h"
#include <string>
#include <vector>

class Graph;
class Edge;
class SplitMem;
class Node;

class Chromosome
{
public:
	Chromosome();
	Chromosome(const ConfigParser &);
	~Chromosome();
	void deleteAllClusters();
	Cluster * getFirstCluster();

	std::size_t getClusterNum(const std::size_t) const;
	std::size_t getNumIntraClustEdges(const std::size_t) const;
	std::size_t getNumClusters() const;
	std::size_t getNumNodes() const;
	std::size_t getNumInitClusters() const;
	std::size_t getNumSplits() const;
	std::size_t getNumMerges() const;
	std::size_t getNumRedis() const;
	std::size_t getNumPointMutations() const;
	std::size_t getNumNodeCrossovers() const;
	std::size_t getNumClusterCrossovers() const;
	std::size_t getCreatedGeneration() const;
	double getFitness() const;
	int getInitType() const;
	std::size_t getLastClustIndex() const;

	void setFirstCluster(Cluster *);
	void setClusterNum(const std::size_t, const std::size_t);
	void setNumIntraClusterEdges(const std::size_t, const std::size_t);
	void setNumClusters(const std::size_t);
	void setNumInitClusters(const std::size_t);
	void setNumNodes(const std::size_t);
	void setNumSplits(const std::size_t);
	void setNumMerges(const std::size_t);
	void setNumRedis(const std::size_t);
	void setNumPointMutations(const std::size_t);
	void setNumNodeCrossovers(const std::size_t);
	void setNumClusterCrossovers(const std::size_t);
	void setGenerationCreated(const std::size_t);
	void setFitness(const double);
	void setInitType(const int);

	bool initialize(Graph &);
	bool initialize(Graph &, std::vector<std::size_t> &);
	void initializeFromClusterNum(Graph &);
	//void labelPropagation(int **labels, int *local_labels, int *label_count, int *max_labels, bool *node_chosen, int *node_ids, Graph *the_graph);

	void randomWalk(Graph &);
	void splitInit(Graph &, SplitMem &, const std::size_t, const std::size_t);
	void addCluster();
	bool deleteCluster(const std::size_t);
	std::size_t findNode(const std::size_t);
	Cluster * getCluster(const std::size_t);
	void calcNumNodes();

	void copy(Chromosome *);
	void copyNonPointers(Chromosome *);
	void reset(const bool, const bool);
	
	bool splitCluster(Graph &, SplitMem &, const std::size_t);
	bool splitClusterV2(Graph &, SplitMem &, const std::size_t);
	bool redistribute(Graph &, SplitMem &, std::size_t, std::size_t);
	bool pointMutation(Graph &, const std::size_t);
	bool canMerge(Graph &, const std::size_t, const std::size_t);
	std::size_t calcNumSharedEdges(Graph &, const std::size_t, const std::size_t);
	bool mergeClusters(const std::size_t, const std::size_t, const bool, const bool);
	bool nodeCrossOver(Graph &, Chromosome *, const std::size_t);
	bool clusterCrossOver(Graph &, Chromosome *, const std::size_t);
	bool mutate(Graph &, SplitMem &);
	bool maintainNumClusters(Graph &, SplitMem &, const std::size_t, const std::size_t);	

	bool moveNode(Graph &, const std::size_t, const std::size_t, const std::size_t);
	bool moveNode(Cluster *, Cluster *, const std::size_t, const std::size_t);
	void recordToFile(std::string &, const char *clust_title);
	int verifyChromosome(Graph &);
	void fixUnconnectedClusters(Graph &, const std::size_t);
	void calculateFitness(Graph &, const int fitness_calc);
	void sortChromosome();
	void replaceNodeNum(const std::vector<std::size_t> &);
	Node *getHeadNode(const std::size_t);
	bool deleteAllEmptyClusters();
	void updateAllClusterNum();
	bool checkClusterNumber();
	void calculateIntraEdges(Graph &);
	void calculateIntraEdgesCore(Graph &, const std::size_t);
	bool isNodeMarginal(Graph &, const std::size_t);
	bool updateLinkedList();

	bool readInCluster(std::string &, Graph &);
	bool appendCluster(const std::string &, Graph &);
	void allocateClusterNum();
	void allocateIntraEdges();

	bool recordOutput(std::string &, const int);
	void readChromosome(FILE *, Graph &, const int);
	void recordObjective(const std::string &, const std::string &) const;

private:
	Cluster *first_cluster;
	std::vector<std::size_t> cluster_num;
	std::vector<std::size_t> num_intra_edges;
	double fitness;
	int init_type;
	std::size_t num_clusters;
	std::size_t num_nodes;

	std::size_t num_init_clusters;
	std::size_t num_splits;
	std::size_t num_merges;
	std::size_t num_redis;
	std::size_t num_point_muts;
	std::size_t num_node_crossovers;
	std::size_t num_cluster_crossovers;
	std::size_t generation_created;	

	double delta_fit;
	std::size_t delta_split;
	std::size_t delta_merges;
	std::size_t delta_redis;
	std::size_t delta_point_muts;
	std::size_t delta_node_cross;
	std::size_t delta_clust_cross;

	const double MUTATION_RATE;
	const double CLUST_SPLIT_RATE;
	const double CLUST_MERGE_RATE;
	const double SPLIT_PROB;
	const double MERGE_PROB;
	const double REDISRIB_PROB;
	const double NODE_MUT_PROB;
	const double MUT_TO_SINGLETON;
};

#endif