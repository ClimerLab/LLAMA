#ifndef CLUSTER_H
#define CLUSTER_H

#include <cstdio>
#include <vector>

class Node;
class Cluster;
class Graph;

class Cluster
{
public:
	Cluster();
	~Cluster();

	Node * getFirstNode();
	Node * getLastNode();
	Cluster * getNextCluster();
	std::size_t getNumNodes() const;
	double getFitness() const;
	bool needFitnessRecalc() const;
	bool needClusterFixed() const;
	std::size_t getDegreeSum() const;
	std::size_t getNumIntraClust() const;

	void setFirstNode(Node *);
	void setLastNode(Node *);
	void setNextCluster(Cluster *);
	void setNumNodes(const std::size_t);
	void setFitness(const double);
	void setRecalcFitness(const bool);
	void setFixCluster(const bool);
	void setDegreeSum(const std::size_t);
	void setNumIntraClust(const std::size_t);

	void addNode(const std::size_t, const std::size_t);
	void addNode(Node *, const std::size_t);
	Node * removeNode(const std::size_t, const std::size_t);
	std::size_t getNodeID(const std::size_t);
	std::size_t getNodeIndex(const std::size_t);
	bool containsNode(const std::size_t);
	void copyNonPointers(Cluster *);
	void calculateDegreeSum(Graph &);
	void sortNodes();
	bool verifyBFS(Graph &, std::vector<std::size_t> &);

private:
	Node *first_node;
	Node *last_node;
	Cluster *next_cluster;
	std::size_t num_nodes;
	double fitness;
	bool recalc_fitness;
	bool fix_cluster;

	std::size_t degree_sum;
	std::size_t intra_clust_edges;
};

#endif