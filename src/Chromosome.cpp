#include "Chromosome.h"
#include "Node.h"
#include "Vertex.h"
#include "Edge.h"
#include "Graph.h"
#include "ConfigParser.h"
#include "UtilityFunction.h"
#include "SplitMem.h"
#include <cstdio>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <stddef.h>
#include <assert.h>

Chromosome::Chromosome() :	num_clusters(0),
							num_nodes(0),
							num_init_clusters(0),
							num_splits(0),
							num_merges(0),
							num_redis(0),
							num_point_muts(0),
							num_node_crossovers(0),
							num_cluster_crossovers(0),
							generation_created(0),
							fitness(0.0),
							init_type(0),
							delta_fit(0.0),
							delta_split(0),
							delta_merges(0),
							delta_redis(0),
							delta_point_muts(0),
							delta_node_cross(0),
							delta_clust_cross(0),
							MUTATION_RATE(0.0),
							CLUST_SPLIT_RATE(0.0),
							CLUST_MERGE_RATE(0.0),
							SPLIT_PROB(0.0),
							MERGE_PROB(0.0),
							REDISRIB_PROB(0.0),
							NODE_MUT_PROB(0.0),
							MUT_TO_SINGLETON(0.0)
{
	first_cluster = NULL;	
}

Chromosome::Chromosome(const ConfigParser &parser) :	num_clusters(0),
														num_nodes(0),
														num_init_clusters(0),
														num_splits(0),
														num_merges(0),
														num_redis(0),
														num_point_muts(0),
														num_node_crossovers(0),
														num_cluster_crossovers(0),
														generation_created(0),
														fitness(0.0),
														init_type(0),
														delta_fit(0.0),
														delta_split(0),
														delta_merges(0),
														delta_redis(0),
														delta_point_muts(0),
														delta_node_cross(0),
														delta_clust_cross(0),
														MUTATION_RATE(parser.getDouble("MUTATION_RATE")),
														CLUST_SPLIT_RATE(parser.getDouble("CLUST_SPLIT_RATE")),
														CLUST_MERGE_RATE(parser.getDouble("CLUST_MERGE_RATE")),
														SPLIT_PROB(parser.getDouble("SPLIT_PROB")),
														MERGE_PROB(parser.getDouble("MERGE_PROB")),
														REDISRIB_PROB(parser.getDouble("REDISRIB_PROB")),
														NODE_MUT_PROB(parser.getDouble("NODE_MUT_PROB")),
														MUT_TO_SINGLETON(parser.getDouble("MUT_TO_SINGLETON"))
{
	first_cluster = NULL;	
}

Chromosome::~Chromosome()
{
	Cluster* tmp = NULL;

	while (first_cluster)
	{
		tmp = first_cluster;
		first_cluster = first_cluster->getNextCluster();
		tmp->setNextCluster(NULL);

		delete tmp;
	}	
}

void Chromosome::deleteAllClusters()
{
	Cluster* tmp = NULL;

	while (first_cluster)
	{
		tmp = first_cluster;
		first_cluster = first_cluster->getNextCluster();
		tmp->setNextCluster(NULL);

		delete tmp;
	}
}
Cluster * Chromosome::getFirstCluster()
{
	return first_cluster;
}

std::size_t Chromosome::getClusterNum(const std::size_t node_id) const
{
	assert(node_id < num_nodes);

	return cluster_num[node_id];
}

std::size_t Chromosome::getNumIntraClustEdges(const std::size_t node_id) const
{
	assert(node_id < num_nodes);

	return num_intra_edges[node_id];
}

std::size_t Chromosome::getNumClusters() const
{
	return num_clusters;
}

std::size_t Chromosome::getNumNodes() const
{
	return num_nodes;
}

std::size_t Chromosome::getNumInitClusters() const
{
	return num_init_clusters;
}

std::size_t Chromosome::getNumSplits() const
{
	return num_splits;
}

std::size_t Chromosome::getNumMerges() const
{
	return num_merges;
}

std::size_t Chromosome::getNumRedis() const
{
	return num_redis;
}

std::size_t Chromosome::getNumPointMutations() const
{
	return num_point_muts;
}

std::size_t Chromosome::getNumNodeCrossovers() const
{
	return num_node_crossovers;
}

std::size_t Chromosome::getNumClusterCrossovers() const
{
	return num_cluster_crossovers;
}

std::size_t Chromosome::getCreatedGeneration() const
{
	return generation_created;
}

double Chromosome::getFitness() const
{
	return fitness;
}

int Chromosome::getInitType() const
{
	return init_type;
}

std::size_t Chromosome::getLastClustIndex() const
{
	if (num_clusters == 0)
		return 0;
	else
		return num_clusters - 1;
}

void Chromosome::setFirstCluster(Cluster * clust_to_set)
{
	first_cluster = clust_to_set;
}

void Chromosome::setClusterNum(const std::size_t node_id, const std::size_t clust_num)
{
	assert(node_id < num_nodes);
	assert(clust_num < num_clusters);
	
	cluster_num[node_id] = clust_num;
}

void Chromosome::setNumIntraClusterEdges(const std::size_t node_id, const std::size_t num_intra)
{
	assert(node_id < num_nodes);

	num_intra_edges[node_id] = num_intra;
}

void Chromosome::setNumClusters(const std::size_t num_clust)
{
	num_clusters = num_clust;
}

void Chromosome::setNumInitClusters(const std::size_t num_init_clust)
{
	num_init_clusters = num_init_clust;
}

void Chromosome::setNumNodes(const std::size_t _num_nodes)
{
	num_nodes = _num_nodes;
}

void Chromosome::setNumSplits(const std::size_t _num_splits)
{
	num_splits = _num_splits;
}

void Chromosome::setNumMerges(const std::size_t _num_merge)
{
	num_merges = _num_merge;
}

void Chromosome::setNumRedis(const std::size_t _num_redis)
{
	num_redis = _num_redis;
}

void Chromosome::setNumPointMutations(const std::size_t _num_mut)
{
	num_point_muts = _num_mut;
}

void Chromosome::setNumNodeCrossovers(const std::size_t _num_cross)
{
	num_node_crossovers = _num_cross;
}

void Chromosome::setNumClusterCrossovers(const std::size_t _num_cross)
{
	num_cluster_crossovers = _num_cross;
}

void Chromosome::setGenerationCreated(const std::size_t gen)
{
	generation_created = gen;
}

void Chromosome::setFitness(const double _fitness)
{
	fitness = _fitness;
}

void Chromosome::setInitType(const int _init_type)
{
	init_type = _init_type;
}

bool Chromosome::initialize(Graph &g)
{
	// Check that the chromosome is empty
	if (this->first_cluster != NULL)
	{
		printf("ERROR in Chromosome::Initialize - Cannot initialize non-empty chromosome.\n");
		return false;
	}

	// Check that the number of nodes to add is valid
	if (g.getNumNodes() == 0)
	{
		printf("ERROR in Chromosome::Initialize - Cannot initialize chromosome with zero nodes.\n");
		return false;
	}

	// Add first cluster
	this->addCluster();

	// Allocate memory for vectors
	cluster_num.resize(g.getNumNodes(), 0);
	num_intra_edges.resize(g.getNumNodes(), 0);	

	// Loop through all nodes in graph
	for(std::size_t i = 0; i < g.getNumNodes(); ++i)
	{
		this->first_cluster->addNode(i, g.getDegree(i));
		this->setClusterNum(i, 0);
	}

	this->calculateIntraEdges(g);

	// Set number of nodes
	this->setNumNodes(g.getNumNodes());

	return true;
}

bool Chromosome::initialize(Graph &g, std::vector<std::size_t> &_cluster_num)
{
	assert(cluster_num.size() == g.getNumNodes());

	Cluster* cur_clust = NULL;

	// Check that the chromosome is empty
	if (this->first_cluster != NULL)
	{
		printf("ERROR in Chromosome::Initialize - Cannot initialize non-empty chromosome.\n");
		return false;
	}
	// Add first cluster
	this->addCluster();

	// Allocate memory for vectors
	cluster_num.resize(g.getNumNodes(), 0);
	num_intra_edges.resize(g.getNumNodes(), 0);	

	// Loop through all nodes in graph
	for (std::size_t i = 0; i < g.getNumNodes(); ++i)
	{
		// Set cluster number
		this->setClusterNum(i, _cluster_num[i]);

		// Check if more clusters need to be added
		while (this->getLastClustIndex() < cluster_num[i])
			this->addCluster();		

		// Get a pointer to the correct cluster
		cur_clust = this->getCluster(cluster_num[i]);

		cur_clust->addNode(i, g.getDegree(i));
	}

	this->calculateIntraEdges(g);

	// Set number of nodes
	this->setNumNodes(g.getNumNodes());

	return true;
}

void Chromosome::initializeFromClusterNum(Graph &g)
{
	assert(cluster_num.size() == g.getNumNodes());

	// Move all nodes to first cluster
	if (num_clusters == 0)
		addCluster();
	else if (num_clusters > 1)
		reset(true, false);


}

/*
void Chromosome::labelPropagation(int **labels, int *local_labels, int *label_count, int *max_labels, bool *node_chosen, int *node_ids, Graph * the_graph)
{
	int cur_index, node_id, end_index, num_diff_labels_found, tmp_label,i_label, most_label_count, max_label_count, new_index;
	bool parent_gen = 0, label_changed = true, label_found;
	Edge *cur_edge;

	// Initialize the labels
	for (int i = 0; i < the_graph->getNumNodes(); i++)
	{
		labels[parent_gen][i] = i;
		node_ids[i] = i;
	}

	// TEST ONLY
	std::string file_name = "C:\\Users\\Kenny\\Documents\\School\\Research\\GA_Clustering_Results\\karate\\labels.csv";
	std::string file_name2 = "C:\\Users\\Kenny\\Documents\\School\\Research\\GA_Clustering_Results\\karate\\node_ids.csv";
	FILE *output;

	// TEST ONLY
	output = fopen(file_name.c_str(), "a+");
	if (output == NULL)
	{
		printf("ERROR in Chromosome::RandomWalk - Could not open file.\n");
		return;
	}
	for (i_label = 0; i_label < the_graph->getNumNodes(); i_label++)
	{
		if (i_label == the_graph->getNumNodes() - 1)
		{
			fprintf(output, "%i\n", labels[parent_gen][i_label]);
		}
		else
		{
			fprintf(output, "%i,", labels[parent_gen][i_label]);
		}
	}
	fclose(output);


	// Loop until no labels change
	while (label_changed == true)
	{
		// Reset flag
		label_changed = false;

		// Reset end index
		end_index = this->GetNumNodes() - 1;

		// Reset which nodes have been chosen this iteration
		//for (int i = 0; i < the_graph->GetNumNodes(); i++)
		//{
//			node_chosen[i] = false;
		//}

		// Loop through the node_ids array until all nodes have been updated
		while (end_index >= 0)
		{
			// Pick a random number
			cur_index = PickNumber(0, end_index);

			// Get the node ID at the chosen index
			node_id = node_ids[cur_index];

			// Switch the node ID at the cur_index with the end_index
			node_ids[cur_index] = node_ids[end_index];
			node_ids[end_index] = node_id;
			// Decrement the end index
			end_index--;

			// TEST ONLY
			/*
			output = fopen(file_name2.c_str(), "a+");
			if (output == NULL)
			{
				printf("ERROR in Chromosome::RandomWalk - Could not open file.\n");
				return;
			}
			for (i_label = 0; i_label < the_graph->GetNumNodes(); i_label++)
			{
				if (i_label == the_graph->GetNumNodes() - 1)
				{
					fprintf(output, "%i\n", node_ids[i_label]);
				}
				else
				{
					fprintf(output, "%i,", node_ids[i_label]);
				}
			}
			fclose(output);
			

			// Get a pointer to the first edge
			cur_edge = the_graph->GetFirstEdge(node_id);

			// Reset number of different labels found
			num_diff_labels_found = 0;

			// Loop through all edges in list
			while (cur_edge != NULL)
			{
				// Reset flag
				label_found = false;

				// Get the label at the index identified by Target ID from the parent generation
				tmp_label = labels[parent_gen][cur_edge->GetTargetId()];

				// Compare 'tmp_label' to all of the different labels found for this node this iteration
				for (i_label = 0; i_label < num_diff_labels_found; i_label++)
				{
					// Check if the label at index 'i' matches tmp_label
					if (local_labels[i_label] == tmp_label)
					{
						label_found = true;
					}
				}

				// Check if a matching label was found
				if (label_found == true)
				{
					// Update the number of labels found with that particular label
					label_count[i_label]++;
				}
				else
				{
					local_labels[num_diff_labels_found] = tmp_label; // Save the value of the new found label
					label_count[num_diff_labels_found] = 1; // Set the count for this label to 1
					num_diff_labels_found++; // Increment the number of different labels found
				}

				// Go to next edge
				cur_edge = cur_edge->GetNextEdge();
			}

			//// Find the label with the highest count ////
			// Reset count
			most_label_count = 0;

			// Loop through all the differnt labels found and find the highest count of any label
			for (i_label = 0; i_label < num_diff_labels_found; i_label++)
			{
				// Check if the count for this lable is above max
				if (label_count[i_label] > most_label_count)
				{
					most_label_count = label_count[i_label];
				}
			}

			// Reset counter
			max_label_count = 0;
			// Loop through all of the different labels found and get the label(s) with the highest count
			for (i_label = 0; i_label < num_diff_labels_found; i_label++)
			{
				// Check if the label count at 'i_label' equals the max count found
				if (label_count[i_label] == most_label_count)
				{
					// Add the label to the list of labels with the max label count
					max_labels[max_label_count] = local_labels[i_label];
					max_label_count++; // Increment the number of labels found with the max label count
				}
			}

			// Pick at random one of the labels with the highest cound
			new_index = PickNumber(0, max_label_count - 1);
			
			// Check if the label has changed
			if (labels[parent_gen][node_id] != max_labels[new_index])
			{
				// Updae the lavel
				labels[parent_gen][node_id] = max_labels[new_index]; 
				label_changed = true; // Set the flag
			}
			
		}


		// TEST ONLY
		output = fopen(file_name.c_str(), "a+");
		if (output == NULL)
		{
			printf("ERROR in Chromosome::RandomWalk - Could not open file.\n");
			return;
		}
		for (i_label = 0; i_label < the_graph->GetNumNodes(); i_label++)
		{
			if (i_label == the_graph->GetNumNodes() - 1)
			{
				fprintf(output, "%i\n", labels[parent_gen][i_label]);
			}
			else
			{
				fprintf(output, "%i,", labels[parent_gen][i_label]);
			}
		}
		fclose(output);

		// Loop through all of the labels
		/*
		for (i_label = 0; i_label < the_graph->GetNumNodes(); i_label++)
		{
			// Check if the labels changes between generations at index 'i_label'
			if (labels[parent_gen][i_label] != labels[!parent_gen][i_label])
			{
				label_changed = true; // Set label changes flag
				break; // Break for loop (we don't need to check other labels)
			}
		}
		

		// Switch which row is the parent generation
		//parent_gen = !parent_gen;
	}

	int test = 4;


	
}
*/

void Chromosome::randomWalk(Graph &g)
{
	Vertex *cur_vert;
	std::size_t index, node_index, neigh_index;

	//this->checkClusterNumber();

	// Loop through all of the nodes in the graph
	for (std::size_t i_node = 0; i_node < g.getNumNodes(); ++i_node)
	{
		// Get pointer to current vertex
		cur_vert = g.getPtrToVertex(i_node);
		
		// Pick the index of a neighbor
		index = PickSizeT(0, cur_vert->getDegree()-1);		

		// Find the cluster that contains the both node
		node_index = getClusterNum(i_node);
		neigh_index = getClusterNum(cur_vert->getTargetId(index));

		// Check if both nodes are in the initial cluster
		if ((node_index == 0) && (neigh_index == 0))
		{
			this->addCluster(); // Create a new cluster

			// Move both nodes to the  new cluster
			this->moveNode(g, node_index, getLastClustIndex(), i_node);
			this->moveNode(g, neigh_index, getLastClustIndex(), cur_vert->getTargetId(index));
		}
		// Check if current node is in the intiali node
		else if (node_index == 0)
		{
			// Move this node to the same cluster as its neighbor
			this->moveNode(g, node_index, neigh_index, i_node);
		}
		// Check if the neighbor is in the intial cluster
		else if(neigh_index == 0)
		{
			// Move this node to the same cluster as its neighbor
			this->moveNode(g, neigh_index, node_index, cur_vert->getTargetId(index));
		}
		// Check that both nodes are in different clusters
		else if(node_index != neigh_index) // Both nodes are in a cluster
		{
			// Mege the two cluster
			this->mergeClusters(node_index, neigh_index, false, true);
			this->calculateIntraEdgesCore(g, node_index); // Recalculate the number of intra cluster edges
			this->getCluster(neigh_index)->setNumIntraClust(0);
		}
	}

	this->setNumMerges(0); // Reset number of merges
	this->setInitType(3);
	this->deleteAllEmptyClusters();	// Delete all empty clusters
	this->updateAllClusterNum();
}

void Chromosome::splitInit(Graph &g, SplitMem &mem, const std::size_t min_num_clust, const std::size_t max_num_clust)
{
	std::size_t num_clusters, split_index;

	// Pick number of clusters between min and max
	num_clusters = PickSizeT(min_num_clust, max_num_clust);

	// Pick a random cluster and split it until the desired number of clustes is reached
	for (std::size_t i_split = 0; i_split < num_clusters; ++i_split)
	{		
		split_index = PickSizeT(0, getLastClustIndex());
		
		if (getCluster(split_index)->getNumNodes() > 1)
			splitCluster(g, mem, split_index);		
	}

	// Set the initial number of cluster for chromosome and reset the number of splits
	setNumInitClusters(getNumSplits());
	setNumSplits(0);

	setNumClusterCrossovers(0);
}

void Chromosome::addCluster()
{
	// Allocate memory for new cluster
	Cluster * new_clust = new Cluster();

	// Check if chromosome is emty
	if (first_cluster == NULL)
		first_cluster = new_clust;
	else
	{
		// Get pointer to last cluster
		Cluster *last_clust = getCluster(getLastClustIndex());

		// Add new cluster to end of linked list
		last_clust->setNextCluster(new_clust);
	}

	// Increment number of clusters
	setNumClusters(getNumClusters() + 1);
}

bool Chromosome::deleteCluster(const std::size_t clust_to_delete_index)
{
	assert(clust_to_delete_index < num_clusters);

	// Get pointer to first cluster in chromosome and delcate null pointer to cluster
	Cluster *cur_clust = getFirstCluster();
	Cluster *pre_cust = NULL;

	// Loop until index reached
	for (std::size_t i_clust = 0; i_clust < clust_to_delete_index; ++i_clust)
	{
		pre_cust = cur_clust; // Set current cluster as previous cluster
		cur_clust = cur_clust->getNextCluster(); // Go to next cluster
	}

	// Check if current cluster is empty
	if (cur_clust->getNumNodes() > 0)
	{
		printf("ERROR in Chromosome::DeleteCluster - Trying to delete cluster with nodes.\n");
		return false;
	}

	// Check if index was found head of list
	if (pre_cust == NULL)
	{
		// Point head of list to next cluster
		this->setFirstCluster(cur_clust->getNextCluster());
	}
	else
	{
		// Point previous cluster to next cluster
		pre_cust->setNextCluster(cur_clust->getNextCluster());
	}

	// Clear up current cluster and delete
	cur_clust->setNextCluster(NULL);
	delete cur_clust;

	// Decrement number of clusters
	this->setNumClusters(this->getNumClusters() - 1);

	return true;
}

std::size_t Chromosome::findNode(const std::size_t node_id)
{
	assert(node_id < num_nodes);
	std::size_t index = num_clusters;

	// Get pointer to first cluster
	Cluster *cur_clust = getFirstCluster();

	// Loop through all cluster
	for (std::size_t i_clust = 0; i_clust < getNumClusters(); ++i_clust)
	{
		// Check if current cluster contains desired node
		if(cur_clust->containsNode(node_id))
		{
			index = i_clust;
			break;
		}			
		
		cur_clust = cur_clust->getNextCluster(); // Go to next cluster				
	}

	if(index == num_clusters)
	{
		printf("ERROR in Chromosome::findNode - Could not find node %lu in chromosome\n", node_id);
		exit(1);
	}

	return index;
}

Cluster * Chromosome::getCluster(const std::size_t index)
{
	assert(index < num_clusters);
	
	// Get pointer to first cluster
	Cluster *cur_cluster = getFirstCluster();

	// Loop until the desired indes is reached (The first cluster is at index 0)
	for (std::size_t i_clust = 0; i_clust < index; ++i_clust)
		cur_cluster = cur_cluster->getNextCluster();

	return cur_cluster;
}

void Chromosome::calcNumNodes()
{
	// Initialize number of node in chromosome
	num_nodes = 0;

	// Get pointer to first cluster
	Cluster *cur_clust = getFirstCluster();

	// Loop through linked list
	while (cur_clust->getNextCluster() != NULL)
	{
		// Add the number of nodes in current cluster to sum
		num_nodes += cur_clust->getNumNodes();

		// Go to next cluster
		cur_clust = cur_clust->getNextCluster();
	}
}

// This function copies the calling chromosome to the 'in_chrom'
void Chromosome::copy(Chromosome * in_chrom)
{
	assert(getNumClusters() > 0);
	assert(in_chrom->getNumClusters() > 0);
	assert(getNumClusters() == in_chrom->getNumClusters());

	// Declare local variables
	Cluster *cur_clust, *clust_to_copy;
	Node *cur_node, *node_to_copy, *pre_node;
	int clust_num;

	// Move all of the nodes from the input chromosome to one cluster
	in_chrom->reset(false, false); // Do not delete cluster

	// Get pointers to the first cluster in both chromosome
	clust_to_copy = this->getFirstCluster(); // Copy the structure/values in the calling chromosome to the input chromosome
	cur_clust = in_chrom->getFirstCluster();

	// Get a pointer to the first node in the destination cluster
	cur_node = cur_clust->getFirstNode();
	pre_node = NULL;
	clust_num = 0;

	// Loop through all clusters in the calling chromosome
	while (clust_to_copy != NULL)
	{
		// Check if the cluster to copy has any nodes. If not skip
		if (clust_to_copy->getNumNodes() > 0)
		{
			// Get a pointer to the first node in the cluster to cop
			node_to_copy = clust_to_copy->getFirstNode();

			// Loop through all nodes in the cluster to copy
			while (node_to_copy != NULL)
			{
				// Set the values of the current node to the values of the node to copy
				cur_node->setID(node_to_copy->getID());
				// Update the cluster number array
				in_chrom->setClusterNum(cur_node->getID(), clust_num);
				// Update number of intra cluster edges
				in_chrom->setNumIntraClusterEdges(cur_node->getID(), this->getNumIntraClustEdges(cur_node->getID()));

				// Go to next node to copy
				node_to_copy = node_to_copy->getNextNode();

				// Go to next node in the current cluster
				pre_node = cur_node;
				cur_node = cur_node->getNextNode();
			}

			// Set the pointer of the last node in the current cluster
			pre_node->setNextNode(NULL);
			cur_clust->setLastNode(pre_node);

			// Set attributes for the current cluster
			clust_to_copy->copyNonPointers(cur_clust);

			// Check if more clusters are needed
			if (in_chrom->getNumClusters() == clust_num+1)
			{
				in_chrom->addCluster();
			}

			// Go to the next clust in the input chromosome
			cur_clust = cur_clust->getNextCluster();

			// Check that a valid cluster was reached
			if (cur_clust != NULL)
			{
				// Set the current node as the first node in the current cluster
				cur_clust->setFirstNode(cur_node);
			}

			// Increment the counter for non-empty clusters
			++clust_num;
		}

		// Go the next cluster
		clust_to_copy = clust_to_copy->getNextCluster();		
	}
	
	in_chrom->deleteAllEmptyClusters();

	// Copy other chromosome attributes
	this->copyNonPointers(in_chrom);
}

void Chromosome::copyNonPointers(Chromosome * in_chrom)
{
	in_chrom->setNumClusters(this->getNumClusters());
	in_chrom->setNumInitClusters(this->getNumInitClusters());
	in_chrom->setNumSplits(this->getNumSplits());
	in_chrom->setNumMerges(this->getNumMerges());
	in_chrom->setNumRedis(this->getNumRedis());
	in_chrom->setNumPointMutations(this->getNumPointMutations());
	in_chrom->setNumNodeCrossovers(this->getNumNodeCrossovers());
	in_chrom->setNumClusterCrossovers(this->getNumClusterCrossovers());
	in_chrom->setGenerationCreated(this->getCreatedGeneration());
	in_chrom->setFitness(this->getFitness());
}

void Chromosome::reset(const bool delete_clusters, const bool reset_clust_num)
{
	if (delete_clusters)
	{
		while (getNumClusters() > 1)
		{
			mergeClusters(0, 1, delete_clusters, true);
		}

		getFirstCluster()->setNumIntraClust(0);
	}
	else
	{
		std::size_t orig_num_custers = getNumClusters();
		
		for(std::size_t i = 1; i < orig_num_custers; ++i)
			mergeClusters(0, i, delete_clusters, false); // Merge the second cluster into the first cluster
		
		this->getFirstCluster()->setNumIntraClust(0);
	}

	// Loop through all nodes
	if (reset_clust_num)
	{
		for (std::size_t i = 0; i < this->getNumNodes(); ++i)
		{
			setClusterNum(i, 0);
			setNumIntraClusterEdges(i, 0);
		}
	}
	
	setNumSplits(0);
	setNumMerges(0);
	setNumRedis(0);
	setNumPointMutations(0);
	setNumClusterCrossovers(0);
	setNumNodeCrossovers(0);
	setGenerationCreated(0);
}

/*
bool Chromosome::splitCluster(Graph &g, SplitMem &mem, const std::size_t cluster_index)
{
	assert(cluster_index < num_clusters);
	assert(getCluster(cluster_index)->getNumNodes() > 1);

	// Declare local variables
	Cluster *clust_to_split;
	Node *cur_node;
	Edge *cur_edge, *last_edge;
	std::size_t tmp_index, node_id, local_index, nodes_part;
	std::size_t index_in_part1, index_in_part2, nodes_in_part1, nodes_in_part2, target_id;
	std::size_t tmp_index1, tmp;
	bool cur_part, node_found;
	
	// Get a pointer to the cluster to split
	clust_to_split = this->getCluster(cluster_index);

	std::vector<std::size_t> part1(clust_to_split->getNumNodes());
	std::vector<std::size_t> part2(clust_to_split->getNumNodes());
	std::vector<std::size_t> nodes_to_pick(clust_to_split->getNumNodes());
	std::vector<std::size_t> index_in_clust(g.getNumNodes());
	std::vector<std::size_t> end_index(clust_to_split->getNumNodes());
	std::vector<bool> in_clust(g.getNumNodes(), false);
	std::vector<bool> in_part(g.getNumNodes(), false);	
	std::vector<bool> part_num(clust_to_split->getNumNodes());
	std::vector<bool> has_edges(clust_to_split->getNumNodes(), true);

		// Get pointer to first node in cluster
	cur_node = clust_to_split->getFirstNode();

	// Set values based on nodes in array to split
	for (std::size_t i_node = 0; i_node < clust_to_split->getNumNodes(); ++i_node)
	{
		in_clust[cur_node->getID()] = true; // Indicate which nodes are in the cluster to split
		index_in_clust[cur_node->getID()] = i_node; // Indicate the index in the cluster of each node
		nodes_to_pick[i_node] = cur_node->getID(); // Add each node to the list of available nodes
		end_index[i_node] = g.getDegree(cur_node->getID()) - 1; // Save the index of the last edge for edge node

		cur_node = cur_node->getNextNode(); // Go to the next edge
	}
	
	// Initialize indicies for queues
	nodes_in_part1 = 0;
	nodes_in_part2 = 0;

	// Pick a node from the cluster
	tmp_index1 = PickSizeT(0, clust_to_split->getNumNodes() - 1);
	// Add it to partition 1 and increment nodes_in_part1
	part1[nodes_in_part1++] = nodes_to_pick[tmp_index1];
	// Indicate it has been added to partition
	in_part[nodes_to_pick[tmp_index1]] = true;
	part_num[tmp_index1] = 1;

	// Swap the node with the last valid node in the array
	tmp = nodes_to_pick[tmp_index1]; // Save the value at 'tmp_index1' 
	// Override the value at 'tmp_index1' with the last value in the array
	nodes_to_pick[tmp_index1] = nodes_to_pick[clust_to_split->getNumNodes() - 1];
	// Override the last variable in the array with the saved temp variable
	nodes_to_pick[clust_to_split->getNumNodes() - 1] = tmp;

	// Pick a node from the cluster and add it to partition 2 (don't allow the last variable in the array to be picked)
	tmp_index = PickSizeT(0, clust_to_split->getNumNodes() - 2);
	// Add it to partition 2 and increment nodes_in_part2
	part2[nodes_in_part2++] = nodes_to_pick[tmp_index];
	// Indicate it has been added to partition
	in_part[nodes_to_pick[tmp_index]] = true;

	// Indicate to which partition the node at 'tmp_index' was added. This needs to be based on the original
	// order of the nodes.
	// Check if the same index value was picked twice
	if (tmp_index == tmp_index1)
	{
		// Since the node that was picked at 'tmp_index' had to be at the last index origianlly, update
		// the 'part_num' at the last index
		part_num[clust_to_split->getNumNodes() - 1] = 0;
	}
	else
	{
		// Update the 'part_num' at the selected index
		part_num[tmp_index] = 0;
	}

	// Swap the last node and the first node picked back
	tmp = nodes_to_pick[tmp_index1]; // Save the value at 'tmp_index1' 
	// Override the value at 'tmp_index1' with the last value in the array
	nodes_to_pick[tmp_index1] = nodes_to_pick[clust_to_split->getNumNodes() - 1];
	// Override the last variable in the array with the saved temp variable
	nodes_to_pick[clust_to_split->getNumNodes() - 1] = tmp;

	// Initialize the current queue to draw from
	cur_part = 0;
	nodes_part = 2; // Accounts for both nodes added to partitions already
	
	// Loop until both partitions are empty
	while ((nodes_in_part1 > 0) || (nodes_in_part2 > 0))
	{
		// Determine which partition to draw from
		if (cur_part == true)
		{
			// Check that queue1 is not empty
			if (nodes_in_part1 > 0)
			{
				// Pick a node from partition 1 and get the ID
				index_in_part1 = PickSizeT(0, nodes_in_part1 - 1);
				node_id = part1[index_in_part1];

				// Get the index into the local graph corresponding to the node ID
				local_index = index_in_clust[node_id];

				// Loop untill vertex has no more edges
				while (has_edges[local_index])
				{
					node_found = false; // Reset flag

					// Pick edge at random
					tmp_index = PickSizeT(0, end_index[local_index]);

					// Get target ID of selected edge
					target_id = g.getPtrToVertex(node_id)->getTargetId(tmp_index);

					// Check if target ID is in the cluster to split and is not in a partition
					if (in_clust[target_id] && !in_part[target_id])
					{
						// Set flag indicating node has been found to add to partition
						node_found = true;
						// Add node to partition and increment number of nodes in partition 1
						part1[nodes_in_part1++] = target_id;
						// Indicate it has been added to a partition
						in_part[target_id] = true;
						// Increment the number of nodes paritioned
						++nodes_part;
						// Indicate that the node has been added to partition 1
						part_num[index_in_clust[target_id]] = 1;
					}

					// Check if the edge at the last valid index was selected
					if(tmp_index != end_index[local_index])
						g.getPtrToVertex(node_id)->swapEdges(tmp_index, end_index[local_index]);
					
					if(end_index[local_index] == 0)
						has_edges[local_index] = false;
					else					
						--end_index[local_index]; // Decrement the index of the last valid edge

					// Check if node was found to add to partition
					if (node_found)
						break; // break while loop					
				}

				// Check if local vertex is empty now				
				if (has_edges[local_index] == false)
				{
					// Override the value at 'index_in_part1' with the value at the last index in the partition
					part1[index_in_part1] = part1[nodes_in_part1 - 1];
					
					// Derement the number of nodes in partition 1
					--nodes_in_part1;
				}
			}

			// Switch current partition
			cur_part = !cur_part;
		}
		else
		{
			// Check that partition2 is not empty
			if (nodes_in_part2 > 0)
			{
				// Pick a node from partition 2 and get the ID
				index_in_part2 = PickNumber(0, nodes_in_part2 - 1);
				node_id = part2[index_in_part2];

				// Get the index into the local graph corresponding to the node ID
				local_index = index_in_clust[node_id];

				// Loop untill vertex has no more edges
				while(has_edges[local_index])
				{
					node_found = false; // Reset flag

					// Pick edge at random
					tmp_index = PickSizeT(0, end_index[local_index]);

					// Get target ID of selected edge
					target_id = g.getPtrToVertex(node_id)->getTargetId(tmp_index);

					// Check if target ID is in the cluster to split and is not in a partition
					if (in_clust[target_id] && !in_part[target_id])
					{
						// Set flag indicating node has been found to add to partition
						node_found = true;
						// Add node to partition and increment number of nodes in partition 2
						part2[nodes_in_part2++] = target_id;
						// Indicate it has been added to a partition
						in_part[target_id] = true;
						// Increment the number of nodes paritioned
						++nodes_part;
						// Indicate that the node has been added to partition 2
						part_num[index_in_clust[target_id]] = 0;
					}

					// Check if the edge at the last valid index was selected
					if (tmp_index != end_index[local_index])
						g.getPtrToVertex(node_id)->swapEdges(tmp_index, end_index[local_index]);

					if(end_index[local_index] == 0)
						has_edges[local_index] = false;
					else					
						--end_index[local_index]; // Decrement the index of the last valid edge

					// Check if node was found to add to partition
					if (node_found)
						break; // break while loop					
				}

				// Check if local vertex is empty now
				if (has_edges[local_index] == false)
				{
					// Override the value at 'index_in_part1' with the value at the last index in the partition
					part2[index_in_part2] = part2[nodes_in_part2 - 1];

					// Derement the number of nodes in partition 2
					--nodes_in_part2;
				}
			}

			// Switch current partition
			cur_part = !cur_part;
		}
	}

	// Check that all nodes were partitioned
	if (nodes_part != clust_to_split->getNumNodes())
	{
		printf("Nodes in part: %lu\n", nodes_part);
		// Print error
		printf("ERROR in Chromosome::SplitCluster - Not all nodes in partition\n");
		exit(1);
		return false;
	}

	// Add new cluster
	this->addCluster();

	// Save the number of nodes in the cluster to split
	std::size_t orig_num_nodes = clust_to_split->getNumNodes();

	// Loop through all nodes in the cluster to split
	for (std::size_t j_node = 0; j_node < orig_num_nodes; ++j_node)
	{
		// Check if the j_node was in partition 2
		if(part_num[j_node] == 0)
		{
			// Move the node to the new cluster
			this->moveNode(g, cluster_index, getLastClustIndex(), nodes_to_pick[j_node]);
		}
	}

	// Get a pointer to the last cluster
	Cluster *new_clust = getCluster(getLastClustIndex());

	// Set flags to recalculate fitness
	clust_to_split->setRecalcFitness(true);
	new_clust->setRecalcFitness(true);

	// Increment the number of splits
	setNumSplits(getNumSplits() + 1);	

	// Return success
	return true;
}
*/

bool Chromosome::splitCluster(Graph &g, SplitMem &mem, const std::size_t cluster_index)
{
	assert(cluster_index < num_clusters);
	assert(getCluster(cluster_index)->getNumNodes() > 1);

	// Declare local variables
	Cluster *clust_to_split;
	Node *cur_node;
	Edge *cur_edge, *last_edge;
	std::size_t tmp_index, node_id, local_index, nodes_part;
	std::size_t index_in_part1, index_in_part2, nodes_in_part1, nodes_in_part2, target_id;
	std::size_t tmp_index1, tmp;
	bool cur_part, node_found;
		
	// Get a pointer to the cluster to split
	clust_to_split = this->getCluster(cluster_index);

	// Initialize arrays
	for (int i = 0; i < g.getNumNodes(); i++)
	{
		mem.in_clust[i] = false;
		mem.in_part[i] = false;
		mem.has_edges[i] = true;
	}

	// Get pointer to first node in cluster
	cur_node = clust_to_split->getFirstNode();

	// Set values based on nodes in array to split
	for (std::size_t i_node = 0; i_node < clust_to_split->getNumNodes(); ++i_node)
	{
		mem.in_clust[cur_node->getID()] = true; // Indicate which nodes are in the cluster to split
		mem.index_in_clust[cur_node->getID()] = i_node; // Indicate the index in the cluster of each node
		mem.nodes_to_pick[i_node] = cur_node->getID(); // Add each node to the list of available nodes
		mem.end_index[i_node] = g.getDegree(cur_node->getID()) - 1; // Save the index of the last edge for edge node

		cur_node = cur_node->getNextNode(); // Go to the next edge
	}
	
	// Initialize indicies for queues
	nodes_in_part1 = 0;
	nodes_in_part2 = 0;

	// Pick a node from the cluster
	tmp_index1 = PickSizeT(0, clust_to_split->getNumNodes() - 1);
	// Add it to partition 1 and increment nodes_in_part1
	mem.part1[nodes_in_part1++] = mem.nodes_to_pick[tmp_index1];
	// Indicate it has been added to partition
	mem.in_part[mem.nodes_to_pick[tmp_index1]] = true;
	mem.part_num[tmp_index1] = 1;

	// Swap the node with the last valid node in the array
	tmp = mem.nodes_to_pick[tmp_index1]; // Save the value at 'tmp_index1' 
	// Override the value at 'tmp_index1' with the last value in the array
	mem.nodes_to_pick[tmp_index1] = mem.nodes_to_pick[clust_to_split->getNumNodes() - 1];
	// Override the last variable in the array with the saved temp variable
	mem.nodes_to_pick[clust_to_split->getNumNodes() - 1] = tmp;

	// Pick a node from the cluster and add it to partition 2 (don't allow the last variable in the array to be picked)
	tmp_index = PickSizeT(0, clust_to_split->getNumNodes() - 2);
	// Add it to partition 2 and increment nodes_in_part2
	mem.part2[nodes_in_part2++] = mem.nodes_to_pick[tmp_index];
	// Indicate it has been added to partition
	mem.in_part[mem.nodes_to_pick[tmp_index]] = true;

	// Indicate to which partition the node at 'tmp_index' was added. This needs to be based on the original
	// order of the nodes.
	// Check if the same index value was picked twice
	if (tmp_index == tmp_index1)
	{
		// Since the node that was picked at 'tmp_index' had to be at the last index origianlly, update
		// the 'part_num' at the last index
		mem.part_num[clust_to_split->getNumNodes() - 1] = 0;
	}
	else
	{
		// Update the 'part_num' at the selected index
		mem.part_num[tmp_index] = 0;
	}

	// Swap the last node and the first node picked back
	tmp = mem.nodes_to_pick[tmp_index1]; // Save the value at 'tmp_index1' 
	// Override the value at 'tmp_index1' with the last value in the array
	mem.nodes_to_pick[tmp_index1] = mem.nodes_to_pick[clust_to_split->getNumNodes() - 1];
	// Override the last variable in the array with the saved temp variable
	mem.nodes_to_pick[clust_to_split->getNumNodes() - 1] = tmp;

	// Initialize the current queue to draw from
	cur_part = 0;
	nodes_part = 2; // Accounts for both nodes added to partitions already

	std::size_t num_try = 0;
	// Loop until both partitions are empty
	while ((nodes_in_part1 > 0) || (nodes_in_part2 > 0))
	{
		// Determine which partition to draw from
		if (cur_part == true)
		{
			// Check that queue1 is not empty
			if (nodes_in_part1 > 0)
			{
				// Pick a node from partition 1 and get the ID
				index_in_part1 = PickSizeT(0, nodes_in_part1 - 1);
				node_id = mem.part1[index_in_part1];

				// Get the index into the local graph corresponding to the node ID
				local_index = mem.index_in_clust[node_id];

				// Loop untill vertex has no more edges
				while (mem.has_edges[local_index])
				{
					node_found = false; // Reset flag

					// Pick edge at random
					tmp_index = PickSizeT(0, mem.end_index[local_index]);

					// Get target ID of selected edge
					target_id = g.getPtrToVertex(node_id)->getTargetId(tmp_index);

					// Check if target ID is in the cluster to split and is not in a partition
					if (mem.in_clust[target_id] && !mem.in_part[target_id])
					{
						// Set flag indicating node has been found to add to partition
						node_found = true;
						// Add node to partition and increment number of nodes in partition 1
						mem.part1[nodes_in_part1++] = target_id;
						// Indicate it has been added to a partition
						mem.in_part[target_id] = true;
						// Increment the number of nodes paritioned
						++nodes_part;
						// Indicate that the node has been added to partition 1
						mem.part_num[mem.index_in_clust[target_id]] = 1;
					}

					// Check if the edge at the last valid index was selected
					if(tmp_index != mem.end_index[local_index])
						g.getPtrToVertex(node_id)->swapEdges(tmp_index, mem.end_index[local_index]);
					
					if(mem.end_index[local_index] == 0)
						mem.has_edges[local_index] = false;
					else					
						--mem.end_index[local_index]; // Decrement the index of the last valid edge

					// Check if node was found to add to partition
					if (node_found)
						break; // break while loop					
				}

				// Check if local vertex is empty now				
				if (mem.has_edges[local_index] == false)
				{
					// Override the value at 'index_in_part1' with the value at the last index in the partition
					mem.part1[index_in_part1] = mem.part1[nodes_in_part1 - 1];
					
					// Derement the number of nodes in partition 1
					--nodes_in_part1;
				}
			}

			// Switch current partition
			cur_part = !cur_part;
		}
		else
		{
			// Check that partition2 is not empty
			if (nodes_in_part2 > 0)
			{
				// Pick a node from partition 2 and get the ID
				index_in_part2 = PickNumber(0, nodes_in_part2 - 1);
				node_id = mem.part2[index_in_part2];

				// Get the index into the local graph corresponding to the node ID
				local_index = mem.index_in_clust[node_id];

				// Loop untill vertex has no more edges
				while(mem.has_edges[local_index])
				{
					node_found = false; // Reset flag

					// Pick edge at random
					tmp_index = PickSizeT(0, mem.end_index[local_index]);

					// Get target ID of selected edge
					target_id = g.getPtrToVertex(node_id)->getTargetId(tmp_index);

					// Check if target ID is in the cluster to split and is not in a partition
					if (mem.in_clust[target_id] && !mem.in_part[target_id])
					{
						// Set flag indicating node has been found to add to partition
						node_found = true;
						// Add node to partition and increment number of nodes in partition 2
						mem.part2[nodes_in_part2++] = target_id;
						// Indicate it has been added to a partition
						mem.in_part[target_id] = true;
						// Increment the number of nodes paritioned
						++nodes_part;
						// Indicate that the node has been added to partition 2
						mem.part_num[mem.index_in_clust[target_id]] = 0;
					}

					// Check if the edge at the last valid index was selected
					if (tmp_index != mem.end_index[local_index])
						g.getPtrToVertex(node_id)->swapEdges(tmp_index, mem.end_index[local_index]);

					if(mem.end_index[local_index] == 0)
						mem.has_edges[local_index] = false;
					else					
						--mem.end_index[local_index]; // Decrement the index of the last valid edge

					// Check if node was found to add to partition
					if (node_found)
						break; // break while loop					
				}

				// Check if local vertex is empty now
				if (mem.has_edges[local_index] == false)
				{
					// Override the value at 'index_in_part1' with the value at the last index in the partition
					mem.part2[index_in_part2] = mem.part2[nodes_in_part2 - 1];

					// Derement the number of nodes in partition 2
					--nodes_in_part2;
				}
			}

			// Switch current partition
			cur_part = !cur_part;
		}
	}

	// Check that all nodes were partitioned
	if (nodes_part != clust_to_split->getNumNodes())
	{
		// Print error
		printf("ERROR in Chromosome::SplitCluster - Not all nodes in partition");
		return false;
	}

	// Add new cluster
	this->addCluster();

	// Save the number of nodes in the cluster to split
	std::size_t orig_num_nodes = clust_to_split->getNumNodes();

	// Loop through all nodes in the cluster to split
	for (std::size_t j_node = 0; j_node < orig_num_nodes; ++j_node)
	{
		// Check if the j_node was in partition 2
		if(mem.part_num[j_node] == 0)
		{
			// Move the node to the new cluster
			this->moveNode(g, cluster_index, getLastClustIndex(), mem.nodes_to_pick[j_node]);
		}
	}

	// Get a pointer to the last cluster
	Cluster *new_clust = getCluster(getLastClustIndex());

	// Set flags to recalculate fitness
	clust_to_split->setRecalcFitness(true);
	new_clust->setRecalcFitness(true);

	// Increment the number of splits
	setNumSplits(getNumSplits() + 1);	

	// Return success
	return true;
}



bool Chromosome::redistribute(Graph &g, SplitMem &mem, std::size_t clust_index1, std::size_t clust_index2)
{
	// Declare local varialbes
	bool result;

	// In order to make book keeping easier, make sure that index1 <= index2
	if (clust_index1 > clust_index2)
	{
		std::size_t tmp = clust_index1;
		clust_index1 = clust_index2;
		clust_index2 = tmp;
	}

	// Merge two cluster and save the result
	result = this->mergeClusters(clust_index1, clust_index2, false, true);
	calculateIntraEdgesCore(g, clust_index1);
	getCluster(clust_index2)->setNumIntraClust(0);

	// Decrement the number of merges
	--num_merges;

	// Check if the merge was successful
	if (result)
	{
		// Split the newly merges cluster
		result = this->splitCluster(g, mem, clust_index1);

		// Decrement the number of splits
		--num_splits;

		// Increment the number of redistributions
		++num_redis;
	}

	return result;
}

bool Chromosome::pointMutation(Graph &g, const std::size_t node_to_mutate)
{
	assert(node_to_mutate < num_nodes);

	std::size_t orig_index, new_index, adj_index, adj_id;
	Cluster *orig_clust;
	Vertex *tmp_vert;
	double prob;

	// Get the index of the cluster that contains node_to_mutate
	orig_index = getClusterNum(node_to_mutate);
	
	// Get a pointer to the original cluster containing node_to_mutate
	orig_clust = getCluster(orig_index);

	// Check if the node will mutate to the same cluster as an adjacent node or create a singletone
	prob = PickProbabilty();
	if (prob < MUT_TO_SINGLETON)//(PickProbabilty() < MUT_TO_SINGLETON)
	{
		// Check if the node is already a singleton
		if (orig_clust->getNumNodes() == 1)
			return true; // There is no need to move a singleton to another singleton cluster.
		else
		{
			// Create a new cluster for the singleton
			addCluster();
			// Get the index of the last cluster
			new_index = getLastClustIndex();
		}
	}
	else
	{
		// Get a pointer to the vertex with 'node_to_mutate' as the source
		tmp_vert = g.getPtrToVertex(node_to_mutate);

		// Pick an edge at random from the vertex
		adj_index = PickSizeT(0, tmp_vert->getDegree() - 1);

		// Get target ID of the edge
		adj_id = tmp_vert->getTargetId(adj_index);		

		// Find the cluster in the calling chromosome that contains the node with ID matching adj_id
		new_index = getClusterNum(adj_id);		
	}

	// Move the 'node_to_mutate' from the original cluster to the new cluster
	if (orig_index != new_index)
	{
		Cluster *new_clust = getCluster(new_index);

		moveNode(g,orig_index, new_index, node_to_mutate);

		// Check if the original cluster still has multiple nodes
		if (orig_clust->getNumNodes() > 1)
			orig_clust->setFixCluster(true); // Set flag to check cluster for unconnected parts
		
		// Increment the counter
		++num_point_muts;

		// Set recalculate flag
		orig_clust->setRecalcFitness(true);
		new_clust->setRecalcFitness(true);
	}

	return true;
}
bool Chromosome::canMerge(Graph &g, const std::size_t cluster_index1, const std::size_t cluster_index2)
{
	Vertex *cur_vert;
	Node *cur_node;

	// Check that index 1 and index 2 refer to different clusters
	if(cluster_index1 == cluster_index2)
		return false;

	Cluster* clust1 = getCluster(cluster_index1);
	Cluster* clust2 = getCluster(cluster_index2);

	// Check that both clusters have nodes
	if((clust1->getNumNodes() == 0) || (clust2->getNumNodes() == 0))
		return false;

	cur_node = clust1->getFirstNode(); // Get pointer to first node in cluster 1

	// Loop through all nodes in cluster 1
	while(cur_node != NULL)
	{
		// Get pointer to vertex in graph
		cur_vert = g.getPtrToVertex(cur_node->getID());

		// Loop through all edges
		for(std::size_t i = 0; i < cur_vert->getDegree(); ++i)
		{
			// Check if the target of the i-th edge is in cluster 2
			if(cluster_num[cur_vert->getTargetId(i)] == cluster_index2)
				return true;
		}

		cur_node = cur_node->getNextNode();
	}

	return false;
}

std::size_t Chromosome::calcNumSharedEdges(Graph &g, const std::size_t cluster_index1, const std::size_t cluster_index2)
{
	Vertex *cur_vert;
	Node *cur_node;
	std::size_t num_edges = 0;

	// Check that index 1 and index 2 refer to different clusters
	if(cluster_index1 == cluster_index2)
		return false;

	Cluster* clust1 = getCluster(cluster_index1);
	Cluster* clust2 = getCluster(cluster_index2);

	// Check that both clusters have nodes
	if((clust1->getNumNodes() == 0) || (clust2->getNumNodes() == 0))
		return false;

	cur_node = clust1->getFirstNode(); // Get pointer to first node in cluster 1

	// Loop through all nodes in cluster 1
	while(cur_node != NULL)
	{
		// Get pointer to vertex in graph
		cur_vert = g.getPtrToVertex(cur_node->getID());

		// Loop through all edges
		for(std::size_t i = 0; i < cur_vert->getDegree(); ++i)
		{
			// Check if the target of the i-th edge is in cluster 2
			if(cluster_num[cur_vert->getTargetId(i)] == cluster_index2)
				++num_edges;
		}

		cur_node = cur_node->getNextNode();
	}
	return num_edges;
}

bool Chromosome::mergeClusters(const std::size_t clust_index1, const std::size_t clust_index2, const bool delete_clust2, const bool update_clust_num)
{
	// Check if cluster indicies are equal
	if (clust_index1 == clust_index2)
	{
		printf("ERROR in Chromosome::MergeClusters - Cannot merge cluster with itselft. (cluster %lu to cluster %lu)\n", clust_index1, clust_index2);
		return false;
	}

	// Get pointer to the clusters to merge
	Cluster *clust1 = getCluster(clust_index1);
	Cluster *clust2 = getCluster(clust_index2);

	// Check that valid 'clust1' was returned
	if (clust1 == NULL)
	{
		printf("ERROR in Chromosome::MergeClusters - Could not get valid cluster at index %lu.\n", clust_index1);
		return false;
	}

	// Check that valid 'clust2' was returned
	if (clust2 == NULL)
	{
		printf("ERROR in Chromosome::MergeClusters - Could not get valid cluster at index %lu.\n", clust_index2);
		return false;
	}

	// Check if the cluster number needs to be updated
	if (update_clust_num)
	{
		// Get a pointer to the first node in the second cluster
		Node *cur_node = clust2->getFirstNode();
		// Loop through all nodes in the cluster and update the cluster number
		while (cur_node != NULL)
		{
			this->setClusterNum(cur_node->getID(), clust_index1);
			cur_node = cur_node->getNextNode();
		}
	}

	// Point the last node of the first cluster to the head node of the second cluster
	clust1->getLastNode()->setNextNode(clust2->getFirstNode());

	// Set the last node of the first cluster to the last node of the second cluster
	clust1->setLastNode(clust2->getLastNode());

	// Update node count
	clust1->setNumNodes(clust1->getNumNodes() + clust2->getNumNodes());
	clust2->setNumNodes(0);

	// Set flag to recalculate fitness
	clust1->setRecalcFitness(true);
	clust2->setRecalcFitness(true);

	// Update degree sums for each cluster
	clust1->setDegreeSum(clust1->getDegreeSum() + clust2->getDegreeSum());
	clust2->setDegreeSum(0);

	// Update node pointers in second cluster
	clust2->setFirstNode(NULL);
	clust2->setLastNode(NULL);

	// Check flag
	if (delete_clust2)
	{
		// Delete second cluster
		if (deleteCluster(clust_index2) == false)
		{
			printf("Chromosome::MergeClusters - Could not delete the cluster at index %lu.\n", clust_index2);
			return false;
		}		
	}
		
	++num_merges;
	return true;
}

bool Chromosome::nodeCrossOver(Graph &g, Chromosome * ref_chrom, const std::size_t node_to_cross)
{
	std::size_t target_id, ref_clust_index, old_clust_index, new_clust_index, index_in_ref_clust;
	Cluster *ref_clust, *old_clust;
	bool node_found = false;

	// Get the index in the calling chromosome that contains the node_to_cross
	old_clust_index = getClusterNum(node_to_cross);
	old_clust = getCluster(old_clust_index);	

	// Get the cluster in the reference chromosome that contains the node_to_cross
	ref_clust_index = ref_chrom->getClusterNum(node_to_cross);
	ref_clust = ref_chrom->getCluster(ref_clust_index);	

	// Check if the node_to_cross is a singleton in the reference cluser
	if (ref_clust->getNumNodes() == 1)
	{
		// Check that the node_to_cross in not a singleton in the calling chromosome
		if (old_clust->getNumNodes() > 1)
		{
			this->addCluster(); // Create a new cluster

			// Get the index of the new cluster
			new_clust_index = this->getLastClustIndex();
		}
		else
		{
			new_clust_index = old_clust_index; // No need to make a new cluster and delete the old one
		}
	}
	else
	{
		// *****
		// This section finds a node adjacent to 'node_to_cross' in the reference cluster
		// *****

		// Get list of nodes adjacent to 'node_to_cross'
		std::vector<std::size_t> avail_nodes = g.getPtrToVertex(node_to_cross)->getIds();
		
		// Initialize the end index in the array
		std::size_t tmp_index, end_index = avail_nodes.size() - 1;		
		while (end_index >= 0)
		{
			// Pick an index at random
			tmp_index = PickSizeT(0, end_index);

			// Pick an adjacent node to the nodet_to_cross
			target_id = avail_nodes[tmp_index];

			// Check if 'target_id' is in the cluster
			node_found = ref_clust->containsNode(target_id);

			if(node_found)
			{
				break;
			}
			else
			{
				// Override the node at 'tmp_index' with the node at 'end_index'
				avail_nodes[tmp_index] = avail_nodes[end_index];
				--end_index;
			}			
		}		

		if (!node_found)
		{
			printf("ERROR in Chromosome::NodeCrossOver - No adjacent nodes to node %lu found in reference cluster.\n", node_to_cross);
			return false;
		}

		// Get the index of the cluster, in the calling chromosome, that contains the target_id
		new_clust_index = this->getClusterNum(target_id);		
	}

	// Check if the node_to_cross will move clusters
	if (old_clust_index != new_clust_index)
	{
		Cluster *new_clust = this->getCluster(new_clust_index);

		// Move the node_to_cross to the new cluster
		this->moveNode(g, old_clust_index, new_clust_index, node_to_cross);
		
		// Check if the original cluster has multiple nodes
		if (old_clust->getNumNodes() > 1)
		{
			// Set flag to check cluster for unconnected parts
			old_clust->setFixCluster(true);
		}
	
		// Inrement number of node crossovers for chromosome
		setNumNodeCrossovers(getNumNodeCrossovers() + 1);

		// Set recalculate flag
		old_clust->setRecalcFitness(true);
		new_clust->setRecalcFitness(true);
	}

	// Return successful node cross over
	return true;
}

bool Chromosome::clusterCrossOver(Graph &g, Chromosome * ref_chrom, const std::size_t clust_to_cross_index)
{
	Cluster *new_clust, *orig_clust;
	Node *cur_node;
	bool node_moved;
	std::vector<bool> cluster_changed;
	int clust_index;
	int orig_num_clusters = num_clusters;

	// Allocate memory for array
	cluster_changed.resize(num_clusters, false);

	// Create new cluster in calling chromosome and get a pointer to it
	addCluster();
	new_clust = getCluster(getLastClustIndex());

	// Get a pointer to the first node in the cluster to cross in the reference chromosome
	cur_node = ref_chrom->getCluster(clust_to_cross_index)->getFirstNode();

	// Loop through all nodes in clust_to_cross
	while (cur_node != NULL)
	{
		// Get a pointer to the cluster in the calling chromosome that contains the current node
		clust_index = getClusterNum(cur_node->getID());
		orig_clust = getCluster(clust_index);

		// Mark cluster as changed
		cluster_changed[clust_index] = true;

		// Move node to new cluster
		node_moved = moveNode(g,clust_index, this->getLastClustIndex(), cur_node->getID());

		// Check if move failed
		if (node_moved == false)
		{
			// Print error
			printf("ERROR in Chromosome::ClusterCrossOver - Could not move node %lu.\n", cur_node->getID());
			// Return status
			return false;
		}

		// Go to next node
		cur_node = cur_node->getNextNode();
	}

	// Get a pointer to the first cluster in the calling chromosome
	Cluster *cur_clust = first_cluster;

	// Check all clusters that were changed for unconnected parts
	for (std::size_t i_clust = 0; i_clust < orig_num_clusters; ++i_clust)
	{
		// Check if cluster was changed
		if (cluster_changed[i_clust] == true)
		{
			// Check that the current cluster is not empty
			if (cur_clust->getNumNodes() > 0)
			{
				// Fix any unconnected parts in cluster i_clust
				this->fixUnconnectedClusters(g, i_clust);
			}
		}

		// Go to next cluster
		cur_clust = cur_clust->getNextCluster();
	}

	// Increment the number of crossovers
	++num_cluster_crossovers;
	return true;
}

bool Chromosome::mutate(Graph &g, SplitMem &mem)
{
	double prob;
	std::size_t node;

	// Pick the number of nodes to attempt to mutate
	int num_nodes_mut = int(getNumNodes()*MUTATION_RATE);

	// Calculate the number of cluster to attempt splits
	int num_clust_split = int(getNumNodes()*CLUST_SPLIT_RATE);

	// Calculate the number of cluster to attempt splits
	int num_clust_merge = int(getNumNodes()*CLUST_MERGE_RATE);

	// Loop through the desired number of clusters and attempt mutation
	//for (int i_clust = 0; i_clust < num_clust_split; i_clust++)
	//{
	// Check if a cluster in the chromosome should split
	prob = PickProbabilty();
	if (prob < SPLIT_PROB)//(PickProbabilty() < SPLIT_PROB)
	{
		// Pick a random cluster to split
		std::size_t clust_split_index = PickSizeT(0, this->getLastClustIndex());

		// Get a pointer to the cluster
		Cluster *clust_to_split = getCluster(clust_split_index);

		// Check if the cluster has enough nodes to split
		if (clust_to_split->getNumNodes() > 1)
		{
			splitCluster(g, mem, clust_split_index);
		}
	}
	//}

	// Loop through the desired number of clusters and attempt mutation
	//for (int i_clust = 0; i_clust < num_clust_merge; i_clust++)
	//{
		// Check if two clusters in chromosome should merge
	prob = PickProbabilty();
	if (prob < MERGE_PROB)//(PickProbabilty() < MERGE_PROB)
	{
		// Check if chromosome has enough cluster to merge
		if (getNumClusters() > 1)
		{
			std::size_t clust_merge_index1, clust_merge_index2;

			// ***************
			// Patch
			// ***************
			if (getNumClusters() == 2)
			{
				clust_merge_index1 = 0;
				clust_merge_index2 = 1;
			}
			else
			{
				// Pick first cluster to merge
				clust_merge_index1 = PickSizeT(0, getLastClustIndex());

				// Initialize second cluster to merge to the first cluster
				clust_merge_index2 = clust_merge_index1;

				// Loop until a second cluster is picked with a different value than the first
				while (clust_merge_index2 == clust_merge_index1)
					clust_merge_index2 = PickSizeT(0, getLastClustIndex());
				
			}

			// Check if clusters can merge
			if(canMerge(g,clust_merge_index1, clust_merge_index2))
			{		
				mergeClusters(clust_merge_index1, clust_merge_index2, false, true);
				calculateIntraEdgesCore(g, clust_merge_index1);
				getCluster(clust_merge_index2)->setNumIntraClust(0);
			}
		}
	}
	//}

	
	prob = PickProbabilty();
	if (prob < REDISRIB_PROB)
	{
		// Check if chromosome has enough cluster to merge
		if (getNumClusters() > 1)
		{
			std::size_t clust_merge_index1, clust_merge_index2;

			// ***************
			// Patch
			// ***************
			if (getNumClusters() == 2)
			{
				clust_merge_index1 = 0;
				clust_merge_index2 = 1;
			}
			else
			{
				// Pick first cluster to merge
				clust_merge_index1 = PickSizeT(0, getLastClustIndex());

				// Initialize second cluster to merge to the first cluster
				clust_merge_index2 = clust_merge_index1;

				// Loop until a second cluster is picked with a different value than the first
				while (clust_merge_index2 == clust_merge_index1)
					clust_merge_index2 = PickSizeT(0, getLastClustIndex());
			}
			
			if (canMerge(g, clust_merge_index1, clust_merge_index2))
				redistribute(g, mem, clust_merge_index1, clust_merge_index2);			
		}
	}

	// Loop through each node and attempt mutation
	for (std::size_t i_node = 0; i_node < num_nodes_mut; ++i_node)
	{
		prob = PickProbabilty();
		if (prob < NODE_MUT_PROB)//(PickProbabilty() < NODE_MUT_PROB)
		{
			// Pick node at random from network
			node = PickSizeT(0, getNumNodes()-1);
			pointMutation(g,node);
		}
	}

	
	// Save the number of clusters
	std::size_t tmp_num_cust = getNumClusters();

	// Get a pointer to the first cluster
	Cluster *cur_clust = getFirstCluster();
	// Loop through all the clusters
	for (int i_clust = 0; i_clust < tmp_num_cust; i_clust++)
	{
		if (cur_clust->needClusterFixed())
			fixUnconnectedClusters(g, i_clust);
		
		cur_clust = cur_clust->getNextCluster(); // Go to next cluster
	}
	
	return true;
}

bool Chromosome::maintainNumClusters(Graph &g, SplitMem &mem, const std::size_t min_num_clust, const std::size_t max_num_clust)
{
	bool num_clust_changed = false;

	// Check if the chromosome has more than the maximum allowed number of clusters
	while (getNumClusters() > max_num_clust)
	{
		// Pick first cluster to merge
		std::size_t clust_merge_index1 = PickSizeT(0, getLastClustIndex());

		// Initialize second cluster to merge to the first cluster
		int clust_merge_index2 = clust_merge_index1;

		// Loop until a second cluster is picked with a different value than the first
		while (clust_merge_index2 == clust_merge_index1)
			clust_merge_index2 = PickSizeT(0, getLastClustIndex());
		
		
		// Check if clusters can merge
		if (canMerge(g, clust_merge_index1, clust_merge_index2))
		{
			mergeClusters(clust_merge_index1, clust_merge_index2, false, true);
			calculateIntraEdgesCore(g, clust_merge_index1);
			getCluster(clust_merge_index2)->setNumIntraClust(0);

			num_clust_changed = true;
		}
	}

	while (getNumClusters() < min_num_clust)
	{
		// Pick a random cluster to split
		std::size_t clust_split_index = PickSizeT(0, getLastClustIndex());

		// Get a pointer to the cluster
		Cluster* clust_to_split = this->getCluster(clust_split_index);

		// Check if the cluster has enough nodes to split
		if (clust_to_split->getNumNodes() > 1)
		{
			splitCluster(g, mem, clust_split_index);
			num_clust_changed = true;
		}
	}

	return num_clust_changed;
}

bool Chromosome::moveNode(Graph &g, const std::size_t old_clust_index, const std::size_t new_clust_index, const std::size_t node_id)
{
	assert(old_clust_index < num_clusters);
	assert(new_clust_index < num_clusters);

	// Check if old and new cluster indicies are the same
	if (old_clust_index == new_clust_index)
	{
		printf("Warning in Chromosome::MoveNode - Moving node from cluster %lu to cluster %lu.\n", old_clust_index, new_clust_index);
		return true;
	}

	// Get pointer to old cluster and new cluser
	Cluster *old_clust = getCluster(old_clust_index);
	Cluster *new_clust = getCluster(new_clust_index);

	// Save the degree of the node to move
	std::size_t degree = g.getDegree(node_id);

	// Remove node from old cluster
	Node *node_to_move = old_clust->removeNode(node_id, degree);

	// Check that a valid node was returned
	if (node_to_move == NULL)
	{
		printf("ERROR in Chromosome::MoveNode - Node %lu not removed from cluster %lu.\n", node_id, old_clust_index);
		return false;
	}

	// Add node to new cluster
	new_clust->addNode(node_to_move, degree);

	// Update cluster number
	setClusterNum(node_id, new_clust_index);

	// Update the intra cluster edge values
	std::size_t old_intra = 0, new_intra = 0, target_id;
	// Get a pointer to vertex 'node_id' in graph
	Vertex *cur_vert = g.getPtrToVertex(node_id);

	// Loop through all edges in the list
	for(std::size_t i = 0; i < cur_vert->getDegree(); ++i)
	{
		// Save the target ID
		target_id = cur_vert->getTargetId(i);
		
		// Check if the neighbor is in the old cluster
		if (getClusterNum(target_id) == old_clust_index)
		{
			++old_intra;
			setNumIntraClusterEdges(target_id, getNumIntraClustEdges(target_id)-1);
		}

		// Check if the neighbot is in the new cluster
		else if (getClusterNum(target_id) == new_clust_index)
		{
			++new_intra;
			setNumIntraClusterEdges(target_id, getNumIntraClustEdges(target_id) + 1);
		}		
	}

	old_clust->setNumIntraClust(old_clust->getNumIntraClust() - old_intra);
	new_clust->setNumIntraClust(new_clust->getNumIntraClust() + new_intra);
	setNumIntraClusterEdges(node_id, new_intra);

	return true;
}

bool Chromosome::moveNode(Cluster * old_clust, Cluster * new_clust, const std::size_t node_id, const std::size_t degree)
{
	// Check that old_clust is a valid pointer
	if (old_clust == NULL)
	{
		printf("ERROR in Chromosome::MoveNode - old_clust is a null pointer.\n");
		return false;
	}

	// Check that new_clust is a valid pointer
	if (new_clust == NULL)
	{
		printf("ERROR in Chromosome::MoveNode - new_clust is a null pointer.\n");
		return false;
	}

	// Remove node from old cluster
	Node *node_to_move = old_clust->removeNode(node_id, degree);

	// Check that a valid node was returned
	if (node_to_move == NULL)
	{
		printf("ERROR in Chromosome::MoveNode - Node %lu not removed from cluster. \n", node_id);
		return false;
	}

	// Add node to new cluster
	new_clust->addNode(node_to_move, degree);

	return true;
}

void Chromosome::recordToFile(std::string &filename, const char * clust_title)
{
	// Declare local variables
	FILE *output;
	Cluster *cur_clust;
	Node *cur_node;

	// Try to open output file
	if((output = fopen(filename.c_str(), "a+")) == NULL)
	{
		printf("ERROR in Chromosome::RecordToFile - Could not open file: %s.\n", filename.c_str());
		exit(1);
	}

	fprintf(output, "Fitness: %0.3f (%f)\n", fitness, fitness);
	fprintf(output, "Initial number of clusters: %lu\n", num_init_clusters);
	fprintf(output, "Number of splits: %lu\n", num_splits);
	fprintf(output, "Number of merges: %lu\n", num_merges);
	fprintf(output, "Number of redistributions %lu\n", num_redis);
	fprintf(output, "Number of point mutations: %lu\n", num_point_muts);
	fprintf(output, "Number of node crossovers: %lu\n", num_node_crossovers);
	fprintf(output, "Number of cluster crossovers: %lu\n", num_cluster_crossovers);
	fprintf(output, "Chromosome created in generation %lu\n", generation_created);

	// Write header
	//fprintf(output, "Cluster title: %s\n", clust_title);
	fprintf(output, "\n");
	fprintf(output, "%lu nodes expected\n", getNumNodes());
	fprintf(output, "Cluster # (# nodes)\t\t|\tNodes\n");

	// Get pointer to first clust
	cur_clust = getFirstCluster();

	// Loop through all clusters in chromosome
	for (std::size_t i_clust = 0; i_clust < getNumClusters(); ++i_clust)
	{
		// Print cluster number and number of nodes in current cluster
		fprintf(output, "%lu (%lu)\t\t\t\t|\t", i_clust, cur_clust->getNumNodes());

		// Get pointer to first node in current cluster
		cur_node = cur_clust->getFirstNode();

		// Loop through all nodes in cluster
		for (std::size_t i_node = 0; i_node < cur_clust->getNumNodes(); ++i_node)
		{
			// Check if this is the last node in the cluster
			if (i_node == cur_clust->getNumNodes() - 1)
				fprintf(output, "%lu\n", cur_node->getID());
			else
				fprintf(output, "%lu, ", cur_node->getID());

			// Go to next node
			cur_node = cur_node->getNextNode();
		}

		// Go to next cluster
		cur_clust = cur_clust->getNextCluster();
	}

	// Print extra new line in file
	fprintf(output, "\n");

	// Close file
	fclose(output);
}

// 0 - Chromosome has no errors
// 1 - The number of reported clustes does not match the calculated number of clusters
// 2 - The number of reported nodes in the chromosome does not match the calculated number of nodes
// 3 - The cluster is claiming to have zero nodes, but as a non null pointer to the first node
// 4 - The cluster is claiming to have zero nodes, but as a non null pointer to the last node
// 5 - The cluster is claiming to have nodes, but has a null pointer to the first node
// 6 - The cluster is claiming to have nodes, but has a null pointer to the last node
// 7 - Number of nodes listed in cluster does not match the number of nodes calculated in cluster
// 8 - The ID of the last node found does not matche the ID of the last node listed
// 9 - Cluster contains unconnected nodes
int Chromosome::verifyChromosome(Graph &g)
{
	// Declare local variables
	std::size_t num_nodes_total = 0, num_nodes_clust = 0, num_clust = 0;
	Cluster *cur_clust;
	Node *cur_node, *last_node, *last_node_found;

	// Get pointer to first cluster
	cur_clust = getFirstCluster();

	// Check if pointer to first cluster is valid
	if (cur_clust == NULL)
	{
		// Check if the number of reported cluster matches
		if (this->getNumClusters() != 0)
			return 1; // Cluster is empty, but claims to have clusters

		// Check if the number of nodes is correct
		if (this->getNumNodes() != 0)
			return 2; // Cluster is empty, but claims to have nodes
	}
	else
	{
		// Loop through all clusters in the linked list
		while (cur_clust != NULL)
		{
			// Get a pointer to the first node
			cur_node = cur_clust->getFirstNode();

			// Check if the cluste claims to have nodes
			if (cur_clust->getNumNodes() == 0)
			{
				if (cur_node != NULL)
					return 3;

				if (cur_clust->getLastNode() != NULL)
					return 4;
			}
			else
			{
				// Get a pointer to the last node
				last_node = cur_clust->getLastNode();

				// Check if the current cluste has a valid first node pointer
				if (cur_node == NULL)
					return 5;

				// Check if the current cluste has a valid first node pointer
				if (last_node == NULL)
					return 6;

				// Reset local counter
				num_nodes_clust = 0;

				// Initialize the last_node_found
				last_node_found = NULL;

				// Loop through all nodes in cluster
				while (cur_node != NULL)
				{
					last_node_found = cur_node; // Update the last_node_found
					num_nodes_clust++; // Increment the number of nodes counted for this cluster
					cur_node = cur_node->getNextNode(); // Go to the next node
				}
				// Add the number of nodes found in the current cluster to the numbe of nodes found in the chromosome
				num_nodes_total += num_nodes_clust;

				// Check that the number of nodes calculated matches the number reported
				if (num_nodes_clust != cur_clust->getNumNodes())
					return 7;
				
				// Check if the last_node_found matches the node pointed to as the last node of the cluster
				if (last_node_found->getID() != last_node->getID())
					return 8;
				
				// Check if all nodes in cluste are connected via BFS
				if (!cur_clust->verifyBFS(g, this->cluster_num))
					return 9;
				
				// Check if the recorded degree sum equals the calculated value
				std::size_t degree_sum = cur_clust->getDegreeSum();
				cur_clust->calculateDegreeSum(g);
				if (degree_sum != cur_clust->getDegreeSum())
					return 10;
				

				std::size_t tmp_intra_sum = cur_clust->getNumIntraClust();
				calculateIntraEdgesCore(g, num_clust);
				if (tmp_intra_sum != cur_clust->getNumIntraClust())
					return 11;
				
			}

			// Increment the number of calculated cluster
			++num_clust;

			// Go to next cluster
			cur_clust = cur_clust->getNextCluster();
		}

		// Check if the number of clusters calculated matches the reported value
		if (num_clust != this->getNumClusters())
			return 1;

		// Check if the number of nodes calculate matched the reported value
		if (num_nodes_total != this->getNumNodes())
			return 2;
	}

	return 0;
}

void Chromosome::fixUnconnectedClusters(Graph &g, const std::size_t clust_to_fix_index)
{
	Cluster *cur_clust;
	Node *cur_node;
	Vertex  *cur_vert;
	std::size_t start_index, end_index, node_id, cur_comp, cur_clust_num;	
	
	// Initialize indicies for queue
	start_index = 0;
	end_index = 0;

	// Get pointer to cluster to fix
	cur_clust = this->getCluster(clust_to_fix_index);

	// If the cluster is empty of a singleton there is no need to fix it
	if (cur_clust->getNumNodes() < 2)
		return;
	
	// Allocate memory for arrays
	std::vector<int> queue(cur_clust->getNumNodes());
	std::vector<std::size_t> node_id_ref(cur_clust->getNumNodes());
	std::vector<int> component(g.getNumNodes());
	std::vector<bool> visisted(g.getNumNodes(), false);	

	// Get a pointer to the first node in the cluster to fix
	cur_node = cur_clust->getFirstNode();

	// Indicate which nodes are in the cluster
	for (int j = 0; j < cur_clust->getNumNodes(); j++)
	{
		node_id_ref[j] = cur_node->getID();
		cur_node = cur_node->getNextNode();
	}

	// Reset the pointer to the first node in the cluster
	cur_node = cur_clust->getFirstNode();

	// Initialize starting component number
	cur_comp = 0;
	cur_clust_num = getClusterNum(cur_node->getID());

	// Loop through all nodes in cluster
	for (std::size_t i_node = 0; i_node < cur_clust->getNumNodes(); ++i_node)
	{
		// Check if node has been visited
		if (!visisted[cur_node->getID()])
		{
			// Add node to queue and increment end_index
			queue[end_index++] = cur_node->getID();
			// Set as visited
			visisted[cur_node->getID()] = true;
			// Set component value
			component[cur_node->getID()] = cur_comp;

			// Loop until queue is empty
			while (start_index < end_index)
			{
				// Get node ID and dequeue
				node_id = queue[start_index++];

				// Get a pointer to the 'node_id' vertex
				cur_vert = g.getPtrToVertex(node_id);

				// Loop through all edges in adjacency list
				for(std::size_t i_edge = 0; i_edge < cur_vert->getDegree(); ++i_edge)
				{
					// Check if the target is in the cluster and unvisited
					if ((!visisted[cur_vert->getTargetId(i_edge)]) &&
						 (cur_clust_num == cluster_num[cur_vert->getTargetId(i_edge)]))
					{
						// Add to queue
						queue[end_index++] = cur_vert->getTargetId(i_edge);
						// Mark as visited
						visisted[cur_vert->getTargetId(i_edge)] = true;
						// Set component
						component[cur_vert->getTargetId(i_edge)] = cur_comp;
					}					
				}
			}

			// Increment current component number
			++cur_comp;
		}

		// Go to next node
		cur_node = cur_node->getNextNode();
	}

	// Check if unconnected parts exist
	if (cur_comp > 1)
	{
		std::size_t orig_num_nodes = cur_clust->getNumNodes(); // This is the number of nodes in the cluster_to_fix prior to moving nodes
		std::size_t ref_index = getLastClustIndex(); // This is the index of the last cluster prior to adding more clusters

		// Add new clusters to chromosome
		for (std::size_t i_clust = 0; i_clust < cur_comp - 1; ++i_clust)
			addCluster();
		
		// Loop through all nodes in cluster, starting at the end
		for (std::size_t j_node = 0; j_node < orig_num_nodes; ++j_node)
		{
			// Check if node is in a different cluster
			if (component[node_id_ref[j_node]] != 0)
			{
				moveNode(g, clust_to_fix_index, ref_index + component[node_id_ref[j_node]], node_id_ref[j_node]);
			}
		}
	}

	// Set Fix Cluster flag
	cur_clust->setFixCluster(false);

	return;
}

void Chromosome::calculateFitness(Graph &g, const int fitness_calc)
{
	// Reset fitness
	fitness = 0;

	// Get pointer to first cluster
	Cluster *cur_clust = getFirstCluster();

	double num_nodes = (double)g.getNumNodes();
	double tmp = num_nodes * (num_nodes - 1);

	// Determine which fitnes function to call based on input
	switch (fitness_calc)
	{
	case (0):
		// Loop through all clusters in chromosome
		while (cur_clust != NULL)
		{
			// Check if the fitness needs to be recalculated for the current cluster
			if (cur_clust->needFitnessRecalc())
			{
				// Calculate fraction of intra-cluster edges in calling cluster. (The edges were counted twice, so divide by 2|E|)
				double intra_fract = (double)cur_clust->getNumIntraClust() / ((double)g.getTotalNumEdges());
				
				// Calculate fraction of expected edges				
				double expect_fract = (double)(cur_clust->getDegreeSum()) / (2*(double)g.getTotalNumEdges());
				
				// Save the fitness value for this cluster
				cur_clust->setFitness(intra_fract - expect_fract * expect_fract);

				// Set flag to indicate the fitness does not need to be recalculated
				cur_clust->setRecalcFitness(false);
			}
			// Add the fitness value for this cluster to the sum
			this->fitness += cur_clust->getFitness();

			// Go to next cluster
			cur_clust = cur_clust->getNextCluster();
		}
		break;

	case(1):
		// Loop through all clusters in chromosome
		while (cur_clust != NULL)
		{
			// Check if the fitness needs to be recalculated for the current cluster
			if (cur_clust->needFitnessRecalc())
			{
				// Calculate fraction of intra-cluster edges in calling cluster. (The edges were counted twice, so divide by 2|E|)
				double intra_fract = (double)cur_clust->getNumIntraClust() / ((double)g.getNumEdges());
				
				// Calculate fraction of expected edges				
				double expect_fract = (double)(cur_clust->getNumNodes() * (cur_clust->getNumNodes() -1));
				expect_fract /= tmp;
				
				// Save the fitness value for this cluster
				cur_clust->setFitness(intra_fract - expect_fract);

				// Set flag to indicate the fitness does not need to be recalculated
				cur_clust->setRecalcFitness(false);
			}
			// Add the fitness value for this cluster to the sum
			fitness += cur_clust->getFitness();

			// Go to next cluster
			cur_clust = cur_clust->getNextCluster();
		}
		fitness =  fitness * g.getNumEdges() / g.getTotalNumEdges();
		break;
	
	default:
		printf("ERROR in Chromosome::CalculateFitness - Unknown fitness function.\n");
		break;
	}
}

void Chromosome::sortChromosome()
{
	Cluster *cur_clust, *pre_clust, *tmp_clust;
	std::size_t i, j;

	// Delete any empty clusters. If at least one cluster was deleted, update the clust_num variable
	if (deleteAllEmptyClusters())
		updateAllClusterNum();
			   
	cur_clust = getFirstCluster();

	// Loop through each cluster
	while(cur_clust != NULL)
	{
		// Rearrange the nodes in current cluster from highest to lowest
		cur_clust->sortNodes();
		
		// Go to next cluster
		cur_clust = cur_clust->getNextCluster();
	}

	// Initialize 
	i = 1;
	
	while (i < getNumClusters())
	{
		j = i;
		while ((j > 0) && (getHeadNode(j - 1)->getID() > getHeadNode(j)->getID()))
		{
			if (j == 1)
			{
				pre_clust = getFirstCluster();
				cur_clust = pre_clust->getNextCluster();

				// Point the previous cluster at the next cluster
				pre_clust->setNextCluster(cur_clust->getNextCluster());

				cur_clust->setNextCluster(pre_clust);

				// Set the current cluster as the head cluster
				setFirstCluster(cur_clust);
			}
			else
			{
				tmp_clust = getCluster(j - 2);
				pre_clust = tmp_clust->getNextCluster();
				cur_clust = pre_clust->getNextCluster();

				// Point the previous cluster at the next cluster
				pre_clust->setNextCluster(cur_clust->getNextCluster());

				// Point the tmp cluster at the current cluster
				tmp_clust->setNextCluster(cur_clust);

				// Point the current cluster at the previuos cluster
				cur_clust->setNextCluster(pre_clust);
			}

			j--;
		}
		i++;
	}	

	return;
}

void Chromosome::replaceNodeNum(const std::vector<std::size_t> &nn_vec)
{
	// Declare local variables
	Cluster *cur_clust;
	Node *cur_node;

	// Get a pointer to the first cluster
	cur_clust = this->getFirstCluster();

	// Loop through all clusters
	while (cur_clust != NULL)
	{
		// Get a pointer to the first node in the current cluster
		cur_node = cur_clust->getFirstNode();

		// Loop through all nodes in the current cluster
		while (cur_node != NULL)
		{
			// Set new node ID value
			cur_node->setID(nn_vec[cur_node->getID()]);

			// Go to next node
			cur_node = cur_node->getNextNode();
		}
		// Go to next cluster
		cur_clust = cur_clust->getNextCluster();
	}
}

Node * Chromosome::getHeadNode(const std::size_t cluster_index)
{
	assert(cluster_index < num_clusters);

	// Get a pointer to the cluster at the specified index
	Cluster *cur_clust = getCluster(cluster_index);

	return cur_clust->getFirstNode();
}

bool Chromosome::deleteAllEmptyClusters()
{
	Cluster *cur_clust, *pre_clust;
	bool deleted_cluster = false;

	// Get a pointer to the first node
	cur_clust = getFirstCluster();
	pre_clust = NULL;

	// Loop through all cluster
	while (cur_clust != NULL)
	{
		// Check if the cluster is empty
		if (cur_clust->getNumNodes() == 0)
		{
			// Check if this is the head clust of the chromosome
			if (pre_clust == NULL)
			{
				// Point the chromosome at the next cluster
				this->setFirstCluster(cur_clust->getNextCluster());

				if (cur_clust->getFirstNode() != NULL)
				{
					printf("ERROR - trying to delete non empty cluster (head node).\n");
				}
				if (cur_clust->getLastNode() != NULL)
				{
					printf("ERROR - trying to delete non empty cluster (last node).\n");
				}

				// Delete the empty cluster
				cur_clust->setNextCluster(NULL);
				delete cur_clust;
				setNumClusters(getNumClusters() - 1);
				deleted_cluster = true;

				// Get the new head cluster
				cur_clust = getFirstCluster();
			}
			else
			{
				// Point the previous cluster at the next cluster
				pre_clust->setNextCluster(cur_clust->getNextCluster());

				if (cur_clust->getFirstNode() != NULL)
				{
					printf("ERROR - trying to delete non empty cluster (head node).\n");
				}
				if (cur_clust->getLastNode() != NULL)
				{
					printf("ERROR - trying to delete non empty cluster (last node).\n");
				}

				// Delete the empty cluster
				cur_clust->setNextCluster(NULL);
				delete cur_clust;
				this->setNumClusters(getNumClusters() - 1);
				deleted_cluster = true;

				// Get the next cluster
				cur_clust = pre_clust->getNextCluster();
			}
		}
		else
		{
			// Go to next cluster
			pre_clust = cur_clust;
			cur_clust = cur_clust->getNextCluster();
		}
	}

	return deleted_cluster;
}

void Chromosome::updateAllClusterNum()
{
	Cluster *cur_clust;
	Node *cur_node;
	std::size_t clust_num = 0;

	// Get a pointer to the first cluster
	cur_clust = getFirstCluster();

	// Loop through all clusters
	while (cur_clust != NULL)
	{	
		// Get a pointer to the first node in the cluster
		cur_node = cur_clust->getFirstNode();
		
		// Loop through all nodes in the cluster
		while (cur_node != NULL)
		{
			setClusterNum(cur_node->getID(), clust_num); // Update cluster number
			cur_node = cur_node->getNextNode(); // Go to next node			
		}
		
		cur_clust = cur_clust->getNextCluster(); // Go to next cluster
		++clust_num; // Increment cluster number
	}
}

bool Chromosome::checkClusterNumber()
{
	Cluster *cur_clust;
	Node *cur_node;
	std::size_t clust_num = 0;

	// Get a pointer to the first cluster
	cur_clust = getFirstCluster();
	
	// Loop through all clusters
	while (cur_clust != NULL)
	{
		// Get a pointer to the first node in the cluster
		cur_node = cur_clust->getFirstNode();

		// Loop through all nodes in the cluster
		while (cur_node != NULL)
		{
			// Check if the record cluster number matches the actual cluster number
			if (clust_num != cluster_num[cur_node->getID()])
			{
				printf("Cluster number is not correct for node %lu\n", cur_node->getID());
				return false;
			}
			
			cur_node = cur_node->getNextNode(); // Go to next node
		}

		cur_clust = cur_clust->getNextCluster(); // Go to next cluster
		clust_num++; // Increment cluster number
	}

	return true;
}

void Chromosome::calculateIntraEdges(Graph &g)
{
	// Get a pointer to the first cluster
	Cluster *cur_clust = getFirstCluster();

	for (std::size_t i_clust = 0; i_clust < getNumClusters(); ++i_clust)
		calculateIntraEdgesCore(g, i_clust);	
}

void Chromosome::calculateIntraEdgesCore(Graph &g, const std::size_t clust_index)
{
	Cluster *cur_clust;
	Node *cur_node;
	Vertex *cur_vert;
	std::size_t num_edges;

	// Get a pointer to the current cluster
	cur_clust = getCluster(clust_index);

	// Reset the number of intra cluster edges
	cur_clust->setNumIntraClust(0);

	// Get a pointer to the first node in the current cluster
	cur_node = cur_clust->getFirstNode();

	// Loop through all nodes in the cluster
	while (cur_node != NULL)
	{
		// Get a pointer to the vertex in the graph
		cur_vert = g.getPtrToVertex(cur_node->getID());

		// Reset local counter
		num_edges = 0;

		// Loop through all edges
		for(std::size_t i = 0; i < cur_vert->getDegree(); ++i)		
		{
			// Check if the target of the egde is in the same cluster as the source
			if (getClusterNum(cur_vert->getTargetId(i)) == clust_index)
				++num_edges; // Increment the number of edges from this node						
		}

		// Record the number of intra cluster egdes for the current node
		this->setNumIntraClusterEdges(num_edges, cur_node->getID());

		// Update the number of intra cluster edges for the current cluster
		cur_clust->setNumIntraClust(cur_clust->getNumIntraClust() + num_edges);

		// Go to next node
		cur_node = cur_node->getNextNode();
	}

	// Account for double counting
	cur_clust->setNumIntraClust(cur_clust->getNumIntraClust()/2);
}

bool Chromosome::isNodeMarginal(Graph &g, const std::size_t node_id)
{
	// Check if node is a boundary node
	if (g.getDegree(node_id) == getNumIntraClustEdges(node_id))
		return false; // If all of the edges are intra cluster, the node is not on the boundary between mulitple clusters
		
	Node *cur_node;
	Cluster *cur_clust;
	int cur_clust_num;
	bool clust_connected;
	// Save the cluster number that the node is currently in
	cur_clust_num = getClusterNum(node_id);

	// Get a pointer to the current cluster
	cur_clust = getCluster(cur_clust_num);

	// Check if node is a singleton
	if (cur_clust->getNumNodes() == 1)
		return true;
	
	// Remove the node from the current cluster
	cur_node = cur_clust->removeNode(node_id, g.getDegree(node_id));
	// Update the cluster number
	setClusterNum(node_id, getNumClusters());

	// Check if the cluster is connected with the node removed
	clust_connected = cur_clust->verifyBFS(g, cluster_num);

	cur_clust->addNode(cur_node, g.getDegree(node_id)); // Add Node back
	setClusterNum(node_id, cur_clust_num); 

	return clust_connected; // If the cluster is still connected after removing the node, the node is marginal
}

bool Chromosome::readInCluster(std::string &file_name, Graph &g)
{
	// Declare local variables
	Cluster *new_clust = NULL, *cur_clust = NULL;
	std::string tmp_str;
	std::istringstream iss;
	std::ifstream input;
	std::size_t clust_num, tmp_num_nodes, node_id;
	char junk;

	// Open the file
	input.open(file_name.c_str());

	// Check if file opened 
	if (!input)
	{
		std::cerr << "ERROR in Chromosome::ReadInCluster - Could not open file: " << file_name.c_str() << std::endl;
		return false;
	}

	// Read the first 3 lines from the file
	std::getline(input, tmp_str);
	std::getline(input, tmp_str);
	std::getline(input, tmp_str);

	while (true)
	{
		char temp = input.peek();

		// Check if all clustes have been read in
		if (temp == '\n')
		{
			break;
		}

		// Read in first cluster from file
		std::getline(input, tmp_str, '|');

		// Clear then update the istringstream
		iss.clear();
		iss.str(tmp_str);

		iss >> clust_num;
		iss >> junk;
		iss >> tmp_num_nodes;

		// Allocate memory for new cluster
		new_clust = new Cluster();

		if (clust_num == 0)
		{
			first_cluster = new_clust;
			cur_clust = new_clust;
		}
		else
		{
			cur_clust->setNextCluster(new_clust);
			cur_clust = cur_clust->getNextCluster();
		}
		++num_clusters; // Increment number of clusters 
		num_nodes += tmp_num_nodes; // Add the number of nodes for the current cluster to the total for the chromosome

		std::getline(input, tmp_str);
		// Clear then update the istringstream
		iss.clear();
		iss.str(tmp_str);

		int i_node = 0;
		while (i_node < tmp_num_nodes)
		{
			iss >> node_id; // Get node ID
			iss >> junk; // Get comma

			cur_clust->addNode(node_id, g.getDegree(node_id));
			i_node++;
		}

	}

	// Close file
	input.close();

	return true;
}

bool Chromosome::appendCluster(const std::string &file_name, Graph &g)
{
	// Declare local variables
	Cluster *new_clust = NULL, *cur_clust = NULL;
	std::string tmp_str;
	std::istringstream iss;
	std::ifstream input;
	std::size_t clust_num, num_nodes_comp, node_id;
	char junk;
	
	// Open file
	input.open(file_name.c_str());

	// Check if file opened
	if (!input)
	{
		std::cerr << "ERROR in Chromosome::AppendCluster - Could not open file: " << file_name.c_str() << std::endl;
		exit(1);
	}

		// Move through file until clusters are reached
	while (true)
	{
		std::getline(input, tmp_str);
		
		if (input.eof())
		{
			std::cerr << "ERROR in Chromosome::AppendCluster - File not formatted correctly: " << file_name.c_str() << std::endl;
			exit(1);
		}
		

		// Check if the clusters have been reached
		if (tmp_str.find("Cluster #") != std::string::npos)
			break;		
	}
	
	// Get a pointer to the last cluster
	cur_clust = getCluster(getLastClustIndex());

	// Read in clusters
	while (true)
	{
		// Peak at the next char
		junk = input.peek();

		// Check if the end of clusters has been reached
		if ((junk == '\n') || (input.eof()))
		{
			break;
		}

		// Read in the cluster information
		std::getline(input, tmp_str, '|');

		// Clear the update the istringstream
		iss.clear();
		iss.str(tmp_str);

		iss >> clust_num;
		iss >> junk;
		iss >> num_nodes_comp;

		// Allocate memory for new cluster
		new_clust = new Cluster();

		if (cur_clust == NULL)
		{
			first_cluster = new_clust;
			cur_clust = new_clust;
		}
		else
		{
			cur_clust->setNextCluster(new_clust);
			cur_clust = cur_clust->getNextCluster();
		}

		++num_clusters;
		
		std::getline(input, tmp_str);
		// Clear then update the istringstream
		iss.clear();
		iss.str(tmp_str);

		std::size_t i_node = 0;
		while (i_node < num_nodes_comp)
		{
			iss >> node_id;
			iss >> junk;

			cur_clust->addNode(node_id, g.getDegree(node_id));
			++i_node;
		}
		setNumNodes(getNumNodes() + num_nodes_comp);
	}

	// Close file
	input.close();
	return true;
}

void Chromosome::allocateClusterNum()
{
	cluster_num.resize(getNumNodes());	
}

void Chromosome::allocateIntraEdges()
{
	num_intra_edges.resize(getNumNodes());	
}

bool Chromosome::recordOutput(std::string &file_name, const int num_edges)
{
	// Declare local variables
	FILE* output;
	
	// Try to open output file
	if((output = fopen(file_name.c_str(), "w")) == NULL)
	{
		printf("ERROR in Chromosome::RecordOutput - Could not open file(%s).\n", file_name.c_str());
		return false;
	}

	// Count number of non-singleton clusters
	std::size_t num_nsc = 0;
	Cluster *cur_clust = getFirstCluster();
	while (cur_clust != NULL)
	{
		//if(cur_clust->getNumNodes() > 1)
			num_nsc++;
		cur_clust = cur_clust->getNextCluster();
	}
	
	// Print clustering info
	fprintf(output, "%lu nodes %lu clusters %i edges\n", getNumNodes(), num_nsc, num_edges);
	// Loop through all nodes (except last node)
	for (std::size_t i_node = 0; i_node < getNumNodes()-1; ++i_node)
		fprintf(output, "%lu ", getClusterNum(i_node)+1);
	// Record the cluster number for the last node in the array
	fprintf(output, "%lu", getClusterNum(getNumNodes() - 1)+1);
		
	// Close File
	fclose(output);

	return true;
}

void Chromosome::readChromosome(FILE *input, Graph &g, const int fit_eq)
{
	// Try to read in number of initial clusters
	if (fscanf(input, "%lu", &num_init_clusters) == 0)
	{
		printf("ERROR in ClusteringResults::ReadPopulation - Population file not formatted correctly.\n");
		exit(1);
	}
	
	// Try to read in number of splits
	if (fscanf(input, "%lu", &num_splits) == 0)
	{
		printf("ERROR in ClusteringResults::ReadPopulation - Population file not formatted correctly.\n");
		exit(1);
	}

	if (fscanf(input, "%lu", &num_merges) == 0)
	{	
		printf("ERROR in ClusteringResults::ReadPopulation - Population file not formatted correctly.\n");
		exit(1);
	}

	if (fscanf(input, "%lu", &num_redis) == 0)
	{
		printf("ERROR in ClusteringResults::ReadPopulation - Population file not formatted correctly.\n");
		exit(1);
	}

	if (fscanf(input, "%lu", &num_point_muts) == 0)
	{
		printf("ERROR in ClusteringResults::ReadPopulation - Population file not formatted correctly.\n");
		exit(1);
	}

	if (fscanf(input, "%lu", &num_node_crossovers) == 0)
	{
		printf("ERROR in ClusteringResults::ReadPopulation - Population file not formatted correctly.\n");
		exit(1);
	}

	if (fscanf(input, "%lu", &num_cluster_crossovers) == 0)
	{
		printf("ERROR in ClusteringResults::ReadPopulation - Population file not formatted correctly.\n");
		exit(1);
	}

	if (fscanf(input, "%lu", &generation_created) == 0)
	{
		printf("ERROR in ClusteringResults::ReadPopulation - Population file not formatted correctly.\n");
		exit(1);
	}

	if (fscanf(input, "%i", &init_type) == 0)
	{
		printf("ERROR in ClusteringResults::ReadPopulation - Population file not formatted correctly.\n");
		exit(1);
	}

	//TEST
	std::size_t expected_num_clust;
	if (fscanf(input, "%lu", &expected_num_clust) == 0)
	{
		printf("ERROR in ClusteringResults::ReadPopulation - Population file not formatted correctly.\n");
		exit(1);
	}
	
	// Read in the next chromosome from the population file
	cluster_num.resize(g.getNumNodes());
	for (std::size_t i_node = 0; i_node < g.getNumNodes(); ++i_node)
	{
		if (fscanf(input, "%lu", &cluster_num[i_node]) == 0)
		{
			printf("ERROR in ClusteringResults::ReadPopulation - Population file not formatted correctly.\n");
			exit(1);
		}
	}

	// Rearrange linked list based on cluster number
	initialize(g, cluster_num);	

	// TEST
	if (expected_num_clust != getNumClusters())
	{
		printf("ERROR in ClusteringResults::ReadPopulation - Number of clustes not calculated correctly");
		exit(1);
	}

	// Calculate the fitness
	calculateFitness(g, fit_eq);	
}

void Chromosome::recordObjective(const std::string &file_name, const std::string &comp_num) const
{
	FILE *output;

	if ((output = fopen(file_name.c_str(), "a")) == NULL)
	{
		printf("ERROR in Chromosome::recordObjective - Could not open file (%s)\n", file_name.c_str());
		exit(1);
	}
	fprintf(output, "%s,%f\n", comp_num.c_str(), fitness);
	fclose(output);
}