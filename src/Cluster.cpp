#include "Cluster.h"
#include "Vertex.h"
#include "Graph.h"
#include "Node.h"
#include "UtilityFunction.h"
#include <queue>
#include <cmath>
#include "assert.h"
#include "algorithm"

Cluster::Cluster() :	num_nodes(0),
						fitness(0.0),
						recalc_fitness(true),
						fix_cluster(false),
						degree_sum(0),
						intra_clust_edges(0)
{
	first_node = NULL;
	last_node = NULL;
	next_cluster = NULL;	
}


Cluster::~Cluster()
{
	if (first_node)
		delete first_node;
	if (next_cluster)
		delete next_cluster;
}

Node * Cluster::getFirstNode()
{
	return first_node;
}

Node * Cluster::getLastNode()
{
	return last_node;
}

Cluster * Cluster::getNextCluster()
{
	return next_cluster;
}

std::size_t Cluster::getNumNodes() const
{
	return num_nodes;
}

double Cluster::getFitness() const
{
	return fitness;
}

bool Cluster::needFitnessRecalc() const
{
	return recalc_fitness;
}

bool Cluster::needClusterFixed() const
{
	return fix_cluster;
}

std::size_t Cluster::getDegreeSum() const
{
	return degree_sum;
}

std::size_t Cluster::getNumIntraClust() const
{
	return intra_clust_edges;
}

void Cluster::setFirstNode(Node * node_to_set)
{
	first_node = node_to_set;
}

void Cluster::setLastNode(Node * node_to_set)
{
	last_node = node_to_set;
}

void Cluster::setNextCluster(Cluster * clust_to_set)
{
	next_cluster = clust_to_set;
}

void Cluster::setNumNodes(const std::size_t _num_nodes)
{
	num_nodes = _num_nodes;
}

void Cluster::setFitness(const double _fitness)
{
	fitness = _fitness;
}

void Cluster::setRecalcFitness(const bool need_recalc)
{
	recalc_fitness = need_recalc;
}

void Cluster::setFixCluster(const bool verifiy_clust)
{
	fix_cluster = verifiy_clust;
}

void Cluster::setDegreeSum(const std::size_t _degree_sum)
{
	degree_sum = _degree_sum;
}

void Cluster::setNumIntraClust(const std::size_t num_edges)
{
	intra_clust_edges = num_edges;
}

void Cluster::addNode(const std::size_t node_id, const std::size_t node_degree)
{
	// Allocate memory for new node and set ID
	Node * new_node = new Node(node_id);

	// Check if cluster has head node
	if (this->first_node == NULL) // No head node
	{
		this->setFirstNode(new_node); // Set new node as head node
		this->setLastNode(new_node); // Set new node as last node
	}
	else
	{
		this->last_node->setNextNode(new_node); // Point existing last node to new node
		this->setLastNode(new_node); // Set new node as the last node
	}

	++num_nodes; // Increment number of nodes
	degree_sum += node_degree; // Update degree sum	
	recalc_fitness = true; // Set flag to recalculate fitness
}

void Cluster::addNode(Node * node_to_add, const std::size_t node_degree)
{
	// Check if cluster has head node
	if (this->getFirstNode() == NULL) // No head node
	{
		this->setFirstNode(node_to_add); // Set new node as head node
		this->setLastNode(node_to_add); // Set new node as last node
	}
	else
	{
		this->getLastNode()->setNextNode(node_to_add); // Point existing last node to new node
		this->setLastNode(node_to_add); // Set new node as the last node
	}

	++num_nodes; // Increment number of nodes
	degree_sum += node_degree; // Update the degree sum
	recalc_fitness = true; // Set flag to recalculate fitness	
}

// This function searches the cluster for the node with an ID that matches 'node_id'. If such a node is found, it is removed
// and returned. If no such node is found a NULL pointer is returned.
Node * Cluster::removeNode(const std::size_t node_id, const std::size_t node_degree)
{
	// Get pointer to first node in the cluster
	Node *node_to_move = this->getFirstNode();
	Node *pre_node = NULL;

	bool node_found = false;

	// Loop until end of linked list is reached or node is found
	while (node_to_move != NULL)
	{
		// Check if the ID of the current node_to_move matches 'node_id'
		if (node_to_move->getID() == node_id) // Node ID matches
		{
			node_found = true; // Set flag indicating a node was found with desired ID
			break; // Break from while loop
		}
		else // Node ID does not match
		{
			pre_node = node_to_move; // Set current node_to_move as previous node
			node_to_move = node_to_move->getNextNode(); // Go to next node in cluster
		}
	}

	// Check if node was found in the cluster
	// This covers the case where the cluster has no nodes and the case where it has nodes, but none with an ID
	// that matches 'node_id'
	if (!node_found)
	{
		printf("ERROR in Cluster::RemoveNode - Could not find node ID %lu in cluster.\n", node_id);
		return NULL;
	}

	// Check if 'node_to_move' is the only node in the cluster
	if (this->getNumNodes() == 1)
	{
		this->setFirstNode(NULL); // Set pointer to first node to NULL
		this->setLastNode(NULL); // Set pointer to last node to NULL
	}
	// Check if the old cluster has multiple nodes and 'node_to_move' is the first
	else if (pre_node == NULL)
	{
		this->setFirstNode(node_to_move->getNextNode()); // Set the next node as the head of linked list
	}
	// Check if the old cluster has multiple nodes and 'node_to_move' is last
	else if (this->getLastNode()->getID() == node_to_move->getID())
	{
		pre_node->setNextNode(NULL); // Update pointer of previuos node
		this->setLastNode(pre_node); // Set the previous node as the last node
	}
	// Old cluster has multiple nodes and 'node_to_move' is neither first nor last
	else
	{
		pre_node->setNextNode(node_to_move->getNextNode());
	}

	// Update the next node for node_to_move
	node_to_move->setNextNode(NULL);

	--num_nodes; // Decrement the number of nodes in the old cluster
	degree_sum -= node_degree; // Update degree sum
	recalc_fitness = true; // Set flag to recalculate fitness

	return node_to_move;
}

std::size_t Cluster::getNodeID(const std::size_t node_index)
{
	assert(first_node != NULL); // Check if cluster is empty
	assert(node_index < num_nodes);

	// Get pointer to first node
	Node *cur_node = getFirstNode();

	// Loop until the node_index is reached
	for (std::size_t i_node = 0; i_node < node_index; ++i_node)
		cur_node = cur_node->getNextNode();

	return cur_node->getID();
}

std::size_t Cluster::getNodeIndex(const std::size_t node_id)
{
	// Initialize local variables
	std::size_t index = 0;
	bool node_found = false;
	Node *cur_node = getFirstNode();

	// Loop until end of list of until node_id is found
	while (cur_node != NULL)
	{
		// Check if the ID of the current node matches the desired ID
		if (cur_node->getID() == node_id)
		{
			node_found = true; // Set flag
			break; // Break while loop
		}
		else
		{
			cur_node = cur_node->getNextNode(); // Go to next node
			++index; // Increment index in linked list of current node
		}
	}

	// Check if the desired node id was found
	if(!node_found)
	{
		printf("ERROR in Cluster::getNodeIndex - No node with ID %lu found\n", node_id);
		exit(1);
	}

	return index; // Return ID of current node
}

bool Cluster::containsNode(const std::size_t node_id)
{
	// Get pointer to first node
	Node *cur_node = getFirstNode();

	// Loop until end of list of until node_id is found
	while (cur_node != NULL)
	{
		// Check if the ID of the current node matches the desired ID
		if (cur_node->getID() == node_id)
			return true;
		
		cur_node = cur_node->getNextNode(); // Go to next node		
	}

	return false; // Node was not found
}

void Cluster::copyNonPointers(Cluster *in_clust)
{
	in_clust->setNumNodes(this->getNumNodes());
	in_clust->setFitness(this->getFitness());
	in_clust->setRecalcFitness(this->needFitnessRecalc());
	in_clust->setFixCluster(this->needClusterFixed());
	in_clust->setDegreeSum(this->getDegreeSum());
	in_clust->setNumIntraClust(this->getNumIntraClust());
}

void Cluster::calculateDegreeSum(Graph &g)
{
	degree_sum = 0; // Reset degree sum

	// Get a pointer to the first node in the cluster
	Node *cur_node = getFirstNode();	

	// Loop through all nodes in the cluster
	while (cur_node != NULL)
	{
		degree_sum += g.getDegree(cur_node->getID());
		cur_node = cur_node->getNextNode();
	}
}

void Cluster::sortNodes()
{
	Node *cur_node = getFirstNode();
	std::vector<std::size_t> nodes(getNumNodes());

	for(std::size_t i = 0; i < getNumNodes(); ++i)
	{
		nodes[i] = cur_node->getID();
		cur_node = cur_node->getNextNode();
	}

	std::sort(nodes.begin(), nodes.end());

	cur_node = getFirstNode();
	for(std::size_t i = 0; i < getNumNodes(); ++i)
	{
		cur_node->setID(nodes[i]);
		cur_node = cur_node->getNextNode();
	}
}

bool Cluster::verifyBFS(Graph &g, std::vector<std::size_t> &clust_num)
{
	Vertex *cur_vert;
	std::size_t start_index, end_index, node_id, node_count, target_id, cur_clust_num;	

	// Initialize indicies for queue
	start_index = 0;
	end_index = 0;

	// Allocate memory for arrays
	std::vector<bool> visisted(g.getNumNodes(), false);
	std::vector<std::size_t> queue(this->getNumNodes());
		

	// Add the first node in the cluster to the queue
	node_id = getFirstNode()->getID();
	queue[end_index++] = node_id;
	visisted[node_id] = true;
	node_count = 1;

	// Save current cluster
	cur_clust_num = clust_num[node_id];

	// Loop until the queue is empty
	while (start_index < end_index)
	{
		// Get node id and dequeue
		node_id = queue[start_index++];

		// Get a pointer to the vertex
		cur_vert = g.getPtrToVertex(node_id);
		
		// Loop through all edges
		for(std::size_t i = 0; i < cur_vert->getDegree(); ++i)
		{
			// Save target id
			target_id = cur_vert->getTargetId(i);

			// Check if the target is unvisited and in current cluster
			if (!visisted[target_id] && (clust_num[target_id] == cur_clust_num))
			{
				// Add to queue
				queue[end_index++] = target_id;
				// Mark as visited
				visisted[target_id] = true;
				// Inrement the number of nodes visisted
				node_count++;
			}
		}		
	}
	
	// Check if all nodes in cluster were reached
	if (node_count == getNumNodes())
	{
		return true; // All nodes reached with BFS
	}
	else
	{
		return false; // Some nodes not reahed with BFS
	}
}