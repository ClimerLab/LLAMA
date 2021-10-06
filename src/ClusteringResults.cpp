#include "ClusteringResults.h"
#include "IslandStat.h"
#include "Cluster.h"
//#include "PopulationStat.h"
#include "UtilityFunction.h"
#include "Graph.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <assert.h>


ClusteringResults::ClusteringResults(ConfigParser &_parser,
									 const std::size_t num_nodes) :	parser(_parser),
									 								POP_SIZE(parser.getSizeT("POP_SIZE")),
																	MAX_GEN(parser.getSizeT("MAX_GEN")),
																	FIT_EQ(parser.getInt("FIT_EQ")),
																	INIT_MODE(parser.getInt("INIT_MODE")),
																	PERC_RND_WALK(parser.getDouble("PERC_RND_WALK")),
																	TOUR_SIZE(parser.getSizeT("TOUR_SIZE")),
																	NUM_FROM_TOUR(parser.getSizeT("NUM_FROM_TOUR")),
																	NUM_ELITE(parser.getSizeT("NUM_ELITE")),
																	PERC_CHILD_CROSS(parser.getDouble("PERC_CHILD_CROSS")),
																	PERC_NODE_CROSS(parser.getDouble("PERC_NODE_CROSS")),
																	PERC_CHILD_CROSS_AND_MUT(parser.getDouble("PERC_CHILD_CROSS_AND_MUT")),
																	SMALL_GRAGH_SIZE(parser.getSizeT("SMALL_GRAGH_SIZE")),
																	SMALL_GEN(parser.getSizeT("SMALL_GEN")),
																	MED_GRAPH_SIZE(parser.getSizeT("MED_GRAPH_SIZE")),
																	MED_GEN(parser.getSizeT("MED_GEN")),
																	MIGRATION_GEN(parser.getSizeT("MIGRATION_GEN")),
																	START_CROSS_TYPE(parser.getBool("START_CROSS_TYPE")),
																	CROSS_TYPE_CHANGE_PERC(parser.getDouble("CROSS_TYPE_CHANGE_PERC")),
																	STATUS_UPDATES(parser.getBool("STATUS_UPDATES")),
																	parent_gen(0),
																	generation(0),
																	min_num_clusters(1),
																	max_num_clusters(num_nodes),
																	cross_type(1),																	
																	using_migration(false),
																	max_gen_reached(false)

{
	fitness.resize(POP_SIZE);
	num_clusters.resize(POP_SIZE);
	chrom_index.resize(POP_SIZE);
	Chromosome new_chrom(parser);
	std::vector<Chromosome> pop_row;
	for(std::size_t i = 0 ; i < POP_SIZE; ++i)
		pop_row.push_back(new_chrom);
	pop.push_back(pop_row);
	pop.push_back(pop_row);	

	// Initialize values
	my_pop_stats.initialize(MAX_GEN + 1);
	mem.initialize(num_nodes);
}


ClusteringResults::~ClusteringResults()
{
}

void ClusteringResults::initializePopulation(Graph &g)
{
	// Add the nodes from the graph to all chromosomes
	for (std::size_t i_chrom = 0; i_chrom < POP_SIZE; i_chrom++)
	{
		pop[0][i_chrom].initialize(g);
		pop[1][i_chrom].initialize(g);
	}

	double max_fit, avg_fit;
	std::size_t max_clust, min_clust, avg_clust;
	int *node_ids, *local_labels, *label_count, *max_labels, **labels;
	bool *node_chosen;
	std::size_t pop_init_change_index;
	//int ver_res;

	// Initialize the chromosomes based on selected mode
	switch (INIT_MODE)
	{
	case (unweighted_split): // Split random clusters until the desired number of clusters is reached
		for (int i_chrom = 0; i_chrom < POP_SIZE; i_chrom++)
		{
			// Check the number of clustes in the chromosome
			if (pop[parent_gen][i_chrom].getNumClusters() > 1)
			{
				pop[parent_gen][i_chrom].reset(true, true); // Merge all nodes into one cluster and delete all extra clusters
			}

			pop[parent_gen][i_chrom].splitInit(g, mem, min_num_clusters, max_num_clusters);

			// Calculate the fitness
			pop[parent_gen][i_chrom].calculateFitness(g, FIT_EQ);
			fitness[i_chrom] = pop[parent_gen][i_chrom].getFitness();
			num_clusters[i_chrom] = pop[parent_gen][i_chrom].getNumClusters();
			chrom_index[i_chrom] = i_chrom;
		}
		break;

		/*
	case (weighted_split):// Split random clusters by weight until the desired number of clusters is reached
		for (int i_chrom = 0; i_chrom < POP_SIZE; i_chrom++)
		{
			// Check the number of clustes in the chromosome
			if (pop_[parent_gen_][i_chrom].GetNumClusters() > 1)
			{
				pop_[parent_gen_][i_chrom].Reset(true, true); // Merge all nodes into one cluster and delete all extra clusters
			}

			int num_clusters, split_index;

			// Pick number of clusters between min and max
			num_clusters = PickNumber(this->min_num_clusters_, this->max_num_clusters_);

			// Pick a random cluster and split it until the desired number of clustes is reached
			for (int i_split = 0; i_split < num_clusters; i_split++)
			{
				split_index = PickNumber(0, pop_[parent_gen_][i_chrom].GetLastClustIndex());
				if (pop_[parent_gen_][i_chrom].GetCluster(split_index)->GetNumNodes() > 1)
					pop_[parent_gen_][i_chrom].SplitClusterByWeight(split_index, the_graph);
			}

			// Set the initial number of cluster for chromosome and reset the number of splits
			pop_[parent_gen_][i_chrom].SetNumInitClusters(pop_[parent_gen_][i_chrom].GetNumSplits());
			pop_[parent_gen_][i_chrom].SetNumSplits(0);

			// Calculate the fitness
			pop_[parent_gen_][i_chrom].CalculateFitness(this->GetFitnessEquation(), the_graph);
			fitness_[i_chrom] = pop_[parent_gen_][i_chrom].GetFitness();
			num_clusters_[i_chrom] = pop_[parent_gen_][i_chrom].GetNumClusters();
			chrom_index_[i_chrom] = i_chrom;
		}
		break;
			*/

	/*
	case(label_propagation): // Label Propagation

		// Allocate memory
		node_ids = new int[the_graph->GetNumNodes()]();
		local_labels = new int[the_graph->GetHighestDegree()]();
		label_count = new int[the_graph->GetHighestDegree()]();
		max_labels = new int[the_graph->GetHighestDegree()]();
		node_chosen = new bool[the_graph->GetNumNodes()]();
		labels = new int*[2];
		for (int i = 0; i < 2; i++)
		{
			labels[i] = new int[the_graph->GetNumNodes()]();
		}

		for (int i_node = 0; i_node < the_graph->GetNumNodes(); i_node++)
		{
			node_ids[i_node] = i_node;
		}

		for (int i_chrom = 0; i_chrom < POP_SIZE; i_chrom++)
		{
			
			pop[!parent_gen][i_chrom].labelPropagation(labels,local_labels,label_count, max_labels, node_chosen, node_ids, the_graph);

			// Calculate the fitness
			pop_[parent_gen_][i_chrom].CalculateFitness(this->GetFitnessEquation(), the_graph);
			fitness_[i_chrom] = pop_[parent_gen_][i_chrom].GetFitness();
			num_clusters_[i_chrom] = pop_[parent_gen_][i_chrom].GetNumClusters();
			chrom_index_[i_chrom] = i_chrom;
		}

		delete[] node_ids;
		delete[] local_labels;
		delete[] label_count;
		delete[] max_labels;
		delete[] node_chosen;
		for (int i = 0; i < 2; i++)
		{
			delete[] labels[i];
		}
		delete[] labels;

		break;
	*/
	case (random_walk):
		for (std::size_t i_chrom = 0; i_chrom < POP_SIZE; ++i_chrom)
		{
			// Check the number of clustes in the chromosome
			if (pop[parent_gen][i_chrom].getNumClusters() > 1)
				pop[parent_gen][i_chrom].reset(true, true); // Merge all nodes into one cluster and delete all extra clusters
			
			pop[parent_gen][i_chrom].randomWalk(g);

			// Set the initial number of cluster for chromosome and reset the number of splits
			pop[parent_gen][i_chrom].setNumInitClusters(pop[parent_gen][i_chrom].getNumClusters());

			// Calculate the fitness
			pop[parent_gen][i_chrom].calculateFitness(g, FIT_EQ);
			fitness[i_chrom] = pop[parent_gen][i_chrom].getFitness();
			num_clusters[i_chrom] = pop[parent_gen][i_chrom].getNumClusters();
			chrom_index[i_chrom] = i_chrom;
		}
		break;

	case(combo):
		pop_init_change_index = (std::size_t)floor(POP_SIZE * PERC_RND_WALK);
		for (std::size_t i_chrom = 0; i_chrom < pop_init_change_index; ++i_chrom)
		{
			// Check the number of clustes in the chromosome
			if (pop[parent_gen][i_chrom].getNumClusters() > 1)
				pop[parent_gen][i_chrom].reset(true, true); // Merge all nodes into one cluster and delete all extra clusters
			
			pop[parent_gen][i_chrom].randomWalk(g);

			// Set the initial number of cluster for chromosome and reset the number of splits
			pop[parent_gen][i_chrom].setNumInitClusters(pop[parent_gen][i_chrom].getNumClusters());

			// Calculate the fitness
			pop[parent_gen][i_chrom].calculateFitness(g, FIT_EQ);
			fitness[i_chrom] = pop[parent_gen][i_chrom].getFitness();
			num_clusters[i_chrom] = pop[parent_gen][i_chrom].getNumClusters();
			chrom_index[i_chrom] = i_chrom;
		}
		
		for (std::size_t i_chrom = pop_init_change_index; i_chrom < POP_SIZE; ++i_chrom)
		{
			// Check the number of clustes in the chromosome
			if (pop[parent_gen][i_chrom].getNumClusters() > 1)
				pop[parent_gen][i_chrom].reset(true, true); // Merge all nodes into one cluster and delete all extra clusters				
			

			pop[parent_gen][i_chrom].splitInit(g, mem, min_num_clusters, max_num_clusters);

			// Calculate the fitness
			pop[parent_gen][i_chrom].calculateFitness(g, FIT_EQ);
			fitness[i_chrom] = pop[parent_gen][i_chrom].getFitness();
			num_clusters[i_chrom] = pop[parent_gen][i_chrom].getNumClusters();
			chrom_index[i_chrom] = i_chrom;
		}
		break;

	default:
		break;
	}

	// Calculate statistics
	CalculateMaxAndAvg(fitness, max_fit, avg_fit);
	CalculateMinMaxAvg(num_clusters,  min_clust, max_clust, avg_clust);
	// Save values 
	my_pop_stats.setFitness(0, max_fit, avg_fit);
	my_pop_stats.setClusters(0, min_clust, max_clust, avg_clust);
	// Sort fitness values
	InsertionSort(fitness, chrom_index);
}

void ClusteringResults::runGA_Core(Graph &g, const std::string &clust_filename)
{
	// Initialize parameters
	std::size_t cur_chrom, node_to_cross, clust_cross_mut_count, clust_to_cross_index,max_gen, final_gen, cross_offset;
	double best_fitness;
	double max_fit, avg_fit;
	std::size_t max_clust, min_clust, avg_clust, best_gen = 0;
	bool num_cluts_changed = false, stead_state = false;
	std::vector<std::size_t> tournament_pool(TOUR_SIZE);

	// Approximatly 80% of the children, after the elites are copied should be created by crossover. This value should be an even number.
	int NUM_CHROM_CROSS = (int)((POP_SIZE - NUM_ELITE)*PERC_CHILD_CROSS);
	// Check if the value is odd
	if (NUM_CHROM_CROSS % 2 == 1)
		NUM_CHROM_CROSS--;
	int NUM_NODE_TO_CROSS = (int)(g.getNumNodes()*PERC_NODE_CROSS);
	int NUM_CLUST_TO_CROSS = 1;
	// Calculate the number of child chromosome to create through crossover then mutation
	int NUM_CROSS_AND_MUT = (int)(NUM_CHROM_CROSS * PERC_CHILD_CROSS_AND_MUT);

	// Calculate the number of child chromosomes to create by mutation only
	int NUM_CHROM_MUT = POP_SIZE - NUM_ELITE - NUM_CHROM_CROSS;

	//RecordFitness(fitness_output,generation, max, avg);
		
	if (g.getNumNodes() <= SMALL_GRAGH_SIZE )
	{
		max_gen = SMALL_GEN;
		final_gen = SMALL_GEN;
	}
	else if (g.getNumNodes() <= MED_GRAPH_SIZE)
	{
		max_gen = MED_GEN;
		final_gen = MED_GEN;
	}
	else
	{
		max_gen = MAX_GEN;
		final_gen = MAX_GEN;
	}
	
	// Check if migration is being used
	if (using_migration)
	{
		if (generation < final_gen / 2)
		{
			max_gen = final_gen / 2;
		}
		// Check if migration limit will be reached before he max_gen
		else if ((generation + MIGRATION_GEN) < max_gen)
		{
			max_gen = generation + MIGRATION_GEN; // Stop evolution when number of migration generations has been rached
		}
		else
		{
			max_gen_reached = true;
		}
	}

	best_fitness = pop[parent_gen][chrom_index[0]].getFitness();
	
	// Loop through evolution code until desired number of generations is reached
	//while (generation < max_gen)

	// Loop through the evolutiond code until the best objective value reaches steady-state
	while(!stead_state)
	{
		// Increment the generation counter
		++generation;

		// Check if status updates should be reported
		if (STATUS_UPDATES)
		{
			// Check if current generation is multiple of 100. If so print
			if (generation % 100 == 0)
				printf("Starting generation %lu (%f)\n", generation, getSortedChromosome(parent_gen, 0)->getFitness());
			
		}

		// If the desired percentage of generations have been reached and the crossover type is the same as the starting value, change the type
		if ((getCrossType() == START_CROSS_TYPE) && ((double)generation / final_gen >= CROSS_TYPE_CHANGE_PERC))
		{
			setCrossType(!START_CROSS_TYPE);			
			printf("Changing to crossover type %i at generation %lu.\n", getCrossType(), generation);
		}

		// Reset chromosome index
		cur_chrom = 0;

		// Save top 'NUM_ELITE' parents in child generation
		for (cur_chrom; cur_chrom < NUM_ELITE; cur_chrom++)
		{
			// Copy top fitness parents to the next generation
			pop[parent_gen][chrom_index[cur_chrom]].copy(&pop[!parent_gen][cur_chrom]);

			// Update the fitness array and chrom_index array
			fitness[cur_chrom] = pop[!parent_gen][cur_chrom].getFitness();
			num_clusters[cur_chrom] = pop[!parent_gen][cur_chrom].getNumClusters();
			chrom_index[cur_chrom] = cur_chrom;
		}

		// Initialy populate the remaining children with the winners of tournaments		
		std::size_t i_chrom = NUM_ELITE; // Initialize this index to the number of elites chosen
		// Loop until all remaining chromosomes are chosen from tournaments
		while (i_chrom < POP_SIZE)
		{
			// Run a tournament (returns indicies in 'tournament_pool' based on sorting their fitness)
			RunTournament(fitness, tournament_pool);

			// Loop through the desired number of individuals to take from each tournament
			for (std::size_t i_winner = 0; i_winner < NUM_FROM_TOUR; ++i_winner)
			{
				// Copy the 'i_winner' chromosome to the 'i_chrom' spot in the next generation
				pop[parent_gen][tournament_pool[i_winner]].copy(&pop[!parent_gen][i_chrom]);

				++i_chrom; // Increment index

				// Check if population size has been reached
				if (i_chrom == POP_SIZE)
					break; // Break from loop
			}
		}

		// Create 20% of the remianing children by mutation
		for (cur_chrom; cur_chrom < NUM_ELITE + NUM_CHROM_MUT; cur_chrom++)
		{
			// Mutate each child
			pop[!parent_gen][cur_chrom].mutate(g, mem);

			// Delete all empty clusters and indicate if any clusters were deleted
			if(pop[!parent_gen][cur_chrom].deleteAllEmptyClusters())
				pop[!parent_gen][cur_chrom].updateAllClusterNum(); // Update cluster numbers afte delete
			
			// Check that the number of clusters in the chromosome is still within the allowed range
			if(pop[!parent_gen][cur_chrom].maintainNumClusters(g, mem, min_num_clusters, max_num_clusters))
			{
				// Delete all empty clusters and indicate if any clusters were deleted
				if (pop[!parent_gen][cur_chrom].deleteAllEmptyClusters())
					pop[!parent_gen][cur_chrom].updateAllClusterNum();
			}

			// Calculate the fitness of the chromosome
			pop[!parent_gen][cur_chrom].calculateFitness(g, FIT_EQ);
			// Add the fitness of the child to the fitnees array
			fitness[cur_chrom] = pop[!parent_gen][cur_chrom].getFitness();
			num_clusters[cur_chrom] = pop[!parent_gen][cur_chrom].getNumClusters();
			// Add in index of the child in the population to the index array
			chrom_index[cur_chrom] = cur_chrom;
		}

		// Reset counter indicating the number of children created through crossover then mutation
		clust_cross_mut_count = 0;

		// Create the remaining children by crossover
		for (cur_chrom; cur_chrom < POP_SIZE; cur_chrom++)
		{
			// Check if cur_chrom is odd
			if (cur_chrom % 2 == 1)
				cross_offset = -1;
			else
				cross_offset = 1;

			// Create children via crossover based on type of crossover selected
			if (cross_type == 0) // NODE_CROSS
			{
				// Loop through the desired number of nodes
				for (std::size_t i_node = 0; i_node < NUM_NODE_TO_CROSS; ++i_node)
				{
					// Pick a node at random
					node_to_cross = PickSizeT(0, pop[!parent_gen][cur_chrom].getNumNodes()-1);					
					pop[!parent_gen][cur_chrom].nodeCrossOver(g, &pop[!parent_gen][cur_chrom + cross_offset], node_to_cross);
				}

				// Save the number of cluster before fixed unconnected parts
				std::size_t num_clust_prior = pop[!parent_gen][cur_chrom].getNumClusters();

				// Get a pointer to the first cluster
				Cluster *cur_clust = pop[!parent_gen][cur_chrom].getFirstCluster();

				// Loop through all clusters
				for (std::size_t i_clust = 0; i_clust < num_clust_prior; ++i_clust)
				{
					// Check if the chromosome needs to be verified
					if (cur_clust->needClusterFixed())
						pop[!parent_gen][cur_chrom].fixUnconnectedClusters(g, i_clust);
					
					cur_clust = cur_clust->getNextCluster(); // Go to next cluster
				}
			}
			else
			{
				for (std::size_t i_clust = 0; i_clust < NUM_CLUST_TO_CROSS; ++i_clust)
				{
					// Pick a cluster at random from the other parent
					clust_to_cross_index = PickSizeT(0, pop[!parent_gen][cur_chrom+cross_offset].getLastClustIndex());

					// Cross over selected cluster
					pop[!parent_gen][cur_chrom].clusterCrossOver(g, &pop[!parent_gen][cur_chrom+cross_offset], clust_to_cross_index);
				}
			}

			// Check how many chromosomes have been created this generation through crossover and mutatution
			if (clust_cross_mut_count <= NUM_CROSS_AND_MUT)
			{
				// Mutate the chromosome
				pop[!parent_gen][cur_chrom].mutate(g, mem);
				++clust_cross_mut_count;
			}

			// Delete all empty clusters and indicate if any clusters were deleted
			if (pop[!parent_gen][cur_chrom].deleteAllEmptyClusters())
				pop[!parent_gen][cur_chrom].updateAllClusterNum();

			// Check that the number of clusters in the chromosome is still within the allowed range
			if (pop[!parent_gen][cur_chrom].maintainNumClusters(g, mem, min_num_clusters, max_num_clusters))
			{
				// Delete all empty clusters and indicate if any clusters were deleted
				if (pop[!parent_gen][cur_chrom].deleteAllEmptyClusters())
					pop[!parent_gen][cur_chrom].updateAllClusterNum();
			}

			// Calculate the fitness of the chromosome
			pop[!parent_gen][cur_chrom].calculateFitness(g, FIT_EQ);
			// Add the fitness of the child to the fitnees array
			fitness[cur_chrom] = pop[!parent_gen][cur_chrom].getFitness();
			num_clusters[cur_chrom] = pop[!parent_gen][cur_chrom].getNumClusters();
			// Add in index of the child in the population to the index array
			chrom_index[cur_chrom] = cur_chrom;
		}

		CalculateMaxAndAvg(fitness, max_fit, avg_fit);
		CalculateMinMaxAvg(num_clusters, min_clust, max_clust, avg_clust);
		// Save values 
		my_pop_stats.setFitness(generation, max_fit, avg_fit);
		my_pop_stats.setClusters(generation, min_clust, max_clust, avg_clust);

		InsertionSort(fitness, chrom_index);

		// Check if the fitness has improved
		Chromosome *tmp = getSortedChromosome(!parent_gen, 0);
		double tmp_fit = tmp->getFitness();
		if (tmp_fit > best_fitness + FITNESS_TOL)
		{
			best_fitness = tmp_fit;
			tmp->setGenerationCreated(generation);
			best_gen = generation;
		}

		if (generation - best_gen > STEADY_STATE_DURATION)
		{
			stead_state = true;
		}

		// Change the parent generation
		parent_gen = !parent_gen;
	}	
}

Chromosome * ClusteringResults::getSortedChromosome(const bool gen, const std::size_t index)
{
	assert(gen < 2);
	assert(index < pop[gen].size());
	return &this->pop[gen][chrom_index[index]];
}

bool ClusteringResults::readPopulation(Graph &g, const std::string &pop_filename)
{
	// Declare local variables
	FILE* input;	
	
	// Try to open file
	if ((input = fopen(pop_filename.c_str(), "r")) == NULL)
	{
		printf("ERROR in ClusteringResults::ReadPopulation - Could not open output file(%s).\n", pop_filename.c_str());
		exit(1);
	}

	// Read in number of generations
	if (fscanf(input, "%lu\n", &generation) == 0)
	{
		printf("ERROR in ClusteringResults::ReadPopulation - Population file not formatted correctly.\n");		
		return false;
	}

	// Loop through all individuals in population
	for (int i_chrom = 0; i_chrom < POP_SIZE; i_chrom++)
	{
		pop[parent_gen][i_chrom].readChromosome(input, g, FIT_EQ);
		
		fitness[i_chrom] = pop[parent_gen][i_chrom].getFitness();
		num_clusters[i_chrom] = pop[parent_gen][i_chrom].getNumClusters();
		chrom_index[i_chrom] = i_chrom;

		// Add the nodes from the graph to the non-parent chromosome
		pop[!parent_gen][i_chrom].initialize(g);
	}
	
	// Close file
	fclose(input);	
	return true;
}

bool ClusteringResults::recordPopulation(const std::string &pop_filename)
{
	// Declare local variables
	FILE* output;

	// Try to open file
	if ((output = fopen(pop_filename.c_str(), "w")) == NULL)
	{
		printf("ERROR in ClusteringResults::RecordPopulation - Could not open output file(%s).\n", pop_filename.c_str());
		exit(1);
	}
	// Record the generation
	fprintf(output, "%lu\n", generation);
	
	// Save the number of nodes
	std::size_t num_nodes = pop[0][0].getNumNodes();
	
	// Loop through all individuals in population
	for (std::size_t i = 0; i < chrom_index.size(); ++i)
	{
		// Record non-cluster number info
		fprintf(output,"%lu ", pop[parent_gen][chrom_index[i]].getNumInitClusters());
		fprintf(output,"%lu ", pop[parent_gen][chrom_index[i]].getNumSplits());
		fprintf(output,"%lu ", pop[parent_gen][chrom_index[i]].getNumMerges());
		fprintf(output,"%lu ", pop[parent_gen][chrom_index[i]].getNumRedis());
		fprintf(output,"%lu ", pop[parent_gen][chrom_index[i]].getNumPointMutations());
		fprintf(output,"%lu ", pop[parent_gen][chrom_index[i]].getNumNodeCrossovers());
		fprintf(output,"%lu ", pop[parent_gen][chrom_index[i]].getNumClusterCrossovers());
		fprintf(output,"%lu ", pop[parent_gen][chrom_index[i]].getCreatedGeneration());		
		fprintf(output,"%i ", pop[parent_gen][chrom_index[i]].getInitType());
		fprintf(output,"%lu ", pop[parent_gen][chrom_index[i]].getNumClusters());

		// Record the cluster number for all nodes, except the last
		for (std::size_t i_node = 0; i_node < num_nodes-1; ++i_node)
			fprintf(output,"%lu ", pop[parent_gen][chrom_index[i]].getClusterNum(i_node));
		
		// Record the cluster number for the last node in the chromosome
		fprintf(output,"%lu\n", pop[parent_gen][chrom_index[i]].getClusterNum(num_nodes-1));
	}

	// Close file
	fclose(output);

	return true;
}

void ClusteringResults::readMigrant(Graph &g, const std::string &mig_filename)
{
	// Declare local variables
	FILE* input;		
	bool chrom_read = false;

	// Try to open file
	if ((input = fopen(mig_filename.c_str(), "r")) == NULL)
	{
		printf("ERROR in ClusteringResults::ReadMigrant - Could not open output file(%s).\n", mig_filename.c_str());
		exit(1);
	}
	
	// Read and then discard population generation
	int tmp=0;
	if (fscanf(input, "%i", &tmp) == 0)
	{
		printf("ERROR in ClusteringResults::ReadMigrant - Population file not formatted correctly.\n");		
		exit(1);
	}
	
	Chromosome tmp_chrom(parser);
	tmp_chrom.readChromosome(input, g, FIT_EQ);

	// Pick index of individual to overwrite with migrant
	if (NUM_ELITE >= POP_SIZE - 1)
	{
		printf("ERROR in ClusteringResults::ReadMigrant - The number of elite chromosomes in population is greater than or equal to the population size.\n");		
		exit(1);
	}
	std::size_t ind_to_over = PickSizeT(NUM_ELITE, POP_SIZE - 1); // Do not overwrite the elite chromosome
	
	// Copy new chromosome to population
	tmp_chrom.copy(&pop[parent_gen][ind_to_over]);

	fitness[ind_to_over] = pop[parent_gen][ind_to_over].getFitness();
	num_clusters[ind_to_over] = pop[parent_gen][ind_to_over].getNumClusters();
	chrom_index[ind_to_over] = ind_to_over;

	// Check if migrated chromosome is now the elite
	if (tmp_chrom.getFitness() > pop[0][0].getFitness())
	{
		// TEST
		double diff = tmp_chrom.getFitness() - pop[0][0].getFitness();
		printf("Migrant has higher fitness than elite. Diff: %f.\n", diff);
		InsertionSort(fitness, chrom_index);
	}
}

void ClusteringResults::setMinNumClusters(const std::size_t min_num_clust)
{
	min_num_clusters = min_num_clust;
}

void ClusteringResults::setMaxNumClusters(const std::size_t max_num_clust)
{
	max_num_clusters = max_num_clust;
}

void ClusteringResults::setCrossType(const bool type)
{
	cross_type = type;
}

void ClusteringResults::setUsingMigration(const bool using_mig)
{
	using_migration = using_mig;
}

void ClusteringResults::setMaxGenReached(const bool _max_gen_reached)
{
	max_gen_reached = _max_gen_reached;
}

std::size_t ClusteringResults::getMinNumClusters() const
{
	return min_num_clusters;
}

std::size_t ClusteringResults::getMaxNumCusters() const
{
	return max_num_clusters;
}

bool ClusteringResults::getCrossType() const
{
	return cross_type;
}

bool ClusteringResults::getUsingMigration() const
{
	return using_migration;
}

bool ClusteringResults::getMaxGenReached() const
{
	return max_gen_reached;
}

std::size_t ClusteringResults::getSmallGraphSize() const
{
	return SMALL_GRAGH_SIZE;
}