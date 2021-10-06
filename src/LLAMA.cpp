#include "LLAMA.h"
#include "Graph.h"
#include "ClusteringResults.h"
#include "IslandStat.h"
#include "UtilityFunction.h"
#include "Chromosome.h"
#include <sstream>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <stdio.h>

int main(int argc, char **argv)
{
	// Check inputs
	if (argc < 9)
	{
		printf("ERROR: Incorrect Inputs.\n");
		printf("Usage GA_Clustering input.gml cfg_file total_num_nodes total_num_edges output_dir output_tag island_num stats_dir total_num_isl(opt)\n");
		exit(1);
	}

	// Initialize parameters
	Graph my_graph;
	clock_t gml_start, gml_finish, ga_start, ga_mid, ga_finish;
	Chromosome* best_chrom;
	IslandStat my_island_stat;
	bool input_file_read, nn_file_read, pop_file_exists = false;
	double running_fitness = 0;
	double best_fitness_comp = 0;
	std::ostringstream ss, island_num_txt;
	std::string input_file, node_num_file, cfg_file, output_dir, output_tag, pop_file, stats_dir, island_dir, cur_isl_txt, obj_file;
	int island_num, total_num_island = 0, total_num_nodes, total_num_edges;
	std::vector<std::size_t> nn_vec;

	// Read in user input to local variables
	input_file = argv[1];
	cfg_file = argv[2];
	total_num_nodes = atoi(argv[3]);
	total_num_edges = atoi(argv[4]);
	output_dir = argv[5];
	output_tag = argv[6];
	cur_isl_txt = argv[7];
	island_num = atoi(argv[7]);
	stats_dir = argv[8];
			
	// Check of optional inputs were passed
	if (argc == 10)
	{
		//pop_file = argv[7]; // Read file name for population
		total_num_island = atoi(argv[9]);
	}
	// Check some of the user settings
	ConfigParser parser(cfg_file);
	bool error_in_settings = checkUserSettings(parser);	

	island_dir = output_dir + "I_" + cur_isl_txt + "/";
	pop_file = island_dir + "popfile.txt";
	obj_file = island_dir + "objective.csv";
	
	// Create output files based on user input
	island_dir += output_tag;
	std::string comp_num("0");

	// Find location of last "_comp" in the string
	std::size_t loc = input_file.rfind("_comp");
	std::size_t end_loc;

	// Check that "_comp" was found
	if (loc != std::string::npos)
	{
		// Find position of ".gml" in string
		end_loc = input_file.find(".gml");

		// Check if valid position was returned
		if (end_loc == std::string::npos)
		{
			// Find position of ".graph" in string
			end_loc = input_file.find(".graph");

			// Check if valid position was returned
			if (end_loc == std::string::npos)
			{
				// Print error and exit if input file type not valid
				printf("ERROR in GA_Clustering - input file extension not recongnized.\n");
				exit(1);
			}
		}

		// Determine component number from file name
		comp_num = input_file.substr(loc + 5, end_loc-(loc+5));
	}

	std::string graph_output = island_dir + "_graph.txt";
	std::string cluster_output = island_dir +"_";
	std::string fitness_output = island_dir + "_fitness_log_";
	std::string stat_output = stats_dir + "stats";
	std::string comp_cluster_output, comp_fitness_output, comp_stat_output, temp_cluster_output, tmp_clust_assign_output, temp_fitness_output;

	int batch_mode = 0;

	// Read input in graph
	gml_start = clock();
	my_graph.readFile(input_file);
	my_graph.setTotalNumNodes(total_num_nodes);
	my_graph.setTotalNumEdges(total_num_edges);

	if(parser.getBool("STATUS_UPDATES"))
	{
		printf("Input file read in.\n");
	}
	gml_finish = clock();

	// Allocate memory for the population based on the number of nodes in the graph
	ClusteringResults results(parser, my_graph.getNumNodes());	

	// Check if component is a singleton
	if (my_graph.getNumNodes() == 1)
	{
		printf("Component is a singleton.\n");
		return 0;
	}
	//else if (my_graph.getNumNodes() <= results.getSmallGraphSize())
	//{
		// Enumerate solutions

	//	return 0;
	//}
	

	// Initialize seed		
	srand(int(std::time(0)+island_num));

	 // Create output file names based on components
	comp_cluster_output = cluster_output + "comp" + comp_num + "_";
	comp_fitness_output = fitness_output + "comp" + comp_num + "_";
	comp_stat_output = stat_output + "_comp" + comp_num + ".csv";	

	int MIN_CLUSTERS, MAX_CLUSTERS;// , step_size;
	MIN_CLUSTERS = 1;
	switch (parser.getInt("MAX_CLUST_SIZE"))
	{
	case (0):
		MAX_CLUSTERS = my_graph.getNumNodes();
		break;
	case (1):
		MAX_CLUSTERS = (int)ceil(sqrt(my_graph.getNumNodes()));
		break;
	default:
		MAX_CLUSTERS = my_graph.getNumNodes();
		break;
	}
			
	//step_size = (MAX_CLUSTERS - MIN_CLUSTERS) / num_islands;

	// Call clustering core for the desired number of islands
	//for (int i_run = 0; i_run < num_islands; i_run++)
	//{
	printf("Starting island %i.\n", island_num);
	int i_run = island_num - 1;
	best_fitness_comp = -1;

	switch (batch_mode)
	{
	case (0):
		my_island_stat.updateMinMax(MIN_CLUSTERS, MAX_CLUSTERS);
		break;
	/*case(1):
		my_island_stat.SetMinNumClust(i_run, MIN_CLUSTERS + i_run * step_size);

		if (i_run < num_islands - 1)
		{
			my_island_stat.SetMaxNumClust(i_run, MIN_CLUSTERS + step_size * (i_run + 2));
		}
		else
		{
			my_island_stat.SetMaxNumClust(i_run, MAX_CLUSTERS);
		}
		break;
	*/
	default:
		break;
	}

	results.setMinNumClusters(my_island_stat.getMinNumClust());
	results.setMaxNumClusters(my_island_stat.getMaxNumClust());

	// Start the clock for this run
	ga_start = clock();

	// Add run information to output file names
	ss.str(std::string()); // Clear stream contents
	ss << i_run + 1; // Inster run number, plus 1 to stream
	temp_cluster_output = comp_cluster_output + ss.str() + ".clust";
	tmp_clust_assign_output = comp_cluster_output + ss.str() + ".out";
	temp_fitness_output = comp_fitness_output + ss.str() + ".csv";

	// Check if population file was passed to program
	if (argc == 10)
	{
		results.setUsingMigration(true);
		printf("Setting Migration\n");

		//Try to read population file
		FILE* pop;
		pop = fopen(pop_file.c_str(), "r");
		// Check if file opened
		if (pop != NULL)
		{
			pop_file_exists = true;
			fclose(pop);
		}
	}
	
	// Check if miragtion is being used
	if(pop_file_exists)
	{
		// Read in population file
		results.readPopulation(my_graph, pop_file);
				
		int island_to_mig = island_num;
		std::string mig_filename;
		
		// Choose an island at random to migrate information from
		while(island_to_mig == island_num)
			island_to_mig = PickNumber(1, total_num_island);

		island_num_txt.str(std::string()); // Clear stream contents
		island_num_txt << island_to_mig; // Inster run number, plus 1 to stream

		mig_filename = output_dir + "I_" + island_num_txt.str() + "/popfile.txt";

		// Read in the top individual from this island
		results.readMigrant(my_graph, mig_filename);
	}
	else
	{
		// Initialize Population
		if (parser.getBool("STATUS_UPDATES"))
			printf("Initializing Population...");
		results.initializePopulation(my_graph);

		std::string init_pop_file = island_dir + "_init_popfile.txt";		
		results.recordPopulation(init_pop_file);
	}

	// Check if status updates are being given
	if (parser.getBool("STATUS_UPDATES"))
		printf("complete\n");
		
	ga_mid = clock();

	// Initialize crossover type
	results.setCrossType(parser.getBool("START_CROSS_TYPE"));

	// Run the Cluster GA
	results.runGA_Core(my_graph, temp_cluster_output);
	
	// Get time after function has completed
	ga_finish = clock();

	// Check if migration is being used, but the max generations has not been reached
	if (results.getUsingMigration() && (!results.getMaxGenReached()))
	{
		// Record population
		results.recordPopulation(pop_file);
	}
	else
	{
		// Record the chromosome with the highest fitness
		printf("Getting best chrome\n");
		best_chrom = results.getSortedChromosome(results.parent_gen, 0);
		//best_chrom->recordToFile(temp_cluster_output, "Chromosome with the highest fitness after sort");
		printf("Recording Objective\n");
		best_chrom->recordObjective(obj_file, comp_num);
		printf("Recording Cluster Assignment\n");
		best_chrom->recordOutput(tmp_clust_assign_output, my_graph.getNumEdges());

		// TEST
		//printf("File name: %s", temp_cluster_output.c_str());
		
		// KEEP
		//results.my_pop_stats_.RecordStats(temp_fitness_output, 0, MAX_GEN, true);

		// Save the times and best chromosome information for this run
		my_island_stat.updateTimes(ga_start, ga_mid, ga_finish);
		my_island_stat.updateFromChromosome(best_chrom);

		// Record statistics
		bool print_headers = true;

		// Try to open the file (file has to exist)
		FILE* test_stat = fopen(comp_stat_output.c_str(), "r");

		// Check if file was opened
		if (test_stat != NULL)
		{
			print_headers = false;
			fclose(test_stat);
		}

		my_island_stat.recordStats(comp_stat_output, island_num, print_headers);
	}
		
	return 0;
}

bool checkUserSettings(ConfigParser &parser)
{
	// Initialize output
	bool error = false;

	// Check if the number of member taken from a tournament is larger than the tournament size
	if (parser.getSizeT("NUM_FROM_TOUR") > parser.getSizeT("TOUR_SIZE"))
	{
		printf("ERROR in Settings.h - NUM_FROM_TOUR cannot be larger than TOUR_SIZE.\n");
		error = true;
	}

	if ((parser.getInt("MAX_CLUST_SIZE") != 0) && (parser.getInt("MAX_CLUST_SIZE") != 1))
	{
		printf("ERROR in Settings.h - Unkown value of MAX_CLUST_SIZE.\n");
		error = true;
	}

	return error;
}

void readNN_File(const std::string &file_name, std::vector<std::size_t> &nn_vec)
{
	// Declare local variables
	FILE *input;
	std::size_t tmp, i_nn = 0;

	// Try to open file
	if((input = fopen(file_name.c_str(), "r")) == NULL)
	{
		printf("ERROR in ReadNN_File - Could not open file (%s).\n", file_name.c_str());
		exit(1);
	}
	
	// Read in values and store them in array until EOF
	while (true)
	{
		// Check if end-of-file was reached
		if (feof(input))
			break; // Break from while loop		

		if (fscanf(input, "%lu\n", &tmp) > 0) // Read in next string
		{			
			// Check if less nodes have been read in than are allocated
			if (i_nn < nn_vec.size())
			{
				nn_vec[i_nn] = tmp; // Convert string to int
				++i_nn; // Increment the index
			}
		}
	}

	// Check if the expected number of node number were read in
	if (i_nn != nn_vec.size())
	{
		printf("ERROR in ReadNN_File - Unexpected number of nodes found in file.\n");
		printf("\tExpcted %lu nodes, but found %lu.\n", nn_vec.size(), i_nn);
		exit(1);
	}	
}