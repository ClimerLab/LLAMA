#include "Graph.h"
#include "Chromosome.h"

#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include "assert.h"

void readNodeNumber(const std::string &, std::vector<std::size_t> &);
void readOutput(const std::string &, std::vector<int> &, std::size_t &, std::size_t &);
void updateClusterAssingment(std::vector<int> &, std::vector<int> &, std::vector<std::size_t> &, const std::size_t);
double readFitness(const std::string &, const std::size_t);
void convertToClusterNumber(std::vector<int> &, std::size_t &, std::vector<std::size_t> &);

int main(int argc, char **argv)
{
	// Check the input format
	if (argc != 7)
	{
		printf("ERROR: Incorrect Inputs.\n");
		printf("Usage: CombComp output_dir output_tag graph_filename singleton_filename nn_dir num_comp\n");
		exit(1);
	}

	// Save inputs to local variables
	std::string out_dir(argv[1]);
	std::string out_tag(argv[2]);
	std::string graph_filename(argv[3]);
	std::string singleton_filename(argv[4]);
	std::string nn_dir(argv[5]);
	int num_comp = atoi(argv[6]);
	
	std::string out_name, nn_name, fitness_filename;
	Chromosome my_chrom;
	Graph my_graph;
	std::vector<int> clust_assignment, local_clust;
	std::vector<std::size_t> nn;
	std::size_t num_local_clusts, num_local_edges, num_clusts = 0, num_edges = 0;
	
	printf("\nCombining Components...\n");
	my_graph.readFile(graph_filename); // Read original graph

	clust_assignment.resize(my_graph.getNumNodes());

	// Read in any singletons
	readNodeNumber(singleton_filename, nn);
	local_clust.resize(nn.size(), -1);
	updateClusterAssingment(clust_assignment, local_clust, nn, num_clusts);

	// Loop through all the components
	for (int i_comp = 1; i_comp <= num_comp; i_comp++)
	{
		if (i_comp % 1000 == 0)
			printf("Comp %i\n",i_comp);
		
		// Update file name based on component and island
		out_name = out_dir + out_tag + "_comp" + std::to_string(i_comp) + "_best.out";
		nn_name = nn_dir + out_tag + "_comp" + std::to_string(i_comp) + ".nn";
		
		readNodeNumber(nn_name, nn);
		readOutput(out_name, local_clust, num_local_clusts, num_local_edges);
		updateClusterAssingment(clust_assignment, local_clust, nn, num_clusts);
		num_clusts += num_local_clusts;
		num_edges += num_local_edges;
	}

	// Read in fitness values
	fitness_filename = out_dir + "objective.csv";
	double fitness = readFitness(fitness_filename, num_comp);
	printf("Fitness from reading: %lf\n", fitness);

	std::vector<std::size_t> cluster_num;
	convertToClusterNumber(clust_assignment, num_clusts, cluster_num);

	my_chrom.initialize(my_graph, cluster_num);
	my_chrom.setFitness(fitness);
	
	// Record the sorted chromosome
	printf("Recording Files\n");
	std::string combined_filename = out_dir + out_tag + "_best.out";
	my_chrom.recordOutput(combined_filename, num_edges);

	combined_filename = out_dir + out_tag + "_final.clust";
	my_chrom.recordToFile(combined_filename, "Final Clusters");
	
	printf("complete\n");
	return 0;
}

void readNodeNumber(const std::string &input_file, std::vector<std::size_t> &nn)
{
	FILE* input;
	std::size_t tmp;
	nn.resize(0);

	if((input = fopen(input_file.c_str(), "r")) == NULL)
	{
		printf("ERROR in CombineComponents - Could not open input file (%s)\n", input_file.c_str());
		exit(1);
	}

	// Read file
	while(true)
	{
		if(fscanf(input, "%lu", &tmp) > 0)
			nn.push_back(tmp);

		if(feof(input))
			break;		
	}

	fclose(input);
}

void readOutput(const std::string &input_file, std::vector<int> &clust_assign, std::size_t &num_clusters, std::size_t &num_edges)
{
	FILE *input;
	std::size_t num_nodes;
	int tmp;
	if((input = fopen(input_file.c_str(), "r")) == NULL)
	{
		printf("ERROR in CombineComponents - Could not open .out file (%s)\n", input_file.c_str());
		exit(1);
	}

	fscanf(input, "%lu", &num_nodes);
	fscanf(input, "%*s");
	fscanf(input, "%lu", &num_clusters);
	fscanf(input, "%*s");
	fscanf(input, "%lu", &num_edges);
	fscanf(input, "%*s");

	clust_assign.resize(num_nodes);

	for(std::size_t i = 0; i < num_nodes; ++i)
	{
		if(fscanf(input, "%d", &tmp) < 1)
		{
			printf("ERROR in CombineComponents - Could not read cluster assignment in (%s)\n", input_file.c_str());
			exit(1);
		}

		clust_assign[i] = tmp;
	}

	fclose(input);
}

void updateClusterAssingment(std::vector<int> &cluster_assignment, std::vector<int> &local_clust, std::vector<std::size_t> &nn, const std::size_t num_clusts)
{
	assert(local_clust.size() == nn.size());

	// Loop through all nodes in local_clust
	for(std::size_t i = 0; i < local_clust.size(); ++i)
	{
		if(local_clust[i] == -1)
			cluster_assignment[nn[i]] = -1;
		else
			cluster_assignment[nn[i]] = local_clust[i] + num_clusts;
	}		
}

double readFitness(const std::string &input_file, const std::size_t num_comps)
{
	double fitness = 0, tmp_fit;
	std::size_t tmp_cmp;
	FILE *input;

	if((input = fopen(input_file.c_str(), "r")) == NULL)
	{
		printf("ERROR in CombineComponents - Could not open file (%s)\n", input_file.c_str());
		exit(1);
	}

	for(std::size_t i = 1; i <= num_comps; ++i)
	{
		if(fscanf(input, "%lu,%lf", &tmp_cmp, &tmp_fit) < 1)
		{
			printf("ERROR in CombineComponents - Fitness file not formatted correctly.\n");
			exit(1);
		}

		if(i != tmp_cmp)
		{
			printf("ERROR in CombineComponents - Unexpected component number.\n");
			exit(1);
		}

		fitness += tmp_fit;
	}

	fclose(input);

	return fitness;
}

void convertToClusterNumber(std::vector<int> &clust_assign, std::size_t &num_clust, std::vector<std::size_t> &clust_num)
{
	clust_num.resize(clust_assign.size());

	for(std::size_t i = 0; i < clust_assign.size(); ++i)
	{
		if(clust_assign[i] == -1)
		{
			clust_num[i] = num_clust;
			num_clust++;
		}
		else
			clust_num[i] = (std::size_t)clust_assign[i] - 1;
	}
}