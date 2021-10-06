#include <string>
#include <sstream>
#include <vector>

std::size_t ReadClusterOutput(std::string &, std::vector<int>&, int &, int &);
std::size_t ReadNodeNumber(std::string &, std::vector<std::size_t>&);
void ConvertClusterAssignment(std::vector<int> &, std::vector<std::size_t>&, const std::size_t, const int, std::vector<int>&);
void RecordResults(std::string &, std::vector<int>&, const int, const int);

int main(int argc, char **argv)
{
    // Check the input format
	if (argc < 6)
	{
		printf("ERROR: Incorrect Inputs.\n");
		printf("Usage CombClustOut output_dir output_tag num_comps total_num_nodes obj_name\n");
		exit(1);
	}

    // Save inputs to local variables
	std::string out_dir(argv[1]);
	std::string out_tag(argv[2]);
    int num_comps = atoi(argv[3]);
    int total_num_nodes = atoi(argv[4]);
    std::string obj_name(argv[5]);
	
    std::vector<int> clust(total_num_nodes, -1);
    std::vector<int> clust_comp(total_num_nodes);
    std::vector<std::size_t> nn_comp(total_num_nodes);
    std::size_t num_out, num_nn;
    int num_clusts, num_edges, total_num_clust = 0, total_num_edges = 0;
    std::size_t num_comps_split = 0;

    // Loop through all components   
    for(std::size_t comp = 0; comp <= num_comps; ++comp)
    {
        std::string comp_filename = out_dir + out_tag + "_comp" + std::to_string(comp); 
        std::string out_filename = comp_filename + "_" + obj_name + ".out";
        std::string nn_filename = comp_filename + ".nn";

        // Read .out file
        num_out = ReadClusterOutput(out_filename, clust_comp, num_clusts, num_edges);

        if(num_clusts > 1)
        {
            printf("Comp %lu split\n", comp);
            ++num_comps_split;
        }

        // Read .nn file
        num_nn = ReadNodeNumber(nn_filename, nn_comp);

        // Check that both function read in the same number of nodes
        if(num_out != num_nn)
        {
                printf("ERROR - Comp %lu .out and .nn files did not have the same number of elements\n", comp);
                printf("Num out: %lu\n", num_out);
                printf("Num nn: %lu\n", num_nn);
                exit(1);
        }

        // Set .out results based on .nn data
        ConvertClusterAssignment(clust_comp, nn_comp, num_out, total_num_clust, clust);
        total_num_clust += num_clusts;
        total_num_edges += num_edges;

        //for(std::size_t i = 0; i < total_num_nodes; ++i)
        //    printf("%i ",clust[i]);
        //printf("\n");
    }
    printf("Number of components split: %lu\n", num_comps_split);

    // Record Results
    std::string output_file = out_dir + out_tag + "_" + obj_name +".out";
    RecordResults(output_file, clust, total_num_clust, total_num_edges);

    return 0;
}

std::size_t ReadClusterOutput(std::string &clust_file, std::vector<int> &clust, int &num_clusters, int &num_edges)
{
    FILE* input;
	char tmp_str[50];
    std::size_t i = 0;
	tmp_str[0] = '\0';
	
	if((input = fopen(clust_file.c_str(), "r")) == NULL)
	{
		printf("ERROR - Cound not open clust.out file (%s)\n.", clust_file.c_str());
		exit(1);
	}
    
    fscanf(input, "%s", tmp_str); // number of nodes
    fscanf(input, "%s", tmp_str); // string
    fscanf(input, "%s", tmp_str); // number of clusters
    num_clusters = atoi(tmp_str);
    fscanf(input, "%s", tmp_str); // string
    fscanf(input, "%s", tmp_str); // number of edges
    num_edges = atoi(tmp_str);
    fscanf(input, "%s", tmp_str); // strings

    while(1)
	{
        // Check if end-of-file was reached
		if (feof(input))
		{
			break; // Break while loop
		}

		if(fscanf(input, "%s", tmp_str) > 0)
        {
            clust[i] = atoi(tmp_str);
            ++i;
        }
	}

	fclose(input);
    return i;
}


std::size_t ReadNodeNumber(std::string & nn_file, std::vector<std::size_t> & node_num)
{
    FILE* input;
	char tmp_str[50];
    std::size_t i = 0;
	tmp_str[0] = '\0';
	
	if((input = fopen(nn_file.c_str(), "r")) == NULL)
	{
		printf("ERROR - Cound not open nn file.");
		exit(1);
	}

    while(1)
	{
        // Check if end-of-file was reached
		if (feof(input))
		{
			break; // Break while loop
		}

		if(fscanf(input, "%s", tmp_str) > 0)
        {
            node_num[i] = atoi(tmp_str);
            ++i;
        }
	}

	fclose(input);
    return i;
}

void ConvertClusterAssignment(std::vector<int> &clust_comp, std::vector<std::size_t>&nn_comp,
                             const std::size_t num_nodes, const int num_existing_clust,
                             std::vector<int> &clust_out)
{
    for(std::size_t i = 0; i < num_nodes; ++i)
    {
        if( clust_comp[i] == -1)
            clust_out[nn_comp[i]] = -1;
        else
            clust_out[nn_comp[i]] = clust_comp[i] + num_existing_clust;
	}    
}

void RecordResults(std::string &output_file, std::vector<int>& clust, const int num_clusts, const int num_edges)
{
    FILE* output;
    if((output = fopen(output_file.c_str(), "w")) == NULL)
    {
        printf("ERROR - Could not open file (%s)\n", output_file.c_str());        
    }

    fprintf(output, "%lu nodes %i clusters %i edges\n", clust.size(), num_clusts, num_edges);
    for(std::size_t i = 0; i < clust.size(); ++i)
        fprintf(output, "%i ", clust[i]);
    fclose(output);
}