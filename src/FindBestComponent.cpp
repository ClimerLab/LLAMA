#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>

#if defined(_WIN32)
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

int main(int argc, char** argv)
{
	// Check the input format
	if (argc < 5)
	{
		printf("ERROR: Incorrect Inputs.\n");
		printf("Usage FindBestComp out_dir out_tag num_comp num_islands \n");
		exit(1);
	}

	// Save inputs to local variables
	std::string out_dir(argv[1]);
	std::string out_tag(argv[2]);
	int num_comp = atoi(argv[3]);
	int num_islands = atoi(argv[4]);

	std::ostringstream ss_isl, ss_i_comp;
	std::string obj_file_name, best_file_name="",new_file_name;
	double tmp_fit;
	int comp_num;
	FILE* obj_file;
	std::vector<double> best_fitness(num_comp, -1);
	std::vector<int> best_island(num_comp, 1);
		
	#if defined(_WIN32)
	std::string new_dir = out_dir + "Best\\";
	//dir_res=_mkdir(new_dir.c_str()); // Create new directory
	#else
	std::string new_dir = out_dir + "Best/";
	//dir_res=mkdir(new_dir.c_str(), 0777);
	#endif
	
	// Initialize fitness values for island 1
	obj_file_name = out_dir + "I_1/objective.csv";
			
	// Try to open file
	if((obj_file = fopen(obj_file_name.c_str(), "r")) == NULL)			
	{
		printf("ERROR in FindBestComponent - Could not open objective file: %s.\n", obj_file_name.c_str());
		exit(1);
	}

	// Read through file untill all components are read
	while (true)
	{
		// Check if end-of-file was reached
		if (feof(obj_file))
			break;

		if (fscanf(obj_file, "%i,%lf\n", &comp_num, &tmp_fit) < 2)
		{
			printf("ERROR in FindBestCompoents - Objective file not formatted correctly.\n");
			exit(1);
		}

		best_fitness[comp_num-1] = tmp_fit;
	}
	// Check that the number of read components matches expected value
	if(comp_num != num_comp)
	{
		printf("ERROR in FindBestCompoents - Objective file not formatted correctly. Number of components does not match expected\n");
		exit(1);
	}
	
	fclose(obj_file);


	// Loop through remaining islands
	for (int i_isl = 2; i_isl <= num_islands; i_isl++)
	{
		ss_isl.str(std::string()); // Clear stream contents
		ss_isl << i_isl; // Send current island to stream

		// Create temp file name
		obj_file_name = out_dir + "I_" + ss_isl.str() + "/objective.csv";
			
		// Try to open file
		if((obj_file = fopen(obj_file_name.c_str(), "r")) == NULL)			
		{
			printf("ERROR in FindBestComponent - Could not open objective file: %s.\n", obj_file_name.c_str());
			exit(1);
		}
		
		// Read through file untill desired component is found
		while (true)
		{
			// Check if end-of-file was reached
			if (feof(obj_file))
				break;

			if (fscanf(obj_file, "%i,%lf\n", &comp_num, &tmp_fit) < 2)
			{
				printf("ERROR in FindBestCompoents - Objective file not formatted correctly.\n");
				exit(1);
			}

			if(tmp_fit > best_fitness[comp_num-1])
			{
				best_fitness[comp_num-1] = tmp_fit;
				best_island[comp_num-1] = i_isl;
			}

		}

		// Close file
		fclose(obj_file);
	}
	
	// Loop through all components
	std::string best_obj_file = new_dir + "objective.csv";
	FILE *best_obj;
	if((best_obj = fopen(best_obj_file.c_str(), "w")) == NULL)
	{
		printf("ERROR in FindBestComponent - Could not open output file (%s)\n", best_obj_file.c_str());
		exit(1);		
	}

	for(int i_comp = 1; i_comp <= num_comp; ++i_comp)
	{
		ss_i_comp.str(std::string()); // Clear stream contents
		ss_i_comp << i_comp; // Send current island to stream

		// Move best file to output directory
		best_file_name = out_dir + "I_" + std::to_string(best_island[i_comp-1]) + "/" + out_tag + "_comp" + ss_i_comp.str() + "_"
		 + std::to_string(best_island[i_comp-1]) + ".out";
		new_file_name = new_dir + out_tag + "_comp" + ss_i_comp.str() + "_best.out"; // Create new file name
		
		if (rename(best_file_name.c_str(), new_file_name.c_str()) != 0)
		{
			printf("ERROR in FindBestComponent - Could not rename file: %s\n", best_file_name.c_str());
		}

		fprintf(best_obj, "%i,%lf\n", i_comp, best_fitness[i_comp-1]);
	}
	fclose(best_obj);

}