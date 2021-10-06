#include "Graph.h"
#include "string.h"

int main(int argc, char **argv)
{
    if(argc != 3)
    {
        printf("Usage: ConGraphToCA input_file output_file\n");
        exit(1);
    }

    std::string input_file(argv[1]);
    std::string output_file(argv[2]);

    Graph my_graph;
    my_graph.readFile(input_file);
    my_graph.saveAsClusterAssignment(output_file);

    return 0;
}