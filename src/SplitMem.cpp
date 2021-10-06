#include "SplitMem.h"

SplitMem::SplitMem()
{

}

SplitMem::SplitMem(const std::size_t num_nodes)
{
    part1.resize(num_nodes);
	part2.resize(num_nodes);
	nodes_to_pick.resize(num_nodes);
    index_in_clust.resize(num_nodes);
    end_index.resize(num_nodes);
	in_clust.resize(num_nodes);
    in_part.resize(num_nodes);
	part_num.resize(num_nodes);
    has_edges.resize(num_nodes);
}

SplitMem::~SplitMem()
{
}

void SplitMem::initialize(const std::size_t num_nodes)
{
    part1.resize(num_nodes);
	part2.resize(num_nodes);
	nodes_to_pick.resize(num_nodes);
    index_in_clust.resize(num_nodes);
    end_index.resize(num_nodes);
	in_clust.resize(num_nodes);
    in_part.resize(num_nodes);
	part_num.resize(num_nodes);
    has_edges.resize(num_nodes);
}