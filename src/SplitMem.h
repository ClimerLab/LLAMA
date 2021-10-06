#ifndef SPLIT_MEM_H
#define SPLIT_MEM_H

#include <vector>

class SplitMem
{
private:
    std::vector<std::size_t> part1;
	std::vector<std::size_t> part2;
	std::vector<std::size_t> nodes_to_pick;
    std::vector<std::size_t> index_in_clust;
    std::vector<std::size_t> end_index;
	std::vector<bool> in_clust;
    std::vector<bool> in_part;	
	std::vector<bool> part_num;
    std::vector<bool> has_edges;

    friend class Chromosome;
	
public:
    SplitMem();
    SplitMem(const std::size_t);
    ~SplitMem();

    void initialize(const std::size_t);
};

#endif