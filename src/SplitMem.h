#ifndef SPLIT_MEM_H
#define SPLIT_MEM_H

#include <vector>

class SplitMem
{
private:
  std::size_t num_nodes;
  std::size_t max_degree;
  std::vector<std::vector<std::size_t>> part_pool;
  std::size_t pool_size[2];
  std::vector<std::size_t> partion_number;
  std::vector<std::vector<std::size_t>> vertices;
  std::vector<std::size_t> num_edges;

  friend class Chromosome;

public:
  SplitMem();
  ~SplitMem();
  void initialize(const std::size_t num_nodes, const std::size_t max_degree);
  void reset();

  void addEdge(const std::size_t src, const std::size_t tgt);
  void removeEdgeByIndex(const std::size_t src, const std::size_t node_index);
  void addNodeToPartition(const std::size_t node_id, const std::size_t partition_id);
  void removeNodeFromPartitionPool(const std::size_t node_id, const std::size_t partition_id);

  std::size_t getPoolValue(const std::size_t paritition_id, const std::size_t index) const;
  std::size_t getPoolSize(const std::size_t partition_id) const;
  std::size_t getNumEdges(const std::size_t node_id) const;
  std::size_t getPartNumber(const std::size_t node_id) const;

  void printPartNumber();
  void printPartPool();
  void printGraph();
};

#endif