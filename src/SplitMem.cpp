#include "assert.h"
#include <cstdio>
#include "SplitMem.h"

SplitMem::SplitMem() : num_nodes(0),
                       max_degree(0) {}

SplitMem::~SplitMem() {}

void SplitMem::initialize(const std::size_t num_nodes, const std::size_t max_degree)
{
  this->num_nodes = num_nodes;
  this->max_degree = max_degree;

  std::vector<std::size_t> part_pool_row(num_nodes-1);
  part_pool.push_back(part_pool_row);
  part_pool.push_back(part_pool_row);

  pool_size[0] = 0;
  pool_size[1] = 0;

  partion_number.resize(num_nodes, 2);

  std::vector<std::size_t> edge_list(max_degree, 0);
  for (std::size_t i = 0; i < num_nodes; ++i) {
    vertices.push_back(edge_list);
  }

  num_edges.resize(num_nodes, 0);
}

void SplitMem::reset() {
  pool_size[0] = 0;
  pool_size[1] = 0;

  for (std::size_t i = 0; i < num_nodes; ++i) {
    partion_number[i] = 2;
    num_edges[i] = 0;
  }  
}

void SplitMem::addEdge(const std::size_t src, const std::size_t tgt) {
  assert(src < num_nodes);
  assert(tgt < num_nodes);

  vertices[src][num_edges[src]] = tgt;
  ++num_edges[src];
}

void SplitMem::removeEdgeByIndex(const std::size_t src, const std::size_t node_index) {
  assert(src < num_nodes);
  assert(node_index < max_degree);

  vertices[src][node_index] = vertices[src][num_edges[src]-1];  
  --num_edges[src];
}

void SplitMem::addNodeToPartition(const std::size_t node_id, const std::size_t partition_id) {
  assert(node_id < num_nodes);

  partion_number[node_id] = partition_id;

  part_pool[partition_id][pool_size[partition_id]] = node_id;
  ++pool_size[partition_id];

  //fprintf(stderr, "Added node %lu to partition %lu\n", node_id, partition_id);
}

void SplitMem::removeNodeFromPartitionPool(const std::size_t node_id, const std::size_t partition_id) {
  assert(node_id < num_nodes);

  bool node_found = false;
  for (std::size_t i = 0; i < pool_size[partition_id]; ++i) {
    if (part_pool[partition_id][i] == node_id) {
      node_found = true;
      part_pool[partition_id][i] = part_pool[partition_id][pool_size[partition_id]-1];
      --pool_size[partition_id];
      return;
    }
  }

  if (!node_found) {
    fprintf(stderr, "Could not find node ID %lu in partition pool %lu\n", node_id, partition_id);
  }
}

std::size_t SplitMem::getPoolValue(const std::size_t paritition_id, const std::size_t index) const {
  assert(index < pool_size[paritition_id]);
  return part_pool[paritition_id][index];
}

std::size_t SplitMem::getPoolSize(const std::size_t partition_id) const {
  assert(partition_id < 2);
  return pool_size[partition_id];
}

std::size_t SplitMem::getNumEdges(const std::size_t node_id) const {
  assert(node_id < num_nodes);
  return num_edges[node_id];
}

std::size_t SplitMem::getPartNumber(const std::size_t node_id) const {
  assert(node_id < num_nodes);
  return partion_number[node_id];
}

void SplitMem::printPartNumber() {
  fprintf(stderr, "Partition Number:");
  for (std::size_t i = 0; i < num_nodes; ++i) {
    fprintf(stderr, " %lu", partion_number[i]);
  }
  fprintf(stderr, "\n");
}

void SplitMem::printPartPool() {
  fprintf(stderr, "Parition Pools\n");
  fprintf(stderr, "Pool 0: (%lu):", pool_size[0]);
  for (std::size_t i = 0; i < part_pool[0].size(); ++i) {
    fprintf(stderr, " %lu", part_pool[0][i]);
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "Pool 1: (%lu):", pool_size[1]);
  for (std::size_t i = 0; i < part_pool[1].size(); ++i) {
    fprintf(stderr, " %lu", part_pool[1][i]);
  }
  fprintf(stderr, "\n");
}

void SplitMem::printGraph() {
  for (std::size_t i = 0; i < num_nodes; ++i) {
    fprintf(stderr, "%lu:", i);
    for (std::size_t j = 0; j < num_edges[i]; ++j) {
      fprintf(stderr, " %lu", vertices[i][j]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}