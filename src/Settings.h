#ifndef SETTINGS_H
#define SETTINGS_H

#include <stddef.h>

// Graph settings
//const bool graph_directed = false;

// Chromosome settings
const double MUT_TO_SINGLETON = 0.01;
const double SPLIT_PROB = 0.15;
const double MERGE_PROB = 0.15;
const double REDISRIB_PROB = .10;
const double NODE_MUT_PROB = 0.05;
const double MUTATION_RATE = .05;
const double CLUST_SPLIT_RATE = 0.01;
const double CLUST_MERGE_RATE = 0.01;

// Clustering Restults settings
const std::size_t POP_SIZE = 100;
const std::size_t NUM_ELITE = 1;
const std::size_t TOUR_SIZE = 8;
const std::size_t NUM_FROM_TOUR = 1;
const std::size_t MAX_GEN = 20000;
const std::size_t MED_GEN = 1000;
const std::size_t SMALL_GEN = 100;
const std::size_t MIGRATION_GEN = 500;
const std::size_t SMALL_GRAGH_SIZE = 10;
const std::size_t MED_GRAPH_SIZE = 100;
const double CROSS_TYPE_CHANGE_PERC = .75;
const bool START_CROSS_TYPE = 1;  // Cluster crossover
const double NODE_CROSS_PROB = 0.5;

const double PERC_CHILD_CROSS = 0.6; // The percentage of children, after the elite chromosomes are saved, that are created from crossover
const double PERC_CHILD_CROSS_AND_MUT = .2; // The percentage of crossover children that are also mutated

//const bool CROSS_TYPE = 1; // 0=NODE_CROSS, 1=CLUSTER_CROSS
const double PERC_NODE_CROSS = 0.05; // The percentage of the nodes to crossover if using NODE_CROSS
const double PERC_CLUST_CROSS = 0.05; // The percentage of the nodes to crossover if using CLUST_CROSS

// General settings
const bool status_updates = true;
const bool record_details = false;
const int INIT_MODE = 4;
const double PERC_RND_WALK = 0.5;
const int FIT_EQ = 0;
const int MAX_CLUST_SIZE = 0;

#endif
