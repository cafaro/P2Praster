/** @file*/

/**
* @author Antonio Mariani
* @date 05/10/19
*/


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <cstring>
#include "boost/functional/hash.hpp"
#include <unordered_map>
#include <chrono>
#include <unordered_set>

using namespace std;


#ifndef P2PRASTER_RASTER_H
#define P2PRASTER_RASTER_H

/**
 *
 * @param [in] seed
 * @param [in] value
 */
void hash_combine(std::size_t& seed, std::size_t value);

/**
 * Struct used to obtain hash for data structure with key: array<int, n>
 */
struct container_hasher {
    template<class T>
    std::size_t operator()(const T& c) const {
        std::size_t seed = 0;
        for(const auto& elem : c) {
            hash_combine(seed, std::hash<typename T::value_type>()(elem));
        }
        return seed;
    }
};

typedef unordered_map<array<int, 2>, double, container_hasher> hashmap; // contains projection result, key is tile, value is its cardinality
typedef vector<unordered_set<array<int, 2>, container_hasher>> vectorSet2; // vector contains clusters. Each cluster (unordered_set) contains its tiles
typedef vector<unordered_set<array<double , 2>, container_hasher>> vectorSet2D; // vector contains clusters. Each cluster (unordered_set) contains its points
typedef vector<unordered_set<array<int, 3>, container_hasher>> vectorSet3; // vector contains clusters. Each cluster (unordered_set) contains its tiles with respective cardinality
typedef unordered_set<array<double, 2>, container_hasher> unSet2D; // cluster that contains its points
typedef unordered_set<array<int, 2>, container_hasher> unSet2; // Cluster that contains tiles
typedef unordered_set<array<int, 3>, container_hasher> unSet3; // Cluster that contains tiles and respective cardinality
typedef unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> hashmapUnset; // contains mapping tiles and the points associated


/**
 * THis function use template in order to use an only function to print both type of cluster
 * that can be used in the two different type of algorithm in the main function. This function
 * can take as argument both type of cluster and print on terminal the number of clusters,
 * the number of tiles clustered and print on a csv each cluster with correspondent tiles.
 *
 * @tparam T - Type of array, can be: array<int, 3> or array<int, 2>
 * @param [in] clusters - Data structure which contains clusters to print
 * @param [in] peerID - Peer id that want to write on file
 * @return 0 in case of success, -7 in case of open file error, -3 in case of bad data structures
 */
template <typename T>
int printClusters(vector<unordered_set<T, container_hasher>> &clusters, int peerID);

/**
 * THis function use template in order to use an only function to print both type of clusters
 * that can be used in the two different type of algorithm in the main function. This function
 * can take as argument both type of cluster and print on file named clustered.csv all dataset
 * point with its clusters. Any point not clustered is assigned at cluster 0.
 *
 * @tparam T - Type of array, can be: array<int, 3> or array<int, 2>
 * @param [in] clusters - Data structure which contains clusters to print
 * @param [in] all_points - Data structure which contains the mapping tiles-points
 * @return 0 in case of success, -7 in case of open file error, -3 in case of bad structure
 */
template <typename T>
int printAllPointsClustered(vector<unordered_set<T, container_hasher>> &clusters, unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points, string outputfile);

template <typename T>
int getAllPointsClustered(vector<unordered_set<T, container_hasher>> &clusters, unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points, vectorSet2D &all_pointsClusters);


/**
 * This function read the dimension of the csv file pass as argument
 *
 * @param [in] name_file - The input file name
 * @param [in,out] row - Number of csv rows
 * @param [in,out] column - Number of csv columns
 * @return Return 0 if success, -7 in case of read file error
 */
int getDim(string name_file, int &row, int &column);

/**
 * This function read the csv file pass as argument and load the dataset
 * into the matrix m
 *
 * @param [in,out] m - Matrix where load dataset
 * @param [in] name_file - Name of dataset file
 * @param [in] column - Number of file column
 * @return Return 0 if success, -7 in case of read file error, -4 in case of read NaN
 */
int loadData(double **m, string name_file, int column);

/**
 * This functions uses a simple for cycle in order to obtains the key of the neighbors of coordinate in input,
 * then search this neighbors into projection and add them into result
 *
 * @param [in] coordinate - Coordinate of tile of which we want to obtain its neighbors
 * @param [in,out] projection - Data structure contains all tiles where search the tile's neighbors
 * @param [in,out] result - Contains all tile's neighbors found into projection
 * @return 0 if success, -2 in case of insert error
 */
int getNeighbors(array<int, 2> coordinate, hashmap &projection, unSet2 &result);

/**
 * This functions uses a simple for cycle in order to obtains the key of the neighbors of coordinate in input,
 * then search this neighbors into projection and squareProjection and add them to result
 *
 * @param [in] coordinate - Coordinate of tile of which we want to obtain its neighbors and its cardinality
 * @param [in,out] squareProjection - Data structure contains all peer's tiles where search the tile's neighbors
 * @param [in,out] projection - Data structure contains all remaining tiles where search the tile's neighbors
 * @param [in,out] result - Contains all tile's neighbors found into projection or squareProjection
 * @return 0 if success, -2 in case of insert error
 */

int getNeighbors(array<int, 3> coordinate, hashmap &squareProjection, hashmap &projection, unSet3 &result);

/**
 * This function performs the first step of raster projection, it simple multiply
 * the points for a scalar and maintains only the integer part of the product,
 * then save into projection the tiles and count its occurrences
 * Note: this function don't apply the threshold to the tiles
 *
 * @param [in] m - Dataset
 * @param [in] precision - Parameter for raster algorithm
 * @param [in,out] projection - Data structure where save the projection result of raster algorithm
 * @param [in] start - First point in the dataset of peer
 * @param [in] end - Last point in the dataset of peer
 * @return 0 in case of success, -2 in case of insert error
 */
int mapToTiles(double **m, double precision, hashmap &projection, int start, int end);

/**
 * This function implements the raster' variant of the algorithm.
 * During the projection phase, using an hashmap to keeps track of the
 * associations tile-point using the tile as key and as value a list
 * with relative associated points
 *
 * @param [in] m - Dataset
 * @param [in] precision - Parameter for raster algorithm (it is the value for the approximation)
 * @param [in] threshold - Parameter for raster algorithm (all tiles with cardinality < threshold are delete)
 * @param [in] n - Number of points in Dataset
 * @param [in,out] projection - Data structure where save the projection result of raster algorithm
 * @param [in,out] all_points - Data structure where save the mapping tiles-points
 * @return 0 in case of success, -10 in case of find error, -2 in case of insert error
 */
int mapToTilesPrime(double **m, double precision, int threshold, int n, hashmap &projection, hashmapUnset &all_points);

/**
 * This function remove from projection structure all tiles with
 * cardinality < threshold
 *
 * @param [in,out] projection - Data structure with all tiles
 * @param [in] threshold - Parameter for raster algorithm (all tiles with cardinality < threshold are delete)
 * @return 0 if success, -3 in case of bad projection structure
 */
int projectionThreshold(hashmap &projection, int threshold);


/**
 * This function take one by one each tile in projection, try to find its neighbors
 * (in projection) and form a cluster with these adjacent tiles.
 * If the cluster has a number of tiles < min_size it is discarded or
 * add it into clusters structure
 *
 *
 * @param [in,out] projection - Data structure with tiles
 * @param [in] min_size - Parameter for raster algorithm (it is the minimum number of tiles in order to form a cluster)
 * @param [in,out] clusters - Data structure where save the clusters found
 * @return 0 in case of success, -3 in case of bad projection structure, -2 in case of insert error
 */
int clusteringTiles(hashmap &projection, int min_size, vectorSet2 &clusters);

/**
 * This function take one by one each tile in squareProjection, try to find its neighbor
 * (in squareProjection and projection) and form a cluster with these adjacent tiles.
 * If the cluster has a number of tiles < min_size it is discarded or
 * add it into clusters structure
 *
 * @param [in,out] squareProjection - Data structure with peer's tiles
 * @param [in,out] projection - Data structure with remaining tiles
 * @param [in] min_size - Parameter for raster algorithm (it is the minimum number of tiles in order to form a cluster)
 * @param [in,out] clusters - Data structure where save the clusters found
 * @return 0 in case of success, -3 if can't complete clustering, -2 in case of insert error
 */
int clusteringTiles(hashmap &squareProjection, hashmap &projection, int min_size, vectorSet3 &clusters);

/**
 * This function compute the mean Shannon entropy and mean Density for each cluster
 * and print the results on terminal
 *
 * @param [in] clusters - Data structure which contains clusters
 * @param [in] all_points - Data structure which contains mapping tile-points for raster prime' variant
 * @param [in] precision - Parameter for raster algorithm
 * @return 0 if success, -1 in case of error on allocation memory
 */
int analyzeClusters(vectorSet2 &clusters, hashmapUnset &all_points, double precision);// raster'



#endif //P2PRASTER_RASTER_H