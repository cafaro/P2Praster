#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include "raster.h"
#include "error.h"


using namespace std;
using namespace std::chrono;



void hash_combine(std::size_t& seed, std::size_t value) {
    seed ^= value + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

int getDim(string name_file, int &row, int &column) {
    // Initialization
    row = 0;
    column = 0;
    ifstream inputFile(name_file);
    if(inputFile.is_open()){

        string s;
        while (getline(inputFile, s)) {

            if (row == 0) {
                istringstream ss(s);
                string line;
                while (getline(ss, line, ',')) {
                    column++;
                }
            }

            row++;

        }

        inputFile.close();

    }
    else{
        cerr << "the file " << name_file << " is not open!" << endl;
        return fileError(__FUNCTION__);
    }

    return 0;
}

int loadData(double **m, string name_file, int column) {

    if(!m){
        cerr << "The argument " << m << " is not a valid pointer" << endl;
        return argumentError(__FUNCTION__);
    }

    ifstream inputFile(name_file);
    if(inputFile.is_open()){

        int row = 0;
        string s;
        while (getline(inputFile, s)) {

            if (s[0] != '#') {
                istringstream ss(s);

                for (int i = 0; i < column; i++) {
                    string line;
                    getline(ss, line, ',');
                    try {
                        m[row][i] = stod(line);
                    }
                    catch (const std::invalid_argument &e) {
                        cerr << "NaN found in file " << name_file << " line " << row
                             << endl;
                        inputFile.close();
                        return arithmeticError(__FUNCTION__);
                    }
                }

            }

            row++;

        }

        inputFile.close();

    }
    else {
        return fileError(__FUNCTION__);
    }

    return 0;
}

int getNeighbors(array<int, 2> coordinate, hashmap &projection, unSet2 &result) {
    int x = coordinate[0];
    int y = coordinate[1];
    int radius = 1; /**!< Hierarchical Level of neighbors involved, 1 means the eight nearest neighbors, 2 also includes the neighbors of neighbors and so on */
    hashmap::iterator it;
    unSet2 neighbors;
    unSet2::iterator it_neighbor;

    /// generation of neighbor coordinate
    for (int i = -radius; i <= radius ; i++) {
        for (int j = -radius; j <= radius ; j++) {
            if (i == 0 && j == 0)
                continue;
            neighbors.insert({x + i, y + j});

        }
    }

    /// if a neighbor is in the projection, remove it from there and add it in result
    it_neighbor = neighbors.begin();
    for (int i = 0; i < neighbors.size(); i++) {
        it = projection.find(*it_neighbor);
        it_neighbor++;
        if (it != projection.end()) {
            auto check = result.insert(it -> first);
            if (!check.second) {
                return insertError(__FUNCTION__);
            }
            projection.erase(it++); // delete tile from projection and increase iterator (after erase (it) is not more valid)
        }
    }

    return 0;
}

int getNeighbors(array<int, 3> coordinate, hashmap &squareProjection, hashmap &projection, unSet3 &result) {
    int x = coordinate[0];
    int y = coordinate[1];
    int radius = 1; /**!< Hierarchical Level of neighbors involved, 1 means the eight nearest neighbors, 2 also includes the neighbors of neighbors and so on */
    hashmap::iterator it;
    unSet2 neighbors;
    unSet2::iterator it_neighbor;

    /// generation of neighbor coordinate
    for (int i = -radius; i <= radius ; i++) {
        for (int j = -radius; j <= radius ; j++) {
            if (i == 0 && j == 0)
                continue;
            neighbors.insert({x + i, y + j});
        }
    }

    /// if a neighbor is in the projection or in the squareProjection, remove it from there and add it in result
    it_neighbor = neighbors.begin();
    for (int i = 0; i < neighbors.size(); i++) {
        it = squareProjection.find(*it_neighbor);
        if (it != squareProjection.end()) {
            auto check = result.insert({(it->first)[0], (it->first)[1],(int) it->second});
            if (!check.second) {
                return insertError(__FUNCTION__);
            }
            squareProjection.erase(it++); // delete tile from squareProjection and increase iterator (after erase (it) is not more valid)
        } else {
            it = projection.find(*it_neighbor);
            if (it != projection.end()) {
                auto check = result.insert({(it->first)[0], (it->first)[1],(int) it->second});
                if (!check.second) {
                    return insertError(__FUNCTION__);
                }
                projection.erase(it++); // delete tile from projection and increase iterator (after erase (it) is not more valid)
            }
        }
        it_neighbor++;
    }

    return 0;
}

int mapToTiles(double **m, double precision, hashmap &projection, int start, int end) {
    double scalar;
    int x, y;
    hashmap::iterator it;
    scalar = pow(10, precision);

    for (int i = start; i <= end; i++) {
        x =(int) (m[i][0] * scalar);
        y =(int) (m[i][1] * scalar);
        array<int, 2> tile = {x,y};
        it = projection.find(tile);
        if (it != projection.end()) {
            it->second++;
        } else {
            auto a = projection.insert({tile, 1.0});
            if (!(a.second)) {
                return insertError(__FUNCTION__);
            }
        }
    }

    return 0;
}

int mapToTilesPrime(double **m, double precision, int threshold, int n, hashmap &projection, hashmapUnset &all_points) {
    double scalar;
    int x, y;
    hashmap::iterator it;
    hashmapUnset::iterator it_map_all_points;

    scalar = pow(10, precision);

    for (int i = 0; i < n; i++) {
        x =(int) (m[i][0] * scalar);
        y =(int) (m[i][1] * scalar);
        array<int, 2> tile = {x,y};
        it = projection.find(tile);
        if (it != projection.end()) {
            it->second++;
            /// if the tile is in the projection then it must be present into all_points
            it_map_all_points = all_points.find(tile);
            if (it_map_all_points == all_points.end()) {
                return findError(__FUNCTION__);
            }

            /// update the tile-points mapping, adding the new point
            auto a = (it_map_all_points -> second).insert({m[i][0],m[i][1]});
            if (!(a.second)) {
                return insertError(__FUNCTION__);
            }
        } else {
            auto a = projection.insert({tile, 1.0});
            if (!(a.second)) {
                return insertError(__FUNCTION__);
            }

            /// create tile-point mapping
            unordered_set<array<double , 2>, container_hasher> point;
            auto b = point.insert({m[i][0], m[i][1]});
            if (!(b.second)) {
                return insertError(__FUNCTION__);
            }
            auto c = all_points.insert({tile, point});
            if (!(c.second)) {
                return insertError(__FUNCTION__);
            }
        }

    }
    /// remove tiles with cardinality < threshold
    projectionThreshold(projection, threshold);

    return 0;

}

int projectionThreshold(hashmap &projection, int threshold) {
    if (projection.empty()) {
        cerr << "Bad projection data structure" << endl;
        return dataError(__FUNCTION__);
    }
    hashmap::iterator it;
    it = projection.begin();
    while(it != projection.end()) {
        if (it -> second < threshold) {
            projection.erase(it++);
        } else {
            it++;
        }
    }
    return 0;

}

int clusteringTiles(hashmap &projection, int min_size, vectorSet2 &clusters) {
    if (projection.empty()) {
        cerr << "Bad projection data structure" << endl;
        return dataError(__FUNCTION__);
    }
    hashmap::iterator iterator;

    /// read and remove all tiles one by one from projection
    while ((iterator = projection.begin()) != projection.end()) {
        // read and remove first element recursively
        array<int, 2> x = iterator->first;
        projection.erase(iterator++);

        /// candidate cluster
        unSet2 visited;
        auto a = visited.insert(x);
        if (!(a.second)) {
            return insertError(__FUNCTION__);
        }

        /// get neighbors of tile in exam
        unSet2 to_check;
        getNeighbors(x, projection, to_check);


        /// for each neighbor, try to find his neighbors recursively in order to add they to a single cluster
        while (!to_check.empty()) {
            array<int, 2> value = *to_check.begin() ;
            to_check.erase((to_check.begin()));
            auto b = visited.insert(value);
            if (!(b.second)) {
                return insertError(__FUNCTION__);
            }

            unSet2 temp;
            getNeighbors(value, projection, temp);
            while (!temp.empty()) {
                auto c = to_check.insert(*temp.begin());
                if (!(c.second)) {
                    return insertError(__FUNCTION__);
                }
                temp.erase((temp.begin()));
            }
        }
        /// validate visited as cluster if is size is >= min size
        if (visited.size() >= min_size) {
            clusters.push_back(visited);
        }
    }

    return 0;
} // unordered

int clusteringTiles(hashmap &squareProjection, hashmap &projection, int min_size, vectorSet3 &clusters) {
    if (squareProjection.empty()) {
        return -3;
    }
    hashmap::iterator iterator;

    /// read and remove all tiles one by one from projection
    while ((iterator = squareProjection.begin()) != squareProjection.end()) {
        // read and remove first element recursively
        array<int, 3> x{};
        x = {(iterator->first)[0], (iterator->first)[1],(int) iterator->second};
        squareProjection.erase(iterator++);

        /// candidate cluster
        unSet3 visited;
        auto a = visited.insert(x);
        if (!(a.second)) {
            return insertError(__FUNCTION__);
        }

        /// get neighbors of tile in exam
        unSet3 to_check;
        int return_value = getNeighbors(x, squareProjection, projection, to_check);
        if(return_value != 0){
            cerr << "getNeighbors() failed" << endl;
            return functionError(__FUNCTION__);
        }

        /// for each neighbor, try to find his neighbors recursively in order to add them to a single cluster
        while (!to_check.empty()) {
            array<int, 3> value = *to_check.begin() ;
            to_check.erase((to_check.begin()));
            auto b = visited.insert(value);
            if (!(b.second)) {
                return insertError(__FUNCTION__);
            }

            unSet3 temp;
            getNeighbors(value, squareProjection, projection, temp);
            while (!temp.empty()) {
                auto c = to_check.insert(*temp.begin());
                if (!(c.second)) {
                    return insertError(__FUNCTION__);
                }
                temp.erase((temp.begin()));
            }
        }
        /// visited is a valid cluster if is size is >= min size
        if (visited.size() >= min_size) {
            clusters.push_back(visited);
        }
    }

    if (clusters.empty())
        return -3;
    else
        return 0;
}

// under the function there are the two variants of T type
template <typename T>
int printAllPointsClustered(vector<unordered_set<T, container_hasher>> &clusters, unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points, string outputfile){
    if (clusters.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }

    if (all_points.empty()) {
        cerr << "Bad all_points data structure" << endl;
        return dataError(__FUNCTION__);
    }

    ofstream outfile("./log/clusters_" + outputfile + ".csv");
    if (!outfile.is_open()) {
        cerr << "can not open file clusters_" + outputfile + ".csv" + "for writing" << endl;
        return fileError(__FUNCTION__);
    }

    int count_not_clustered = 0;
    int count_tiles = 0;
    int count_points = 0;
    double outputOnFile = !outputfile.empty();

    hashmapUnset::iterator it_map_all_points;
    unordered_set<array<double , 2>, container_hasher>::iterator it_set_all_points;
    typename unordered_set<T, container_hasher>::iterator it_tiles;


    /************ for each cluster in clusters ************/
    for (int j = 0; j < clusters.size(); j++) {
        //cout << "Cluster n° " << j << " with size " << cluster.at(j).size() << ": " << endl;
        if (clusters.at(j).empty()) {
            cerr << "Bad cluster structure" << endl;
            outfile.clear();
            outfile.close();
            return dataError(__FUNCTION__);
        }
        it_tiles = clusters.at(j).begin(); // pointer to start of j-th cluster in clusters (cluster = list of tiles, clusters = list of clusters)
        /************ for each tile in cluster j-th ************/
        for (int i = 0; i < clusters.at(j).size(); i++) {
            it_map_all_points = all_points.find({(*it_tiles)[0], (*it_tiles)[1]}); // try to find in all_points the tile (with its list of points) from cluster
            if (it_map_all_points != all_points.end()) {
                if ((it_map_all_points -> second).empty()) {
                    cerr << "Bad all_points structure" << endl;
                    outfile.clear();
                    outfile.close();
                    return dataError(__FUNCTION__);

                }
                it_set_all_points = (it_map_all_points -> second).begin(); // pointer to the first element in the list of points associated to the found tile
                /************ for each point in the tile ************/
                for (int k = 0; k < (it_map_all_points -> second).size(); k++) {
                    outfile << (*it_set_all_points)[0] << ",";
                    outfile << (*it_set_all_points)[1] << ",";
                    outfile << j + 1 <<endl;
                    it_set_all_points++;
                    count_points++;
                }
                all_points.erase(it_map_all_points++);
            }
            it_tiles++;
            count_tiles++;
        }
    }
    // in case of points/tiles not clustered
    if (!all_points.empty()) {
        it_map_all_points = all_points.begin(); // first tile remaining in all_points
        /************ for each tile that is not in the clusters ************/
        for (int i = 0; i < all_points.size(); i++) {
            if (it_map_all_points != all_points.end()) {
                if ((it_map_all_points -> second).empty()) {
                    cerr << "Bad all_points structure" << endl;
                    outfile.clear();
                    outfile.close();
                    return dataError(__FUNCTION__);
                }
                it_set_all_points = (it_map_all_points -> second).begin(); // pointer to the first element in the list of points associated to the found tile
                /************ for each point in the tiles that is not in the clusters ************/
                for (int k = 0; k < (it_map_all_points -> second).size(); k++) {
                    outfile << (*it_set_all_points)[0] << ",";
                    outfile << (*it_set_all_points)[1] << ",";
                    outfile << 0 <<endl;
                    it_set_all_points++;
                    count_not_clustered++;
                }
                it_map_all_points++;
            }
        }
    }

    outfile.close();
    if (!outputOnFile) {
        cout << "\n\n";
        cout << "Points not clustered: " << count_not_clustered << endl;
        cout << "Points clustered: " << count_points << endl;
        cout << "Total points analyzed " << count_points + count_not_clustered << endl;
        cout << "Tile not clustered: " << all_points.size() << endl;
        cout << "Tiles clustered: " << count_tiles << endl;
        cout << "Clusters: " << clusters.size() << endl;
    } else {
        outfile.open("./log/metrics_" + outputfile + ".txt", ofstream::app);

        if (outfile.is_open()) {
            outfile << "\n\n";
            outfile << "Points not clustered: " << count_not_clustered << endl;
            outfile << "Points clustered: " << count_points << endl;
            outfile << "Total points analyzed " << count_points + count_not_clustered << endl;
            outfile << "Tile not clustered: " << all_points.size() << endl;
            outfile << "Tiles clustered: " << count_tiles << endl;
            outfile << "Clusters: " << clusters.size() << endl;
        }

        outfile.close();
    }


    return 0;
}
/**
 * These declarations are necessary for compiler in order to correctly link
 * the header of a function with its implementation in presence of a template
 */
template int printAllPointsClustered<array<int, 3>>(vectorSet3 &clusters, hashmapUnset &all_points, string outputfile);
template int printAllPointsClustered<array<int, 2>>(vectorSet2 &clusters, hashmapUnset &all_points, string outputfile);


// under the function there are the two variants of T type
template <typename T>
int printClusters(vector<unordered_set<T, container_hasher>> &clusters, int peerID) {
    if (clusters.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }

    ofstream outfile(to_string(peerID) +".csv");
    if (!outfile.is_open())
        return fileError(__FUNCTION__);

    cout <<  "n° cluster: " << clusters.size() << endl;
    typename unordered_set<T, container_hasher>::iterator it;

    int count_tiles = 0;
    for (int j = 0; j < clusters.size(); j++) {
        outfile << "Cluster n° " << j << " with size " << clusters.at(j).size() << ": " << endl;
        if (clusters.at(j).empty()) {
            cerr << "Bad cluster structure" << endl;
            outfile.clear();
            outfile.close();
            return dataError(__FUNCTION__);
        }
        it = clusters.at(j).begin(); // pointer to start of j-th cluster (cluster = list of tiles)
        for (int i = 0; i < clusters.at(j).size(); i++) {
            count_tiles++; // count the total number of clustered tiles
            outfile << (*it)[0] << ",";
            outfile << (*it)[1] << ",";
            outfile << j << endl;
            it++; // next tile of the current cluster
        }
    }
    outfile.close();
    cout << "Tiles clustered: " << count_tiles << endl;

    return 0;

}
/**
 * These declarations are necessary for compiler in order to correctly link
 * the header of a function with its implementation in presence of a template
 */
template int printClusters<array<int, 3>>(vectorSet3 &clusters, int peerID);
template int printClusters<array<int, 2>>(vectorSet2 &clusters, int peerID);

int analyzeClusters(vectorSet2 &clusters, hashmapUnset &all_points, double precision) {
    int returnValue = -1;
    double scalar = pow(10, precision);
    double area = scalar * scalar;
    double mean_shannon;
    double mean_density;
    hashmapUnset::iterator it_map_all_points;
    unordered_set<array<double , 2>, container_hasher>::iterator it_set_all_points;
    unSet2::iterator it_tiles;

    double* shannon = nullptr;
    double* density = nullptr;
    double* size_cluster = nullptr;
    double** tile_points = nullptr;

    shannon = new (nothrow) double[clusters.size()];
    if(!shannon)
        return memoryError(__FUNCTION__);


    density = new (nothrow) double[clusters.size()];
    if(!density){
        returnValue = memoryError(__FUNCTION__);
        goto ON_EXIT;
    }


    size_cluster= new (nothrow) double[clusters.size() * sizeof(double)];
    if(!size_cluster){
        returnValue = memoryError(__FUNCTION__);
        goto ON_EXIT;
    }
    double pji; // probability that a point is in the i-th tile of the j-th cluster



    tile_points = new (nothrow) double*[clusters.size()];
    if(!tile_points){
        returnValue = memoryError(__FUNCTION__);
        goto ON_EXIT;
    }
    // allocating vectors that will contain n° of points for each tile i of respective cluster j
    for (int j = 0; j < clusters.size(); j++) {
        //cout << "Cluster n° " << j << " with size " << clusters.at(j).size() << ": " << endl;
        tile_points[j] = new (nothrow) double[clusters.at(j).size()];
        if (!tile_points[j]) {
            cerr << "Not enough memory allocating the " << j << "-th tile_points array" << endl;
            returnValue = memoryError(__FUNCTION__);
            goto ON_EXIT;
        }
    }

    /************ for each cluster in clusters ************/
    for (int j = 0; j < clusters.size(); j++) {
        //cout << "Cluster n° " << j << " with size " << cluster.at(j).size() << ": " << endl;
        it_tiles = clusters.at(j).begin(); // pointer to start of j-th cluster in clusters (cluster = list of tiles, clusters = list of clusters)
        /************ for each tile in cluster j-th ************/
        for (int i = 0; i < clusters.at(j).size(); i++) {  // clusters.at(j).size() represent the number of tiles contained in the j-th cluster
            it_map_all_points = all_points.find((*it_tiles)); // try to find in all_points the tile (with its list of points) from cluster
            size_cluster[j] += (it_map_all_points -> second).size();
            tile_points[j][i] = (it_map_all_points -> second).size();
            it_tiles++;
        }
    }
    /************ for each cluster in clusters ************/
    for (int j = 0; j < clusters.size(); j++) {
        /************ for each tile in cluster j-th ************/
        for (int i = 0; i < clusters.at(j).size(); i++) {  // clusters.at(j).size() represent the number of tiles contained in the j-th cluster
            pji = tile_points[j][i]/size_cluster[j];
            shannon[j] += -(pji * log2(pji));
        }
    }

    // calculating the density for each cluster
    /************ for each cluster in clusters ************/
    for (int j = 0; j < clusters.size(); j++) {
        density[j] = size_cluster[j]/(clusters.at(j).size() * area);
    }

    mean_shannon = 0;
    mean_density = 0;
    for (int j = 0; j < clusters.size(); j++) {
        mean_shannon += shannon[j];
        mean_density += density[j];
    }

    mean_shannon = mean_shannon/clusters.size();
    mean_density = mean_density/clusters.size();

    cout << "Shannon mean: " << mean_shannon << endl;
    cout << "Density mean: " << mean_density << endl;

    returnValue = 0;

    ON_EXIT:

    if(shannon != nullptr)
        delete[] shannon, shannon = nullptr;

    if(density != nullptr)
        delete[] density, density = nullptr;

    if(size_cluster != nullptr)
        delete[] size_cluster, size_cluster = nullptr;

    for(int x = 0; x < clusters.size(); x++)
        if(tile_points[x] != nullptr)
            delete[] tile_points[x], tile_points[x] = nullptr;

    if(tile_points != nullptr)
        delete[] tile_points, tile_points = nullptr;

    return returnValue;
}

template<typename T>
int getAllPointsClustered(vector<unordered_set<T, container_hasher>> &clusters,
                          unordered_map<array<int, 2>, unordered_set<array<double, 2>, container_hasher>, container_hasher> &all_points,
                          vectorSet2D &all_pointsClusters) {

    int count_points = 0;
    if (clusters.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }

    if (all_points.empty()) {
        cerr << "Bad all_points data structure" << endl;
        return dataError(__FUNCTION__);
    }

    hashmapUnset::iterator it_map_all_points;
    unordered_set<array<double , 2>, container_hasher>::iterator it_set_all_points;
    typename unordered_set<T, container_hasher>::iterator it_tiles;


    /************ for each cluster in clusters ************/
    for (int j = 0; j < clusters.size(); j++) {
        unordered_set<array<double , 2>, container_hasher> s;
        //cout << "Cluster n° " << j << " with size " << cluster.at(j).size() << ": " << endl;
        if (clusters.at(j).empty()) {
            cerr << "Bad cluster structure" << endl;
            return dataError(__FUNCTION__);
        }
        it_tiles = clusters.at(j).begin(); // pointer to start of j-th cluster in clusters (cluster = list of tiles, clusters = list of clusters)
        /************ for each tile in cluster j-th ************/
        for (int i = 0; i < clusters.at(j).size(); i++) {

            it_map_all_points = all_points.find({(*it_tiles)[0], (*it_tiles)[1]}); // try to find in all_points the tile (with its list of points) from cluster
            if (it_map_all_points != all_points.end()) {
                if ((it_map_all_points -> second).empty()) {
                    cerr << "Bad all_points structure" << endl;
                    return dataError(__FUNCTION__);

                }
                it_set_all_points = (it_map_all_points -> second).begin(); // pointer to the first element in the list of points associated to the found tile
                /************ for each point in the tile ************/
                for (int k = 0; k < (it_map_all_points -> second).size(); k++) {
                    s.insert({ (*it_set_all_points)[0], (*it_set_all_points)[1]});
                    it_set_all_points++;
                    count_points++;
                }
                it_map_all_points++;
            }

            it_tiles++;
        }
        all_pointsClusters.push_back(s);
    }


    return 0;
}
template int getAllPointsClustered<array<int, 3>>(vectorSet3 &clusters, hashmapUnset &all_points, vectorSet2D &all_pointsClusters);
template int getAllPointsClustered<array<int, 2>>(vectorSet2 &clusters, hashmapUnset &all_points, vectorSet2D &all_pointsClusters);