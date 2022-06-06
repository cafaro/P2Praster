/** @file*/

/**
* @author Antonio Mariani
* @date 27/11/19
*/

#include <algorithm>
#include "statistic.h"
#include "error.h"
using namespace std;

int getPrecision(vectorSet2D &C, vectorSet2D &T, double &Precision) {
    if (T.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }

    if (C.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }
    int returnValue;

    int *n = nullptr; /// n[i] is the size of i-th cluster
    int *max_nij = nullptr; /// nij is the cardinality of intersection between the clusters C[i] and T[i]
    double *_precision = nullptr;

    n = new (nothrow) int[C.size()];
    if (!n)
        return memoryError(__FUNCTION__);

    max_nij = new (nothrow) int[C.size()];
    if (!max_nij) {
        returnValue = memoryError(__FUNCTION__);
        goto ON_EXIT;
    }

    _precision = new (nothrow) double[C.size()];
    if(!_precision) {
        returnValue = memoryError(__FUNCTION__);
        goto ON_EXIT;
    }

    /// Get dimension of each cluster
    for (int i = 0; i < C.size(); ++i) {
        n[i] = C.at(i).size();
        if (!n[i]) {
            cerr << "Cluster can't be empty" << endl;
            returnValue = dataError(__FUNCTION__);
            goto ON_EXIT;
        }
    }

    /// Get max_nij
    for (int i = 0; i < C.size(); ++i) {
        returnValue = getMax_nij(T, C.at(i), max_nij[i]);
        if (returnValue < 0 ) {
            cerr << "Can't compute max_nij" << endl;
            goto ON_EXIT;
        }
    }

    /// compute precision for each clusters
    for (int i = 0; i < C.size(); ++i) {
        _precision[i] = (double) max_nij[i]/n[i];
    }

    /// compute average precision
    Precision = 0;
    for (int j = 0; j < C.size(); ++j) {
        Precision += _precision[j];
    }
    Precision /= C.size();

    returnValue = 0;

    ON_EXIT:

    if (n != nullptr)
        delete[] n, n = nullptr;

    if (max_nij != nullptr)
        delete[] max_nij, max_nij = nullptr;

    if (_precision != nullptr)
        delete[] _precision, _precision = nullptr;

    return returnValue;
}

int getNij(unSet2D &Ti, unSet2D &Ci) {
    if (Ti.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }

    if (Ci.empty()) {
        cerr << "Bad cluster data structure" << endl;
        return dataError(__FUNCTION__);
    }

    unSet2D intersection;

    /// Cluster intersection
    for (auto j : Ci) {
        if (Ti.find(j) != Ti.end()){
            auto check = intersection.insert(j);
            if (!(check.second)) {
                return insertError(__FUNCTION__);
            }
        }
    }
    return intersection.size();
}



int getMax_nij(vectorSet2D &T, unSet2D &Ci, int &max_nij) {
    if (T.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }

    if (Ci.empty()) {
        cerr << "Bad cluster data structure" << endl;
        return dataError(__FUNCTION__);
    }

    unSet2D intersection;
    unSet2D::iterator it;

    /// clusters intersection
    max_nij = 0;
    for (auto & i : T) {
        if (i.empty()) {
            cerr << "Bad cluster data structure" << endl;
            return dataError(__FUNCTION__);
        }

        for (auto j : Ci) {
            if (i.find(j) != i.end()){
                auto check = intersection.insert(j);
                if (!(check.second)) {
                    return insertError(__FUNCTION__);
                }
            }
        }

        if (intersection.size() > max_nij)
            max_nij = intersection.size();
        intersection.clear();
    }

    return 0;
}

int getMax_nij_mji(vectorSet2D &T, unSet2D &Ci, int &max_nij, int &mji) {
    if (T.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }

    if (Ci.empty()) {
        cerr << "Bad cluster data structure" << endl;
        return dataError(__FUNCTION__);
    }

    unSet2D intersection;
    unSet2D::iterator it;

    /// clusters intersection
    max_nij = 0;
    for (auto & i : T) {
        if (i.empty()) {
            cerr << "Bad cluster data structure" << endl;
            return dataError(__FUNCTION__);
        }

        for (auto j : Ci) {
            if (i.find(j) != i.end()) {
                auto check = intersection.insert(j);
                if (!(check.second)) {
                    return insertError(__FUNCTION__);
                }
            }
        }

        if (intersection.size() > max_nij) {
            max_nij = intersection.size(); // intersection set with highest cardinality
            mji = i.size(); // cardinality of the majority partition Tj_i
        }

        intersection.clear();
    }

    return 0;
}

int getRecall(vectorSet2D &C, vectorSet2D &T, double &Recall) {
    if (T.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }

    if (C.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }

    int returnValue;

    int *mji = nullptr; /// is the cardinality of major partition Tj_i
    int *max_nij = nullptr; /// nij is the cardinality of intersection between the clusters C[i] and T[j]
    double *_recall = nullptr;

    mji = new (nothrow) int[C.size()];
    if (!mji)
        return memoryError(__FUNCTION__);

    max_nij = new (nothrow) int[C.size()];
    if (!max_nij) {
        returnValue = memoryError(__FUNCTION__);
        goto ON_EXIT;
    }

    _recall = new (nothrow) double[C.size()]();
    if(!_recall) {
        returnValue = memoryError(__FUNCTION__);
        goto ON_EXIT;
    }

    /// Get max_nij and mji
    for (int i = 0; i < C.size(); ++i) {
        returnValue = getMax_nij_mji(T, C[i], max_nij[i], mji[i]);
        if (returnValue) {
            cerr << "Can't compute max_nij/mji" << endl;
            goto ON_EXIT;
        }
    }

    /// compute recall for each clusters
    for (int i = 0; i < C.size(); ++i) {
        if (max_nij[i] && mji[i]) {
            _recall[i] = (double) max_nij[i]/mji[i];
        }

    }

    /// compute average Recall
    Recall = 0;
    for (int i = 0; i < C.size(); ++i) {
        Recall += _recall[i];
    }
    Recall /= C.size();

    returnValue = 0;

    ON_EXIT:

    if (mji != nullptr)
        delete[] mji, mji = nullptr;

    if (max_nij != nullptr)
        delete[] max_nij, max_nij = nullptr;

    if (_recall != nullptr)
        delete[] _recall, _recall = nullptr;

    return returnValue;
}

double getFMeasure(double Recall, double Precision) {
    return 2.0*(Recall*Precision)/(Recall+Precision);
}

int getConditionalEntropy(vectorSet2D &C, vectorSet2D &T, int n, double &Entropy) {
    if (T.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }

    if (C.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }
    if (n <= 0) {
        cerr <<"Dataset points can't be <= 0 " << endl;
        return argumentError(__FUNCTION__);
    }

    int returnValue;
    double *pCi = nullptr; /// pCi[i] is the probability of cluster Ci: n_i/n
    double **pij = nullptr; /// pij[i][j] is the probability that a point in cluster i also belongs to partition j: nij/n
    double *pij_storage = nullptr;

    pCi = new (nothrow) double[C.size()];
    if (!pCi)
        return memoryError(__FUNCTION__);

    pij_storage = new (nothrow) double[C.size() * T.size()]();
    if (!pij_storage) {
        returnValue = memoryError(__FUNCTION__);
        goto ON_EXIT;
    }

    pij = new (nothrow) double*[C.size()];
    if (!pij) {
        returnValue = memoryError(__FUNCTION__);
        goto ON_EXIT;
    }

    for (int i = 0; i < C.size(); ++i) {
        pij[i] = &pij_storage[i*T.size()];
        if (!pij[i]) {
            returnValue = memoryError(__FUNCTION__);
            goto ON_EXIT;
        }
    }

    /// compute pCi
    returnValue = getPCi(C, n, pCi);
    if(returnValue) {
        cerr << "Can't compute pCi" << endl;
        goto ON_EXIT;
    }

    /// compute pij
    returnValue = getPij(C, T, pij, n);
    if(returnValue) {
        cerr << "Can't compute pij" << endl;
        goto ON_EXIT;
    }

    /// compute H(T|C)
    Entropy = 0;
    for (int i = 0; i < C.size(); ++i) {
        for (int j = 0; j < T.size(); ++j) {
            if (pij[i][j]) {
                Entropy -= pij[i][j] * log2(pij[i][j] / pCi[i]);
            }
        }
    }

    returnValue = 0;

    ON_EXIT:

    if (pCi != nullptr)
        delete[] pCi, pCi = nullptr;

    if (pij != nullptr)
        delete[] pij, pij = nullptr;

    if (pij_storage != nullptr)
        delete[] pij_storage, pij_storage = nullptr;

    return returnValue;
}

int getPij(vectorSet2D &C, vectorSet2D &T, double **pij, int n) {
    if (T.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }

    if (C.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }

    if (n <= 0) {
        cerr <<"Dataset points can't be <= 0 " << endl;
        return argumentError(__FUNCTION__);
    }

    for (int i = 0; i < C.size(); ++i) {
        for (int j = 0; j < T.size(); ++j) {
            int nij = getNij(T.at(j), C.at(i));
            if (nij < 0)
                return functionError(__FUNCTION__);
            pij[i][j] = (double)nij/n;
        }
    }
    return 0;
}

int getPCi(vectorSet2D &C, int n, double *pCi) {
    if (C.empty()) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }
    if (n <= 0) {
        cerr <<"Dataset points can't be <= 0 " << endl;
        return argumentError(__FUNCTION__);
    }

    for (int i = 0; i < C.size(); ++i) {
        pCi[i] = (double) C.at(i).size()/n;
    }
    return 0;
}
