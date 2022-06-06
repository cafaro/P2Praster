/** @file*/

/**
* @author Antonio Mariani
* @date 27/11/19
*/

#include "raster.h"

#ifndef P2PRASTER_STATISTIC_H
#define P2PRASTER_STATISTIC_H

/**
 * This function computes the average precision of the clusters.
 * Given cluster Ci , let j_i denote the partition that contains the maxi-
 * mum number of points from Ci , i.e., j_i = max {nij } on k clusters.
 * nij is the cardinality of the intersection set between Ci and Tj
 * The precision of a cluster Ci is the same as its purity.
 * It measures the fraction of points in Ci from the majority partition Tj_i
 *
 *
 * @param [in] C - Distributed clustering
 * @param [in] T - Ground-truth clustering
 * @param [in,out] Precision - Where precision is stored
 * @return 0 if success, -1 in case of memory error, -3 in case of data error, -2 in case of function error
 */
int getPrecision(vectorSet2D &C, vectorSet2D &T, double &Precision);


/**
 * This function computes the probability that a point
 * in cluster i also belongs to partition j
 *
 * @param [in] C - Distributed clustering
 * @param [in] T - Ground-truth clustering
 * @param [in,out] pij - Variables where result is stored
 * @param [in] n - Number of dataset points
 * @return 0 if success, -13 in case of function error, -99 in case of argument error
 */
int getPij(vectorSet2D &C, vectorSet2D &T, double **pij, int n);


/**
 * This function computes the probability of cluster Ci: n_i/n
 *
 * @param [in] C - Distributed clustering
 * @param [in] n - Number of dataset points
 * @param [in,out] pCi - Variable where result is stored
 * @return 0 if success, -3 in case of data error, -99 in case of argument error
 */
int getPCi(vectorSet2D &C, int n, double *pCi);


/**
 * This function computes the intersection set between
 * the two clusters Ti and Ci and return its cardinality
 *
 * @param [in] Ti - Cluster from ground-truth clustering
 * @param [in] Ci - Cluster from distributed clustering
 * @return Set size (>=0) if success, -2 in case of insert error, -3 in case of data error
 */
int getNij(unSet2D &Ti, unSet2D &Ci);


/**
 * Given cluster Ci , let j_i denote the partition that contains the maxi-
 * mum number of points from Ci , i.e., j_i = max {nij } on k clusters.
 * This function returns nij_i, which is the intersection set between Ci
 * and Tj (j=1,...k) with the higher cardinality
 *
 * @param [in] T - Ground-truth clustering
 * @param [in] Ci - Cluster from distributed clustering
 * @param [in,out] max_nij
 * @return 0 if success, -3 in case of data error, -2 in case of insert error
 */
int getMax_nij(vectorSet2D &T, unSet2D &Ci, int &max_nij);


/**
 * Given cluster Ci , let j_i denote the partition that contains the maxi-
 * mum number of points from Ci , i.e., j_i = max {nij } on k clusters.
 * This function returns nij_i, which is the intersection set between Ci
 * and Tj (j=1,...k) with the higher cardinality and the Tj_i cardinality
 *
 * @param [in] T - Ground-truth clustering
 * @param [in] Ci - Cluster from distributed clustering
 * @param [in,out] max_nij - Where nij_i is stored
 * @param [in,out] mji - Where Tj_i cardinality is stored
 * @return 0 if success, -3 in case of data error, -2 in case of insert error
 */
int getMax_nij_mji(vectorSet2D &T, unSet2D &Ci, int &max_nij, int &mji);


/**
 * This function computes the average recall of the clusters.
 * Given cluster Ci , let j_i denote the partition that contains the maxi-
 * mum number of points from Ci , i.e., j_i = max {nij } on k clusters.
 * nij is the cardinality of the intersection set between Ci and T_j
 * It measures the fraction of point in partition Tji shared in common
 * with cluster Ci.
 *
 *
 * @param [in] C - Distributed clustering
 * @param [in] T - Ground-truth clustering
 * @param [in,out] Recall - Where recall is stored
 * @return 0 if success, -1 in case of memory error, -3 in case of data error, -2 in case of function error
 */
int getRecall(vectorSet2D &C, vectorSet2D &T, double &Recall);


/**
 * This function computes the F1 measure.
 * The F-measure is the harmonic mean of the precision and recall values for each cluster.
 *
 * @param [in] Recall - Variable that contains recall value
 * @param [in] Precision - Variable that contains precision value
 * @return F1 measure
 */
double getFMeasure(double Recall, double Precision);


/**
 * This function computes the conditional entropy of a clustering C.
 * The more a cluster’s members are split into different partitions, the higher
 * the conditional entropy. For a perfect clustering, the conditional entropy value
 * is zero, whereas the worst possible conditional entropy value is log k
 * (k is the number of clusters contained in T)
 *
 * @param [in] C - Distributed clustering
 * @param [in] T - Ground-truth clustering
 * @param [in] n - Number of dataset points
 * @param [in,out] Entropy - Variable where result is stored
 * @return 0 if success, -3 in case of data error, -13 on case of function error, -99 in case of argument error
 */
int getConditionalEntropy(vectorSet2D &C, vectorSet2D &T, int n, double &Entropy);

#endif //P2PRASTER_STATISTIC_H
