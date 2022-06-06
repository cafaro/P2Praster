/** @file*/

/**
* @author Antonio Mariani
* @date 16/11/19
*/

#include <iostream>
#include "error.h"
using namespace std;

/**
 *
 * @param [in] nameFunction - Name of the calling function
 * @return -1
 */
int memoryError(string nameFunction) {
    cerr << nameFunction << ": Not enough memory" << endl;
    return -1;
}

/**
 *
 * @param [in] nameFunction - Name of the calling function
 * @return -2
 */
int insertError(string nameFunction) {
    cerr << nameFunction << ": Insert Error" << endl;
    return -2;
}

/**
 *
 * @param [in] nameFunction - Name of the calling function
 * @return -3
 */
int dataError(string nameFunction) {
    cerr << nameFunction << ": Corrupt data Error" << endl;
    return -3;
}

/**
 *
 * @param [in] nameFunction - Name of the calling function
 * @return -4
 */
int arithmeticError(string nameFunction) {
    cerr << nameFunction << ": Arithmetic Error" << endl;
    return -4;
}


/**
 *
 * @param [in] nameFunction - Name of the calling function
 * @return -7
 */
int fileError(string nameFunction) {
    cerr << nameFunction << ": Can't open file: " << endl;
    return -7;
}

/**
 *
 * @param [in] nameFunction - Name of the calling function
 * @return -8
 */
int readDatasetError(string nameFunction) {
    cerr << nameFunction << ": Can't load dataset: " << "clustered.csv" << endl;
    return -8;
}

/**
 *
 * @param [in] nameFunction - Name of the calling function
 * @return -9
 */
int partitionError(string nameFunction) {
    cerr << nameFunction << ": Error partition dataset" << endl;
    return -9;
}

/**
 *
 * @param [in] nameFunction - Name of the calling function
 * @return -10
 */
int findError(string nameFunction) {
    cerr << nameFunction << ": Find Error" << endl;
    return -10;
}

/**
 *
 * @param [in] nameFunction - Name of the calling function
 * @return -10
 */
int mergeError(string nameFunction) {
    cerr << nameFunction << ": Find Error" << endl;
    return -11;
}

/**
 *
 * @param [in] nameFunction - Name of the calling function
 * @return -12
 */
int graphError(string nameFunction) {
    cerr << nameFunction << ": Graph Error" << endl;
    return -12;
}

/**
 *
 * @param [in] nameFunction - Name of the calling function
 * @return -13
 */
int functionError(string nameFunction) {
    cerr << nameFunction << ": One of the called functions failed" << endl;
    return -13;
}

/**
 *
 * @param [in] nameFunction - Name of the calling function
 * @return -99
 */
int argumentError(string nameFunction) {
    cerr << nameFunction << ": One of the arguments is not valid" << endl;
    return -99;
}




