/**
 * @addtogroup error_pca
 * @ingroup pca
 * @{
 * @file error_pca.h
 *
 * @brief function to manage error types
 */


#ifndef ERROR_PCA_H
#define ERROR_PCA_H

/**
 * print a specific lfmm error message
 *
 * @param msg   the string to recognize the error type
                It can be "cmd", "option", "missing", "basic", 
                "specific" or internal otherwise.
 * @param file  the name of a file (for "option" and "specific" error)
 */
void print_error_pca(char* msg, char* file);

#endif // ERROR_PCA_H

/** @} */
