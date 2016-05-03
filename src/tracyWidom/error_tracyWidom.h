/**
 * @addtogroup error_tracyWidom
 * @ingroup tracyWidom
 * @{
 * @file error_tracyWidom.h
 *
 * @brief function to manage error types
 */

#ifndef ERROR_TRACYWIDOM_H
#define ERROR_TRACYWIDOM_H

/**
 * print a specific tracyWidom error message
 *
 * @param msg   the string to recognize the error type. 
		It can be "cmd", "option", "missing", "basic", 
		"specific" or internal otherwise.
 * @param file  the name of a file (for "option" and "specific" error)
 */
void print_error_tracyWidom(char *msg, char *file);

#endif                          // ERROR_TRACYWIDOM_H

/** @} */
