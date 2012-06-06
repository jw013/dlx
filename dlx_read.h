#ifndef DLX_READ_H
#define DLX_READ_H

#include "dlx.h"

/**
 * Read sparse binary matrix from text stream, and construct a DLX matrix.
 * See below for details on the input format.
 *
 * @param dlx		output dlx matrix.  The function will perform the
 * 			necessary memory allocation for every pointer field in
 * 			dlx.  The caller is responsible for free()ing the same
 * 			pointer fields.  The row_off member contains offsets to
 * 			dlx_node structures in the nodes buffer.  The row id
 * 			for each node is set to point to its row_off entry.
 * @param stream	text stream representation of the input matrix.
 * 			The only characters allowed are '0', '1', and newline.
 * 			Each newline terminates the preceding row.
 *
 * @return 0 for success, -1 for memory allocation error, -2 for malformed
 * 		input, -3 for I/O error.
 */
int dlx_read_matrix(struct dlx_matrix *dlx, FILE *stream);

#endif
