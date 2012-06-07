#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dlx_read.h"

/*
 * Overview
 * ========
 * Common exit codes
 * Dynamic array code.
 * Binary CSR matrix from text stream.
 * Building DLX matrx from CSR matrix.
 * Glue function: DLX from text stream
 */

enum exit_code_t {
	DLXR_ESUCCESS,
	DLXR_EALLOC,		/* memory allocation error */
	DLXR_EDATAERR,	/* malformed input */
	DLXR_EIOERR
};

/* dynamic array @{ */

/**
 * A dynamic array manager module for managing memory allocation on a pointer
 * that automatically grows as needed.  This could be in a separate file
 * ("translation unit") except it is only needed here, so there's no need to
 * complicate the compilation process.
 *
 * Only supports append and read operations.
 * Treat all members as private, except export allows access to *data.  Also,
 * the memory allocated for data is independent of the containing dynamic
 * array.  Make sure to export the pointer value before destroying the dynamic
 * array, as destroying a dynamic array does NOT free the contained data.
 */
struct size_t_darray {
	size_t sz_max;
	size_t sz_cur;
	size_t i_next;
	size_t *data;	/* export returns a copy of this */
};

/**
 * Attempt to allocate a new size_t_darray.
 *
 * @param sz_initial	number of elements to allocate space for
 * @return NULL on failure, otherwise pointer to new darray on success.
 */
static struct size_t_darray *
size_t_darray_create(size_t sz_initial)
{
	struct size_t_darray *array;

	if ((array = malloc(sizeof(*array))) == NULL)
		return NULL;
	/* rely on calloc to catch multiplication integer overflow */
	if ((array->data = calloc(sz_initial, sizeof(size_t))) == NULL) {
		free(array);
		return NULL;
	}
	array->sz_max = sz_initial;
	array->sz_cur = 0;
	array->i_next = 0;
	return array;
}

/**
 * Destroy a size_t_darray object.  Does NOT free the data field, so make
 * sure to export before calling destroy.
 */
static void
size_t_darray_destroy(struct size_t_darray *array)
{
	free(array);
}

/**
 * Private.
 * Expand allocated storage for size_t_darray.
 * @return -1 for memory allocation error, 0 for success.
 */
static int
size_t_darray_grow(struct size_t_darray *array)
{
	size_t new_size;
	void *new_ptr;

	/*
	 * new_size	= 2			(case: sz_max == 1)
	 * 		= 1.5 * array->sz_max 	(case: sz_max > 1)
	 * 		= SIZE_MAX		(case: case 2 overflows)
	 */
	new_size = array->sz_max / 2;
	if (new_size == 0) new_size++;
	if (new_size > SIZE_MAX - array->sz_max)	/* overflow test   */
		new_size = SIZE_MAX;
	else
		new_size += array->sz_max;
	if (new_size <= array->sz_cur)		/* failed to increase size */
		return -1;
	if ((new_ptr = realloc(array->data, new_size)) == NULL)
		return -1;

	array->sz_max = new_size;
	array->data = new_ptr;
	return 0;
}

/**
 * Append value to the end of the array, growing the array automatically, if
 * necessary.
 *
 * @param val		value to append
 * @return 0 on success, -1 on failure.
 */
static int
size_t_darray_append(struct size_t_darray *array, size_t val)
{
	if (array->sz_cur > array->sz_max - sizeof(val)) {
		if (size_t_darray_grow(array) != 0)
			return -1;
	}
	array->data[array->i_next++] = val;
	array->sz_cur += sizeof(val);
	return 0;
}

/**
 * Shrink allocated storage size to only the currently used portion.
 * @return 0 on success, -1 on failure.
 */
static int
size_t_darray_trim(struct size_t_darray *array)
{
	void *new_ptr;

	if ((new_ptr = realloc(array->data, array->sz_cur)) == NULL &&
			array->sz_cur != 0)
		return -1;
	array->sz_max = array->sz_cur;
	array->data = new_ptr;
	return 0;
}

/**
 * @param length	if not NULL, number of elements in array is placed here
 * @return data pointer
 */
static size_t *
size_t_darray_export(struct size_t_darray *array, size_t *length)
{
	if (length) *length = array->i_next;
	return array->data;
}

/**
 * @return number of elements in array
 */
static size_t
size_t_darray_size(struct size_t_darray *array)
{
	return array->i_next;
}

/* @} END dynamic array code */

/* csr matrix @{ */

/**
 * Compressed row storage representation of a binary sparse matrix, a.k.a
 * compressed sparse row.  Since binary means the only possible non-zero value
 * is one, the array of values is unnecessary.  See the Wikipedia article on
 * sparse matrices.[1]  The implementation limits the size of the sparse
 * representation to the range of the size_t type.
 *
 * All members are public.
 *
 * [1]: https://en.wikipedia.org/wiki/Sparse_matrix
 */
struct binary_csr_matrix {
	size_t  row_ptr_size;	/* number of rows, m + 1 */
	size_t *col_ind;	/* column indices of each non-zero value */
	size_t *row_ptr;	/* row pointers; last entry is size of col_ind */
};

/**
 * Read sparse binary matrix from text stream, and store it in CSR format.
 * See below for details on the input format.
 *
 * @param csr		output location.  The function will perform the
 * 			necessary malloc calls for the col_ind and row_ptr
 * 			fields.  The caller is responsible for free()ing the
 * 			two pointer fields.
 * @param n		if not NULL, the number of columns in the widest row
 * 			encountered will be stored here.
 * @param stream	text stream representation of the input matrix.
 * 			The only characters allowed are '0', '1', and newline.
 * 			Each newline begins a new row.
 *
 * @return 0 for success, -1 for memory allocation error, -2 for malformed
 * 		input, -3 for I/O error.
 */
static int
read_bcsr(struct binary_csr_matrix *csr, size_t *n, FILE *stream)
{
	struct size_t_darray *col_ind_darray;
	struct size_t_darray *row_ptr_darray;
	int c;
	int ret = DLXR_ESUCCESS;
	size_t max_cols = 0;	/* width of widest row so far */
	size_t col = 0;		/* number of columns in the current row so far */
	int fl_newline = 1;	/* flag: last character was a newline or
				   alternatively, this next character is the
				   first in a new line */

	/* initial sizes rather arbitrarily chosen */
	if ((col_ind_darray = size_t_darray_create(512)) == NULL ||
			(row_ptr_darray = size_t_darray_create(256)) == NULL)
		return -DLXR_EALLOC;

	/* first row always starts at the first index (0) of col_ind */
	size_t_darray_append(row_ptr_darray, 0);
	while ((c = getc(stream)) != EOF) {
		if (c == '1') {
			size_t_darray_append(col_ind_darray, col++);
			fl_newline = 0;
		} else if (c == '0') {
			col++;
			fl_newline = 0;
		} else if (c == '\n') { /* EOL */
			/* mark the end of the current row */
			size_t_darray_append(row_ptr_darray,
					size_t_darray_size(col_ind_darray));
			if (max_cols < col)
				max_cols = col;
			col = 0;
			fl_newline = 1;
		} else { /* invalid char */
			ret = -DLXR_EDATAERR;
			break;
		}
	}

	if (c == EOF) {
		if (ferror(stream)) {
			ret = -DLXR_EIOERR;
		} else if (feof(stream) && !fl_newline) {
			size_t_darray_append(row_ptr_darray,
					size_t_darray_size(col_ind_darray));
			if (max_cols < col)
				max_cols = col;
		}
		/* ignore EOF immediately following newline */
	}

	/* trim both dynamic arrays */
	if (ret == DLXR_ESUCCESS &&
			(size_t_darray_trim(col_ind_darray) != 0 ||
			 size_t_darray_trim(row_ptr_darray) != 0   )) {
		ret = -DLXR_EALLOC;
	}
	/* wrap up, free resources */
	if (ret == DLXR_ESUCCESS) {
		if (n) *n = max_cols;
		csr->col_ind = size_t_darray_export(col_ind_darray, NULL);
		csr->row_ptr = size_t_darray_export(row_ptr_darray,
				&csr->row_ptr_size);
	} else { /* error return */
		csr->col_ind = NULL;
		csr->row_ptr = NULL;
		csr->row_ptr_size = 0;
		free(size_t_darray_export(col_ind_darray, NULL));
		free(size_t_darray_export(row_ptr_darray, NULL));
	}
	size_t_darray_destroy(col_ind_darray);
	size_t_darray_destroy(row_ptr_darray);
	return ret;
}

/* @} END csr matrix code */

/* CSR to DLX @{ */

/**
 * Build a dlx_matrix from the provided binary sparse CSR matrix.
 *
 * @param csr	input binary matrix
 * @param dlx	output dlx matrix.  The function will perform the necessary
 * 		memory allocation for every pointer field in dlx.  The caller
 * 		is responsible for free()ing the same pointer fields.  The
 * 		row_off member contains offsets to dlx_node structures in the
 * 		nodes buffer.  The row id for each node is set to point to its
 * 		row_off entry.
 * @param ncols	number of columns (since csr does not contain this information)
 * @return 0 on success, -1 on failure.
 */
static int
bcsr_to_dlx(struct binary_csr_matrix *csr, struct dlx_matrix *dlx, size_t ncols)
{
	struct dlx_hnode **pheaders;
	size_t i, j;
	size_t ri_start, ri_end;

	if ((dlx->headers = calloc(ncols, sizeof *dlx->headers)) &&
			(dlx->nodes = calloc(csr->row_ptr[csr->row_ptr_size - 1],
					     sizeof *dlx->nodes)) &&
			(dlx->row_off = calloc(csr->row_ptr_size,
					       sizeof *dlx->row_off)) &&
			(pheaders = calloc(ncols, sizeof *pheaders))) {
		/* all allocations succeeded */
	} else {
		/* an allocation failed, abort */
		free(dlx->headers);
		free(dlx->nodes);
		free(dlx->row_off);
		free(pheaders);
		return -DLXR_EALLOC;
	}

	dlx_make_header_row(&dlx->root, dlx->headers, ncols);

	/* process nodes row-wise */
	dlx->row_off[0] = 0; /* first row must start at index 0 of dlx->nodes */
	for (i = 0; i < csr->row_ptr_size - 1; i++) {
		ri_start = csr->row_ptr[i];
		ri_end = csr->row_ptr[i + 1];
		for (j = ri_start; j < ri_end; j++)
			pheaders[j - ri_start] = dlx->headers + csr->col_ind[j];
		dlx_make_row(dlx->nodes + ri_start, dlx->row_off + i,
				ri_end - ri_start);
		dlx_add_row(dlx->nodes + ri_start, pheaders, ri_end - ri_start);
		dlx->row_off[i + 1] = ri_end;
	}

	dlx->n_col = ncols;
	dlx->n_row = csr->row_ptr_size - 1;

	free(pheaders);
	return DLXR_ESUCCESS;
}

/* @} END CSR to DLX */

/* see dlx_read.h */
int
dlx_read_matrix(struct dlx_matrix *dlx, FILE *stream)
{
	int ret;
	size_t ncols;
	struct binary_csr_matrix csr;

	csr.col_ind = NULL;
	csr.row_ptr = NULL;
	if ((ret = read_bcsr(&csr, &ncols, stream)) == 0) {
		ret = bcsr_to_dlx(&csr, dlx, ncols);
	}

	free(csr.col_ind);
	free(csr.row_ptr);
	return ret;
}
