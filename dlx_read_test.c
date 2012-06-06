#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "dlx_read.h"

static void
die(int e, char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	exit(e);
}

/**
 * Read matrix from stdin, and print its dimensions and the 0-indexed row
 * numbers of the dlx solution, if any.
 */
int
main(void)
{
	struct dlx_matrix dlx;
	struct dlx_srow *solutions;
	size_t n;
	size_t i = 1;
	ptrdiff_t rid;

	switch (dlx_read_matrix(&dlx, stdin)) {
	case -1: die(-1, "memory allocation error\n"); break;
	case -2: die(-1, "invalid input\n"); break;
	case -3: die(-1, "I/O error\n"); break;
	default: break;
	}

	printf("Dimensions: [%zu, %zu]\n", dlx.n_row, dlx.n_col);

	if ((solutions = calloc(dlx.n_row, sizeof(*solutions))) == NULL)
		die(-1, "memory allocation\n");

	if ((n = dlx_exact_cover(solutions, &dlx.root, 0, &i)) == 0 &&
			dlx.n_col != 0)
		die(-1, "no solution found\n");

	/* print row indices */
	for (i = 0; i < n; i++) {
		rid = (size_t *) solutions[i].row_node->row_id - dlx.row_off;
		printf("%ld", (long) rid);
		if (i + 1 != n) fputs(",", stdout);
	}
	puts("");

	return 0;
}
