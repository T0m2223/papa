#include <algorithm>
#include <omp.h>
#include "helpers.hpp"

unsigned long SequenceInfo::gpsa_sequential(float** S) {
    unsigned long visited = 0;

	// Boundary
    for (unsigned int i = 1; i < rows; i++) {
        S[i][0] = i * gap_penalty;
		visited++;
	}

    for (unsigned int j = 0; j < cols; j++) {
        S[0][j] = j * gap_penalty;
		visited++;
	}
	
	// Main part
	for (unsigned int i = 1; i < rows; i++) {
		for (unsigned int j = 1; j < cols; j++) {
			float match = S[i - 1][j - 1] + (X[i - 1] == Y[j - 1] ? match_score : mismatch_score);
			float del = S[i - 1][j] + gap_penalty;
			float insert = S[i][j - 1] + gap_penalty;
			S[i][j] = std::max({match, del, insert});
		
			visited++;
		}
	}

    return visited;
}

unsigned long SequenceInfo::gpsa_taskloop(float** S, long grain_size=1, int block_size_x=1, int block_size_y=1)
{
	unsigned long visited = 0;

	// copy everything local
	const long  r  = rows;
	const long  c  = cols;
	const float gp = gap_penalty;
	const float ms = match_score;
	const float mm = mismatch_score;
	const char* xp = X.data();
	const char* yp = Y.data();
	const long  bx = block_size_x;
	const long  by = block_size_y;

	const long num_blocks_x   = (r - 1 + bx - 1) / bx;
	const long num_blocks_y   = (c - 1 + by - 1) / by;
	const long num_anti_diags = num_blocks_x + num_blocks_y - 1;

	// configuration vars are firstprivate so every task works on a private copy
	#pragma omp parallel shared(S, visited, xp, yp) \
		firstprivate(r, c, gp, ms, mm, bx, by, num_blocks_x, num_blocks_y, num_anti_diags)
	#pragma omp single
	{
		// Boundary
		#pragma omp taskloop shared(S) firstprivate(r, gp) reduction(+:visited) grainsize(grain_size)
		for (long i = 1; i < r; i++)
		{
			S[i][0] = i * gp;
			visited++;
		}

		#pragma omp taskloop shared(S) firstprivate(c, gp) reduction(+:visited) grainsize(grain_size)
		for (long j = 0; j < c; j++)
		{
			S[0][j] = j * gp;
			visited++;
		}

		// Main part: one taskloop per anti-diagonal
		for (long anti_diag = 0; anti_diag < num_anti_diags; ++anti_diag)
		{
			const long first_block = std::max(0l, anti_diag - num_blocks_y + 1);
			const long last_block  = std::min(anti_diag, num_blocks_x - 1);

			#pragma omp taskloop shared(S, xp, yp) \
				firstprivate(r, c, gp, ms, mm, bx, by, anti_diag, first_block, last_block) \
				reduction(+:visited) grainsize(grain_size)
			for (long block = first_block; block <= last_block; ++block)
			{
				const long i0 = 1 + block * bx;
				const long i1 = std::min(i0 + bx, r);
				const long j0 = 1 + (anti_diag - block) * by;
				const long j1 = std::min(j0 + by, c);
				for (long i = i0; i < i1; ++i)
				{
					for (long j = j0; j < j1; ++j)
					{
						float match  = S[i - 1][j - 1] + (xp[i - 1] == yp[j - 1] ? ms : mm);
						float del    = S[i - 1][j] + gp;
						float insert = S[i][j - 1] + gp;
						S[i][j] = std::max({match, del, insert});
						visited++;
					}
				}
			}
		}
	}

    return visited;
}

unsigned long SequenceInfo::gpsa_tasks(float** S, long grain_size=1, int block_size_x=1, int block_size_y=1) {
	unsigned long visited = 0;

	const long  r  = rows;
	const long  c  = cols;
	const float gp = gap_penalty;
	const float ms = match_score;
	const float mm = mismatch_score;
	const char* xp = X.data();
	const char* yp = Y.data();
	const long  bx = block_size_x;
	const long  by = block_size_y;

	const long num_blocks_x = (r - 1 + bx - 1) / bx;
	const long num_blocks_y = (c - 1 + by - 1) / by;

	#pragma omp parallel shared(S, visited, xp, yp) \
		firstprivate(r, c, gp, ms, mm, bx, by, num_blocks_x, num_blocks_y)
	#pragma omp single
	{
		// Boundary
		#pragma omp taskgroup task_reduction(+:visited)
		{
			#pragma omp task shared(S) firstprivate(r, gp) in_reduction(+:visited)
			for (long i = 1; i < r; ++i)
			{
				S[i][0] = i * gp;
				visited++;
			}

			#pragma omp task shared(S) firstprivate(c, gp) in_reduction(+:visited)
			for (long j = 0; j < c; ++j)
			{
				S[0][j] = j * gp;
				visited++;
			}
		}

		#pragma omp taskgroup task_reduction(+:visited)
		{
			for (long p = 0; p < num_blocks_x; ++p)
			{
				const long i0 = 1 + p * bx;
				const long i1 = std::min(i0 + bx, r);

				for (long q = 0; q < num_blocks_y; ++q)
				{
					const long j0 = 1 + q * by;
					const long j1 = std::min(j0 + by, c);

					// missing neighbours use S[0][0], which no block task writes,
					// so it is a harmless read-only token.
					const long up_i   = (p > 0) ? i0 - bx : 0;
					const long up_j   = (p > 0) ? j0      : 0;
					const long left_i = (q > 0) ? i0      : 0;
					const long left_j = (q > 0) ? j0 - by : 0;

					// depend on the block above and on the block to the left
					#pragma omp task shared(S, xp, yp) \
						firstprivate(r, c, gp, ms, mm, bx, by, i0, i1, j0, j1) \
						depend(inout: S[i0][j0]) \
						depend(in: S[up_i][up_j], S[left_i][left_j]) \
						in_reduction(+:visited)
					{
						for (long i = i0; i < i1; ++i)
						{
							for (long j = j0; j < j1; ++j)
							{
								float match  = S[i - 1][j - 1] + (xp[i - 1] == yp[j - 1] ? ms : mm);
								float del    = S[i - 1][j] + gp;
								float insert = S[i][j - 1] + gp;
								S[i][j] = std::max({match, del, insert});
								visited++;
							}
						}
					}
				}
			}
		}
	}

    return visited;
}
