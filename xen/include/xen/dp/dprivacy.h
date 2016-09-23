#ifndef _DPRIVACY_H
#define _DPRIVACY_H

#include <xen/float/softfloat.h>

#define MAX_INDEX_LEN 14 

/*
 * For this specific application, the original data vector contains positive data, and 
 * the values increase monotonically. So the data type for "base" and "diff" are 
 * unsigned int. 
*/

struct dp_struct {
	uint32_t base; /* X[b], b is the largest power of 2 less than index i */
	float32_t base_n; /* X'[b] */	
	/* diff[k] = X[j] - X[b], the kth least significant bit of j is 1 */
	uint32_t diff[MAX_INDEX_LEN]; 
	float32_t diff_n[MAX_INDEX_LEN]; /* diff_n[k] = diff[k] + lap */
	uint32_t index;
	uint32_t bit_len; /* the length of the binary representation of index */
	float32_t e_inverse;
};

void dp_initialize(struct dp_struct *dp, float32_t e);
int32_t dp_add_noise(struct dp_struct *dp, uint32_t value);

#endif /* _DPRIVACY_H */
