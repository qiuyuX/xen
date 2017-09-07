/*
 * This file implements the d*-private algorithm which adds noise to the 
 * sensitive data to defeat side-channel attack. The scale of the noise 
 * increases as the number of query increases. The recurrence for adding 
 * noise is defined as following:
 * 
 * X'[i] = X'[G(i)] + (X[i] - X[G(i)]) + r_i
 *
 * G(i) = 0                               : if i = 1
 *      = i / 2 	                  : if i = D(i) and i >= 2
 *	= i - D(i)                        : if i > D(i)
 * 
 * r_i  = Lap(1 / e)                      : if i = D(i)
 *      = Lap(floor(log2(i)) / e)         : otherwise
 *
 * D(i) is defined as the largest power of 2 which devides i. For example, 
 * let i be 10110 (binary representation), then D(i) is 10.
 *
*/

#include <xen/dp/dprivacy.h>
#include <xen/dp/laplace.h>
#include <xen/random.h>

int rdtsc_count;
int times;

void dp_initialize(struct dp_struct *dp, float32_t e)
{
	float32_t one_f;
	
	one_f = ui32_to_f32(1);
	dp->base = 0;
	dp->base_n = ui32_to_f64(0);
	dp->index = 0;
	dp->e_inverse = f32_div(one_f, e);
	dp->previous = 0;
	dp->tsc_count = 0;
	
	rdtsc_count = 0;
	times = 0;
}

void dp_refresh(struct dp_struct *dp)
{
	dp->base = 0;
	dp->base_n = ui32_to_f64(0);
	dp->index = 0;
	dp->previous = 0;
	
	rdtsc_count = 0;
	times++;
}

// update the new bit length of the index. otherwise return -1.
static int is_new_length(uint32_t index)
{
        int bit_len = -1;
        int j = 1;
        int i;

        for (i = 0; i < sizeof(uint32_t) * 8; i++) {
                if (j & index) {
                        if (bit_len == -1) 
				bit_len = i + 1;
                        else 
				return -1;
                }
                j = j << 1;
        }

        return bit_len;
}

uint64_t dp_add_noise(struct dp_struct *dp, uint64_t value)
{
	int bit_len;
	int64_t result;
	float64_t value_n; /* output after adding noise */
	float64_t tmp_f64;
	float32_t tmp_f32, zero_f;

	zero_f = ui32_to_f32(0);

        /*
	if (rdtsc_count >= 1800150)
		dp_refresh(dp);
        */

	dp->index++;

	bit_len = is_new_length(dp->index);
	
	if (bit_len > 0) {
		dp->bit_len = bit_len;
		tmp_f64 = ui64_to_f64(value - dp->base); 
		dp->base_n = f64_add(dp->base_n, tmp_f64);
		tmp_f32 = get_fast_laplace(zero_f, dp->e_inverse);
		tmp_f64 = f32_to_f64(tmp_f32);
		dp->base_n = f64_add(dp->base_n, tmp_f64);
		dp->base = value;
		value_n.v = dp->base_n.v;
	}
	else {
		int i, j;
		uint64_t sum;
		/* lsbth right most bit of index is the first 1 appears in right */  
		int lsb = -1; 
		sum = 0;
		j = 1;
		
		// calculate the new difference for the difference vector
		for (i = 0; i < dp->bit_len - 1; i++) {
			if (dp->index & j) {	
				if (lsb == -1)
					lsb = i;
				else 
					sum += dp->diff[i];
			}		
			j = j << 1;
		}

		dp->diff[lsb] = value - dp->base - sum;

		// sample a laplace noise and calculate a noised difference vector	
		tmp_f32 = ui32_to_f32(dp->bit_len - 1);
		tmp_f32 = f32_mul(tmp_f32, dp->e_inverse);
		tmp_f32 = get_fast_laplace(zero_f, tmp_f32);
		dp->diff_n[lsb] = f32_to_f64(tmp_f32);
		tmp_f64 = ui64_to_f64(dp->diff[lsb]);
		dp->diff_n[lsb] = f64_add(dp->diff_n[lsb], tmp_f64);
		
		// calculate the noised output	
		value_n.v = dp->base_n.v;	
		j = 1;
		
		for (i = 0; i < dp->bit_len - 1; i++) {
			if (dp->index & j)
				value_n = f64_add(value_n, dp->diff_n[i]);
			j = j << 1;
		}
	}
	
	result = f64_to_i64_r_minMag(value_n, true); 	

	if (result < dp->previous)
		return dp->previous;
	else {
		dp->previous = result;
		return dp->previous;
	}
}

uint64_t uniform_add_noise(struct dp_struct *dp, uint64_t value)
{
	float32_t tmp1_f32, tmp2_f32;
	uint64_t result;
	uint32_t scale = 2000;
	uint32_t half_scale = 1000;
	uint32_t max = 0xFFFFFFFF;
	uint32_t ran = get_random();

	tmp1_f32 = ui32_to_f32(ran);
	tmp2_f32 = ui32_to_f32(max);
	tmp1_f32 = f32_div(tmp1_f32, tmp2_f32);
	tmp2_f32 = ui32_to_f32(scale);
	tmp1_f32 = f32_mul(tmp1_f32, tmp2_f32);

	result = f32_to_ui64_r_minMag(tmp1_f32, true);
	result += value;
	result -= half_scale;

	if (result < dp->previous)
		return dp->previous;
	else {
		dp->previous = result;
		return dp->previous;
	}
}

uint64_t dlone_add_noise(struct dp_struct *dp, uint64_t value)
{
	int64_t result;
	float64_t tmp_f64, result_f64;
	float32_t zero_f, tmp_f32;

	zero_f = ui32_to_f32(0);			
//	tmp_f32 = ui32_to_f32(1200);			
	tmp_f32 = ui32_to_f32(4000);			
	tmp_f32 = get_fast_laplace(zero_f, tmp_f32);
	tmp_f64 = f32_to_f64(tmp_f32);
	result_f64 = ui64_to_f64(value); 
	result_f64 = f64_add(result_f64, tmp_f64);
	result = f64_to_i64_r_minMag(result_f64, true); 	
	
	if (result < dp->previous)
		return dp->previous;
	else {
		dp->previous = (uint64_t)result;
		return dp->previous;
	}
}

uint64_t resolution_add_noise(struct dp_struct *dp, uint64_t value)
{
	int64_t result;
	float64_t val_f64, tmp_f64;
	
	val_f64 = ui64_to_f64(value);
	tmp_f64 = ui64_to_f64(5000);
	tmp_f64 = f64_div(val_f64, tmp_f64);
	result = f64_to_i64_r_minMag(tmp_f64, true); 	
	result = result * 5000;

	return (uint64_t)result;
}
