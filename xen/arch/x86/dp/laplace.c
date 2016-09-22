#include <xen/dp/laplace.h>
#include <xen/random.h>

#define FLOAT_TO_UINT32(f_v, u_v) \
        float_data.f = f_v; \
        u_v = float_data.u

#define UINT32_TO_FLOAT(u_v, f_v) \
        float_data.u = u_v; \
        f_v = float_data.f

union { 
        float f;
        uint32_t u;
} float_data;

uint32_t DP_MAX = 0xFFFFFFFF;

static float32_t fastlog(uint32_t x)
{
	float32_t y_f32, mx_f32, tmp_f32;
	
	y_f32 = ui32_to_f32(x);
	mx_f32.v = (x & 0x007FFFFF) | 0x3F000000;
	
	// y *= 1.1920928955078125e-7f;
	FLOAT_TO_UINT32(1.1920928955078125e-7f, tmp_f32.v);
	y_f32 = f32_mul(y_f32, tmp_f32);

	// y -= 124.22551499f;
	FLOAT_TO_UINT32(124.22551499f, tmp_f32.v);
	y_f32 = f32_sub(y_f32, tmp_f32);

	// y -= 1.498030302f * mx.f
	FLOAT_TO_UINT32(1.498030302f, tmp_f32.v);
	tmp_f32 = f32_mul(mx_f32, tmp_f32);
	y_f32 = f32_sub(y_f32, tmp_f32);
	
	// y -= 1.72587999f / (0.3520887068f + mx.f)
	FLOAT_TO_UINT32(0.3520887068f, tmp_f32.v);
	mx_f32 = f32_add(tmp_f32, mx_f32);
	FLOAT_TO_UINT32(1.72587999f, tmp_f32.v);
	tmp_f32 = f32_div(tmp_f32, mx_f32);
	y_f32 = f32_sub(y_f32, tmp_f32);

	// y *= 0.69314718f 
	FLOAT_TO_UINT32(0.69314718f, tmp_f32.v);
	y_f32 = f32_mul(y_f32, tmp_f32);

	return y_f32;
}

float32_t get_fast_laplace(float32_t u_f32, float32_t b_f32)
{
	/*
	 * generating the laplace variables. 
	 * pdf: f(x) = (1/2b) * exp(-|x - u|/b) 
	 * x_x
	 */
	uint32_t ran = 0;
	float32_t e_f32, tmp_f32, zero_f32;
	FLOAT_TO_UINT32(0, zero_f32.v);
		
	while(ran == 0 || ran == DP_MAX) {
		ran = get_random();
	}

	e_f32 = ui32_to_f32(ran);
	tmp_f32 = ui32_to_f32(DP_MAX);
	e_f32 = f32_div(e_f32, tmp_f32);
	FLOAT_TO_UINT32(0.5, tmp_f32.v);
	e_f32 = f32_sub(e_f32, tmp_f32);

	/* X = u - (b * sgn(U) * ln(1 - 2|U|)) */
	if(f32_lt(zero_f32, e_f32)) { // e > 0
		// result = u - b * fastlog(1 - 2 * e)
		FLOAT_TO_UINT32(2, tmp_f32.v);
		e_f32 = f32_mul(tmp_f32, e_f32);
		FLOAT_TO_UINT32(1, tmp_f32.v);
		e_f32 = f32_sub(tmp_f32, e_f32);

		e_f32 = fastlog(e_f32.v);
		e_f32 = f32_mul(b_f32, e_f32);
		e_f32 = f32_sub(u_f32, e_f32);
	}
	else if(f32_eq(zero_f32, e_f32)) { // e == 0
		e_f32.v = u_f32.v;
	}
	else {
		// result = u + b * fastlog(1 + 2 * e)
		FLOAT_TO_UINT32(2, tmp_f32.v);
		e_f32 = f32_mul(tmp_f32, e_f32);
		FLOAT_TO_UINT32(1, tmp_f32.v);
		e_f32 = f32_add(tmp_f32, e_f32);

		e_f32 = fastlog(e_f32.v);
		e_f32 = f32_mul(b_f32, e_f32);
		e_f32 = f32_add(u_f32, e_f32);
	}

	return e_f32;
}
