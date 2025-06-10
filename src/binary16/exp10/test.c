#include <stdint.h>
#include <stdio.h>

typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;

_Float16 cr_exp10f16(_Float16 x){
	uint16_t x0 = 0xc786; // binary representation of x0 in order to compare uint16_t rather than flt16
 	_Float16 x1 = 0x1.34p2; // largest _Float16 such that 10^x1 <= MAX_FLOAT16 < 10^x1+
	static const float tb[] = // tabulate value of 10^(i/2^6) for i in [-2^5*log(2)/log(10), 2^5*log(2)/log(10)], size(tb) = 39
		{0x8.13b01p-4, 0x8.5f6eep-4, 0x8.adf41p-4, 0x8.ff59ap-4,
		 0x9.53ba8p-4, 0x9.ab32bp-4, 0xa.05df3p-4, 0xa.63dep-4,  
		 0xa.c54e5p-4, 0xb.2a506p-4, 0xb.9305cp-4, 0xb.ff911p-4,  
		 0xc.70165p-4, 0xc.e4bacp-4, 0xd.5da52p-4, 0xd.dafd7p-4,  
		 0xe.5ced3p-4, 0xe.e39f8p-4, 0xf.6f41p-4, 0x1p+0,  
		 0x1.0960c6p+0, 0x1.13198p+0, 0x1.1d2d64p+0, 0x1.279fcap+0,  
		 0x1.32742ap+0, 0x1.3dae18p+0, 0x1.49515p+0, 0x1.5561aap+0,  
		 0x1.61e326p+0, 0x1.6ed9eap+0, 0x1.7c4a4p+0, 0x1.8a38ap+0,  
		 0x1.98a9a4p+0, 0x1.a7a218p+0, 0x1.b726fp+0, 0x1.c73d52p+0,  
		 0x1.d7ea92p+0, 0x1.e93436p+0, 0x1.fb1ffcp+0};
	b16u16_u v = {.f = x};
	if (v.u > x0) return (_Float16) 0x1p-25f;
	else if (x > x1) return (_Float16) 0x1.ffcp15f + 0x1p4f; 
	else {
		float log10_on_log2 = 0x1.a934fp1f;
		float minus_log2_on_log10 = -0x1.344136p-2f;
		float k = __builtin_roundevenf((float) x * log10_on_log2); 
		float xp = __builtin_fmaf(minus_log2_on_log10, k, (float) x);
		int i = 0x40 * xp;
		printf("i : %d\n", i);
		float xpp = __builtin_fmaf((float) i, -0x1p-6, xp); // x = klog(2)/log(10) + i/2^6 + xpp
																		   									// So, 10^x = 2^k * 10^(i/2^6) * 10^xpp

		// result
		xpp = 1.0 + xpp * (0x1.26bb1cp1 + xpp * (0x1.53524ep1 + xpp * (0x1.046efap1 + xpp * 0x1.2b9e52p0)));
		b32u32_u w = {.f = xpp * tb[i + 19]};
    w.u += (int) k << 23;
  	return w.f;
	}
}

int main(void) {
	printf("ret : %a\n", (float) cr_exp10f16(0x1p+2));
	return 0;
}
