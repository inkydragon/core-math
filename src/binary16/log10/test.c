#include <stdint.h>
#include <stdio.h>
#include <fenv.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;
b16u16_u neginf = {.u = 0xfc00};

_Float16 cr_log10f16(_Float16 x){
	b16u16_u t = {.f = x};
	if (t.u == 0) return neginf.f;
	else if (t.u >> 10 >= 0x1f) {
		if (t.u == 0x8000) return neginf.f;
		else if (t.u >> 15) return 0.0 / 0.0;
		else return x + x;
	}
	float log10_2 = 0x1.344136p-2;
	static const float tb[] = 
		{0x0p+0, 0x3.6bd21p-8, 0x6.bd7e48p-8, 0x9.f688dp-8,  
		 0xd.1854fp-8, 0x1.02428cp-4, 0x1.31b306p-4, 0x1.5fe804p-4,  
		 0x1.8cf178p-4, 0x1.b8de4ep-4, 0x1.e3bc1ap-4, 0x2.0d97acp-4,  
		 0x2.367ce4p-4, 0x2.5e76d4p-4, 0x2.858fdp-4, 0x2.abd18cp-4,  
		 0x2.d1451p-4, 0x2.f5f2e8p-4, 0x3.19e2f2p-4, 0x3.3d1cf8p-4,  
		 0x3.5fa7c8p-4, 0x3.818a28p-4, 0x3.a2ca6p-4, 0x3.c36e68p-4,  
		 0x3.e37bdcp-4, 0x4.02f818p-4, 0x4.21e82p-4, 0x4.4050c8p-4,  
		 0x4.5e3698p-4, 0x4.7b9dep-4, 0x4.988acp-4, 0x4.b5013p-4};
	static const float tl[] =
		{0x1p+0, 0xf.83e1p-4, 0xf.0f0f1p-4, 0xe.a0ea1p-4,  
		 0xe.38e39p-4, 0xd.d67c9p-4, 0xd.79436p-4, 0xd.20cfc8p-4,
		 0xc.ccccdp-4, 0xc.7ce0cp-4, 0xc.30c31p-4, 0xb.e82fap-4,  
		 0xb.a2cbc8p-4, 0xb.60b61p-4, 0xb.21643p-4, 0xa.e4c41p-4,  
		 0xa.aaaabp-4, 0xa.72ed36p-4, 0xa.3d70ap-4, 0xa.0a0a1p-4,  
		 0x9.d89d9p-4, 0x9.a90e8p-4, 0x9.7b426p-4, 0x9.4f209p-4,  
		 0x9.24925p-4, 0x8.fb824p-4, 0x8.d3dcbp-4, 0x8.ad8f3p-4,  
		 0x8.88938ap-4, 0x8.64bfd8p-4, 0x8.4218b8p-4, 0x8.208136p-4};
	b32u32_u xf = {.f = x};
	int expo = (xf.u >> 23) - 127; // used float instead of flaot16 to avoid working with subnormalized
	int i = (xf.u & 0x007c0000) >> 18;
	printf("i : %d\n", i);
	xf.f = 0x1p-23 * (xf.u & 0x0003ffff) * tl[i];
	// We have, x = 2^expo * (1 + i2⁻⁵ + xf.f)
	// Thus, log(x) = expo log(2) + log(1 + i2⁻⁵) + log(1 + xf.f / (1 + i2⁻⁵))
	xf.f *= __builtin_fmaf(__builtin_fmaf(0x1.555554p-2, xf.f, -0.5f), xf.f, 1.0f) * 0x1.bcb7b2p-2;
	return __builtin_fmaf(log10_2, (float) expo, xf.f + tb[i]); 
}

int main(void) {
	printf("%a\n", (float) cr_log10f16(0x1.f84p+6));
	printf("expected : %a\n", 0x1.0ccp+1);
	fesetround(FE_TOWARDZERO);
	//printf("%a\n", (float) cr_log10f16(0x1.394p-1));
	//printf("expected : %a\n", -0x1.b4cp-3);
	fesetround(FE_UPWARD);
	printf("%a\n", (float) cr_log10f16(0x1.f4p+9));
	printf("expected : %a\n", 0x1.8p+1);
	//printf("%a\n", (float) cr_log10f16(0x1.388p+13));
	//printf("expected : %a\n", 0x1p+2);
}
