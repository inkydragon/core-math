/* Correctly rounded exponential function for binary64 values.

Copyright (c) 2022 Alexei Sibidanov.

This file is part of the CORE-MATH project
(https://core-math.gitlabpages.inria.fr/).

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdint.h>
#include <x86intrin.h>

typedef uint64_t u64;
typedef union {double f; u64 u;} b64u64_u;

static inline double muldd(double xh, double xl, double ch, double cl, double *l){
  double ahlh = ch*xl, alhh = cl*xh, ahhh = ch*xh, ahhl = __builtin_fma(ch, xh, -ahhh);
  ahhl += alhh + ahlh;
  ch = ahhh + ahhl;
  *l = (ahhh - ch) + ahhl;
  return ch;
}

static inline double opolydd(double xh, double xl, int n, const double c[][2], double *l){
  int i = n-1;
  double ch = c[i][0], cl = c[i][1];
  while(--i>=0){
    ch = muldd(xh,xl,ch,cl,&cl);
    double th = ch + c[i][0], tl = (c[i][0] - th) + ch;
    ch = th;
    cl += tl + c[i][1];
  }
  *l = cl;
  return ch;
}

static double __attribute__((noinline)) as_exp_accurate(double x, double t, double th, double tl, double *l){
  static const double ch[][2] = 
    {{0x1p+0, 0x1.6c16bd194535dp-94}, {0x1p-1, -0x1.8259d904fd34fp-93},
     {0x1.5555555555555p-3, 0x1.53e93e9f26e62p-57}, {0x1.5555555555555p-5, 0x1.58e26c116d8a7p-59},
     {0x1.11111113e93e9p-7, -0x1.2b8b0942cd7a5p-61}, {0x1.6c16c169400a7p-10, -0x1.7e96f34fa93d8p-65}};
  const double l2h = 0x1.62e42ffp-13, l2l = 0x1.718432a1b0e26p-47, l2ll = 0x1.9ff0342542fc3p-102;
  double dx = x - l2h*t, dxl = l2l*t, dxll = l2ll*t + __builtin_fma(l2l,t,-dxl);
  double dxh = dx + dxl; dxl = (dx - dxh) + dxl + dxll;
  double fl, fh = opolydd(dxh,dxl,6,ch, &fl);
  fh = muldd(dxh,dxl,fh,fl,&fl);
  fh = muldd(th,tl,fh,fl,&fl);
  double zh = th + fh, zl = (th-zh) + fh;
  double uh = zh + tl, ul = ((zh-uh) + tl) + zl;
  double vh = uh + fl, vl = ((uh-vh) + fl) + ul;
  *l = vl;
  return vh;
}

static __attribute__((noinline)) double as_exp_database(double x, double f){
  static const double db[][3] = {
    {-0x1.981587ad4542fp+7, 0x1.8c0d4140c77b7p+0, 0x1.8p-105}, 
    {-0x1.59f038076039cp+6, 0x1.2c0fa76a0e15fp+0, 0x1.8p-54}, 
    {-0x1.8f80e06f3a04cp+4, 0x1.f80aafa92b498p+0, -0x1p-106}, 
    {-0x1.78ba0840d79e4p+4, 0x1.05960be55d55ap+0, 0x1p-105}, 
    {-0x1.e8bdbfcd9144ep+3, 0x1.f3e558cf4de54p+0, 0x1p-107}, 
    {-0x1.02393d5976769p+1, 0x1.1064b2c103ddbp+0, -0x1.8p-54},
    {-0x1.cc37ef7de7501p+0, 0x1.534d4de870713p+0, 0x1p-106}, 
    {-0x1.2a9cad9998262p+0, 0x1.3ef1e9b3a81c8p+0, -0x1p-106}, 
    {-0x1.0a54d87783d6fp+0, 0x1.69cef05657108p+0, 0x1p-105}, 
    {-0x1.bdc7955d1482cp-1, 0x1.acb8cf13bc769p+0, -0x1p-54},
    {-0x1.22e24fa3d5cf9p-1, 0x1.2217147b85eaap+0, -0x1.8p-54}, 
    {-0x1.ea16274b0109bp-3, 0x1.9309142b73ea6p+0, -0x1p-105}, 
    {-0x1.d3f3799439415p-3, 0x1.976a4c9985f5bp+0, 0x1.8p-54}, 
    {-0x1.8aeb636f3ce35p-3, 0x1.a634ae87df6aep+0, 0x1.8p-109}, 
    {-0x1.290ea09e36479p-3, 0x1.baded30cbf1c4p+0, -0x1p-110}, 
    {-0x1.daf693d64fadap-4, 0x1.c7f14af0a08ebp+0, -0x1p-108}, 
    {-0x1.bd44fdaed819fp-4, 0x1.cb4287f11060ap+0, +0x1p-54},
    {-0x1.a8f783d749a8fp-4, 0x1.cd8abd4de5c33p+0, 0x1.8p-54}, 
    {-0x1.a4187f2ca71f9p-6, 0x1.f309f46111221p+0, -0x1.8p-54}, 
    {-0x1.8c56ff5326197p-6, 0x1.f3c35328f1d5dp+0, -0x1p-54},
    {-0x1.5c5ed0ec83666p-6, 0x1.f53a751d7db49p+0, 0x1.8p-105}, 
    {-0x1.54511e930898cp-7, 0x1.fab5c6e464e0dp+0, 0x1.8p-54}, 
    {-0x1.381126525f9d9p-7, 0x1.fb25a83c4532p+0, 0x1p-105}, 
    {-0x1.1ff9b8e8b38bep-7, 0x1.fb85251a3f26fp+0, 0x1p-108}, 
    {-0x1.7fb235d76cce7p-8, 0x1.fd02d98c24bbbp+0, -0x1.8p-54}, 
    {-0x1.d792b60084f92p-11, 0x1.ff8a28e429f33p+0, -0x1p-105}, 
    {-0x1.ceff32831e2c2p-12, 0x1.ffc6235eee28dp+0, -0x1.8p-54}, 
    {-0x1.ce3f638d0c742p-12, 0x1.ffc63b5617d47p+0, 0x1.8p-105}, 
    {-0x1.09a285e42a59bp-12, 0x1.ffdeccc2df65cp+0,  0x1p-55},
    {-0x1.a2fefefd580dfp-13, 0x1.ffe5d0bb7eabfp+0, 0x1.8p-111}, 
    {-0x1.664716b68a409p-14, 0x1.fff4cde6a0bfbp+0, -0x1.8p-54},
    {-0x1.0401ae48409b5p-28, 0x1.ffffffdf7fca3p-1, 0x1.8p-55}, 
    {-0x1.ff171507f8ba5p-30, 0x1.fffffff007475p-1, 0x1.8p-55}, 
    {-0x1.85e60704a3a9cp-30, 0x1.fffffff3d0cfdp-1, -0x1.8p-55}, 
    {-0x1.4bd46601ae1efp-31, 0x1.fffffffad0ae7p-1, -0x1.8p-55}, 
    {-0x1.000000000002p-46, 0x1.fffffffffff8p-1, 0x1.8p-139}, 
    {-0x1.c000000000031p-47, 0x1.fffffffffff9p-1, 0x1p-139}, 
    {-0x1.8000000000024p-47, 0x1.fffffffffffap-1, 0x1p-140}, 
    {-0x1.4000000000019p-47, 0x1.fffffffffffbp-1, 0x1.8p-141}, 
    {-0x1.000000000001p-47, 0x1.fffffffffffcp-1, 0x1.8p-142}, 
    {-0x1.8000000000012p-48, 0x1.fffffffffffdp-1, 0x1p-143}, 
    {-0x1.0000000000008p-48, 0x1.fffffffffffep-1, 0x1.8p-145}, 
    {-0x1.8000000000009p-49, 0x1.fffffffffffe8p-1, 0x1p-146}, 
    {-0x1.0000000000004p-49, 0x1.ffffffffffffp-1, 0x1.8p-148}, 
    {-0x1.0000000000002p-50, 0x1.ffffffffffff8p-1, 0x1.8p-151}, 
    {-0x1.0000000000001p-51, 0x1.ffffffffffffcp-1, 0x1.8p-154}, 
    {0x1p-53, 0x1.0000000000001p+0, -0x1.8p-54}, 
    {0x1.fffffffffffffp-53, 0x1.0000000000001p+0, -0x1.8p-158}, 
    {0x1.ffffffffffffep-52, 0x1.0000000000002p+0, -0x1.8p-155}, 
    {0x1.ffffffffffffcp-51, 0x1.0000000000004p+0, -0x1.8p-152}, 
    {0x1.ffffffffffff8p-50, 0x1.0000000000008p+0, -0x1.8p-149}, 
    {0x1.7fffffffffff7p-49, 0x1.000000000000cp+0, -0x1p-147}, 
    {0x1.ffffffffffffp-49, 0x1.000000000001p+0, -0x1.8p-146}, 
    {0x1.7ffffffffffeep-48, 0x1.0000000000018p+0, -0x1p-144}, 
    {0x1.fffffffffffep-48, 0x1.000000000002p+0, -0x1.8p-143}, 
    {0x1.3ffffffffffe7p-47, 0x1.0000000000028p+0, -0x1.8p-142}, 
    {0x1.7ffffffffffdcp-47, 0x1.000000000003p+0, -0x1p-141}, 
    {0x1.bffffffffffcfp-47, 0x1.0000000000038p+0, -0x1p-140}, 
    {0x1.fffffffffffcp-47, 0x1.000000000004p+0, -0x1.8p-140}, 
    {0x1.ffff7fffe0001p-36, 0x1.000000001ffffp+0, 0x1.8p-54}, 
    {0x1.8cd99fffb319cp-34, 0x1.0000000063367p+0, -0x1.8p-54}, 
    {0x1.413f33fe6ce07p-31, 0x1.00000002827e7p+0, -0x1.8p-54}, 
    {0x1.f8e165f8388f7p-30, 0x1.00000007e3859p+0, 0x1.8p-54}, 
    {0x1.282aae3169ee7p-20, 0x1.00001282ab8e7p+0, -0x1.8p-54}, 
    {0x1.e5ba92aa9b1efp-20, 0x1.00001e5baaf77p+0, 0x1.8p-54}, 
    {0x1.7486e56dbca44p-14, 0x1.0005d22c869a6p+0, 0x1.8p-107}, 
    {0x1.ba07d73250de7p-14, 0x1.0006e83736f8dp+0, -0x1.8p-54}, 
    {0x1.0fd88e5260142p-13, 0x1.00087ee8887fcp+0, 0x1.8p-106}, 
    {0x1.a02a866632735p-12, 0x1.001a03fab6f63p+0, 0x1.8p-54}, 
    {0x1.458f7365fd894p-8, 0x1.01465ece08736p+0, -0x1p-106}, 
    {0x1.6a4d1af9cc989p-8, 0x1.016b4df3299d7p+0, 0x1.8p-54}, 
    {0x1.5a75293a5dcdap-6, 0x1.05789640bc8adp+0, -0x1.8p-54}, 
    {0x1.42ea46949b3c7p-5, 0x1.0a4ae9718080cp+0, -0x1p-54},
    {0x1.777d71396a75ep-5, 0x1.0c01d404fae34p+0, 0x1.8p-105}, 
    {0x1.7c8bb0cf5d16p-5, 0x1.0c2c2efbe696p+0, 0x1.8p-105}, 
    {0x1.46370d915991bp-4, 0x1.1538ea18a4585p+0, 0x1p-106}, 
    {0x1.be2caeebfc83bp-4, 0x1.1d761d8637563p+0, 0x1.8p-106}, 
    {0x1.0948d39a41695p-3, 0x1.23677186be25p+0, 0x1p-105}, 
    {0x1.a065fefae814fp-3, 0x1.39b8021bc065dp+0, -0x1.8p-54}, 
    {0x1.f6e4c3ced7c72p-3, 0x1.47408cb9583cep+0, 0x1.8p-107}, 
    {0x1.1a0408712e00ap-2, 0x1.512b3126454f3p+0, 0x1.8p-106}, 
    {0x1.bcab27d05abdep-2, 0x1.8b367381d82f5p+0, -0x1.8p-54}, 
    {0x1.005ae04256babp-1, 0x1.a65d89abf3d1fp+0, -0x1.8p-54}, 
    {0x1.62f71c4656b61p-1, 0x1.000976581ce4ep+0, 0x1p-106}, 
    {0x1.273c188aa7b14p+2, 0x1.93295a96ec6ebp+0, -0x1.8p-54}, 
    {0x1.83d4bcdebb3f4p+2, 0x1.ac50b409c8aeep+0, 0x1p-111}, 
    {0x1.08f51434652c3p+4, 0x1.daac459b157e5p+0, 0x1p-106}, 
    {0x1.1d5c2daebe367p+4, 0x1.a8c02e974c315p+0, -0x1p-107}, 
    {0x1.c44ce0d716a1ap+4, 0x1.b890ca8637ae2p+0, -0x1.8p-108}, 
    {0x1.e07e71bfcf06fp+5, 0x1.91ec4412c344fp+0, 0x1p-109}, 
    {0x1.6474c604cc0d7p+6, 0x1.7a8f65ad009bdp+0, -0x1p-107}, 
    {0x1.54cd1fea7663ap+7, 0x1.c90810d354618p+0, 0x1p-109}, 
    {0x1.2da9e5e6af0bp+8, 0x1.27d6fe867d6f6p+0, 0x1p-105}, 
    {0x1.d6479eba7c971p+8, 0x1.62a88613629b6p+0, -0x1p-109}, 
    {0x1.eb9914d4ac1c8p+8, 0x1.2b67eff65dce8p+0, -0x1p-105}, 
  };
  int a = 0, b = sizeof(db)/sizeof(db[0]) - 1, m = (a+b)/2;
  while (a <= b) {
    if (db[m][0] < x)
      a = m + 1;
    else if (db[m][0] == x) {
      f = db[m][1] + db[m][2];
      break;
    } else
      b = m - 1;
    m = (a + b)/2;
  }
  return f;
}

double cr_exp(double x){
  static const double t0[][2] = {
    {0x0p+0, 0x1p+0}, {-0x1.19083535b085ep-56, 0x1.02c9a3e778061p+0},
    {0x1.d73e2a475b466p-55, 0x1.059b0d3158574p+0}, {0x1.186be4bb285p-57, 0x1.0874518759bc8p+0},
    {0x1.8a62e4adc610ap-54, 0x1.0b5586cf9890fp+0}, {0x1.03a1727c57b52p-59, 0x1.0e3ec32d3d1a2p+0},
    {-0x1.6c51039449b3ap-54, 0x1.11301d0125b51p+0}, {-0x1.32fbf9af1369ep-54, 0x1.1429aaea92dep+0},
    {-0x1.19041b9d78a76p-55, 0x1.172b83c7d517bp+0}, {0x1.e5b4c7b4968e4p-55, 0x1.1a35beb6fcb75p+0},
    {0x1.e016e00a2643cp-54, 0x1.1d4873168b9aap+0}, {0x1.dc775814a8494p-55, 0x1.2063b88628cd6p+0},
    {0x1.9b07eb6c70572p-54, 0x1.2387a6e756238p+0}, {0x1.2bd339940e9dap-55, 0x1.26b4565e27cddp+0},
    {0x1.612e8afad1256p-55, 0x1.29e9df51fdee1p+0}, {0x1.0024754db41d4p-54, 0x1.2d285a6e4030bp+0},
    {0x1.6f46ad23182e4p-55, 0x1.306fe0a31b715p+0}, {0x1.32721843659a6p-54, 0x1.33c08b26416ffp+0},
    {-0x1.63aeabf42eae2p-54, 0x1.371a7373aa9cbp+0}, {-0x1.5e436d661f5e2p-56, 0x1.3a7db34e59ff7p+0},
    {0x1.ada0911f09ebcp-55, 0x1.3dea64c123422p+0}, {-0x1.ef3691c309278p-58, 0x1.4160a21f72e2ap+0},
    {0x1.89b7a04ef80dp-59, 0x1.44e086061892dp+0}, {0x1.3c1a3b69062fp-56, 0x1.486a2b5c13cdp+0},
    {0x1.d4397afec42e2p-56, 0x1.4bfdad5362a27p+0}, {-0x1.4b309d25957e4p-54, 0x1.4f9b2769d2ca7p+0},
    {-0x1.07abe1db13cacp-55, 0x1.5342b569d4f82p+0}, {0x1.9bb2c011d93acp-54, 0x1.56f4736b527dap+0},
    {0x1.6324c054647acp-54, 0x1.5ab07dd485429p+0}, {0x1.ba6f93080e65ep-54, 0x1.5e76f15ad2148p+0},
    {-0x1.383c17e40b496p-54, 0x1.6247eb03a5585p+0}, {-0x1.bb60987591c34p-54, 0x1.6623882552225p+0},
    {-0x1.bdd3413b26456p-54, 0x1.6a09e667f3bcdp+0}, {-0x1.bbe3a683c88aap-57, 0x1.6dfb23c651a2fp+0},
    {-0x1.16e4786887a9ap-55, 0x1.71f75e8ec5f74p+0}, {-0x1.0245957316dd4p-54, 0x1.75feb564267c9p+0},
    {-0x1.41577ee04993p-55, 0x1.7a11473eb0187p+0}, {0x1.05d02ba15797ep-56, 0x1.7e2f336cf4e62p+0},
    {-0x1.d4c1dd41532d8p-54, 0x1.82589994cce13p+0}, {-0x1.fc6f89bd4f6bap-54, 0x1.868d99b4492edp+0},
    {0x1.6e9f156864b26p-54, 0x1.8ace5422aa0dbp+0}, {0x1.5cc13a2e3976cp-55, 0x1.8f1ae99157736p+0},
    {-0x1.75fc781b57ebcp-57, 0x1.93737b0cdc5e5p+0}, {-0x1.d185b7c1b85dp-54, 0x1.97d829fde4e5p+0},
    {0x1.c7c46b071f2bep-56, 0x1.9c49182a3f09p+0}, {-0x1.359495d1cd532p-54, 0x1.a0c667b5de565p+0},
    {-0x1.d2f6edb8d41e2p-54, 0x1.a5503b23e255dp+0}, {0x1.0fac90ef7fd32p-54, 0x1.a9e6b5579fdbfp+0},
    {0x1.7a1cd345dcc82p-54, 0x1.ae89f995ad3adp+0}, {-0x1.2805e3084d708p-57, 0x1.b33a2b84f15fbp+0},
    {-0x1.5584f7e54ac3ap-56, 0x1.b7f76f2fb5e47p+0}, {0x1.23dd07a2d9e84p-55, 0x1.bcc1e904bc1d2p+0},
    {0x1.11065895048dep-55, 0x1.c199bdd85529cp+0}, {0x1.2884dff483cacp-54, 0x1.c67f12e57d14bp+0},
    {0x1.503cbd1e949dcp-56, 0x1.cb720dcef9069p+0}, {-0x1.cbc3743797a9cp-54, 0x1.d072d4a07897cp+0},
    {0x1.2ed02d75b3706p-55, 0x1.d5818dcfba487p+0}, {0x1.c2300696db532p-54, 0x1.da9e603db3285p+0},
    {-0x1.1a5cd4f184b5cp-54, 0x1.dfc97337b9b5fp+0}, {0x1.39e8980a9cc9p-55, 0x1.e502ee78b3ff6p+0},
    {-0x1.e9c23179c2894p-54, 0x1.ea4afa2a490dap+0}, {0x1.dc7f486a4b6bp-54, 0x1.efa1bee615a27p+0},
    {0x1.9d3e12dd8a18ap-54, 0x1.f50765b6e454p+0}, {0x1.74853f3a5931ep-55, 0x1.fa7c1819e90d8p+0}
  };
  static const double t1[][2] = {
    {0x0p+0, 0x1p+0}, {0x1.ae8e38c59c72ap-54, 0x1.000b175effdc7p+0},
    {-0x1.7b5d0d58ea8f4p-58, 0x1.00162f3904052p+0}, {0x1.4115cb6b16a8ep-54, 0x1.0021478e11ce6p+0},
    {-0x1.d7c96f201bb2ep-55, 0x1.002c605e2e8cfp+0}, {0x1.84711d4c35eap-54, 0x1.003779a95f959p+0},
    {-0x1.0484245243778p-55, 0x1.0042936faa3d8p+0}, {-0x1.4b237da2025fap-54, 0x1.004dadb113dap+0},
    {-0x1.5e00e62d6b30ep-56, 0x1.0058c86da1c0ap+0}, {0x1.a1d6cedbb948p-54, 0x1.0063e3a559473p+0},
    {-0x1.4acf197a00142p-54, 0x1.006eff583fc3dp+0}, {-0x1.eaf2ea42391a6p-57, 0x1.007a1b865a8cap+0},
    {0x1.da93f90835f76p-56, 0x1.0085382faef83p+0}, {-0x1.6a79084ab093cp-55, 0x1.00905554425d4p+0},
    {0x1.86364f8fbe8f8p-54, 0x1.009b72f41a12bp+0}, {-0x1.82e8e14e3110ep-55, 0x1.00a6910f3b6fdp+0},
    {-0x1.4f6b2a7609f72p-55, 0x1.00b1afa5abcbfp+0}, {-0x1.e1a258ea8f71ap-56, 0x1.00bcceb7707ecp+0},
    {0x1.4362ca5bc26f2p-56, 0x1.00c7ee448ee02p+0}, {0x1.095a56c919d02p-54, 0x1.00d30e4d0c483p+0},
    {-0x1.406ac4e81a646p-57, 0x1.00de2ed0ee0f5p+0}, {0x1.b5a6902767e08p-54, 0x1.00e94fd0398ep+0},
    {-0x1.91b206085932p-54, 0x1.00f4714af41d3p+0}, {0x1.427068ab22306p-55, 0x1.00ff93412315cp+0},
    {0x1.c1d0660524e08p-54, 0x1.010ab5b2cbd11p+0}, {-0x1.e7bdfb3204be8p-54, 0x1.0115d89ff3a8bp+0},
    {0x1.843aa8b9cbbc6p-55, 0x1.0120fc089ff63p+0}, {-0x1.34104ee7edae8p-56, 0x1.012c1fecd613bp+0},
    {-0x1.2b6aeb6176892p-56, 0x1.0137444c9b5b5p+0}, {0x1.a8cd33b8a1bb2p-56, 0x1.01426927f5278p+0},
    {0x1.2edc08e5da99ap-56, 0x1.014d8e7ee8d2fp+0}, {0x1.57ba2dc7e0c72p-55, 0x1.0158b4517bb88p+0},
    {0x1.b61299ab8cdb8p-54, 0x1.0163da9fb3335p+0}, {-0x1.90565902c5f44p-54, 0x1.016f0169949edp+0},
    {0x1.70fc41c5c2d54p-55, 0x1.017a28af25567p+0}, {0x1.4b9a6e145d76cp-54, 0x1.018550706ab62p+0},
    {-0x1.008eff5142bfap-56, 0x1.019078ad6a19fp+0}, {-0x1.77669f033c7dep-54, 0x1.019ba16628de2p+0},
    {-0x1.09bb78eeead0ap-54, 0x1.01a6ca9aac5f3p+0}, {0x1.371231477ece6p-54, 0x1.01b1f44af9f9ep+0},
    {0x1.5e7626621eb5ap-56, 0x1.01bd1e77170b4p+0}, {-0x1.bc72b100828a4p-54, 0x1.01c8491f08f08p+0},
    {-0x1.ce39cbbab8bbep-57, 0x1.01d37442d507p+0}, {0x1.16996709da2e2p-55, 0x1.01de9fe280ac8p+0},
    {-0x1.c11f5239bf536p-55, 0x1.01e9cbfe113efp+0}, {0x1.e1d4eb5edc6b4p-55, 0x1.01f4f8958c1c6p+0},
    {-0x1.afb99946ee3fp-54, 0x1.020025a8f6a35p+0}, {-0x1.8f06d8a148a32p-54, 0x1.020b533856324p+0},
    {-0x1.2bf310fc54eb6p-55, 0x1.02168143b0281p+0}, {-0x1.c95a035eb4176p-54, 0x1.0221afcb09e3ep+0},
    {-0x1.491793e46834cp-54, 0x1.022cdece68c4fp+0}, {-0x1.3e8d0d9c4909p-56, 0x1.02380e4dd22adp+0},
    {-0x1.314aa16278aa4p-54, 0x1.02433e494b755p+0}, {0x1.48daf888e965p-55, 0x1.024e6ec0da046p+0},
    {0x1.56dc8046821f4p-55, 0x1.02599fb483385p+0}, {0x1.45b42356b9d46p-54, 0x1.0264d1244c719p+0},
    {-0x1.082ef51b61d7ep-56, 0x1.027003103b10ep+0}, {0x1.2106ed0920a34p-56, 0x1.027b357854772p+0},
    {-0x1.fd4cf26ea5d0ep-54, 0x1.0286685c9e059p+0}, {-0x1.09f8775e78084p-54, 0x1.02919bbd1d1d8p+0},
    {0x1.64cbba902ca28p-58, 0x1.029ccf99d720ap+0}, {0x1.4383ef231d206p-54, 0x1.02a803f2d170dp+0},
    {0x1.4a47a505b3a46p-54, 0x1.02b338c811703p+0}, {0x1.e47120223468p-54, 0x1.02be6e199c811p+0},
  };
  const double s = 0x1.71547652b82fep+12;
  double v0 = __builtin_fma(x, s, 0x1.8000004p+25);
  b64u64_u jt = {.f = v0};
  __m128d v = _mm_set_sd (v0);
  __m128i tt = {~((1<<27)-1l),0};
  v = _mm_and_pd(v,(__m128d)tt);
  double tn = v[0] - 0x1.8p25, t = tn;
  b64u64_u ix = {.f = x};
  u64 aix = ix.u & (~0ul>>1);
  if(__builtin_expect(aix>=0x40862e42fefa39f0ul, 0)){
    if(aix>0x7ff0000000000000ul) return x;
    if(!(ix.u>>63)) return 0x1p1023*2.0;
    if(aix>=0x40874910d52d3052ul) return 0x1.5p-1022*0x1p-55;
    if(aix>0x40874385446d71c3ul) return 0x1.5p-1022*0x1p-53;
  } else if(__builtin_expect(aix<0x3c90000000000000ul, 0)){
    return 1 + x;
  }
  long i1 = (jt.u>>27)&0x3f, i0 = (jt.u>>33)&0x3f, ie = (long)(jt.u<<13)>>52;
  double t0h = t0[i0][1], t0l = t0[i0][0];
  double t1h = t1[i1][1], t1l = t1[i1][0];
  double th = t0h*t1h, tl = t0h*t1l + t1h*t0l + __builtin_fma(t0h,t1h,-th);
  const double l2h = 0x1.62e42ffp-13, l2l = 0x1.718432a1b0e26p-47;
  double dx = (x - l2h*t) + l2l*t, dx2 = dx*dx;
  static const double ch[] = {0x1p+0, 0x1p-1, 0x1.5555555aaaaaep-3, 0x1.55555551c98cp-5};
  double p = dx*((ch[0] + dx*ch[1]) + dx2*(ch[2] + dx*ch[3]));
  double fh = th, fl = tl + th*p;
  if(__builtin_expect(ie<=-1022, 0)){
    double vh = fh + fl; fl += fh - vh;
    ie += 512;
    vh *= 0x1p-512;
    fl *= 0x1p-512;
    const double eps = 0x1.8p-575;
    b64u64_u sd = {.u = (u64)(0x3ffl+ie)<<52}, su = {.u = (u64)(0x3ffl-ie)<<52};
    double vd = vh*sd.f, vdl = (vh - vd*su.f) + fl;
    double fp = __builtin_fma(sd.f, vdl+eps, vd);
    double fm = __builtin_fma(sd.f, vdl-eps, vd);
    if(__builtin_expect( fp != fm, 0)) {
      fh = as_exp_accurate(x,t,th,tl,&fl);
      vh = fh + fl; fl += fh-vh;
      vh *= 0x1p-512;
      fl *= 0x1p-512;
      vd = vh*sd.f, vdl = (vh - vd*su.f) + fl;
      fp = __builtin_fma(sd.f, vdl, vd);
    }
    fh = fp;
  } else {
    double eps = 1.64e-19, fp = fh + (fl + eps), fm = fh + (fl - eps);
    double vh = fh + fl;
    if(__builtin_expect( fp != fm, 0)){
      fh = as_exp_accurate(x,t,th,tl,&fl);
      vh = fh + fl; fl += fh-vh;
      double afl = __builtin_fabs(fl);
      if(__builtin_fabs(afl-0x1p-53) < 0x1p-105 || __builtin_fabs(afl - 0x1p-52)< 0x1p-102 ||
	 afl<0x1p-102 || afl == 0x1p-54) vh = as_exp_database(x, vh);
    }
    fh = vh;
    __m128i sb; sb[0] = ie<<52;
    __m128d r = _mm_set_sd (fh);
    r = (__m128d)_mm_add_epi64((__m128i)r, sb);
    fh = r[0];
  }
  return fh;
}
