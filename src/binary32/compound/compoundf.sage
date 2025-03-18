def inv():
   Zmin = 0
   Zmax = 0
   T = []
   U = []
   maxerr = 0 # maximal absolute error on u[i] ~ -log(r[i])/log(2)
   for i in range(6,29+1):
      a = 1/2+i/32
      b = 1/2+(i+1)/32
      c = (a+b)/2
      p = 53
      while true:
         R = RealField(p)
         if a==1 or b==1: # force r=1 near 1
            r = R(1)
         else:
            r = R(1/c)
         R = r.exact_rational()
         zmin = R*a-1
         zmax = R*b-1
         u = min(RR(a).ulp(),RR(b).ulp())
         u = u*r.ulp()
         # r*x-1 is an integer multiple of u
         kmin = abs(zmin/u)
         kmax = abs(zmax/u)
         if max(kmin,kmax)>2^53:
            p = p-1
            continue
         # print (i, get_hex(r))
         Zmin = min(Zmin,zmin)
         Zmax = max(Zmax,zmax)
         break
      T.append(r)
      u = RR(n(-log(R)/log(2),100))
      err = abs(n(u.exact_rational()+log(R)/log(2),100))
      if err>maxerr:
         maxerr = err
         print (i,err)
      U.append(u)
   print ("  static const double inv[] = {")
   s = "   "
   for r in T:
      s = s + " " + get_hex(r) + ","
      if len(s)+10>=80:
         print (s)
         s = "   "
   if s!="   ":
      print (s)
   print ("  };")
   print ("  static const double log2inv[] = {")
   s = "   "
   for r in U:
      s = s + " " + get_hex(r) + ","
      if len(s)+20>=80:
         print (s)
         s = "   "
   if s!="   ":
      print (s)
   print ("  };")
   print (Zmin, Zmax, err)

def inv2():
   Zmin = 0
   Zmax = 0
   T = []
   U = []
   for i in range(46):
      a = 1/2+(i+13)/64
      b = 1/2+(i+14)/64
      c = (a+b)/2
      p = 53
      while true:
         R = RealField(p)
         if a==1 or b==1: # force r=1 near 1
            r = R(1)
         else:
            r = R(1/c)
         R = r.exact_rational()
         zmin = R*a-1
         zmax = R*b-1
         u = min(RR(a).ulp(),RR(b).ulp())
         u = u*r.ulp()
         # r*x-1 is an integer multiple of u
         kmin = abs(zmin/u)
         kmax = abs(zmax/u)
         if max(kmin,kmax)>2^53:
            p = p-1
            continue
         # print (i, get_hex(r))
         Zmin = min(Zmin,zmin)
         Zmax = max(Zmax,zmax)
         break
      T.append(r)
      h = RR(n(-log(R)/log(2),200))
      l = RR(n(-log(R)/log(2)-h.exact_rational(),200))
      U.append((h,l))
   print ("  static const double inv2[] = {")
   s = "   "
   for r in T:
      s = s + " " + get_hex(r) + ","
      if len(s)+11>=80:
         print (s)
         s = "   "
   if s!="   ":
      print (s)
   print ("  };")
   print ("  static const double log2inv2[][2] = {")
   for h,l in U:
      print ("    {" + get_hex(h) + "," + get_hex(l) + "},")
   print ("  };")
   print (Zmin, Zmax)

# exp2_tables()
# 21 1.0410278456845496999236326786e-16
# 1.0410278456845496999236326786e-16
# thus the maximal absolute error between u[i] and 2^r[i] is
# 1.0410278456845496999236326786e-16 < 2^-53.092
def exp2_tables():
   T = [RR((i-16)/2^5) for i in [0..32]]
   print ("  static const double exp2_T[] = {")
   s = "   "
   for r in T:
      s = s + " " + get_hex(r) + ","
      if len(s)+10>=80:
         print (s)
         s = "   "
   if s!="   ":
      print (s)
   print ("  };")
   U = [RR(n(2^x.exact_rational(),100)) for x in T]
   print ("  static const double exp2_U[] = {")
   s = "   "
   for r in U:
      s = s + " " + get_hex(r) + ","
      if len(s)+20>=80:
         print (s)
         s = "   "
   if s!="   ":
      print (s)
   print ("  };")
   maxerr = 0
   for i in range(len(T)):
      x = T[i]
      X = x.exact_rational()
      err = abs(n(2^X-U[i].exact_rational(),100))
      if err>maxerr:
         print (i,err)
         maxerr = err
   return maxerr

# return true if x-1 is exact
def x_minus_one_exact(x):
   if x != x: # x = NaN
      return false
   if abs(x) == R24("+Inf"):
      return false
   y = x-R24(1)
   return y.exact_rational() == x.exact_rational()-1
   
# return the 'ulp' of the interval x i.e. max(ulp(t)) for t in x
# this internal routine is used below
def RIFulp(x):
   return max(x.lower().ulp(),x.upper().ulp())

def analyze_p1():
   z = RIF(-1/32,17/512)
   z2 = z*z
   err_z2 = RIFulp(z2)
   P = ["0x0p0", "0x1.71547652b82f8p0", "-0x1.71547652b85dep-1",
    "0x1.ec709dc45dfb6p-2", "-0x1.7154764584224p-2", "0x1.2776b6f101774p-2",
    "-0x1.ec718114b01e7p-3", "0x1.a6c377acc0adbp-3", "-0x1.6f7027053a9e5p-3"]
   P = [RIF(RR(x,16)) for x in P]
   c7 = P[8]*z+P[7]
   err_c7 = RIFulp(c7)
   c5 = P[6]*z+P[5]
   err_c5 = RIFulp(c5)
   c3 = P[4]*z+P[3]
   err_c3 = RIFulp(c3)
   c1 = P[2]*z+P[1]
   err_c1 = RIFulp(c1)
   z4 = z2*z2
   err_z4 = RIFulp(z4)+2*z2.abs().upper()*err_z2
   c5 = c7*z2+c5
   err_c5 = err_c7*z2.abs().upper()+c7.abs().upper()*err_z2+err_c5+RIFulp(c5)
   c1 = c3*z2+c1
   err_c1 = err_c3*z2.abs().upper()+c3.abs().upper()*err_z2+err_c1+RIFulp(c1)
   c1 = c5*z4+c1
   err_c1 = err_c5*z4.abs().upper()+c5.abs().upper()*err_z4+err_c1+RIFulp(c1)
   relerr_c1 = err_c1/c1.abs().lower() # bound on relative error on c1
   res = z*c1
   # since c1 has relative error <= relerr_c1, the product z*c1 has
   # relative error < 2^-52, and the Sollya relative error is < 2^-49.809,
   # the total relative error is bounded by
   # (1+relerr_c1)*(1+2^-52)*(1+2^-49.809)-1
   R = RealField(100)
   relerr = (1+R(relerr_c1))*(1+R(2^-52))*(1+R(2^-49.809))-1
   return RR(relerr), res.lower(), res.upper()

# analyze the maximal relative error of _log2p1() for |x| >= 2^-29
# analyze_log2p1()
# e= 0 i= 11 relerr= 4.72359317396951e-15
# thus the maximal relative error is 4.72359317396951e-15 < 2^-47.589
def analyze_log2p1():
   R = RealField(100)
   err0 = R(2^-58.198) # additional relative error for x >= 2^53   
   p = RIF(R("-0x1.773d5d60f4006p-5",16), R("0x1.8eb1333703407p-5",16))
   log2inv = ["-0x1.0c10500d63aa6p-1", "-0x1.d6753e032ea0fp-2", "-0x1.91bba891f1709p-2", "-0x1.49a784bcd1b8bp-2", "-0x1.24407ab0e073ap-2", "-0x1.acf5e2db4ec94p-3", "-0x1.5c01a39fbd688p-3", "-0x1.08c588cda79e4p-3", "-0x1.663f6fac91316p-4", "0x0p+0", "0x0p+0", "0x1.1bb32a600549dp-4", "0x1.e0b1ae8f2fd56p-4", "0x1.22dadc2ab3497p-3", "0x1.8a8980abfbd32p-3", "0x1.dac22d3e441d3p-3", "0x1.169c05363f158p-2", "0x1.32bfee370ee68p-2", "0x1.5dfdcf1eeae0ep-2", "0x1.7b89f02cf2aadp-2", "0x1.a8ff971810a5ep-2", "0x1.c819dc2d45fe4p-2", "0x1.e7df5fe538ab3p-2", "0x1.042bd4b9a7c99p-1"]
   log2inv = [R(x,16) for x in log2inv]
   n = len(log2inv)
   assert n==24, "n==24"
   relerr_p1 = R(2^-49.058)
   abserr_p1 = relerr_p1*p.abs().upper() # bound on absolute error for p1
   maxerr = 0
   for e in [-29..128]:
      for i in range(24):
         t = RIF(log2inv[i]) + p
         u = R(e) + t
         if log2inv[i]==0 and e==0:
            relerr_u = relerr_p1
         else:
            abserr_u = RIFulp(u) + RIFulp(t) + abserr_p1
            relerr_u = abserr_u/u.abs().lower()
         # for x >= 2^53, we have an additional relative error err0
         if e >= 53:
            relerr_u = (1+relerr_u)*(1+err0)-1
         if relerr_u>maxerr:
            print ("e=", e, "i=", i, "relerr=", relerr_u)
            maxerr = relerr_u

# analyze_q1()
# (7.81483493906268e-14, 0.989110713948377, 1.01088928605163)
# thus the absolute error is bounded by 7.81483493906268e-14 < 2^-43.540
def analyze_q1():
   err0 = RR(2^-43.549) # Sollya error on q1
   z = RIF(-2^-6,2^-6)
   Q = ["1.0", "0x1.62e42fef6d01ap-1", "0x1.ebfbdff7feebap-3", "0x1.c6b167e579beep-5", "0x1.3b2b3428d06dep-7"]
   Q = [RIF(RR(x,16)) for x in Q]
   z2 = z*z
   err_z2 = RIFulp(z2)
   c3 = Q[4]*z+Q[3]
   err_c3 = RIFulp(c3)
   c0 = Q[1]*z+Q[0]
   err_c0 = RIFulp(c0)
   c2 = c3*z+Q[2]
   err_c2 = err_c3*z.abs().upper()+RIFulp(c2)
   res = c2*z2+c0
   err_res = err_c2*z2.abs().upper()+c2.abs().upper()*err_z2+err_c0+RIFulp(res)
   err_res += err0 # Sollya error
   return err_res, res.lower(), res.upper()
