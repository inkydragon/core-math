def table1(out):
   f = open(out,"w")
   f.write ("static const float T1[] = {\n")
   s = ""
   for i1 in range(2^8):
      u1 = 118*2^7+i1*2^3
      x1 = R24(asbfloat16(u1))
      y1 = exp(x1)
      v = get_hex(y1)
      t = " " + v + ","
      if len(s)+len(t)>=79:
         f.write (s + "\n")
         s = t
      else:
         s += t
   if s!="":
      f.write (s + "\n")
   f.write ("};\n")
   f.close()

def table2(out):
   f = open(out,"w")
   f.write ("static const float T2[] = {\n")
   s = ""
   for i2 in range(2^7):
      h,l = divmod(i2,2^3)
      u2 = (118+h)*2^7+l
      e = u2>>7
      x2 = R24(2^(e-118-16)*l)
      y2 = exp(x2)
      v = get_hex(y2)
      t = " " + v + ","
      if len(s)+len(t)>=79:
         f.write (s + "\n")
         s = t
      else:
         s += t
   if s!="":
      f.write (s + "\n")
   f.write ("};\n")
   f.close()
      
