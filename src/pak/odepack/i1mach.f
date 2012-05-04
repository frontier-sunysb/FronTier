      integer function i1mach(i)
      implicit none
c
c  i/o unit numbers.
c
c    i1mach( 1) = the standard input unit.
c
c    i1mach( 2) = the standard output unit.
c
c    i1mach( 3) = the standard punch unit.
c
c    i1mach( 4) = the standard error message unit.
c
c  words.
c
c    i1mach( 5) = the number of bits per integer storage unit.
c
c    i1mach( 6) = the number of characters per character storage unit.
c                 for fortran 77, this is always 1.  for fortran 66,
c                 character storage unit = integer storage unit.
c
c  integers.
c
c    assume integers are represented in the s-digit, base-a form
c
c               sign ( x(s-1)*a**(s-1) + ... + x(1)*a + x(0) )
c
c               where 0 .le. x(i) .lt. a for i=0,...,s-1.
c
c    i1mach( 7) = a, the base.
c
c    i1mach( 8) = s, the number of base-a digits.
c
c    i1mach( 9) = a**s - 1, the largest magnitude.
c
c  floating-point numbers.
c
c    assume floating-point numbers are represented in the t-digit,
c    base-b form
c
c               sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
c
c               where 0 .le. x(i) .lt. b for i=1,...,t,
c               0 .lt. x(1), and emin .le. e .le. emax.
c
c    i1mach(10) = b, the base.
c
c  single-precision
c
c    i1mach(11) = t, the number of base-b digits.
c
c    i1mach(12) = emin, the smallest exponent e.
c
c    i1mach(13) = emax, the largest exponent e.
c
c  double-precision
c
c    i1mach(14) = t, the number of base-b digits.
c
c    i1mach(15) = emin, the smallest exponent e.
c
c    i1mach(16) = emax, the largest exponent e.
c
c  to alter this function for a particular environment,
c  the desired set of data statements should be activated by
c  removing the c from column 1.  also, the values of
c  i1mach(1) - i1mach(4) should be checked for consistency
c  with the local operating system.  for fortran 77, you may wish
c  to adjust the data statement so imach(6) is set to 1, and
c  then to comment out the executable test on i .eq. 6 below.
c  on rare machines a static statement may need to be added.
c  (but probably more systems prohibit it than require it.)
c
c  for ieee-arithmetic machines (binary standard), the first
c  set of constants below should be appropriate, except perhaps
c  for imach(1) - imach(4).
c
      integer i
      integer imach(16),output,sanity
c
      equivalence (imach(4),output)
c
c     machine constants for ieee arithmetic machines, such as the at&t
c     3b series, motorola 68000 based machines (e.g. sun 3 and at&t
c     pc 7300), and 8087 based micros (e.g. ibm pc and at&t 6300).
c
       data imach( 1) /    5 /
       data imach( 2) /    6 /
       data imach( 3) /    7 /
       data imach( 4) /    6 /
       data imach( 5) /   32 /
       data imach( 6) /    4 /
       data imach( 7) /    2 /
       data imach( 8) /   31 /
       data imach( 9) / 2147483647 /
       data imach(10) /    2 /
       data imach(11) /   24 /
       data imach(12) / -125 /
       data imach(13) /  128 /
       data imach(14) /   53 /
       data imach(15) / -1021 /
       data imach(16) /  1024 /, sanity/987/
c
c     machine constants for amdahl machines.
c
c      data imach( 1) /   5 /
c      data imach( 2) /   6 /
c      data imach( 3) /   7 /
c      data imach( 4) /   6 /
c      data imach( 5) /  32 /
c      data imach( 6) /   4 /
c      data imach( 7) /   2 /
c      data imach( 8) /  31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /  16 /
c      data imach(11) /   6 /
c      data imach(12) / -64 /
c      data imach(13) /  63 /
c      data imach(14) /  14 /
c      data imach(15) / -64 /
c      data imach(16) /  63 /, sanity/987/
c
c     machine constants for the burroughs 1700 system.
c
c      data imach( 1) /    7 /
c      data imach( 2) /    2 /
c      data imach( 3) /    2 /
c      data imach( 4) /    2 /
c      data imach( 5) /   36 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   33 /
c      data imach( 9) / z1ffffffff /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -256 /
c      data imach(13) /  255 /
c      data imach(14) /   60 /
c      data imach(15) / -256 /
c      data imach(16) /  255 /, sanity/987/
c
c     machine constants for the burroughs 5700 system.
c
c      data imach( 1) /   5 /
c      data imach( 2) /   6 /
c      data imach( 3) /   7 /
c      data imach( 4) /   6 /
c      data imach( 5) /  48 /
c      data imach( 6) /   6 /
c      data imach( 7) /   2 /
c      data imach( 8) /  39 /
c      data imach( 9) / o0007777777777777 /
c      data imach(10) /   8 /
c      data imach(11) /  13 /
c      data imach(12) / -50 /
c      data imach(13) /  76 /
c      data imach(14) /  26 /
c      data imach(15) / -50 /
c      data imach(16) /  76 /, sanity/987/
c
c     machine constants for the burroughs 6700/7700 systems.
c
c      data imach( 1) /   5 /
c      data imach( 2) /   6 /
c      data imach( 3) /   7 /
c      data imach( 4) /   6 /
c      data imach( 5) /  48 /
c      data imach( 6) /   6 /
c      data imach( 7) /   2 /
c      data imach( 8) /  39 /
c      data imach( 9) / o0007777777777777 /
c      data imach(10) /   8 /
c      data imach(11) /  13 /
c      data imach(12) / -50 /
c      data imach(13) /  76 /
c      data imach(14) /  26 /
c      data imach(15) / -32754 /
c      data imach(16) /  32780 /, sanity/987/
c
c     machine constants for ftn4 on the cdc 6000/7000 series.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   60 /
c      data imach( 6) /   10 /
c      data imach( 7) /    2 /
c      data imach( 8) /   48 /
c      data imach( 9) / 00007777777777777777b /
c      data imach(10) /    2 /
c      data imach(11) /   47 /
c      data imach(12) / -929 /
c      data imach(13) / 1070 /
c      data imach(14) /   94 /
c      data imach(15) / -929 /
c      data imach(16) / 1069 /, sanity/987/
c
c     machine constants for ftn5 on the cdc 6000/7000 series.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   60 /
c      data imach( 6) /   10 /
c      data imach( 7) /    2 /
c      data imach( 8) /   48 /
c      data imach( 9) / o"00007777777777777777" /
c      data imach(10) /    2 /
c      data imach(11) /   47 /
c      data imach(12) / -929 /
c      data imach(13) / 1070 /
c      data imach(14) /   94 /
c      data imach(15) / -929 /
c      data imach(16) / 1069 /, sanity/987/
c
c     machine constants for convex c-1.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   32 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -128 /
c      data imach(13) /  127 /
c      data imach(14) /   53 /
c      data imach(15) /-1024 /
c      data imach(16) / 1023 /, sanity/987/
c
c     machine constants for the cray 1, xmp, 2, and 3.
c
c      data imach( 1) /     5 /
c      data imach( 2) /     6 /
c      data imach( 3) /   102 /
c      data imach( 4) /     6 /
c      data imach( 5) /    64 /
c      data imach( 6) /     8 /
c      data imach( 7) /     2 /
c      data imach( 8) /    63 /
c      data imach( 9) /  777777777777777777777b /
c      data imach(10) /     2 /
c      data imach(11) /    47 /
c      data imach(12) / -8189 /
c      data imach(13) /  8190 /
c      data imach(14) /    94 /
c      data imach(15) / -8099 /
c      data imach(16) /  8190 /, sanity/987/
c
c     machine constants for the data general eclipse s/200.
c
c      data imach( 1) /   11 /
c      data imach( 2) /   12 /
c      data imach( 3) /    8 /
c      data imach( 4) /   10 /
c      data imach( 5) /   16 /
c      data imach( 6) /    2 /
c      data imach( 7) /    2 /
c      data imach( 8) /   15 /
c      data imach( 9) /32767 /
c      data imach(10) /   16 /
c      data imach(11) /    6 /
c      data imach(12) /  -64 /
c      data imach(13) /   63 /
c      data imach(14) /   14 /
c      data imach(15) /  -64 /
c      data imach(16) /   63 /, sanity/987/
c
c     machine constants for the harris slash 6 and slash 7.
c
c      data imach( 1) /       5 /
c      data imach( 2) /       6 /
c      data imach( 3) /       0 /
c      data imach( 4) /       6 /
c      data imach( 5) /      24 /
c      data imach( 6) /       3 /
c      data imach( 7) /       2 /
c      data imach( 8) /      23 /
c      data imach( 9) / 8388607 /
c      data imach(10) /       2 /
c      data imach(11) /      23 /
c      data imach(12) /    -127 /
c      data imach(13) /     127 /
c      data imach(14) /      38 /
c      data imach(15) /    -127 /
c      data imach(16) /     127 /, sanity/987/
c
c     machine constants for the honeywell dps 8/70 series.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /   43 /
c      data imach( 4) /    6 /
c      data imach( 5) /   36 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   35 /
c      data imach( 9) / o377777777777 /
c      data imach(10) /    2 /
c      data imach(11) /   27 /
c      data imach(12) / -127 /
c      data imach(13) /  127 /
c      data imach(14) /   63 /
c      data imach(15) / -127 /
c      data imach(16) /  127 /, sanity/987/
c
c     machine constants for the ibm 360/370 series,
c     the xerox sigma 5/7/9 and the sel systems 85/86.
c
c      data imach( 1) /   5 /
c      data imach( 2) /   6 /
c      data imach( 3) /   7 /
c      data imach( 4) /   6 /
c      data imach( 5) /  32 /
c      data imach( 6) /   4 /
c      data imach( 7) /   2 /
c      data imach( 8) /  31 /
c      data imach( 9) / z7fffffff /
c      data imach(10) /  16 /
c      data imach(11) /   6 /
c      data imach(12) / -64 /
c      data imach(13) /  63 /
c      data imach(14) /  14 /
c      data imach(15) / -64 /
c      data imach(16) /  63 /, sanity/987/
c
c     machine constants for the interdata 8/32
c     with the unix system fortran 77 compiler.
c
c     for the interdata fortran vii compiler replace
c     the z's specifying hex constants with y's.
c
c      data imach( 1) /   5 /
c      data imach( 2) /   6 /
c      data imach( 3) /   6 /
c      data imach( 4) /   6 /
c      data imach( 5) /  32 /
c      data imach( 6) /   4 /
c      data imach( 7) /   2 /
c      data imach( 8) /  31 /
c      data imach( 9) / z'7fffffff' /
c      data imach(10) /  16 /
c      data imach(11) /   6 /
c      data imach(12) / -64 /
c      data imach(13) /  62 /
c      data imach(14) /  14 /
c      data imach(15) / -64 /
c      data imach(16) /  62 /, sanity/987/
c
c     machine constants for the pdp-10 (ka processor).
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   36 /
c      data imach( 6) /    5 /
c      data imach( 7) /    2 /
c      data imach( 8) /   35 /
c      data imach( 9) / "377777777777 /
c      data imach(10) /    2 /
c      data imach(11) /   27 /
c      data imach(12) / -128 /
c      data imach(13) /  127 /
c      data imach(14) /   54 /
c      data imach(15) / -101 /
c      data imach(16) /  127 /, sanity/987/
c
c     machine constants for the pdp-10 (ki processor).
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   36 /
c      data imach( 6) /    5 /
c      data imach( 7) /    2 /
c      data imach( 8) /   35 /
c      data imach( 9) / "377777777777 /
c      data imach(10) /    2 /
c      data imach(11) /   27 /
c      data imach(12) / -128 /
c      data imach(13) /  127 /
c      data imach(14) /   62 /
c      data imach(15) / -128 /
c      data imach(16) /  127 /, sanity/987/
c
c     machine constants for pdp-11 fortrans supporting
c     32-bit integer arithmetic.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   32 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -127 /
c      data imach(13) /  127 /
c      data imach(14) /   56 /
c      data imach(15) / -127 /
c      data imach(16) /  127 /, sanity/987/
c
c     machine constants for pdp-11 fortrans supporting
c     16-bit integer arithmetic.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   16 /
c      data imach( 6) /    2 /
c      data imach( 7) /    2 /
c      data imach( 8) /   15 /
c      data imach( 9) / 32767 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -127 /
c      data imach(13) /  127 /
c      data imach(14) /   56 /
c      data imach(15) / -127 /
c      data imach(16) /  127 /, sanity/987/
c
c     machine constants for the prime 50 series systems
c     wtih 32-bit integers and 64v mode instructions,
c     supplied by igor bray.
c
c      data imach( 1) /            1 /
c      data imach( 2) /            1 /
c      data imach( 3) /            2 /
c      data imach( 4) /            1 /
c      data imach( 5) /           32 /
c      data imach( 6) /            4 /
c      data imach( 7) /            2 /
c      data imach( 8) /           31 /
c      data imach( 9) / :17777777777 /
c      data imach(10) /            2 /
c      data imach(11) /           23 /
c      data imach(12) /         -127 /
c      data imach(13) /         +127 /
c      data imach(14) /           47 /
c      data imach(15) /       -32895 /
c      data imach(16) /       +32637 /, sanity/987/
c
c     machine constants for the sequent balance 8000.
c
c      data imach( 1) /     0 /
c      data imach( 2) /     0 /
c      data imach( 3) /     7 /
c      data imach( 4) /     0 /
c      data imach( 5) /    32 /
c      data imach( 6) /     1 /
c      data imach( 7) /     2 /
c      data imach( 8) /    31 /
c      data imach( 9) /  2147483647 /
c      data imach(10) /     2 /
c      data imach(11) /    24 /
c      data imach(12) /  -125 /
c      data imach(13) /   128 /
c      data imach(14) /    53 /
c      data imach(15) / -1021 /
c      data imach(16) /  1024 /, sanity/987/
c
c     machine constants for the univac 1100 series.
c
c     note that the punch unit, i1mach(3), has been set to 7
c     which is appropriate for the univac-for system.
c     if you have the univac-ftn system, set it to 1.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   36 /
c      data imach( 6) /    6 /
c      data imach( 7) /    2 /
c      data imach( 8) /   35 /
c      data imach( 9) / o377777777777 /
c      data imach(10) /    2 /
c      data imach(11) /   27 /
c      data imach(12) / -128 /
c      data imach(13) /  127 /
c      data imach(14) /   60 /
c      data imach(15) /-1024 /
c      data imach(16) / 1023 /, sanity/987/
c
c     machine constants for vax.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   32 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -127 /
c      data imach(13) /  127 /
c      data imach(14) /   56 /
c      data imach(15) / -127 /
c      data imach(16) /  127 /, sanity/987/
c
c  ***  issue stop 777 if all data statements are commented...
      if (sanity .ne. 987) stop 777
      if (i .lt. 1  .or.  i .gt. 16) go to 10
c
      i1mach = imach(i)
c/6s
c/7s
      if(i.eq.6) i1mach=1
c/
      return
   10 write(output,1999) i
 1999 format(' i1mach - i out of bounds',i10)
      stop
      end
