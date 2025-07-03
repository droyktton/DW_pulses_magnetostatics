a=1;b=1;
set print 'tauyampvsk0.dat';
set term png
f(x)=a*(1-exp(-x/b))
do for [k0 in '1.0 0.1 0.01 0.001 0.0001']{
  set out sprintf("fit%s.png",k0)
  fit [0:] f(x) 'monitor_k0_'.k0.'.dat' u ($1-2950):2 via a, b;
  plot [0:] 'monitor_k0_'.k0.'.dat' u ($1-2950):2, f(x);
  print sprintf("%f %f %f", k0+0.0, a, b); pause 1
}
