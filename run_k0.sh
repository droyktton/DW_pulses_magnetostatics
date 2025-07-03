
for h in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
	for k0 in 0.000001 0.00001 0.0001 0.001 0.01
	do
		echo $h $k0
		
  		./cuerda 256 $h $k0 1.0 100.0 0.1

  		mv monitor2.dat "monitor_k0_"$k0"_h_"$h".dat"
  		#mv monitor.dat monitor_k0_$k0.dat
  		#mv monitorstrob.dat monitorstrob_k0_$k0.dat
	done
done


#gnuplot -p -e "\
#set xla 't k_0^{0.9}'; set yla 'x k_0^{0.9}'; tmax=3000; \
#p [:] for[k0 in '0.000001 0.00001 0.0001 0.001 0.01']  'monitor_k0_'.k0.'.dat' u ((\$1-tmax)*k0**0.9):(\$2*k0**0.9) w lp t ''.k0, x*0.1"
