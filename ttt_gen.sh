#!/bin/bash
F=(bhc gol col pws she45 she47 she410 har64 ras20)

#Funcoes Faceis:
#F=(bea bhc boo bra eas gol mat shc ros2 sch2 shu zak2 dej har34 \
#  col perm perm0 pws she45 she47 she410 zak5 har64 sch6 tri6 ras10 zak10 ssq10 tri1) 

# Funcoes dificeis:
#F=(zak5 har64 sch6 tri6 ras10 zak10 ssq10 tri10 zak20 ras20 ros20 ssq20 pow\
#F=(ack30 lev sph gri10 gri20 ros5 ros10 dxp)

declare -A funcDim=([eas]=2 [shc]=2 [bea]=2 [tri6]=6 [tri10]=10 [boo]=2 [she45]=5 [she410]=10 [zak5]=5 [bra]=2 [gol]=2 [ros2]=2 [shu]=2 \
		   	[har34]=3 [she47]=7 [ros5]=5 [bhc]=2 \
           	[mat]=2 [ros2]=2 [shu]=2 [dej]=3 [col]=4 [perm]=4 [ras2]=2 [sch2]=2 [zak2]=2 [perm0]=4 [pws]=4 [ras5]=5 [sch6]=6 [har64]=6 \
           	[ros5]=5 [zak5]=5 [ros10]=10 [gri10]=10 [ras10]=10 [ack10]=10 [zak10]=10 [ssq10]=10 [gri20]=20 [ros20]=20 [zak20]=20 [ras20]=20 [ssq20]=20 \
           	[dxp]=25 [lev]=30 [pow]=24 [ack30]=30 [sph]=30) 
 
declare -A funcFullName=([bea]=BEALE [bhc]=BOHACHEVSKY [boo]=BOOTH [bra]=BRANIN [eas]=EASOM [gol]=GOLDSTEINPRICE [mat]=MATYAS [shc]=HUMP 
						[ros2]=ROSENBROCK [sch2]=SCHWEFEL [shu]=SHUBERT [zak2]=ZAKHAROV [dej]=SPHERE [har34]=HARTMANN \
						[col]=COLVILLE [perm]=PERM [perm0]=PERM0 [pws]=POWERSUM [she45]=SHEKEL [she47]=SHEKEL [she410]=SHEKEL [ros5]=ROSENBROCK \
						[zak5]=ZAKHAROV [har64]=HARTMANN [sch6]=SCHWEFEL [tri6]=TRID \
						[gri10]=GRIEWANK [ras10]=RASTRIGIN [ros10]=ROSENBROCK [ssq10]=SUMSQUARES [tri10]=TRID [zak10]=ZAKHAROV \
						[gri20]=GRIEWANK [ras20]=RASTRIGIN [ros20]=ROSENBROCK [ssq20]=SUMSQUARES [zak20]=ZAKHAROV [pow]=POWELL [dxp]=DIXONPRICE \
						[ack30]=ACKLEY [lev]=LEVY [sph]=SPHERE)




for func in ${F[@]}; do
	echo "Instance ${funcFullName[$func]}"

	# while read config; do
	# 	echo "[$(date)]: Iniciou a exec. para a configuração $config" | tee -a log_v2.txt
	# 	./CGrasp -i ${funcFullName[$func]} --nvar ${funcDim[$func]} $config
	# done < configurations


	instance=""
	if [ "${funcFullName[$func]}" == "SHEKEL" ]; then
		instance="${funcFullName[$func]}4${funcDim[$func]}"
	else
		instance=${funcFullName[$func]}${funcDim[$func]}
	fi 

	echo $instance
	perl tttplots.pl -f ./tttplots_data/c_"$instance"_cfo_file

	# perl tttplots.pl -f ./tttplots_data/dmc_"$instance"_cfo_file
	# gnuplot -e "algorithm='dmc'; function='$instance'" ttt-plot-compare2.gpl

	# perl tttplots.pl -f ./tttplots_data/xdmc_"$instance"_cfo_file
	# gnuplot -e "algorithm='xdmc'; function='$instance'" ttt-plot-compare2.gpl

	# perl tttplots.pl -f ./tttplots_data/mxdmc_"$instance"_cfo_file
	# gnuplot -e "algorithm='mxdmc'; function='$instance'" ttt-plot-compare2.gpl

	# perl tttplots.pl -f ./tttplots_data/rmxdmc_"$instance"_cfo_file
	# gnuplot -e "algorithm='rmxdmc'; function='$instance'" ttt-plot-compare2.gpl

	perl tttplots.pl -f ./tttplots_data/amxdmc_"$instance"_cfo_file
	gnuplot -e "algorithm='amxdmc'; function='$instance'" ttt-plot-compare2.gpl
done


#./script_exp_gap.sh
#perl tttplots.pl -f c_sol_file
#perl tttplots.pl -f rmxdmc_sol_file
#gnuplot ttt-plot-compare.gpl    