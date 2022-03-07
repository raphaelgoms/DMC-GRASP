#!/bin/bash
#F=(ras20)

# Funções Faceis:
#F=(bea boo bra eas mat shc ros2 sch2 shu zak2 dej har34 \
#	zak5  sch6 tri6  ssq10 tri10 gri10 pow dxp ros20\
#	zak10 ssq20 zak20 gri10 pow dxp ros20 \
#	 perm perm0 gri20 ack30 lev sph)

# Funções Medias: (hs = 0.01)
#F=(zak20) RODAR !!!
#F=(bhc gol col pws she45 she47 she410 har64 \
#	ras20 perm0 ros5 perm ros10 \
#	gri10 ros20 gri20 ack30) 
#F=(ras20 perm0 ros5 perm ros10)
 

# Funções Difíceis: 
#F=(gri10 ros20 gri20 ack30)

# Funções Super Difíceis: (hs = 0.001)
F=(ros50) # zak50 ros50 zak100 ros100 ros500)

#F=(ros10 gri10 pow dxp ros20 perm perm0 \
# gri20 ack30 lev sph)

declare -A funcDim=([eas]=2 [shc]=2 [bea]=2 [tri6]=6 [tri10]=10 [boo]=2 [she45]=5 [she410]=10 [zak5]=5 [bra]=2 [gol]=2 [ros2]=2 [shu]=2 \
		   	[har34]=3 [she47]=7 [ros5]=5 [bhc]=2 \
           	[mat]=2 [ros2]=2 [shu]=2 [dej]=3 [col]=4 [perm]=4 [ras2]=2 [sch2]=2 [zak2]=2 [perm0]=4 [pws]=4 [ras5]=5 [sch6]=6 [har64]=6 \
           	[ros5]=5 [zak5]=5 [ros10]=10 [gri10]=10 [ras10]=10 [ack10]=10 [zak10]=10 [ssq10]=10 [gri20]=20 [ros20]=20 [zak20]=20 [ras20]=20 [ssq20]=20 \
           	[dxp]=25 [lev]=30 [pow]=24 [ack30]=30 [sph]=30 [zak50]=50 [ros50]=50 [zak100]=100 [ros100]=100 [ros500]=500)

declare -A funcFullName=([bea]=BEALE [bhc]=BOHACHEVSKY [boo]=BOOTH [bra]=BRANIN [eas]=EASOM [gol]=GOLDSTEINPRICE [mat]=MATYAS [shc]=HUMP 
						[ros2]=ROSENBROCK [sch2]=SCHWEFEL [shu]=SHUBERT [zak2]=ZAKHAROV [dej]=SPHERE [har34]=HARTMANN \
						[col]=COLVILLE [perm]=PERM [perm0]=PERM0 [pws]=POWERSUM [she45]=SHEKEL [she47]=SHEKEL [she410]=SHEKEL [ros5]=ROSENBROCK \
						[zak5]=ZAKHAROV [har64]=HARTMANN [sch6]=SCHWEFEL [tri6]=TRID \
						[gri10]=GRIEWANK [ras10]=RASTRIGIN [ros10]=ROSENBROCK [ssq10]=SUMSQUARES [tri10]=TRID [zak10]=ZAKHAROV \
						[gri20]=GRIEWANK [ras20]=RASTRIGIN [ros20]=ROSENBROCK [ssq20]=SUMSQUARES [zak20]=ZAKHAROV [pow]=POWELL [dxp]=DIXONPRICE \
						[ack30]=ACKLEY [lev]=LEVY [sph]=SPHERE [zak50]=ZAKHAROV [ros50]=ROSENBROCK [zak100]=ZAKHAROV [ros100]=ROSENBROCK [ros500]=ROSENBROCK)
E=(10)



while read config; do
echo "[$(date)]: Iniciou a exec. para a configuração $config" | tee -a log_v2.txt
	for func in ${F[@]}; do
		echo "[$(date)]: >> Função $func" | tee -a log_v2.txt
		./CGrasp -i ${funcFullName[$func]} --nvar ${funcDim[$func]} $config
	done
echo "[$(date)]: Finalizou para a configuração $config" | tee -a log_v2.txt
done < configurations


