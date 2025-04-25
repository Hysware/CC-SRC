# 该程序在运行时，需添加一个参数，代表beta2,目的：分为多个后台运行，减小运行时间
timeStart=$(date +%s)

gfortran -o ./CC-SRC/cc/ccfull ./CC-SRC/cc/ccfull-sc.f  #compile the sorce file



for((beta2=289;beta2<389;beta2+=100))   #beta2*1000,可以将beta2改为一个可以从其他地方读取的参数
do
	for((beta3=$1;beta3<$1+5;beta3+=5))   #beta3*1000
	do
		for((beta4=100;beta4<220;beta4+=10))   # beta4 *1000
		do
			# build files
			num1=$(printf "%03d" $beta2)
			num2=$(printf "%03d" $beta3)
			num3=$(printf "%03d" $beta4)
			dirname=$num1$num2$num3
			mkdir $dirname
			# copy files
			cp ./CC-SRC/cc/ccfull-sc.inp ./$dirname/
			cp ./CC-SRC/cc/ccfull ./$dirname/
			cd $dirname
			# changes 
                        ###### here may be need changes
			b2=$(echo "scale=3;$beta2/1000"|bc)
			b3=$(echo "scale=3;$beta3/1000"|bc)
			b4=$(echo "scale=3;$beta4/1000"|bc)


                        ## also need to be changed with diff reactions.
		#	sed -e 2a\3,$b2,$b3,$b4 ./ccfull-sc.inp > cc.inp.tmp  
			sed -e 3a\0.082,$b2,$b4,5 ./ccfull-sc.inp > cc.inp.tmp
			mv ./cc.inp.tmp ./ccfull-sc.inp
			sed "3d" ./ccfull-sc.inp > cc.inp.tmp
			tr -d '\r' < ./cc.inp.tmp > ./ccfull-sc.inp

			sed -e 4a\1.012,$b3,3,1 ./ccfull-sc.inp > cc.inp.tmp
			mv ./cc.inp.tmp ./ccfull-sc.inp
			sed "4d" ./ccfull-sc.inp > cc.inp.tmp
			tr -d '\r' < ./cc.inp.tmp > ./ccfull-sc.inp
			rm ./cc.inp.tmp
			
			./ccfull < ccfull-sc.inp | tee -a output.txt
			
			# result collect

			cd ..
# rm -rf $dirname
done
done
done

timeEnd=$(date +%s)
time=$((timeEnd-timeStart))
echo "程序运行了 $time s"
