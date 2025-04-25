############### chi^2
## 采用相对路径表示
# dirname0='/media/DATA4T2/ZhangHY/ccfull-7/'
dirname0='../'
dirname1='/qel.dat'
dirname2='origin_Exp.xlsx'
dir_Exp=$dirname0$dirname2
# dir_fig='/media/DATA4T2/ZhangHY/ccfull-7/result/cc7fig/'
dir_fig='./fig/'
fig_type='.png'
fig_barri_type='_barrier.png'
#python3 ./get_data0.py $dir_Exp 1  # 实验的能量值

for((beta2=31;beta2<32;beta2+=2))   #beta2*1000
do
       for((beta3=-50;beta3<60;beta3+=10))   #beta3*1000
        do
                for((beta4=0;beta4<350;bet4+=50))   # beta4 *1000
                do
			# build files
                        num1=$(printf "%03d" $beta2)
                        num2=$(printf "%03d" $beta3)
                        num3=$(printf "%03d" $beta4)
                        dirname=$num1$num2$num3
			
			dir=$dirname0$dirname0$dirname$dirname1
			echo $dirname			
                        cross_fig=$dir_fig$dirname$fig_type
                        barri_fig=$dir_fig$dirname$fig_barri_type
			python3 ./chi2xlsx.py $dir $dir_Exp $cross_fig $barri_fig
                        
                        # 下面将输出的chi2存入一个数组以及对应的文件夹进行存储，以对其进行排序			
			
done
done
done

