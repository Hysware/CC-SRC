# just change potenteil parameters
timeStart=$(date +%s)

gfortran -o ./CC-SRC/cc/ccfull ./CC-SRC/cc/ccfull-sc.f  #compile the sorce file

for((Vw=$1;Vw<$1+5;Vw+=10))   # vw * 10
do
  for((Rw=$2;Rw<$2+5;Rw+=10))   # Rw *100
  do
    for((Aw=10;Aw<200;Aw+=20))   # Aw *1000
    do
      # build files
      num1=$(printf "%03d" $Vw)
      num2=$(printf "%03d" $Rw)
      num3=$(printf "%03d" $Aw)
      dirname=$num1$num2$num3
      mkdir $dirname
      # copy files
      cp ./CC-SRC/cc/ccfull-sc.inp ./$dirname/
      cp ./CC-SRC/cc/ccfull ./$dirname/
      cd $dirname
      # changes 
      ###### here may be need changes
      vw=$(echo "scale=3;$Vw/10  "|bc)
      rw=$(echo "scale=3;$Rw/100 "|bc)
      aw=$(echo "scale=3;$Aw/1000"|bc)



      ## change potenteil parameters, row 8 in input file.
      sed -e 8a\ $vw,$rw,$aw ./ccfull-sc.inp > cc.inp.tmp
      mv ./cc.inp.tmp ./ccfull-sc.inp
      sed "8d" ./ccfull-sc.inp > cc.inp.tmp
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
