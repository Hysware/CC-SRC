#!/bin/bash

data="484180,493278,502924,511940,518232,527591,535624,545566,555582,571089,581156,588649,605817,613154,623293,632508,641867,652070,661195,671461,678628,686628"  # 数据以逗号分隔

# 用 IFS（内部字段分隔符）设置逗号为分隔符，然后循环
IFS=',' read -ra nums <<< "$data"
echo 11
for i in "${nums[@]}"; do
echo 11

    echo "当前 i = $i"
    mkdir $i
    cp ccfull-sc.inp ccfull ./$i/
    cd ./$i/
    e=$(echo "scale=4;$i/10000"|bc)

    sed -e "11a\\$e,$e,2" ./ccfull-sc.inp > cc.inp.tmp    ### change V0
    mv ./cc.inp.tmp ./ccfull-sc.inp
    sed "11d" ./ccfull-sc.inp > cc.inp.tmp
    tr -d '\r' < ./cc.inp.tmp > ./ccfull-sc.inp

    rm ./cc.inp.tmp

    ./ccfull < ccfull-sc.inp | tee -a output.txt

    cd ..

done
