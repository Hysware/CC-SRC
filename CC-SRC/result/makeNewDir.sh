## check whether the file exit,if there is no $dir exit,make a new one
## 2024.1.4
# run in here
targetDir='./fig/'
cd $targetDir
nameArray=("174deg" "168deg" "163deg" "157deg" "152deg")
for dir in ${nameArray[@]}
do
  if [ ! -e $dir ];then  # not exist
    echo "no $dir ..."
    mkdir $dir
    echo "mkdir $dir"
    # exit 1
  else                   # exist
    echo "$dir exist in here"
  fi
done
