source ./functions.sh
for d in */ 
do 
    cd "$d"
    create_files
    cd ../
done
