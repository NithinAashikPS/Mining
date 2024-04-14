sudo apt update

sudo apt upgrade 

sudo apt install libgmp-dev

nvcc main.cu -o out -lgmp
./out
