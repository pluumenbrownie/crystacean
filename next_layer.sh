#!/bin/bash

for d in $(find -name "final_*" -type d); 
    do echo $d; 
    python /home/wessel/Documents/bachelorproject/cli.py from-dft-folder -c 4.5 -m 2 $d next;
done

python /home/wessel/Documents/bachelorproject/cli.py cp2kify --destructive next next final_0000/testlarge.in final_0000/BASIS final_0000/run_cp2k.sh;

# python /home/wessel/Documents/bachelorproject/cli.py cp2kify test_0000/ out_0000 final_0000/testlarge.in final_0000/BASIS final_0000/run_cp2k.sh

#for d in $(find final_* -type d); do echo $d ;cp testlarge.in $d ; cd $d ;rm *.out;sbatch run_cp2k.sh ;cd ../;done
# python ~/Documents/bachelorproject/cli.py from-dft-folders -c 3.5 -t -m 0 -f -d 0.1 $d ../../exports/second_layer_data/
