#!/bin/bash
# run all *.m file in directory and redirect output to file with same name
fstr="\n\nCPU time: %S+%U sec\tMax. resident set size: %M KB\t Elapsed: %e sec."
mfiles=(biohydrogenation crauste daisy_mamil3 daisy_mamil4 fitzhugh_nagumo harmonic hiv lotka_volterra seir vanderpol)
mkdir -p outputs
for file in ${mfiles[@]}; do
    echo $file;
    for i in {0..9}; do
        for j in 1 2 3; do 
            filename=${file}_${i}_${j};
            echo "Running ${filename}.m, output to outputs/${filename}.out";
            /usr/bin/time -f "$fstr" matlab -nodisplay -nosplash -nodesktop -r "run ${filename}/${filename}.m; exit" &>outputs/${filename}.out
            # break
        done
    done
done
