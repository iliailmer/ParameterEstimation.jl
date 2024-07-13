#!/bin/bash
# run all *.m file in directory and redirect output to file with same name
fstr="\n\nCPU time: %S+%U sec\tMax. resident set size: %M KB\t Elapsed: %e sec."
mfiles=$(ls biohydrogenation* crauste* daisy_mamil* fitzhugh_nagumo* harmonic* hiv* lotka_volterra* seir* vanderpol*)
mkdir -p outputs
for file in ${mfiles[@]}; do
    echo "Running $file, output to outputs/${file%.m}.out"
    /usr/bin/time -f "$fstr" matlab -nodisplay -nosplash -nodesktop -r "run $file; exit" &>outputs/${file%.m}.out
    # break
done
