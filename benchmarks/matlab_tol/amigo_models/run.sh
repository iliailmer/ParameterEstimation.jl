# run all *.m file in directory and redirect output to file with same name
fstr="\n\nCPU time: %S+%U sec\tMax. resident set size: %M KB\t Elapsed: %e sec."
mfiles=(crauste.m daisy_mamil3.m FHN.m LotkaVolterra.m simple.m daisy_ex3_local.m daisy_mamil4_local.m HIV.m)
mkdir -p outputs
for file in ${mfiles[@]}; do
    echo "Running $file, output to outputs/${file%.m}.out"
    /usr/bin/time -f "$fstr" matlab -nodisplay -nosplash -nodesktop -r "run $file; exit" &>outputs/${file%.m}.out
    # break
done
