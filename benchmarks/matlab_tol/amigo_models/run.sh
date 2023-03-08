# run all *.m file in directory and redirect output to file with same name
fstr="\n\nCPU time: %S+%U sec\tMax. resident set size: %M KB\t Elapsed: %e sec."
for file in *.m; do
    echo "Running $file, output to ${file%.m}}"
    /usr/bin/time -f $fstr -o ${file%.m} matlab -nodisplay -nosplash -nodesktop -r "run $file; exit" >${file%.m}
done
