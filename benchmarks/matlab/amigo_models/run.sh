/usr/bin/time -f "\n\nCPU time: %S+%U sec\tMax. resident set size: %M KB\t Elapsed: %e sec." matlab -nodisplay -nosplash -nodesktop -r "run(\"./HIV.m\"); exit;"
