# shell script to run all julia scripts in this directory

for each in *.jl; do
    echo "Running $each"
    # if system is macos use gtime
    os_type=$(uname)
    if [ "$os_type" = "Darwin" ]; then
        gtime -f "\n\nCPU time: %S+%U sec\tMax. resident set size: %M KB\t Elapsed: %e sec." julia --project=../../ $each
    else
        /usr/bin/time -f "\n\nCPU time: %S+%U sec\tMax. resident set size: %M KB\t Elapsed: %e sec." julia --project=../../ $each
    fi
done
