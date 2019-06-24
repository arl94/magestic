#!/bin/bash

if [[ "$@" == "" ]]; then
  set -- "all"
fi

snakemake --cluster-config cluster.json \
  --cluster "{cluster.sbatch} \
  --cpus-per-task {cluster.n} --mem {cluster.mem} -t {cluster.time} \
  {cluster.moreoptions}" --jobs 999 --jobname "{rulename}.{jobid}" \
  --timestamp --keep-going \
  --latency-wait 600 --restart-times 2 \
  --use-conda "$@"
