#!/bin/bash
#$ -l highp,h_rt=60:00:00,h_data=8G
#$ -pe shared 5
#$ -N pOrtho2_10perc
#$ -cwd
#$ -m bea
#$ -o ./pOrtho_10.out
#$ -e ./pOrtho_10.err
#$ -M ab08028
#$ -t 1-15

proteinortho5=/u/home/a/ab08028/bin/proteinortho_v5.15/proteinortho5.pl
i=$((SGE_TASK_ID - 1))

$proteinortho5 -synteny -cov=25 -e=1e-06 -blastpath=/u/home/a/ab08028/bin/ncbi-blast-2.5.0+/bin/ -cpus=5 -verbose -step=2 -startat=${i} -stopat=${i} -graph amel.faa btau.faa  cfam.faa  ecab.faa  elut.faa  fcat.faa  hsap.faa  mfur.faa pbra.faa pvam.faa mluc.faa ttru.faa oros.faa umar.faa lwed.faa

