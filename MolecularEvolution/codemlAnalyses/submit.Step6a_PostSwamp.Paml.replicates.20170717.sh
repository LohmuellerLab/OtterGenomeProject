# this uses the correct if statements to determine filter type
# qsub script branch date G/S/N cleandata SwampDate
# otters ancestral branch
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.sh otters 20170717 S 1 20170713

# elut:
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.sh elut 20170717 S 1 20170713

# pbra:
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.sh pbra 20170717 S 1 20170713

# replicates 2 and 3 : (omega 0.2, omega 2)
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep2.sh otters 20170717 S 1 20170713
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep3.sh otters 20170717 S 1 20170713
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep2.sh elut 20170717 S 1 20170713
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep2.sh pbra 20170717 S 1 20170713
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep3.sh elut 20170717 S 1 20170713
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep3.sh pbra 20170717 S 1 20170713
