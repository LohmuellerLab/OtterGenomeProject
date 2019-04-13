# this uses the correct if statements to determine filter type
# qsub script branch date G/S/N cleandata SwampDate
# otters ancestral branch
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep.sh otters 20170922 GMED 0
# sea otter
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep.sh elut 20170922 GMED 0
# giant otter
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep.sh pbra 20170922 GMED 0


# replicates 2 and 3 : (omega 0.2, omega 2)
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep2.sh otters 20170922 GMED 0
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep2.sh elut 20170922 GMED 0
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep2.sh pbra 20170922 GMED 0

qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep3.sh otters 20170922 GMED 0
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep3.sh elut 20170922 GMED 0
qsub Step6a_PAML_postGuidance.wGSN.ClnData.Opts.rep3.sh pbra 20170922 GMED 0
