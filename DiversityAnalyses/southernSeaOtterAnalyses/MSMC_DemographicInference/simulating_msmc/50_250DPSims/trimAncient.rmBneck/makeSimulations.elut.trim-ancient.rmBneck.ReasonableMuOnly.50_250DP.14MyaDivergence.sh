
reasonableMu=8.644114e-09

for mu in $reasonableMu
do
# 30Mb chunks
python simulateOutput.TrimAncient-AndBottleneck_20180125.elut.InputMu_50_250DP.py --mu $mu > macs.Simulation.Elut.MSMC.20180125.Trim-Ancient.rmBneck.50_250DP.30MbChunks.Mu.$mu.sh
# 100kb chunks
python simulateOutput.TrimAncient-AndBottleneck_20180125.elut.InputMu.100kbChunks.py --mu $mu > macs.Simulation.Elut.MSMC.20180125.Trim-Ancient.rmBneck.50_250DP.100kbChunks.Mu.$mu.sh

done
