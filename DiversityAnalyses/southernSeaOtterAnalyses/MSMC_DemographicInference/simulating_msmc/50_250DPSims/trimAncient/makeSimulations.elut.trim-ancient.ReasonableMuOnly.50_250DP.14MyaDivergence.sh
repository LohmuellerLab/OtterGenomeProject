reasonableMu=8.644114e-09


for mu in $reasonableMu
do
python simulateOutput.TrimAncient-20171220.elut.InputMu.50_250DP.py --mu $mu > macs.Simulation.Elut.MSMC.20180125.Trim-Ancient.50_250DP.30MbChunks.Mu.$mu.sh
python simulateOutput.TrimAncient-20171220.elut.InputMu.100kbChunks.50_250DP.py --mu $mu > macs.Simulation.Elut.MSMC.20180125.Trim-Ancient.50_250DP.100kbChunks.Mu.$mu.sh

done
