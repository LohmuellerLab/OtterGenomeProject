reasonableMu=8.644114e-09

for mu in $reasonableMu
do
python simulateOutput.TrimAncient-20171220.pbra.InputMu.50_250DP.py --mu $mu > macs.Simulation.Pbra.MSMC.TrimAncient.50_250DP.Mu.$mu.sh
done
