reasonableMu=8.644114e-09

for mu in $reasonableMu
do
python simulateOutput.20171220.pbra.InputMu_50_250DPFilter.py --mu $mu > macs.Simulation.Pbra.MSMC.Full.50_250DP.Mu.$mu.sh
done
