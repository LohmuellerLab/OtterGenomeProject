reasonableMu=8.644114e-09

for mu in $reasonableMu
do
python simulateOutput.20190216.elut2.InputMu_50_250DPFilter.py --mu $mu > macs.Simulation.elut2.MSMC.Full.50_250DP.Mu.$mu.sh
done
