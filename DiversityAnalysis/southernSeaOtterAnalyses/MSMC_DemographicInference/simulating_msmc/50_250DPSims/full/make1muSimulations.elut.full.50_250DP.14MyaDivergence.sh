reasonableMu=8.644114e-09


for mu in $reasonableMu
do
python simulateOutput.20171220.elut.InputMu.50_250DP.py --mu $mu > macs.Simulation.Elut.MSMC.20171220.Full.50_250DP.Mu.$mu.sh
done
