# Run the simulation that are need for Figures 8 and 9 in the paper
# The effect of the heuristic components on the mean lead time under the offline setting in a
# 9x5 grid with 95% confidence intervals (figure 8) and the run times (figure 9)
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -R 0.1 -L --greedy -H
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 5 -E 1 -L -R 0.1 -M 5 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 5 -I 5 -E 1 -L -R 0.1 -M 5 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 5 -E 1 -L -R 0.1 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 5 -I 5 -E 1 -L -R 0.1 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 20 -E 1 -L -R 0.1 -M 5 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 20 -I 20 -E 1 -L -R 0.1 -M 5 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 20 -E 1 -L -R 0.1 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 20 -I 20 -E 1 -L -R 0.1 --full

python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -R 0.2 -L --greedy
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 5 -E 1 -L -R 0.2 -M 5 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 5 -I 5 -E 1 -L -R 0.2 -M 5 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 5 -E 1 -L -R 0.2 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 5 -I 5 -E 1 -L -R 0.2 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 20 -E 1 -L -R 0.2 -M 5 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 20 -I 20 -E 1 -L -R 0.2 -M 5 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 20 -E 1 -L -R 0.2 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 20 -I 20 -E 1 -L -R 0.2 --full

python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -R 0.4 -L --greedy
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 5 -E 1 -L -R 0.4 -M 5 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 5 -I 5 -E 1 -L -R 0.4 -M 5 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 5 -E 1 -L -R 0.4 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 5 -I 5 -E 1 -L -R 0.4 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 20 -E 1 -L -R 0.4 -M 5 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 20 -I 20 -E 1 -L -R 0.4 -M 5 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 20 -E 1 -L -R 0.4 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 6 -f figure_8_and_9.csv -S 1600 -T 20 -I 20 -E 1 -L -R 0.4 --full

python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -R 0.1 -L --greedy
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 5 -E 1 -L -R 0.1 -M 5 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 5 -I 5 -E 1 -L -R 0.1 -M 5 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 5 -E 1 -L -R 0.1 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 5 -I 5 -E 1 -L -R 0.1 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 20 -E 1 -L -R 0.1 -M 5 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 20 -I 20 -E 1 -L -R 0.1 -M 5 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 20 -E 1 -L -R 0.1 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 20 -I 20 -E 1 -L -R 0.1 --full

python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -R 0.2 -L --greedy
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 5 -E 1 -L -R 0.2 -M 5 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 5 -I 5 -E 1 -L -R 0.2 -M 5 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 5 -E 1 -L -R 0.2 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 5 -I 5 -E 1 -L -R 0.2 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 20 -E 1 -L -R 0.2 -M 5 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 20 -I 20 -E 1 -L -R 0.2 -M 5 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 20 -E 1 -L -R 0.2 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 20 -I 20 -E 1 -L -R 0.2 --full

python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -R 0.4 -L --greedy
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 5 -E 1 -L -R 0.4 -M 5 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 5 -I 5 -E 1 -L -R 0.4 -M 5 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 5 -E 1 -L -R 0.4 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 5 -I 5 -E 1 -L -R 0.4 
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 20 -E 1 -L -R 0.4 -M 5 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 20 -I 20 -E 1 -L -R 0.4 -M 5 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 20 -E 1 -L -R 0.4 --full
python ./EscortFlowSim_v5.py --seed 1 -x 9 -y 5 -O 4 0 -e 4 -f figure_8_and_9.csv -S 1600 -T 20 -I 20 -E 1 -L -R 0.4 --full