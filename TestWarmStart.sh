
python EscortFlowSim_v7.py -L -x 13 -y 7 -O 6 0 -R 0.4 -S 500 --hybrid --hybrid_ratio 1 -f Warmstart.csv  --warmstart greedy
python EscortFlowSim_v7.py -L -x 13 -y 7 -O 6 0 -R 0.4 -S 500 --hybrid --hybrid_ratio 1 -f Warmstart.csv  --warmstart ilp
python EscortFlowSim_v7.py -L -x 13 -y 7 -O 6 0 -R 0.4 -S 500 --hybrid --hybrid_ratio 1 -f Warmstart.csv  --warmstart ilp greedy
python EscortFlowSim_v7.py -L -x 13 -y 7 -O 6 0 -R 0.4 -S 500 --hybrid --hybrid_ratio 1 -f Warmstart.csv