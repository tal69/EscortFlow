# Full factorial on 9x5 Grid with output at (4,0).  4 and 6 escorts, arrival rate 0.2, and 0.4 request/sec
# Factor: Surrogate, Hybrid, PLPR, Limited attention

python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 6 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2 -I  --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.2     --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4 -I  --hybrid # Surrogate and Hybrid
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     --full # None
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     # Surrogate only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     --full --hybrid # Hybrid only
python EscortFlowSim_v8.py --seed 0  -o 0.2 -x 9 -y 5 -O 4 0 -e 4 -M 7 -f FullFactorial9x5-o0.2.csv -S 1000 -R 0.4     --hybrid # Surrogate and Hybrid