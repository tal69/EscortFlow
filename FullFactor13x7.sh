# Full factorial on 13x7 Grid with output at (4,0).  4 and 8 escorts, arrival rate 0.1, 0.2, and 0.4 request/sec
# Factor: Surrogate, Hybrid, Greedy Warmstart

python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 20 -E 1 -T 20   --full # None
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 5  -E 1 -T 5 # Surrogate only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 20 -E 1 -T 20 --full --hybrid # Hybrid only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 20  -E 1 -T 20 --full --warmstart greedy zero # Warmstart only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 5  -E 1 -T 5 --hybrid # Surrogate and Hybrid
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 5  -E 1 -T 5 --warmstart greedy zero # Surrogate and Greedy
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 20  -E 1 -T 20 --full --hybrid --warmstart greedy zero # Hybrid and warmstart
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 5  -E 1 -T 5 --warmstart greedy zero --hybrid # Surrogate and Greedy

python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 20 -E 1 -T 20   --full # None
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 5  -E 1 -T 5 # Surrogate only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 20 -E 1 -T 20 --full --hybrid # Hybrid only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 20  -E 1 -T 20 --full --warmstart greedy zero # Warmstart only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 5  -E 1 -T 5 --hybrid # Surrogate and Hybrid
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 5  -E 1 -T 5 --warmstart greedy zero # Surrogate and Greedy
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 20  -E 1 -T 20 --full --hybrid --warmstart greedy zero # Hybrid and warmstart
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.1 -I 5  -E 1 -T 5 --warmstart greedy zero --hybrid # Surrogate and Greedy

python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 20 -E 1 -T 20   --full # None
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 5  -E 1 -T 5 # Surrogate only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 20 -E 1 -T 20 --full --hybrid # Hybrid only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 20  -E 1 -T 20 --full --warmstart greedy zero # Warmstart only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 5  -E 1 -T 5 --hybrid # Surrogate and Hybrid
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 5  -E 1 -T 5 --warmstart greedy zero # Surrogate and Greedy
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 20  -E 1 -T 20 --full --hybrid --warmstart greedy zero # Hybrid and warmstart
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 5  -E 1 -T 5 --warmstart greedy zero --hybrid # Surrogate and Greedy

python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 20 -E 1 -T 20   --full # None
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 5  -E 1 -T 5 # Surrogate only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 20 -E 1 -T 20 --full --hybrid # Hybrid only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 20  -E 1 -T 20 --full --warmstart greedy zero # Warmstart only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 5  -E 1 -T 5 --hybrid # Surrogate and Hybrid
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 5  -E 1 -T 5 --warmstart greedy zero # Surrogate and Greedy
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 20  -E 1 -T 20 --full --hybrid --warmstart greedy zero # Hybrid and warmstart
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.2 -I 5  -E 1 -T 5 --warmstart greedy zero --hybrid # Surrogate and Greedy

python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 20 -E 1 -T 20   --full # None
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 5  -E 1 -T 5 # Surrogate only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 20 -E 1 -T 20 --full --hybrid # Hybrid only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 20  -E 1 -T 20 --full --warmstart greedy zero # Warmstart only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 5  -E 1 -T 5 --hybrid # Surrogate and Hybrid
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 5  -E 1 -T 5 --warmstart greedy zero # Surrogate and Greedy
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 20  -E 1 -T 20 --full --hybrid --warmstart greedy zero # Hybrid and warmstart
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 8 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 5  -E 1 -T 5 --warmstart greedy zero --hybrid # Surrogate and Greedy

python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 20 -E 1 -T 20   --full # None
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 5  -E 1 -T 5 # Surrogate only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 20 -E 1 -T 20 --full --hybrid # Hybrid only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 20  -E 1 -T 20 --full --warmstart greedy zero # Warmstart only
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 5  -E 1 -T 5 --hybrid # Surrogate and Hybrid
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 5  -E 1 -T 5 --warmstart greedy zero # Surrogate and Greedy
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 20  -E 1 -T 20 --full --hybrid --warmstart greedy zero # Hybrid and warmstart
python EscortFlowSim_v7.py --seed 0 -x 13 -y 7 -O 6 0 -e 16 -f FullFactorial13x7.csv -S 1600 -R 0.4 -I 5  -E 1 -T 5 --warmstart greedy zero --hybrid # Surrogate and Greedy
