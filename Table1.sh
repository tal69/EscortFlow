python EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 3-8 -l 1  -r 1-100 --gurobi -f table1_escortflow.csv  --time_limit=300
python EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 3-8 -l 1  -r 1-100 --gurobi -f table1_escortflow.csv --time_limit=300

python EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 3-8 -l 1  --bnc -r 1-100 --gurobi -f table1_escortflow_bnc.csv   --time_limit=300
python EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 3-8 -l 1 --bnc  -r 1-100 --gurobi -f table1_escortflow_bnc.csv --time_limit=300

python LoadFlowStatic.py -x 10 -y 10 -O 0 0 -e 3-8 -l 1  -r 1-100 --gurobi -f tabel1_loadflow.csv  --time_limit=300
python LoadFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 3-8 -l 1  -r 1-100 --gurobi -f tabel1_loadflow.csv  --time_limit=300