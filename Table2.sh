python EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 8:4:20 -l 4  -r 1-100 --gurobi -f table2_escortflow.csv  --time_limit=300
python EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 8:4:20 -l 4  -r 1-100 --gurobi -f table2_escortflow.csv --time_limit=300
python EscortFlowStatic.py -x 27 -y 10 -O 4 0 13 0 22 0  -e 8:4:20 -l 4  -r 1-100 --gurobi -f table2_escortflow.csv --time_limit=300

python EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 8:4:20 -l 4  --bnc -r 1-100 --gurobi -f table2_escortflow_bnc.csv   --time_limit=300
python EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 8:4:20 -l 4 --bnc  -r 1-100 --gurobi -f table2_escortflow_bnc.csv --time_limit=300
python EscortFlowStatic.py -x 27 -y 10 -O 4 0 13 0 22 0  -e 8:4:20 -l 4 --bnc -r 1-100 --gurobi -f table2_escortflow.csv --time_limit=300

python LoadFlowStatic.py -x 10 -y 10 -O 0 0 -e 8:4:20 -l 4 -r 1-100 --gurobi -f table2_loadflow.csv --time_limit=300
python LoadFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 8:4:20 -l 4  -r 1-100 --gurobi -f table2_loadflow.csv --time_limit=300
python LoadFlowStatic.py -x 27 -y 10 -O 4 0 13 0 22 0  -e 8:4:20 -l 4  -r 1-100 --gurobi -f table2_loadflow.csv --time_limit=300