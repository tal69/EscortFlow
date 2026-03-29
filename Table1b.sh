python EscortFlowStatic.py -x 9 -y 5 -O 4 0  -e 8:4:16 -l 4 --bnc  -r 1-100 --gurobi -f table1b_escortflow_lp.csv --lp
python EscortFlowStatic.py -x 13 -y 7 -O 6 0   -e 8:4:16 -l 4 --bnc -r 1-100 --gurobi -f table1b_escortflow_lp.csv --lp
python EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 8:4:16 -l 4  -r 1-100 --gurobi -f table1b_escortflow_lp.csv  --lp
python EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 8:4:16 -l 4  -r 1-100 --gurobi -f table1b_escortflow_lp.csv --lp
python EscortFlowStatic.py -x 27 -y 10 -O 4 0 13 0 22 0  -e 8:4:16 -l 4  -r 1-100 --gurobi -f table1b_escortflow_lp.csv --lp

python LoadFlowStatic.py -x 9 -y 5 -O 4 0 -e 8:4:16 -l 4 -r 1-100 --gurobi -f table1b_loadflow.csv --lp
python LoadFlowStatic.py -x 13 -y 7 -O 6 0 -e 8:4:16 -l 4 -r 1-100 --gurobi -f table1b_loadflow.csv --lp
python LoadFlowStatic.py -x 10 -y 10 -O 0 0 -e 8:4:16 -l 4 -r 1-100 --gurobi -f table1b_loadflow.csv --lp
python LoadFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 8:4:16 -l 4  -r 1-100 --gurobi -f table1b_loadflow.csv --lp
python LoadFlowStatic.py -x 27 -y 10 -O 4 0 13 0 22 0  -e 8:4:16 -l 4  -r 1-100 --gurobi -f table1b_loadflow.csv --lp