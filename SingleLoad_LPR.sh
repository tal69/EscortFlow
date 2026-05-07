python EscortFlowStatic.py -x 9 -y 5 -O 4 0 -e 3-8 -l 1  -r 1-100 --gurobi -f table1a_escortflowlp.csv  --lp
python EscortFlowStatic.py -x 13 -y 7 -O 6 0 -e 3-8 -l 1  -r 1-100 --gurobi -f table1a_escortflowlp.csv  --lp
python EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 3-8 -l 1  -r 1-100 --gurobi -f table1a_escortflowlp.csv  --lp
python EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 3-8 -l 1  -r 1-100 --gurobi -f table1a_escortflowlp.csv --lp
python EscortFlowStatic.py -x 27 -y 10 -O 4 0 13 0 22 0 -e 3-8 -l 1  -r 1-100 --gurobi -f table1a_escortflowlp.csv --lp

python LoadFlowStatic.py -x 9 -y 5 -O 4 0 -e 3-8 -l 1  -r 1-100 --gurobi -f table1a_loadflowlp.csv  --lp
python LoadFlowStatic.py -x 13 -y 7 -O 6 0 -e 3-8 -l 1  -r 1-100 --gurobi -f table1a_loadflowlp.csv  --lp
python LoadFlowStatic.py -x 10 -y 10 -O 0 0 -e 3-8 -l 1  -r 1-100 --gurobi -f tabel1a_loadflowlp.csv  --lp
python LoadFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 3-8 -l 1  -r 1-100 --gurobi -f tabel1a_loadflowlp.csv  --lp
python LoadFlowStatic.py -x 27 -y 10 -O 4 0 13 0 22 0 -e 3-8 -l 1  -r 1-100 --gurobi -f table1a_loadflowlp.csv --lp