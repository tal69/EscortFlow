python EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 3-8 -l 1  -r 1-100 --gurobi -f table1_escortflow.csv  --work_limit=500
python EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 3-8 -l 1  -r 1-100 --gurobi -f table1_escortflow.csv --work_limit=500

python EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 3-8 -l 1  --bnc -r 1-100 --gurobi -f table1_escortflow_bnc.csv   --work_limit=500
python EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 3-8 -l 1 --bnc  -r 1-100 --gurobi -f table1_escortflow_bnc.csv --work_limit=500

python LoadFlowStatic.py -x 10 -y 10 -O 0 0 -e 3-8 -l 1  -r 1-100 --gurobi -f tabel1_loadflow.csv  --work_limit=500
python LoadFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 3-8 -l 1  -r 1-100 --gurobi -f tabel1_loadflow.csv  --work_limit=500

python EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 8:4:20 -l 4  -r 1-100 --gurobi -f table2_escortflow.csv  --work_limit=500
python EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 8:4:20 -l 4  -r 1-100 --gurobi -f table2_escortflow.csv --work_limit=500
python EscortFlowStatic.py -x 27 -y 10 -O 4 0 13 0 22 0  -e 8:4:20 -l 4  -r 1-100 --gurobi -f table2_escortflow.csv --work_limit=500

python EscortFlowStatic.py -x 10 -y 10 -O 0 0 -e 8:4:20 -l 4  --bnc -r 1-100 --gurobi -f table2_escortflow_bnc.csv   --work_limit=500
python EscortFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 8:4:20 -l 4 --bnc  -r 1-100 --gurobi -f table2_escortflow_bnc.csv --work_limit=500
python EscortFlowStatic.py -x 27 -y 10 -O 4 0 13 0 22 0  -e 8:4:20 -l 4 --bnc -r 1-100 --gurobi -f table2_escortflow.csv --work_limit=500

python LoadFlowStatic.py -x 10 -y 10 -O 0 0 -e 8:4:20 -l 4 -r 1-100 --gurobi -f table2_loadflow.csv --work_limit=500
python LoadFlowStatic.py -x 16 -y 10 -O 4 0 11 0  -e 8:4:20 -l 4  -r 1-100 --gurobi -f table2_loadflow.csv --work_limit=500
python LoadFlowStatic.py -x 27 -y 10 -O 4 0 13 0 22 0  -e 8:4:20 -l 4  -r 1-100 --gurobi -f table2_loadflow.csv --work_limit=500