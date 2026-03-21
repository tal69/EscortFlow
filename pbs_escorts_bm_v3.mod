/*********************************************
 * OPL 22.11 Model
 * Author: Tal Raviv
 * Creation Date: Aug 4, 2020,  Updated March 16, 2023, December 2023, December 2025 (v3)
 * The escort flow model (multi-loads, block movements, stay/leave/continue)
 *********************************************/


float temp;

float time_limit = ...;
int threads = ...;

execute {cplex.tilim = time_limit;
		cplex.threads = threads;
		var before = new Date();
		temp = before.getTime();} 


string file_export = ...;
string file_res = ...;

float beta = ...; // flowtime weight - can be zero
float gamma = ...; // movement weight  (assumed > 0)

tuple Loc {
  key int x;
  key int y;
} 

 tuple Move {
   key int orig_x;
   key int orig_y;
   key int dest_x;
   key int dest_y;
}


int Lx = ...;
int Ly = ...;
int T  = ...;

range Tr = 0..T;
range Xr = 0..Lx-1;
range Yr = 0..Ly-1;


{Loc} locations = ...; // L in the paper

{Move} movesA = ...;
{Move} movesE = ...;

{Loc} E = ...; // set of initial escort  initial locations
{Loc} A = ...;  // Set of retrieved item initial locations
{Loc} O = ...;  // Output points locations
{Loc} not_O = locations diff O;

string retrieval_mode = ...;  // after retrieval load may  'stay', 'leave', or 'continue'

int loads_to_retrieve = card(A diff O);

{Loc} NA[locations] = ...;
{Loc} NE[locations] = ...;


{Move} CellCover[locations] = ...;
{Move} MoveCover[movesA] = ...;



// Commodity 1 is the target loads, and 2 is the escorts 
float supply[l in locations][i in 1..2] = 0;

execute {
  for (var l in A) supply[l][1] = 1;   // supply of retrieved loads   
  for (l in E) supply[l][2] = 1;  // supply of escorts
}


dvar boolean xA[movesA,Tr];  // flow of target loads
dvar boolean xE[movesE,Tr];  // flow of target loads
dvar int+ q[O];  // auxiliary to facilitate cuts


dexpr float NumberOfMovements = (sum(m in movesE, t in Tr )  (abs(m.orig_x-m.dest_x)+ abs(m.orig_y-m.dest_y))*xE[m,t]);
dexpr float FlowTime = sum(l2 in O)  sum(l1 in NA[l2] : l2 != l1) sum(t in 1..T) (t+1)* xA[<l1.x, l1.y, l2.x, l2.y>,t];  
dexpr float CalcMakespan = max(l2 in O, l1 in NA[l2], t in Tr : l2 != l1) (t+1)* xA[<l1.x, l1.y, l2.x, l2.y>,t];  


minimize gamma * NumberOfMovements + beta * sum(l in O) q[l];
subject to
{
  	 		 
  if (retrieval_mode == "stay") {  // not in the paper
    // Flow conservation at nodes for escorts 
          escort_flow_conservation_stay1: forall ( l1 in locations, t in 1..T)
  	 		sum(l2 in NE[l1]) xE[<l2.x, l2.y, l1.x, l1.y>,t-1] == sum(l2 in NE[l1]) xE[<l1.x, l1.y, l2.x, l2.y>,t];
    
    	 // (2) Flow conservation at nodes for loads 
  		load_flow_conservation_stay2: forall ( l1 in locations, t in 1..T)
  	 		sum( l2 in NA[l1]) xA[<l2.x, l2.y, l1.x, l1.y>,t-1] == sum(l2 in NA[l1]) xA[<l1.x, l1.y, l2.x, l2.y>,t]; 	 	
  }  
  
  if (retrieval_mode == "continue") {  // described in the text at the end of Section 3.1
  
       // Flow conservation at nodes for escorts 
          escort_flow_conservation_cont: forall ( l1 in locations, t in 1..T)     
  	 		sum(l2 in NE[l1]) xE[<l2.x, l2.y, l1.x, l1.y>,t-1] == sum(l2 in NE[l1]) xE[<l1.x, l1.y, l2.x, l2.y>,t];
  		
  	 	flow_conservation_cont2: forall ( l1 in locations diff O, t in 1..T) 
  	 	sum( l2 in NA[l1]) xA[<l2.x, l2.y, l1.x, l1.y>,t-1]  == 
  	 	sum(l2 in NA[l1]) xA[<l1.x, l1.y, l2.x, l2.y>,t];
  }  
  
  if (retrieval_mode == "leave") {  // described in the text at the end of Section 3.1
  		flow_conservation_leave1: forall (t in 1..T, l1 in O) xA[<l1.x, l1.y, l1.x, l1.y>,t] == 
  	 		sum(l2 in NA[l1]: l2 != l1) xA[<l2.x, l2.y, l1.x, l1.y>,t-1] ;
  	 	
  	 	flow_conservation_leave2: forall (l1 in O, t in 1..T)   
  	 	xA[<l1.x, l1.y, l1.x, l1.y>,t-1] + sum( l2 in NE[l1]) xE[<l2.x, l2.y, l1.x, l1.y>,t-1]  ==  
  	 		sum(l2 in NE[l1]) xE[<l1.x, l1.y, l2.x, l2.y>,t];
  	 	
  
  	 	flow_conservation_leave3: forall ( l1 in locations diff O, t in 1..T)   
  	 	sum( l2 in NA[l1]) xA[<l2.x, l2.y, l1.x, l1.y>,t-1]  == 
  	 		sum(l2 in NA[l1]) xA[<l1.x, l1.y, l2.x, l2.y>,t];
  	 		
  	 	escort_flow_conservation_leave4: forall ( l1 in locations diff O, t in 1..T)     
  	 		sum(l2 in NE[l1]) xE[<l2.x, l2.y, l1.x, l1.y>,t-1] == sum(l2 in NE[l1]) xE[<l1.x, l1.y, l2.x, l2.y>,t];
  	 	
  }
  
  // (3) in the paper - target loads stay at the IOs
  target_stay: forall (t in Tr, l1 in O) sum(l2 in NA[l1] : l2 != l1) xA[<l1.x, l1.y, l2.x, l2.y>,t] == 0;
    
  // (4) and (5) in the paper - supply (initial locations)
  supplyA: forall ( l1 in locations) sum(l2 in NA[l1]) xA[<l1.x, l1.y, l2.x, l2.y>,0] == supply[l1,1];
  supplyB: forall ( l1 in locations) sum(l2 in NE[l1]) xE[<l1.x, l1.y, l2.x, l2.y>,0] == supply[l1,2];



  // (6) in the paper - Generalized capacity constraint

  Capacity: forall ( l1 in locations, t in Tr )  
  	sum(l2 in NA[l1]) xA[<l1.x, l1.y, l2.x, l2.y>,t]  + sum(l2 in NE[l1]) xE[<l1.x, l1.y, l2.x, l2.y>,t]  <= 1;
  
  
    // (7) in the paper - constraint - avoid conflicts
  
	avoid_conflicts: forall ( l in locations, t in Tr)  sum(m in CellCover[l]) xE[m,t] <= 1;
	
	
	//  target load movements allowed
	target_load_movements: forall(m in movesA, t in Tr : <m.orig_x, m.orig_y> != <m.dest_x, m.dest_y>) 
		xA[m,t] <= sum(m1 in MoveCover[m])  xE[m1 ,t];	
	
  //  target load movements enforced
  forall(l in locations, t in Tr) 1 - xA[<l.x, l.y, l.x, l.y>, t] >= sum(m in CellCover[l]) xE[m,t];
  	
	
  // (9) in the paper - all target loads arrive at output cells
  all_arrive: sum(l2 in O, l1 in NA[l2] : l2 != l1) sum(t in Tr) xA[<l1.x, l1.y, l2.x, l2.y>,t] == loads_to_retrieve;


	// Integrality cut 
	forall(l2 in O ) sum(t in Tr, l1 in NA[l2]: l1 != l2) (t+1)* xA[<l1.x, l1.y, l2.x, l2.y> ,t] == q[l2]; 

}


main {
   var before = new Date();
   thisOplModel.generate();
   var success = cplex.solve()
   
   var after = new Date();
   var CpuTime = (after.getTime() - before.getTime())/1000;
   var f_res = new IloOplOutputFile(thisOplModel.file_res, true); 	
  
   if (success) {
     if (thisOplModel.file_export != "") {
       
       // Write load movements into file
	   	var f = new IloOplOutputFile(thisOplModel.file_export);
	   	f.write("[")
	  	for (var t = 0; t<= thisOplModel.CalcMakespan; t++) {
			f.write("[");		
			for (var m in thisOplModel.movesE){
				if ((thisOplModel.xE[thisOplModel.movesE.find(m)][t] > 0.99) && (m.orig_x != m.dest_x || m.orig_y != m.dest_y)) {
								
					if (m.dest_x < m.orig_x )  // escort move left, loads right
					    for(var x = m.dest_x; x < m.orig_x; x++) f.write("((", x ,",",m.dest_y ,"), (", x+1 ,",",m.orig_y ,")),");
					    
					if (m.dest_x > m.orig_x )  // escort move right, loads left
					    for(var x = m.orig_x; x < m.dest_x; x++) f.write("((", x+1 ,",",m.dest_y ,"), (", x ,",",m.orig_y ,")),");
					    
					if (m.dest_y < m.orig_y )  // escort move down, loads up
					    for(var y = m.dest_y; y < m.orig_y; y++) f.write("((", m.dest_x ,",",y ,"), (", m.orig_x ,",",y+1 ,")),"); 
					    
					if (m.dest_y > m.orig_y )  // escort move up, loads down
					    for(var y = m.orig_y; y < m.dest_y; y++) f.write("((", m.dest_x ,",",y+1 ,"), (", m.orig_x ,",",y ,")),");  	 	    
			  	} 
 			}		
 			if  (thisOplModel.retrieval_mode == "leave" && t>0) {	 
				for (var m in thisOplModel.movesA){
				  	for (var l in thisOplModel.O) {
				  	  if (thisOplModel.xA[thisOplModel.movesA.find(m)][t-1] > 0.99 
				  	  && m.dest_x == l.x && m.dest_y == l.y && (m.orig_x != l.x || m.orig_y != l.y))  	  
				  	  	f.write("((", m.dest_x ,",",m.dest_y ,"), (None, None)),");
				  	}			  	
				}
 			}							
			f.write("],");
 		}
 		f.write("]");	  	
	
	  
	 	f.close();
		   
	} // if	file_export != ""
	
	
	f_res.write(",", thisOplModel.CalcMakespan);
	
	
	f_res.write( ",", thisOplModel.FlowTime, ",", thisOplModel.NumberOfMovements,",", 
	cplex.getObjValue(), ",",cplex.getBestObjValue(),",",CpuTime, ",,", 0, ",",cplex.getCplexStatus() );

  
  }  
  else { // could not find a feasible solution
  
	  f_res.write(",-,-,-,-,",cplex.getBestObjValue(),",",CpuTime, ",,", 0, ",",cplex.getCplexStatus()  );	
	
  } 
  f_res.close(); 
  
} 
