/*********************************************
 * OPL 22.11 Model
 * Author: Tal Raviv
 * Creation Date: Aug 4, 2020,  Updated March 16, 2023, December 2023
 * The escort flow model (multi-loads, non block movements, stay/leave/continue)
 * Lighter model with neighborhood sets 
 * Now with warmstart for single load (4/12/2023)
 *********************************************/

float temp;

float time_limit = ...;

execute {cplex.tilim = time_limit;
		var before = new Date();
		temp = before.getTime();} 

string file_export = ...;
string file_res = ...;

float alpha = ...; // C_max weight - can be zero
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
   float cost;
}

int Lx = ...;
int Ly = ...;
int T  = ...;

range Tr = 0..T;
range Xr = 0..Lx-1;
range Yr = 0..Ly-1;

{Move} moves = {};

{Loc} locations = {}; // L in the paper

{Loc} E = ...; // set of initial escort  initial locations
{Loc} A = ...;  // Set of retrived item initial locations
{Loc} O = ...;  // Output points locations  - NOT SURE WE NEED THIS
{Loc} not_O = locations diff O;
//{Loc} A_not_O = A diff O;

string retrieval_mode = ...;  // after reterival load may  'stay', 'leave', or 'continue'

int loads_to_retreive = card(A diff O);

execute {
  for (var x in Xr) for (var y in Yr) locations.add(x,y);
}   

{Loc} N0[locations];


execute  {
  for (var x in Xr) for (var y in Yr) {
    if (x< Lx-1) moves.add(x,y,x+1,y, gamma); // right		 
    if (y< Ly-1) moves.add(x,y,x,y+1, gamma);  // up
    if(x>0)moves.add(x,y,x-1,y, gamma);  // left 
    if(y>0) moves.add(x,y,x,y-1, gamma);  //down
    moves.add(x,y,x,y, 0);  // stay still 
  }
}

execute  {
  for (var l in locations) {
	N0[l].add(l);
	if (l.x< Lx-1) N0[l].add(l.x+1,l.y)
	if (l.y< Ly-1)  N0[l].add(l.x,l.y+1)
	if(l.x>0)   N0[l].add(l.x-1,l.y)
   if(l.y>0)   N0[l].add(l.x,l.y-1)
  }
}



// Commodity 1 is target loads, and 2 is escorts 
float supply[l in locations][i in 1..2] = 0;

execute {
  for (var l in A) supply[l][1] = 1; // supply of retrieved loads   
  for (l in E) supply[l][2] = 1;  // supply of escorts 
     
  
}


dvar boolean x[moves,Tr,1..2];  // flow on movement arcs - 1 target loads, 2 - escorts
dvar int+ z;  // makespan if alpha> 0   (it could be float but it better for the solver to have it as int)
dvar int+ q[O];  // auxilary to faciltate integrality cuts 


dexpr float NumberOfMovements = (sum(m in moves, t in Tr )  m.cost*x[m,t,2])/gamma;
dexpr float FlowTime = sum(l2 in O)  sum(l1 in N0[l2] : l2 != l1) sum(t in 1..T) (t+1)* x[<l1.x, l1.y, l2.x, l2.y>,t,1];  
dexpr float CalcMakespan = max(l2 in O, l1 in N0[l2], t in Tr : l2 != l1) (t+1)* x[<l1.x, l1.y, l2.x, l2.y>,t,1];  


// ******  Initialize MipStart (advanced solution from heuristic) ********

int x_val[moves,Tr, 1..2];

tuple tInitSol {
  int x1;
  int y1;
  int x2;
  int y2;
  int t;
  int k;
}

// *****  Warmstart for single retrieval-stay problems
{tInitSol} init_sol = ...;
execute {
    if (Opl.card(A)==1 & Opl.card(init_sol) > 0) { 
	 	for (var a in init_sol)  x_val[moves.get(a.x1,a.y1,a.x2,a.y2)][a.t][a.k] = 1;
	    cplex.addMIPStart(x, x_val);
  }    
}


minimize alpha* z + sum(m in moves, t in Tr)  m.cost  *x[m,t,2] + beta * sum(l in O) q[l];

//minimize alpha* z + sum(m in moves, t in Tr)  m.cost  *x[m,t,2] + beta * 
//   sum(m in moves: <m.orig_x, m.orig_y> in locations diff O ) sum(t in Tr) x[m,t,1];
// I am not sure why but it seems that the first formulation of the OF converges faster - maybe somthing with integrality cuts that the 
// solver generates

subject to
{
  
  //debug: forall( m in moves, t in Tr, k in 1..2) x[m,t,k] == x_val[m,t,k];

 
 
  // Flow conservation at nodes for retrieved loads (both commodities)
  if (retrieval_mode == "stay") {
    	 // (2) in the paper
  		flow_conservation_stay: forall ( l1 in locations, t in 1..T, k in 1..2)   
  	 	sum( l2 in N0[l1]) x[<l2.x, l2.y, l1.x, l1.y>,t-1, k]  == 
  	 	sum(l2 in N0[l1]) x[<l1.x, l1.y, l2.x, l2.y>,t,k];
  	 	
  	 	
}  
  
  if (retrieval_mode == "continue") {  // described in the text at the end of Section 3.1
  		flow_conservation_cont1: forall ( l1 in locations, t in 1..T) 
  	 	sum( l2 in N0[l1]) x[<l2.x, l2.y, l1.x, l1.y>,t-1, 2]  == 
  	 	sum(l2 in N0[l1]) x[<l1.x, l1.y, l2.x, l2.y>,t,2];
  	 	
  	 	flow_conservation_cont2: forall ( l1 in locations diff O, t in 1..T) 
  	 	sum( l2 in N0[l1]) x[<l2.x, l2.y, l1.x, l1.y>,t-1, 1]  == 
  	 	sum(l2 in N0[l1]) x[<l1.x, l1.y, l2.x, l2.y>,t,1];
  }  
  
  if (retrieval_mode == "leave") {  // described in the text at the end of Section 3.1  equations (12)-(14)
    	
    	
    	flow_conservation_leave1: forall (t in 1..T, l1 in O) x[<l1.x, l1.y, l1.x, l1.y>,t,1] == 
  	 		sum(l2 in N0[l1]: l2 != l1) x[<l2.x, l2.y, l1.x, l1.y>,t-1,1] ;
  	 	
  	 	flow_conservation_leave2: forall (l1 in O, t in 1..T)   
  	 	x[<l1.x, l1.y, l1.x, l1.y>,t-1, 1] + sum( l2 in N0[l1]) x[<l2.x, l2.y, l1.x, l1.y>,t-1, 2]  ==  
  	 		sum(l2 in N0[l1]) x[<l1.x, l1.y, l2.x, l2.y>,t, 2];
  	 	
  	 	flow_conservation_leave3: forall ( l1 in locations diff O, t in 1..T, k in 1..2)   
  	 	sum( l2 in N0[l1]) x[<l2.x, l2.y, l1.x, l1.y>,t-1, k]  == 
  	 		sum(l2 in N0[l1]) x[<l1.x, l1.y, l2.x, l2.y>,t,k];
  	 	
  }
  
  // (3) in the paper - target loads never leaves the output cell
  target_load_never_leaves: forall (t in Tr, l1 in O) sum(l2 in N0[l1]: l1 != l2) x[<l1.x, l1.y, l2.x, l2.y>,t,1] == 0;
    
  // (4) and (5) in the paper - supply (initial locations)
  supply_cons: forall ( l1 in locations, k in 1..2) sum(l2 in N0[l1]) x[<l1.x, l1.y, l2.x, l2.y>,0,k] == supply[l1,k];

  // (6) in the paper - Generalized capacity constraint
  capacity: forall ( l1 in locations, t in Tr )  
  	sum(l2 in N0[l1]) sum(i in 1..2) x[<l1.x, l1.y, l2.x, l2.y>,t,i]  <= 1;
  
 
  // (7) in the paper - constraint - avoid conflicts
	avoid_conflicts: forall ( l1 in locations, l2 in N0[l1] : l2!=l1 )  forall(t in Tr)
  	x[<l1.x, l1.y, l2.x, l2.y> ,t,2] + 
  	sum(l3 in N0[l2] )  x[<l2.x, l2.y, l3.x, l3.y>, t,2] <= 1;
	
	
	// (8) in the paper - target load movements
	target_load_movement: forall(m in moves, t in Tr : <m.orig_x, m.orig_y> != <m.dest_x, m.dest_y>) 
		x[m,t,1] <= x[<m.dest_x, m.dest_y, m.orig_x, m.orig_y>,t,2];	
	
	
  // (9) in the paper - all target loads arrive at output cells
  all_arrive: sum(l2 in O, l1 in N0[l2] : l2 != l1)
  sum(t in Tr) x[<l1.x, l1.y, l2.x, l2.y>,t,1] == loads_to_retreive;	
  
    
  // (10) in the paper - assign makespan to z
 
 	
	if (alpha > 0)   
	{
	  if (retrieval_mode == "stay" || card(A)==1) forall(l2 in O, l1 in N0[l2]: l1 != l2 ) sum(t in Tr) (t+1)* x[<l1.x, l1.y, l2.x, l2.y> ,t,1] + 1<= z;
	  else forall(l2 in O, l1 in N0[l2],t in Tr: l1 != l2 )  t* x[<l1.x, l1.y, l2.x, l2.y> ,t,1] + 1<= z;
	}			
	else alpha1:  z== 0; // just to supress the warning (z has never been used)

// Integrality cut 
int1: forall(l2 in O ) sum(t in Tr, l1 in N0[l2]: l1 != l2) (t+1)* x[<l1.x, l1.y, l2.x, l2.y> ,t,1] == q[l2]; 
int2: sum(m in moves: <m.orig_x, m.orig_y> in locations diff O ) sum(t in Tr) x[m,t,1] == sum( l in O) q[l];

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
		var lFirstTimeStep = true;
	  	for (var t = 0; t<= thisOplModel.CalcMakespan; t++) {
	  	  	if (lFirstTimeStep == false) {f.write(","); }
			else { lFirstTimeStep = false; }
			f.write("[");		
			var lFirstInTimeStep = true;
			for (var m in thisOplModel.moves){
				if ((thisOplModel.x[thisOplModel.moves.find(m)][t][2] > 0.99) && (m.orig_x != m.dest_x || m.orig_y != m.dest_y)) {
				
				  	if (lFirstInTimeStep == false) {f.write(","); }
					else { lFirstInTimeStep = false; }
			   		f.write("((", m.dest_x ,",",m.dest_y ,"), (", m.orig_x ,",",m.orig_y ,"))");
			   			    
			  	} 
			  	if (thisOplModel.retrieval_mode == "leave") {
				  	for (var l in thisOplModel.O) {
					  	if (t > 0 && thisOplModel.x[thisOplModel.moves.find(m)][t-1][1] > 0.99 && m.dest_x == l.x && 
					  	m.dest_y == l.y && (m.orig_x != l.x || m.orig_y != l.y))
					   	{
					   	  	  if (lFirstInTimeStep == false) {f.write(","); }
							  else { lFirstInTimeStep = false;  } 
					   		  f.write("((", m.dest_x ,",",m.dest_y ,"), (None, None))");
					   		  
					   	} 
					   	if (t==0 && thisOplModel.x[thisOplModel.moves.find(m)][0][1] > 0.99 && m.orig_x == l.x && m.orig_y == l.y) {
					   	   	  if (lFirstInTimeStep == false) {f.write(","); }
							  else { lFirstInTimeStep = false;  } 
					   		  f.write("((", l.x ,",",l.y ,"), (None, None))");
					   		
					   	}
	  				}
    			}	  							
			}
			f.write("]");
 		}
 		f.write("]");	  	
	
	  
	  f.close();
		   
	} // if	file_export != ""
	
	if (thisOplModel.alpha > 0) { 
		f_res.write(",", thisOplModel.z);
	} 	
	else {
		f_res.write(",", thisOplModel.CalcMakespan);
	}
	
	f_res.write( ",", thisOplModel.FlowTime, ",", thisOplModel.NumberOfMovements,",", 
	cplex.getObjValue(), ",",cplex.getBestObjValue(),",",CpuTime, ",",cplex.getCplexStatus() );

  
  }  
  else { // could not find a feasible solution
  
	  f_res.write(",-,-,-,-,",cplex.getBestObjValue(),",",CpuTime,  ",",cplex.getCplexStatus()  );	
	
  } 
  f_res.close(); 
  //f_output.close();
 }
