/*********************************************
 * OPL 22.11 Model
 * Author: Tal Raviv
 * Creation Date: Aug 4, 2020,  Updated March 16, 2023, December 2023
 * The escort flow model (multi-loads, non block movements, stay/leave/continue)
 * Lighter model with neighborhood sets 
 * Now with warmstart for single load (4/12/2023)
 * 11/2023 Rolling horizon versio with penalty for load location and no arrival constraint
 *********************************************/

float temp;

float time_limit = ...;

execute {cplex.tilim = time_limit;
		var before = new Date();
		temp = before.getTime();} 

string file_export = ...;

float distance_penalty = ...; // Penalty for distance  (large surgate c_max)
float time_penalty = ...; //  Penalty for flowtime (if 1  this is exactly flow time)
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
int T  = ...;   // total number of periods in the model including fractional
int T_int = ...;  //  total number of periods represented by boolean variables
int T_exec = ...; // Execution horizon, the state of the system at the end of this period is the output
// T > T_int > T_exec

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


float p[locations];  // location penalty

execute {  // calculate location penlaty 
  for (var l in locations) p[l] = (Lx+Ly+1)*(distance_penalty+ distance_penalty);
  
  for (var i in O) for (var l in locations) {
    var d = Math.abs(i.x-l.x) + Math.abs(i.y-l.y);
    if ((time_penalty + distance_penalty*d) < p[l] ) p[l] = distance_penalty+ distance_penalty*d;
    if(d==0) {
      if (retrieval_mode == 'leave')  {p[l] = gamma} else {p[l] = 0};  // make it leave
  }     
}
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

dvar float x[moves,Tr,1..2] in 0..1;  // flow on movement arcs - 1 target loads, 2 - escorts
dvar boolean x_int[moves, 0..T_int, 1..2];


dexpr float NumberOfMovements = (sum(m in moves, t in 0..(T_exec-1) )  m.cost*x[m,t,2])/gamma;
dexpr float makespan_rh = max(t in 0..(T_exec-1),l1 in not_O, l2 in N0[l1]) (t+1)*x[<l1.x, l1.y, l2.x, l2.y>,t,1];
dexpr float flowtime_rh = sum(l2 in O, l1 in N0[l2]: l1 != l2) sum( t in 0..(T_exec-1)) (t+1)*x[<l1.x, l1.y, l2.x, l2.y>,t,1];   

// ******  Initialize MipStart (advanced solution from heuristic) ********

int warm_start = ...;

int x_val[moves,0..T_int, 1..2] ;

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
    if (warm_start == 1) { 
	 	for (var a in init_sol)  x_val[moves.get(a.x1,a.y1,a.x2,a.y2)][a.t][a.k] = 1;
	    cplex.addMIPStart(x_int, x_val);	  
  	}    
}


minimize sum(m in moves, t in Tr)  (m.cost   *x[m,t,2] + p[<m.dest_x, m.dest_y>] * x[m,t,1]); 


subject to
{
  
  //debug: forall( m in moves, t in Tr, k in 1..2) x[m,t,k] == x_val[m,t,k];
  
  // Partial LP relaxation
  forall ( m in moves, t in 0..T_int, k in 1..2) x_int[m,t,k] == x[m,t,k];

 
 
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
	
 }




main {
   var before = new Date();
   thisOplModel.generate();
   var success = cplex.solve()
   var after = new Date();
   var CpuTime = (after.getTime() - before.getTime())/1000;
    
   if (success) {
	   var f  = new IloOplOutputFile("end_of_exec_horizon.txt");
	   f.write("[");
	   for (var m in thisOplModel.moves) {
	 		if (thisOplModel.x[thisOplModel.moves.find(m)][thisOplModel.T_exec][1] > 0.99) {
	 			f.write("(", m.orig_x ,",",m.orig_y,"),")
	    	}
	   }    
	   f.writeln("]");
	   
	  
	   f.write("[");
	   for (var m in thisOplModel.moves) {
	 		if (thisOplModel.x[thisOplModel.moves.find(m)][thisOplModel.T_exec][2] > 0.99) {
	 			f.write("(", m.orig_x ,",",m.orig_y,"),")
	    	}
	   }    
	   f.writeln("]");
	   f.writeln(cplex.getCplexStatus());
	   f.writeln(CpuTime);
	   f.writeln(thisOplModel.makespan_rh);
	   f.writeln(thisOplModel.flowtime_rh);
	   f.writeln(thisOplModel.NumberOfMovements);
	  
	   f.writeln(cplex.getObjValue());
	   f.writeln(cplex.getBestObjValue());
	   f.close()
	   	
     
     if (thisOplModel.file_export != "") {  // for animation 
       // Write load movements into file
	   	var f = new IloOplOutputFile(thisOplModel.file_export);
	   	f.write("[")
	  	for (var t = 0; t< thisOplModel.T_exec; t++) {
	 		f.write("[");		
			for (var m in thisOplModel.moves){
				if ((thisOplModel.x[thisOplModel.moves.find(m)][t][2] > 0.99) && (m.orig_x != m.dest_x || m.orig_y != m.dest_y)) {
			   		f.write("((", m.dest_x ,",",m.dest_y ,"), (", m.orig_x ,",",m.orig_y ,")),");
			  	} 
			  	if (thisOplModel.retrieval_mode == "leave") {
				  	for (var l in thisOplModel.O) {
					  	if (t > 0 && thisOplModel.x[thisOplModel.moves.find(m)][t-1][1] > 0.99 && m.dest_x == l.x && 
					  	m.dest_y == l.y && (m.orig_x != l.x || m.orig_y != l.y))
					   	{
					   		  f.write("((", m.dest_x ,",",m.dest_y ,"), (None, None)),"); 
					   	} 
					   	if (t==0 && thisOplModel.x[thisOplModel.moves.find(m)][0][1] > 0.99 && m.orig_x == l.x && m.orig_y == l.y) { 
					   		  f.write("((", l.x ,",",l.y ,"), (None, None)),");
					   		
					   	}
	  				}
    			}	  							
			}
			f.write("],");
 		}
 		f.write("]");	  	
	  f.close();
	} 
	
	if (thisOplModel.warm_start == 1)  {  // for warm starting next horizon
			var f = new IloOplOutputFile("warm_start.txt");
			f.writeln(thisOplModel.x);
			f.close()
	}
	
  }  
  else { // could not find a feasible solution
  	  
	  var f  = new IloOplOutputFile("end_of_exec_horizon.txt");  // write an empty end of horizon file
	  f.writeln("[]");
	  f.writeln("[]");
	  f.writeln(-1);  // in leiu of Cplex staus
	  f.writeln(CpuTime);
	  f.writeln(0);
	  f.writeln(0);
	  f.writeln(0);
	  f.writeln(0);
	  f.writeln(0);
	  f.close();
  } 
}
