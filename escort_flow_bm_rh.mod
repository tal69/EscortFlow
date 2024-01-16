/*********************************************
 * OPL 22.11 Model
 * Author: Tal Raviv
 * Creation Date: Aug 4, 2020,  Updated March 16, 2023, December 2023
 * The escort flow model (multi-loads block movements, stay/leave/continue)
 * Lighter model with neighborhood sets 
 * 
 * 11/2023 Rolling horizon versio with penalty for load location and no arrival constraint
 * 13/12/2023 bm version 
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


// ********************  Construct network and neighborhood sets ***********************

{Loc} NA[locations];
{Loc} NE[locations];
{Move} movesA = {};
{Move} movesE = {};


execute  {
	  for (var l in locations) {
			NA[l].add(l);
			movesA.add(l.x,l.y,l.x,l.y, 0);  // stay still 
			if (l.x< Lx-1) {
			 	NA[l].add(l.x+1,l.y);
			 	movesA.add(l.x,l.y,l.x+1,l.y, gamma); // right		 
			}
			if (l.y< Ly-1)  {
			  	NA[l].add(l.x,l.y+1); 
			  	movesA.add(l.x,l.y,l.x,l.y+1, gamma);  // up
			}
			if(l.x>0){  
			 	NA[l].add(l.x-1,l.y);
			 	movesA.add(l.x,l.y,l.x-1,l.y, gamma);  // left 
			}	 
		    if(l.y>0) { 
		      	NA[l].add(l.x,l.y-1)
		      	movesA.add(l.x,l.y,l.x,l.y-1, gamma);  //down
		   }
		   for (var x in Xr) {
		     movesE.add(l.x,l.y,x,l.y, Math.abs(l.x-x));
			 NE[l].add(x,l.y);
			 			   	
		   }
		   
		   for (var y in Yr) { 
			   	if (y != l.y)  {
				     movesE.add(l.x,l.y,l.x,y, Math.abs(l.y-y));
					 NE[l].add(l.x,y);
	 			}				   	
		   }     
	  }
	  
}

{Move} CellCover[locations];
{Move} MoveCover[movesA];  // Continue with the creation of these sets

execute {
  for (var m in movesE) {
    if (m.orig_y == m.dest_y)  {// horizontal movment including idle
    	for (var x = Math.min(m.orig_x, m.dest_x); x <= Math.max(m.orig_x, m.dest_x); x++)  CellCover[locations.find(x, m.orig_y)].add(m);
    	for (x = m.orig_x; x < m.dest_x; x++) MoveCover[movesA.find(x+1, m.orig_y, x, m.orig_y)].add(m);  // right movement	
    	for (x = m.orig_x; x > m.dest_x; x--) MoveCover[movesA.find(x-1, m.orig_y, x, m.orig_y)].add(m);  // left movement	
    	
  	}    
    if (m.orig_y != m.dest_y)  // vertical movment
    	for (var y = Math.min(m.orig_y, m.dest_y); y <= Math.max(m.orig_y, m.dest_y); y++)  CellCover[locations.find(m.orig_x,y)].add(m);
    	for (y = m.orig_y; y < m.dest_y; y++) MoveCover[movesA.find(m.orig_x, y+1, m.orig_x, y)].add(m);  // up movement	
    	for (y = m.orig_y; y > m.dest_y; y--) MoveCover[movesA.find(m.orig_x, y-1, m.orig_x, y)].add(m);  // down movement	
    	
    
  }    
    
    
  
}


// ****************************************************************



// Commodity 1 is target loads, and 2 is escorts 
float supply[l in locations][i in 1..2] = 0;

execute {
  for (var l in A) supply[l][1] = 1; // supply of retrieved loads   
  for (l in E) supply[l][2] = 1;  // supply of escorts 
     
  
}

dvar float xE[movesE,Tr] in 0..1;  // flow on escort movement arcs 
dvar float xA[movesA,Tr] in 0..1;  // flow on target load movement 

dvar boolean xE_int[movesE, 0..T_int];
dvar boolean xA_int[movesA, 0..T_int];

dexpr float NumberOfMovements = (sum(m in movesE, t in 0..(T_exec-1) )  m.cost*xE[m,t]);

dexpr float makespan_rh = max(t in 0..(T_exec-1),l1 in not_O, l2 in NA[l1]) (t+1)*xA[<l1.x, l1.y, l2.x, l2.y>,t];
dexpr float flowtime_rh = sum(l2 in O, l1 in NA[l2]: l1 != l2) sum( t in 0..(T_exec-1)) (t+1)*xA[<l1.x, l1.y, l2.x, l2.y>,t];   




// ******  Initialize MipStart (advanced solution from heuristic) ********

int warm_start = ...;

//int xA_val[movesA,0..T_int];
int xE_val[movesE,0..T_int];
int xE_val0[movesE,0..T_int];
//int xA_val0[movesA,0..T_int];

tuple tInitSol {
  int x1;
  int y1;
  int x2;
  int y2;
  int t;
}


// *****  Warmstart for single retrieval-stay problems
{tInitSol} init_solA = ...;
{tInitSol} init_solE = ...;


execute {
    if (warm_start == 1) { 
	 	//for (var a in init_solA)  xA_val[movesA.get(a.x1,a.y1,a.x2,a.y2)][a.t] = 1;
	 	for (var a in init_solE)  xE_val[movesE.get(a.x1,a.y1,a.x2,a.y2)][a.t] = 1;
	 	
	 	
	    //cplex.addMIPStart(xA_int, xA_val);
	    cplex.addMIPStart(xE_int, xE_val);	  
  	}    
  	
  	for (e in E) {
  	  for(var t = 0; t<= T_int; t++) {
  	     xE_val0[movesE.get(e.x,e.y,e.x,e.y)][t] = 1;
  	  }
  	}
  	
//  	for (e in A) {
//  	  for(var t = 0; t<= T_int; t++) {
//  	     xA_val0[movesA.get(e.x,e.y,e.x,e.y)][t] = 1;
//  	  }
//  	}
//  	cplex.addMIPStart(xA_int, xA_val0);
	cplex.addMIPStart(xE_int, xE_val0);	  
  	
}
 
 
minimize sum(m in movesE, t in Tr) m.cost*xE[m,t] + sum(m in movesA, t in Tr)  p[<m.dest_x, m.dest_y>] * xA[m,t]; 


subject to
{
  
  //debug: forall( m in moves, t in Tr, k in 1..2) x[m,t,k] == x_val[m,t,k];
  
  // Partial LP relaxation
  forall ( m in movesA, t in 0..T_int) xA_int[m,t] == xA[m,t];
  forall ( m in movesE, t in 0..T_int) xE_int[m,t] == xE[m,t];

 
 
  // Flow conservation at nodes for retrieved loads (both commodities)
  if (retrieval_mode == "stay") {
    	 // (2) in the paper
  		flow_conservation_stay_A: forall ( l1 in locations, t in 1..T)   
  	 	sum( l2 in NA[l1]) xA[<l2.x, l2.y, l1.x, l1.y>,t-1]  == 
  	 	sum(l2 in NA[l1]) xA[<l1.x, l1.y, l2.x, l2.y>,t];
  	 	
  	 	flow_conservation_stay_E: forall ( l1 in locations, t in 1..T)   
  	 	sum( l2 in NE[l1]) xE[<l2.x, l2.y, l1.x, l1.y>,t-1]  == 
  	 	sum(l2 in NE[l1]) xE[<l1.x, l1.y, l2.x, l2.y>,t];
  	 	
  	 	
  	 	
}  
  
  if (retrieval_mode == "continue") {  // described in the text at the end of Section 3.1
  		flow_conservation_cont1: forall ( l1 in locations, t in 1..T) 
  	 	sum( l2 in NE[l1]) xE[<l2.x, l2.y, l1.x, l1.y>,t-1]  == 
  	 	sum(l2 in NE[l1]) xE[<l1.x, l1.y, l2.x, l2.y>,t];
  	 	
  	 	flow_conservation_cont2: forall ( l1 in locations diff O, t in 1..T) 
  	 	sum( l2 in NA[l1]) xA[<l2.x, l2.y, l1.x, l1.y>,t-1]  == 
  	 	sum(l2 in NA[l1]) xA[<l1.x, l1.y, l2.x, l2.y>,t];
  }  
  
  if (retrieval_mode == "leave") {  // described in the text at the end of Section 3.1  equations (12)-(14)
    	
    	
    	flow_conservation_leave1: forall (t in 1..T, l1 in O) xA[<l1.x, l1.y, l1.x, l1.y>,t] == 
  	 		sum(l2 in NA[l1]: l2 != l1) xA[<l2.x, l2.y, l1.x, l1.y>,t-1] ;
  	 	
  	 	flow_conservation_leave2: forall (l1 in O, t in 1..T)   
  	 	xA[<l1.x, l1.y, l1.x, l1.y>,t-1] + sum( l2 in NE[l1]) xE[<l2.x, l2.y, l1.x, l1.y>,t-1]  ==  
  	 		sum(l2 in NE[l1]) xE[<l1.x, l1.y, l2.x, l2.y>,t];
  	 	
  	 	flow_conservation_leave3: forall ( l1 in locations diff O, t in 1..T)   
  	 		sum( l2 in NA[l1]) xA[<l2.x, l2.y, l1.x, l1.y>,t-1]  == sum(l2 in NA[l1]) xA[<l1.x, l1.y, l2.x, l2.y>,t];
  	 	
  	 	flow_conservation_leave4: forall ( l1 in locations diff O, t in 1..T)   
  	 		sum( l2 in NE[l1]) xE[<l2.x, l2.y, l1.x, l1.y>,t-1]  == sum(l2 in NE[l1]) xE[<l1.x, l1.y, l2.x, l2.y>,t];
  	 		
  }
  
  // (3) in the paper - target loads never leaves the output cell
  target_load_never_leaves: forall (t in Tr, l1 in O) sum(l2 in NA[l1]: l1 != l2) xA[<l1.x, l1.y, l2.x, l2.y>,t] == 0;
    
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
	   for (var m in thisOplModel.movesA) {
	 		if (thisOplModel.xA[thisOplModel.movesA.find(m)][thisOplModel.T_exec] > 0.99) {
	 			f.write("(", m.orig_x ,",",m.orig_y,"),")
	    	}
	   }    
	   f.writeln("]");
	   
	  
	   f.write("[");
	   for (var m in thisOplModel.movesE) {
	 		if (thisOplModel.xE[thisOplModel.movesE.find(m)][thisOplModel.T_exec] > 0.99) {
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
	   	
     
     if (thisOplModel.file_export != "") {
       
       // Write load movements into file
	   	var f = new IloOplOutputFile(thisOplModel.file_export);
	   	f.write("[")
	  	for (var t = 0; t< thisOplModel.T_exec; t++) {
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
 			if  (thisOplModel.retrieval_mode == "leave") {	 
				for (var m in thisOplModel.movesA){
				  	for (var l in thisOplModel.O) {
				  	  if (t > 0 && thisOplModel.xA[thisOplModel.movesA.find(m)][t-1] > 0.99 
				  	  && m.dest_x == l.x && m.dest_y == l.y && (m.orig_x != l.x || m.orig_y != l.y))  	  
				  	  	f.write("((", m.dest_x ,",",m.dest_y ,"), (None, None)),");
				  	  if (t==0 && thisOplModel.xA[thisOplModel.movesA.find(m)][0] > 0.99 && m.orig_x == l.x && m.orig_y == l.y)  
					   		  f.write("((", l.x ,",",l.y ,"), (None, None)),");
					}			  	
				}
 			}	
 									
			f.write("],");
 		}
 		f.write("]");	  	
	
	  
	 	f.close();
		   
	} // if	file_export != ""
	
	if (thisOplModel.warm_start == 1)  {  //  **** for warm starting next horizon ****
			f = new IloOplOutputFile("warm_start.dat");
			f.write("init_solA = {");
			for (var t=thisOplModel.T_exec ;t < thisOplModel.T_exec+ thisOplModel.T_int; t++) {
			  	f.writeln("  ");
				for (m in thisOplModel.movesA) if (thisOplModel.xA[thisOplModel.movesA.find(m)][t] > 0.99) {
				  f.write("<",m.orig_x," ",m.orig_y," ",m.dest_x," ",m.dest_y," ",t-thisOplModel.T_exec,"> ");
				}
 			}			 
 			f.writeln("};");
 			
 			
			f.write("init_solE = {");
			for (var t=thisOplModel.T_exec ;t < thisOplModel.T_exec+ thisOplModel.T_int; t++) {
			  	f.writeln("  ");
				for (m in thisOplModel.movesE) if (thisOplModel.xE[thisOplModel.movesE.find(m)][t] > 0.99) {
				  f.write("<",m.orig_x," ",m.orig_y," ",m.dest_x," ",m.dest_y," ",t-thisOplModel.T_exec,"> ");
				}
 			}			 
 			f.writeln("};");
 			
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
