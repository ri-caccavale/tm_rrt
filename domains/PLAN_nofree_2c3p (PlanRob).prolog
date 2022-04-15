%
%PLANNING DOMAIN FOR 2-CARTS TASK_AND_MOTION 
%

getConstraint(X):-
	yield([],_),
	constraint(X,Constraint),
	yield(Constraint,_),!.

getEffect(X):-
	yield([],_),
	effect(X,Effect),
	yield(Effect,_),!.

getTarget(X):-
	yield([],_),
	target(X,Target),
	yield(Target,_),!.

%probably unused
getSchemaFromEffect(X):-
	yield([],_),
	findall(Schema,includeEffect(Schema,X),SchemaList),
	yield(SchemaList,_),!.

includeEffect(Schema,X):-
	result(Schema,Result),
	subset(X,Result),
	ground(Schema).



getTaskList:-
	yield([],_),
	findall(T,task(T),TaskList),
	yield(TaskList,_),!.

getVarList:-
	yield([],_),
	findall(V,variable(V),VarList),
	yield(VarList,_),!.

getPoseList:-
	yield([],_),
	findall(pose(P,X,Y,W),pose(P,X,Y,W),PoseList),
	yield(PoseList,_),!.

getObjectList:-
	yield([],_),
	findall(O,object(O,_),ObjectList),
	yield(ObjectList,_),!.

getObjectShapes(O):-
	yield([],_),
	object(O,ShapesList),
	yield(ShapesList,_),!.


%get goal and initial state statically from here, it is more confortable!
getGoalState:-
	yield([],_),
	yield([free(p3),on(c1,p2),on(c2,p1)],_),!.

getInitialState:-
	yield([],_),
	yield([free(p3),on(c1,p1),on(c2,p2)],_),!.

getMap:-
	yield([],_),
	yield(["/maps/map_simple.png"],_),!.
	%yield(["/maps/map_3carts.png"],_),!.



%list of tasks, objects, poses, variables
task(bwPick(X,Y)):-object(X,_),pose(Y,_,_,_).
task(bwPut(X,Y)):-object(X,_),pose(Y,_,_,_).

target(bwPick(X,_),[X]).
target(bwPut(_,Y),[Y]).

object(c1,[
	shape(-0.625,-0.47,0.0,0.15,0.1,0.01),
  	shape( 0.625,-0.47,0.0,0.15,0.1,0.01),
  	shape(-0.625, 0.47,0.0,0.15,0.1,0.01),
  	shape( 0.625, 0.47,0.0,0.15,0.1,0.01) ]).
object(c2,[
	shape(-0.625,-0.47,0.0,0.15,0.1,0.01),
  	shape( 0.625,-0.47,0.0,0.15,0.1,0.01),
  	shape(-0.625, 0.47,0.0,0.15,0.1,0.01),
  	shape( 0.625, 0.47,0.0,0.15,0.1,0.01) ]).

pose(p1,2,0,0).
pose(p2,-2,0,0).
pose(p3,-2,-2,0).

variable(free(X)):-pose(X,_,_,_).
variable(carry(X)):-object(X,_).
variable(on(X,Y)):-object(X,_),pose(Y,_,_,_).
variable(carrying).

%blockworld (il goal è sempre tra gli effetti)
effect(bwMoveOn(X,Y),[on(X,Y),-free(Y)]).

effect(bwPick(X,Y),[-on(X,Y),free(Y),carry(X),carrying]).

effect(bwPut(X,gnd),[-carry(X),free(X),-carrying]).
effect(bwPut(X,Y),[on(X,Y),-carry(X),-free(Y),-carrying]).

effect(bwFree(X),[free(X)]).
%effect(bwPlan(_),[on(a,b),on(b,c)]).

%probably unused
result(S,R):-goal(S,G),effect(S,E),append(G,E,R).
result(S,R):-goal(S,R).
result(S,R):-effect(S,R).

%blockworld
constraint(bwPick(X,Y),[on(X,Y),-carrying]).

constraint(bwPut(X,gnd),[carry(X)]).
constraint(bwPut(X,Y),[carry(X),free(Y)]).

constraint(bwFree(_),[-carrying]).
constraint(bwLeaves,[carriyng]).

%constraint(reachColor(C),[C.present]).
%constraint(place(Obj,_,ground),[knapsack.Obj]).
