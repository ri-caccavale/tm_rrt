%
% PLANNING DOMAIN WITH 2 CARTS AND 3 POSES
%


% problem formulation
map(["/maps/map_simple.png"]).

goal_state([free(p3),on(c1,p2),on(c2,p1)]).

initial_state([free(p3),on(c1,p1),on(c2,p2)]).


% lists of tasks, objects, poses, variables

task(bwPick(X,Y)):-object(X,_),pose(Y,_,_,_).
task(bwPut(X,Y)):-object(X,_),pose(Y,_,_,_).

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

% lists of targets, effects, constraints

target(bwPick(X,_),[X]).
target(bwPut(_,Y),[Y]).

effect(bwMoveOn(X,Y),[on(X,Y),-free(Y)]).

effect(bwPick(X,Y),[-on(X,Y),free(Y),carry(X),carrying]).

effect(bwPut(X,gnd),[-carry(X),free(X),-carrying]).
effect(bwPut(X,Y),[on(X,Y),-carry(X),-free(Y),-carrying]).

effect(bwFree(X),[free(X)]).

constraint(bwPick(X,Y),[on(X,Y),-carrying]).

constraint(bwPut(X,gnd),[carry(X)]).
constraint(bwPut(X,Y),[carry(X),free(Y)]).

constraint(bwFree(_),[-carrying]).
constraint(bwLeaves,[carriyng]).
