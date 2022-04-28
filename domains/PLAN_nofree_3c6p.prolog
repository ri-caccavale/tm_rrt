%
% PLANNING DOMAIN WITH 3 CARTS AND 6 POSES
%

% problem formulation
map(["/maps/map_simple.png"]).

goal_state([free(p1),free(p2),free(p3),on(c1,p4),on(c2,p5),on(c3,p6)]).

initial_state([free(p4),free(p5),free(p6),on(c1,p1),on(c2,p2),on(c3,p3)]).

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
object(c3,[
	shape(-0.625,-0.47,0.0,0.15,0.1,0.01),
  	shape( 0.625,-0.47,0.0,0.15,0.1,0.01),
  	shape(-0.625, 0.47,0.0,0.15,0.1,0.01),
  	shape( 0.625, 0.47,0.0,0.15,0.1,0.01) ]).

pose(p1,2,0,0).
pose(p2,3,-2,0).
pose(p3,3,2,0).
pose(p4,-2,0,0).
pose(p5,-2,-2,0).
pose(p6,-2,2,0).

variable(free(X)):-pose(X,_,_,_).
variable(carry(X)):-object(X,_).
variable(on(X,Y)):-object(X,_),pose(Y,_,_,_).
variable(carrying).

% lists of targets, effects, constraints

target(bwPick(_,Y),[Y]).
target(bwPut(_,Y),[Y]).

effect(bwPick(X,Y),[-on(X,Y),free(Y),carry(X),carrying]).
effect(bwPut(X,Y),[on(X,Y),-carry(X),-free(Y),-carrying]).

constraint(bwPick(X,Y),[on(X,Y),-carrying]).
constraint(bwPut(X,Y),[carry(X),free(Y)]).

