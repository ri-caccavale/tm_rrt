%
% PLANNING DOMAIN WITH 3 CARTS AND 6 POSES
%

% problem formulation
map(["/maps/map_anomaly2.png"]).

goal_state([free(p2),free(p3),free(p6),on(c1,p1),on(c2,p4),on(c3,p5)]).

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

% x, y, yaw
pose(p1,3,-1,0).
pose(p2,3,-2.5,0).
pose(p3,3,-4,0).
pose(p4,3,2,0).
pose(p5,-3,2,0).
pose(p6,0,-4,0).

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

