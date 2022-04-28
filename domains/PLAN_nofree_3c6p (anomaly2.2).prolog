%
% PLANNING DOMAIN WITH 3 CARTS AND 6 POSES
%

% problem formulation
map(["/maps/map_anomaly2_X.png"]).

goal_state([free(p5),free(p6),free(p2),on(c1,p1),on(c2,p8),on(c3,p3)]).

initial_state([free(p5),free(p6),free(p8),on(c1,p1),on(c2,p2),on(c3,p3)]).

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
 	shape( 0.625, 0.47,0.0,0.15,0.1,0.01) ]). %% c3p3 used as DISTRACTIONS
%object(c4,[
%	shape(-0.625,-0.47,0.0,0.15,0.1,0.01),
%  	shape( 0.625,-0.47,0.0,0.15,0.1,0.01),
% 	shape(-0.625, 0.47,0.0,0.15,0.1,0.01),
% 	shape( 0.625, 0.47,0.0,0.15,0.1,0.01) ]).

% x, y, yaw
pose(p1,3.75,3.5,0).
pose(p2,3.75,1,0).
pose(p3,3.75,-1.5,0).
%pose(p4,3.75,-4,0).
pose(p5,-3,-4,0).
pose(p6,-3,-1.5,0).
%pose(p7,-3,1,0).
pose(p8,-3,3.5,0).

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

