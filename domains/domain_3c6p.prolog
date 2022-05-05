%
% PLANNING DOMAIN WITH 3 CARTS AND 6 POSES
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Occupancy grid in png for motion planning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map(["/maps/map_simple.png"]).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symbolic Goal State
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goal_state([free(p1),free(p2),free(p3),on(c1,p4),on(c2,p5),on(c3,p6)]).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symbolic Initial State 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initial_state([free(p4),free(p5),free(p6),on(c1,p1),on(c2,p2),on(c3,p3)]).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%State Variables with types for arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variable(free(X)):- pose(X,_,_,_).
variable(carry(X)):-object(X,_).
variable(on(X,Y)):-object(X,_),pose(Y,_,_,_).
variable(carrying).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Operators - argument types
operator(bwPick(X,Y)):-object(X,_),pose(Y,_,_,_).
operator(bwPut(X,Y)):-object(X,_),pose(Y,_,_,_).

% Operators - Effects
effect(bwPick(X,Y),[-on(X,Y),free(Y),carry(X),carrying]).
effect(bwPut(X,Y),[on(X,Y),-carry(X),-free(Y),-carrying]).


%Operator - Preconditions
precondition(bwPick(X,Y),[on(X,Y),-carrying]).
precondition(bwPut(X,Y),[carry(X),free(Y)]).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G function
% (G function gets geometric information from poses defined below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
target(bwPick(_,Y),[Y]).
target(bwPut(_,Y),[Y]).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain terms with motion configuration (objects and poses)
% Objects: objects are associated with bounding boxes
% Poses: symbolic poses are associated with coordinates in the configuration space   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Carts 
object(c1,X):-geometric_bounding(c1,X).
object(c2,X):-geometric_bounding(c2,X).
object(c3,X):-geometric_bounding(c3,X).

% Poses
pose(p1,2,0,0).
pose(p2,3,-2,0).
pose(p3,3,2,0).
pose(p4,-2,0,0).
pose(p5,-2,-2,0).
pose(p6,-2,2,0).

% Geometric bounding
geometric_bounding(_,[
	shape(-0.625,-0.47,0.0,0.15,0.1,0.01),
  	shape( 0.625,-0.47,0.0,0.15,0.1,0.01),
  	shape(-0.625, 0.47,0.0,0.15,0.1,0.01),
  	shape( 0.625, 0.47,0.0,0.15,0.1,0.01) ]).


