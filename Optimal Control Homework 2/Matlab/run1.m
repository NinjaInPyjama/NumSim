clear
clc

%% General definitions

F = [2 2 3; 7 5 4; 4 6 5; 7 8 6; 4 4 6; 1 7 8; 8 8 8; 8 8 8]';
F0 = [1 1 2; 6 3 4; 1.5 1.5 4; 0 1 1; 3 3 2; 1 1 4; 1.5 1.5 1.5; 0 0 0]';

alpha = 0.9;

%% Aufgabe 1a) (Value function Iteration)

[V_funIt,U_funIt,~] = valFunIteration(F,F0,alpha);

disp('For the value function iterations:')
disp(' ')
disp('V* = ')
disp(V_funIt)
disp('Optimal feedback = ')
disp(U_funIt)

%% Aufgabe 1f) (Linear Program)

[V_linProg, U_linProg] = valFunItAsLinProg(F,F0,alpha);

disp('For the linear programm:')
disp(' ')
disp('V* = ')
disp(V_funIt)
disp('Optimal feedback = ')
disp(U_funIt)