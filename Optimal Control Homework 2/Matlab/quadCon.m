function [C, Ceq] = quadCon(T,c,z)
    C = z'*T*z - c;
    Ceq = 0;
end