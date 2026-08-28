function [ B ] = BSPLINE( A,t,h )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
t1=A(1);
t2=A(2);
t3=A(3);
t4=A(4);
t5=A(5);
if t<=t1
    B=0;
else if t>t1 & t<=t2
        B=1/(6*h^3)*(t-t1)^3;
    else if t>t2 & t<=t3
            B=1/6+1/(2*h)*(t-t2)+1/(2*h^2)*(t-t2)^2-1/(2*h^3)*(t-t2)^3;
        else if t>t3 & t<=t4
                B=1/6+1/(2*h)*(t4-t)+1/(2*h^2)*(t4-t)^2-1/(2*h^3)*(t4-t)^3;
            else if t>t4 & t<=t5
                    B=1/(6*h^3)*(t5-t)^3;
                else if t>t5
                    B=0;
                    end
                end
            end
        end
    end
end
              
end