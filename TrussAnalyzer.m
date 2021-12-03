function TrussAnalyzer(inputFile, units, Lengths)
load(inputFile, 'C', 'Sx', 'Sy', 'X', 'Y', 'L');
if(nargin < 2)
   units = "oz"; 
end
if(nargin < 3)
   Lengths = 0; 
end
sizex = size(X);
numj = max(sizex);
distances = zeros(numj,numj);
for i = 1:numj
    for j = 1:numj
        distances(i,j) = sqrt(power(X(i) - X(j),2) + power(Y(i) - Y(j),2));
    end
end
[joints, members] = size(C);
xcomp = zeros(joints,members);
ycomp = zeros(joints,members);
for j = 1:members
    for i = 1:joints
        if(C(i,j) == 1)
            l = i;
            break;
        end
    end
    for i = l:joints
        if(C(i,j) == 1)
            m = i;
        end
    end
    if((X(m) - X(l)) ~= 0)
        xcomp(l,j) = (X(m) - X(l))/distances(l,m);
        xcomp(m,j) = (X(l) - X(m))/distances(l,m);
    end
    if((Y(m) - Y(l)) ~= 0)
        ycomp(l,j) = (Y(m) - Y(l))/distances(l,m);
        ycomp(m,j) = (Y(l) - Y(m))/distances(l,m);
    end
end
A = [[xcomp, Sx];[ycomp, Sy]];
T = A\transpose(L);
Td = char(abs(sign(A)));
for i = 1:members
    if(T(i) < 0)
        Td(i) = 'C';
    else
        Td(i) = 'T';
    end
end
fprintf("EK301 Section A3: Ananth Sanjay, Eric Xie, William Nilsen - 11/5/2021\n");
fprintf("Load: %d %s\n", sum(L), units);
fprintf("Member forces in %s:\n", units);
for i = 1:members
    fprintf("m%d: %.3f (%c)\n",i,abs(T(i)),Td(i));
end
fprintf("Reaction forces in %s:\n", units);
fprintf("Sx1: %.3f\n", T(members + 1));
fprintf("Sy1: %.3f\n", T(members + 2));
fprintf("Sy2: %.3f\n", T(members + 3));
if(max(Lengths) > 1)
   cost = sum(Lengths) + 10*joints;
   fprintf("Cost of truss: $%.2f\n", cost);
   R = T(1:members)/sum(L); 
   P = zeros(1,members);
   for i = 1:members
       P(i) = 2570/(Lengths(i)*Lengths(i));
       if(T(i) >= 0)
          R(i) = 0.00000000001; 
       end
   end
   W = P./transpose(R);
   fprintf("Critical member %d buckles at load %.3f, R = %.3f\n", find(W == min(W),1), -1*min(W), R(find(W == min(W),1)));
end
end