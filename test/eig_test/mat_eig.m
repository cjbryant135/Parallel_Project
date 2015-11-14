N = 4;

A = zeros(N,N);
B = zeros(N,N);

for i = 1:N 
   for j = 1:N
       A(i,j) = i^(2)+3*j+i^j;
   end
   %A(i,i) = 1
   if i <= N/2 
       B(i,i) = -1;
   else
       B(i,i) = 1;
   end
end


[V,D] = eig(A,B)
