function Y = admittance_matrix(bus_n,branch_n,Admittance,lined) 

Bus_Number = bus_n;
Branch_Number = branch_n;
Y = zeros(Bus_Number,Bus_Number);
line_val = lined;

% Formation of the Off Diagonal Elements...
 for k = 1:Branch_Number
     Y(line_val(k,1),line_val(k,2)) = Y(line_val(k,1),line_val(k,2)) - Admittance(k);%互导纳的求解
     Y(line_val(k,2),line_val(k,1)) = Y(line_val(k,1),line_val(k,2));%互导纳对称性
 end
 
 % Formation of Diagonal Elements....
 for m = 1:Bus_Number%自导纳的求解（所有和节点有链接的支路的导纳之和）
     for n = 1:Branch_Number
         if line_val(n,1) == m
             Y(m,m) = Y(m,m) + Admittance(n);
         elseif line_val(n,2) == m
             Y(m,m) = Y(m,m) + Admittance(n);
         end
     end
  end