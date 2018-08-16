%% Set a River Domain %%
function M = setRiver(Nx, Ny, p, A, scenario)
    %Probabilities for simulation
    M.p = p;
    M.ProbMoveInsideRiver = p;
    M.ProbMoveOutsideRiver = 0.99;
    
    %Total of states and size of grid
    M.St = Nx * Ny;
    M.Nx = Nx;
    M.Ny = Ny;
    
    %Index of river position
    M.indexRiver=zeros(M.Ny,M.Nx);
    for (i=1:M.Ny*M.Nx)
        M.indexRiver(i)=i;
    end
    
    %Map of river
    M.RiverMap = zeros(M.Ny,M.Nx);
    M.RiverMap(1,1+1:M.Ny-1)=3; %bridge
    M.RiverMap(1+1:M.Ny-1,1+1:M.Ny-1)=2; %river
    M.RiverMap(1:M.Ny,1)=1;  M.RiverMap(1:M.Ny,M.Ny)=1; %land
    M.RiverMap(M.Ny,1+1:M.Ny-1)=-3; %waterfall
    M.transict = zeros(M.Nx*M.Ny, A);
    
    %scenario true is a acumulated positive and false a negative 
    if scenario == true
       M.recomp=zeros(M.Nx*M.Ny,4);
       M.recomp((M.Nx-1)*M.Ny+2,2)=1;
    else if scenario == false
            M.recomp=ones(M.Nx*M.Ny,4)*-1;
            M.recomp((M.Nx-1)*M.Ny+2,2)=0;
         end
    end
    
    M.drawback = zeros(M.St,1); %return states
    sx = 1;
    
    while sx < M.Nx
       num = (sx*M.Ny)+1;
       M.drawback(num)=1;
       sx = sx +1;
    end
    
    M.r=zeros(M.St,1);
    M.r(M.St-Ny+1)=1;
    coord = [0 0];
    for m=1:M.Nx*M.Ny
        m = m-1;
        coord(2) = mod(m,M.Ny);
        coord(1) = mod( floor(m/M.Ny) , M.Nx);
        coord = coord + 1;
        
   %probabilities to move in and outside river
   if ((coord(1) > 1 && coord(1) < M.Nx) && (coord(2) > 1 && coord(2) < M.Ny)) %is in river?
        M.transict(m+1,1) = M.ProbMoveInsideRiver;
        M.transict(m+1,2) = M.ProbMoveInsideRiver;
        M.transict(m+1,3) = M.ProbMoveInsideRiver;
        M.transict(m+1,4) = M.ProbMoveInsideRiver;
   else if (coord(2)==1 && coord(1)>1)
        M.transict(m+1,1) = 0;
        M.transict(m+1,2) = 0;
        M.transict(m+1,3) = 0;
        M.transict(m+1,4) = 0;
  else
        M.transict(m+1,1) = M.ProbMoveOutsideRiver;
        M.transict(m+1,2) = M.ProbMoveOutsideRiver;
        M.transict(m+1,3) = M.ProbMoveOutsideRiver;
        M.transict(m+1,4) = M.ProbMoveOutsideRiver; 
       end
   end
   end    
   M.V = zeros(M.Nx*M.Ny);  
   M.A = A;  
   
   %origin of grid
   M.Px = 1;
   M.Py = 1;
   
   %trasiction Matrix
   M.T = zeros(M.St,A,M.St);
   %expectated mapping S,A
   M.SA = zeros(M.St,A);
   for z=1:M.St 
   if mod(z,M.Ny)==0 
        n1 = z;
   else
        n1 =  z + 1;         % N
   end
    
   if mod(z-1,M.Ny)==0
        n2 = z;
   else 
        n2 = z - 1;          % S
   end
        n3 = z + M.Ny;       % O
        n4 = z - M.Ny;       % L
        
   if (n1 <= 0 || n1 > M.St)
        n1 = z;
   end
   if (n2 <= 0 || n2 > M.St)
        n2 = z; 
   end
   if (n3 <= 0 || n3 > M.St)
        n3 = z; 
   end
   if (n4 <= 0 || n4 > M.St)
        n4 = z; 
   end
          
   if (n1>0)
        if mod(z,M.Ny)==0
            M.T(z,1,n1)=0;
        else
            if (M.transict(z,1)==M.ProbMoveOutsideRiver)
                      M.T(z,1,n1)= M.transict(z,1);
                      M.T(z,1,z)= 1-M.transict(z,1);
              else if  (M.transict(z,1)==M.ProbMoveInsideRiver)
                     M.T(z,1,n1)= M.transict(z,1);
                     M.T(z,1,n2)= 1-M.transict(z,1);
              
              else if (M.transict(z,1)==0)
               M.T(z,1,1) = 1;
              end
              end
          end
       end
   end
   if (n2>0)
        if n2==z
            M.T(z,2,z)=0;   
        else  if (M.transict(z,2)==M.ProbMoveOutsideRiver)
                      M.T(z,2,n2)=M.transict(z,2);
                      M.T(z,2,z)=1-M.transict(z,2);
              else if  (M.transict(z,2)==M.ProbMoveInsideRiver)
                     M.T(z,2,z) =  M.transict(z,2);
                     M.T(z,2,n2) = 1-M.transict(z,2);
              else if (M.transict(z,2)==0)
               M.T(z,2,1) = 1;
              end
              end
       end  
       end
  end
  if (n3>0)
        if n3==z
        	M.T(z,3,z)=0;  
        elseif (M.transict(z,3)==M.ProbMoveOutsideRiver)
                      M.T(z,3,n3)=M.transict(z,3);
                      M.T(z,3,z)=1-M.transict(z,3);
        else if  (M.transict(z,3)==M.ProbMoveInsideRiver)
                     M.T(z,3,n3)= M.transict(z,3);
                     M.T(z,3,n2)= 1-M.transict(z,3);
                else if (M.transict(z,3)==0)
               M.T(z,3,1) = 1;
                end
       end    
       end
   end
          
   if (n4>0)
        if n4==z
            M.T(z,4,z)=0;  
            elseif (M.transict(z,4)==M.ProbMoveOutsideRiver)
                      M.T(z,4,n4)=M.transict(z,4);
                      M.T(z,4,z)=1-M.transict(z,4);
              else if  (M.transict(z,4)==M.ProbMoveInsideRiver)
                     M.T(z,4,n4)= M.transict(z,4);
                     M.T(z,4,n2)= 1-M.transict(z,4);
             
              else if (M.transict(z,4)==0)
               M.T(z,4,1) = 1;
              end
              end  
       end
    end
    M.SA(z,1) = n1;
    M.SA(z,2) = n2;
    M.SA(z,3) = n3;
    M.SA(z,4) = n4;
   end
end
