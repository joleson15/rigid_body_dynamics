function [A, N, n_pts, CoM, areas, masses] = bouy(n, m_total, side_length)
    
    n_pts = 6*n^2;
    
    h = side_length/n;

    A = zeros(n_pts, 3); %locations of points
    N = zeros(n_pts, 3); %normal vectors corresponding to each point
    ct = 0; %counter


    %For simplicity, I will hard code the individual links, but for the
    %future, the MODULO operator is the key to expanding to any 'n'

   %face 1:

   z1 = -side_length/2;
   for y1 = (-side_length/2 + h/2):h:(side_length/2 - h/2)
       for x1 = (-side_length/2 + h/2):h:(side_length/2 - h/2)
           ct = ct + 1;
           A(ct,:) = [x1, y1, z1];
           N(ct, :) = [0,0,-1];
       end
   end

   %face 2:

   z2 = side_length/2;
   for y2 = (-side_length/2 + h/2):h:(side_length/2 - h/2)
       for x2 = (-side_length/2 + h/2):h:(side_length/2 - h/2)
           ct = ct + 1;
           A(ct,:) = [x2, y2, z2];
           N(ct, :) = [0,0,1];
       end
   end
   
   %face 3: 

   y3 = -side_length/2;
   for z3 = (-side_length/2 + h/2):h:(side_length/2 - h/2)
       for x3 = (-side_length/2 + h/2):h:(side_length/2 - h/2)
           ct = ct + 1;
           A(ct,:) = [x3, y3, z3];
           N(ct, :) = [0,-1,0];
       end
   end
   
   %face 4: 

   y4 = side_length/2;
   for z4 = (-side_length/2 + h/2):h:(side_length/2 - h/2)
       for x4 = (-side_length/2 + h/2):h:(side_length/2 - h/2)
           ct = ct + 1;
           A(ct,:) = [x4, y4, z4];
           N(ct, :) = [0,1,0];
       end
   end

   %face 5: 

   x5 = -side_length/2;
   for z5 = (-side_length/2 + h/2):h:(side_length/2 - h/2)
       for y5 = (-side_length/2 + h/2):h:(side_length/2 - h/2)
           ct = ct + 1;
           A(ct,:) = [x5, y5, z5];
           N(ct, :) = [-1,0,0];
       end
   end
   %face 6: 
   
   x6 = side_length/2;
   for z6 = (-side_length/2 + h/2):h:(side_length/2 - h/2)
       for y6 = (-side_length/2 + h/2):h:(side_length/2 - h/2)
           ct = ct + 1;
           A(ct,:) = [x6, y6, z6];
           N(ct, :) = [1,0,0];
       end
   end

   CoM = [mean(A(:,1)), mean(A(:,2)), mean(A(:,3))];

   areas = h^2*ones(n_pts, 1);
   masses = (m_total/n_pts) * ones(n_pts, 1);

   