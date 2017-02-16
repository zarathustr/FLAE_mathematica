%Fast Linear Attitude Estimator from Vector Observations
%Authors: Jin Wu, Zebo Zhou et al.
%Copytight (c) 2016

function [Q,W,iter,time_poly,time_eig] = FLAE( r_base, b_base, weights )
    iter=0;
    
    tic;
    MM=zeros(3,3);
   
    for i=1:length(b_base(1,:))
        
        MM=MM+weights(i)*b_base(:,i)*r_base(:,i)';
    end
    
    Hx1=MM(1,1);    Hx2=MM(1,2);    Hx3=MM(1,3);
    Hy1=MM(2,1);    Hy2=MM(2,2);    Hy3=MM(2,3);
    Hz1=MM(3,1);    Hz2=MM(3,2);    Hz3=MM(3,3);
    
    W=[Hx1+Hy2+Hz3,-Hy3+Hz2,-Hz1+Hx3,-Hx2+Hy1;
       -Hy3+Hz2, Hx1-Hy2-Hz3,Hx2+Hy1,Hx3+Hz1;
       -Hz1+Hx3,Hx2+Hy1,Hy2-Hx1-Hz3,Hy3+Hz2;
       -Hx2+Hy1,Hx3+Hz1,Hy3+Hz2,Hz3-Hy2-Hx1];
   

   c=det(W);
   b=8*Hx3*Hy2*Hz1 - 8*Hx2*Hy3*Hz1 - 8*Hx3*Hy1*Hz2 + 8*Hx1*Hy3*Hz2 + 8*Hx2*Hy1*Hz3 - 8*Hx1*Hy2*Hz3;
   a=-2*Hx1*Hx1 - 2*Hx2*Hx2 - 2*Hx3*Hx3 - 2*Hy1*Hy1 - 2*Hy2*Hy2 - 2*Hy3*Hy3 - 2*Hz1*Hz1 - 2*Hz2*Hz2 - 2*Hz3*Hz3;

   time_poly=toc;
   
   tic;
   T0 = 2*a^3 + 27*b^2 - 72*a*c;
   T1 = (T0 + sqrt(-4*(a^2 + 12*c)^3 + T0^2))^(1/3);
   T2 = sqrt(-4*a + 2^(4/3)*(a^2 + 12*c)/T1 + 2^(2/3)*T1);
   
   lambda1 =   0.20412414523193150818310700622549*(T2 - sqrt(-T2^2 - 12*a - 12 *2.4494897427831780981972840747059*b/T2));
   lambda2 =   0.20412414523193150818310700622549*(T2 + sqrt(-T2^2 - 12*a - 12 *2.4494897427831780981972840747059*b/T2));
   lambda3 =  -0.20412414523193150818310700622549*(T2 + sqrt(-T2^2 - 12*a + 12 *2.4494897427831780981972840747059*b/T2));
   lambda4 =  -0.20412414523193150818310700622549*(T2 - sqrt(-T2^2 - 12*a + 12 *2.4494897427831780981972840747059*b/T2));
   
   lambda=max([real(lambda1) real(lambda2) real(lambda3) real(lambda4)]);
   
   G=W-lambda*eye(4);
   
   
   
   pivot = G(1,1);  
   G(1,:) = G(1,:)/pivot;
   G(2,:) = G(2,:) - G(2, 1)*G(1,:);
   G(3,:) = G(3,:) - G(3, 1)*G(1,:);
   G(4,:) = G(4,:) - G(4, 1)*G(1,:);

   
   pivot = G(2,2);
   G(2,:) = G(2,:)/pivot;
   G(1,:) = G(1,:) - G(1, 2)*G(2,:);
   G(3,:) = G(3,:) - G(3, 2)*G(2,:);
   G(4,:) = G(4,:) - G(4, 2)*G(2,:);
   
   pivot = G(3,3);
   G(3,:) = G(3,:)/pivot;
   G(1,:) = G(1,:) - G(1, 3)*G(3,:);
   G(2,:) = G(2,:) - G(2, 3)*G(3,:);
   G(4,:) = G(4,:) - G(4, 3)*G(3,:);
   
   q=[G(1,4);G(2,4);G(3,4);-1];
 
   Q=q./norm(q);
   
   time_eig=toc;
   
end

