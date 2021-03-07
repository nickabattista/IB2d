%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Performs Half step of third order WENO scheme to solve advection eqn
%
%   Author: Matea Santiago
%   Date: February 2021
%   Institution (when created): UC Merced
%
%   Returns: spatial derivatives in x and y of concentration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Cx,Cy]=WENO_3O(C,Uavg,Vavg,dx,dy,dt,Lx,Ly)
C=C';
Cex=[C(end-2,:); C(end-1,:); C(end,:); C; C(1,:); C(2,:); C(3,:)];
Cex=[Cex(:,end-2) Cex(:,end-1) Cex(:,end) Cex Cex(:,1) Cex(:,2) Cex(:,3)];
[Nx,Ny]=size(C);

XEx = (0-3*dx:dx:Lx+2*dx); 
YEx = (0-3*dx:dx:Ly+2*dx); 
Uavg=Uavg';
Vavg=Vavg';



for i=1:Nx
  for j=1:Ny
DD1_0=(Cex(i+1,j+3)-Cex(i,j+3))/dx;
DD1_1=(Cex(i+2,j+3)-Cex(i+1,j+3))/dx;

DD1_4=(Cex(i+5,j+3)-Cex(i+4,j+3))/dx;
DD1_5=(Cex(i+6,j+3)-Cex(i+5,j+3))/dx;

  %i+3
DD1_2=(Cex(i+3,j+3)-Cex(i+2,j+3))/dx;
DD1_3=(Cex(i+4,j+3)-Cex(i+3,j+3))/dx;

  %i+2 i+1 i
DD2_0=(DD1_1-DD1_0)/(2*dx);
  %i+6 i+5 i+4
DD2_4=(DD1_5-DD1_4)/(2*dx);

  %i+3 i+2 i+1 
DD2_1=(DD1_2-DD1_1)/(2*dx);
  %i+4 i+3 i+2
DD2_2=(DD1_3-DD1_2)/(2*dx);
  %i+5 i+4 i+3 
DD2_3=(DD1_4-DD1_3)/(2*dx);

%i+3 i+2 i+1 i 
DD3_1=(DD2_1-DD2_0)/(3*dx);
%i+4 i+3 i+2 i+1  
DD3_2=(DD2_2-DD2_1)/(3*dx);
%i+5 i+4 i+3 i+2 
DD3_3=(DD2_3-DD2_2)/(3*dx);
%i+6 i+5 i+4 i+3
DD3_4=(DD2_4-DD2_3)/(3*dx);

if Uavg(i,j)>0

              %i+3 i+2 i+1 i 
              S1=DD1_2+DD2_1*(XEx(i+3)-XEx(i+3)+XEx(i+3)-XEx(i+2))+DD3_1*((XEx(i+3)-XEx(i+2))*((XEx(i+3)-XEx(i+1))));
              %i+4 i+3 i+2 i+1
              S2=DD1_2+DD2_1*(XEx(i+3)-XEx(i+3)+XEx(i+3)-XEx(i+2))+DD3_2*((XEx(i+3)-XEx(i+2))*((XEx(i+3)-XEx(i+1))));
               %i+5 i+4 i+3 i+2 
              S3=DD1_2+DD2_2*(XEx(i+3)-XEx(i+3)+XEx(i+3)-XEx(i+2))+DD3_3*((XEx(i+3)-XEx(i+2))*((XEx(i+3)-XEx(i+4))));
              
              gamma1=1/10;
              gamma2=6/10;
              gamma3=3/10;
              
           %i+3 i+2 i+1 i 
              B1=13/12*(Cex(i+1,j+3)-2*Cex(i+2,j+3)+Cex(i+3,j+3))^2+1/4*(Cex(i+1,j+3)-4*Cex(i+2,j+3)+3*Cex(i+3,j+3))^2;
              %i+4 i+3 i+2 i+1  
              B2=13/12*(Cex(i+2,j+3)-2*Cex(i+3,j+3)+Cex(i+4,j+3))^2+1/4*(Cex(i+2,j+3)-Cex(i+4,j+3))^2;
              %i+5 i+3 i+2 i+1
              B3=13/12*(Cex(i+3,j+3)-2*Cex(i+4,j+3)+Cex(i+5,j+3))^2+1/4*(Cex(i+3,j+3)-4*Cex(i+4,j+3)+3*Cex(i+5,j+3))^2;
                
              epsilon=10^-6;
              Wt1=gamma1/(epsilon+B1)^2;
              Wt2=gamma2/(epsilon+B2)^2;
              Wt3=gamma3/(epsilon+B3)^2;
        
        W1=Wt1/(Wt1+Wt2+Wt3);
        W2=Wt2/(Wt1+Wt2+Wt3);
        W3=Wt3/(Wt1+Wt2+Wt3);
else
              % i+1 i+2 i+3 i+4
              S1=DD1_3+DD2_2*(XEx(i+3)-XEx(i+4))+DD3_2*((XEx(i+3)-XEx(i+2))*((XEx(i+3)-XEx(i+4))));
              %i+5 i+4 i+3 i+2 
              S2=DD1_3+DD2_2*(XEx(i+3)-XEx(i+4))+DD3_3*((XEx(i+3)-XEx(i+2))*((XEx(i+3)-XEx(i+4))));
              %i+6 i+5 i+4 i+3
              S3=DD1_3+DD2_3*(XEx(i+3)-XEx(i+4))+DD3_4*((XEx(i+3)-XEx(i+5))*((XEx(i+3)-XEx(i+4))));

              gamma1=3/10;
              gamma2=6/10;
              gamma3=1/10;
               
                %i +4 i+3 i+2 i+1
              B1=13/12*(Cex(i+1,j+3)-2*Cex(i+2,j+3)+Cex(i+3,j+3))^2+1/4*(Cex(i+1,j+3)-4*Cex(i+2,j+3)+3*Cex(i+3,j+3))^2;
              %i+5 i+4 i+3 i+2  
              B2=13/12*(Cex(i+2,j+3)-2*Cex(i+3,j+3)+Cex(i+4,j+3))^2+1/4*(Cex(i+2,j+3)-Cex(i+4,j+3))^2;
              %i+6 i+5 i+3 i+2 
              B3=13/12*(Cex(i+3,j+3)-2*Cex(i+4,j+3)+Cex(i+5,j+3))^2+1/4*(Cex(i+3,j+3)-4*Cex(i+4,j+3)+3*Cex(i+5,j+3))^2;
                
              epsilon=10^-6;
              Wt1=gamma1/(epsilon+B1)^2;
              Wt2=gamma2/(epsilon+B2)^2;
              Wt3=gamma3/(epsilon+B3)^2;
        
        W1=Wt1/(Wt1+Wt2+Wt3);
        W2=Wt2/(Wt1+Wt2+Wt3);
        W3=Wt3/(Wt1+Wt2+Wt3);
end
Cx(i,j)=S1*W1+S2*W2+S3*W3;


DD1_0=(Cex(i+3,j+1)-Cex(i+3,j))/dy;
DD1_1=(Cex(i+3,j+2)-Cex(i+3,j+1))/dy;

DD1_4=(Cex(i+3,j+5)-Cex(i+3,j+4))/dy;
DD1_5=(Cex(i+3,j+6)-Cex(i+3,j+5))/dy;

  %i+3
DD1_2=(Cex(i+3,j+3)-Cex(i+3,j+2))/dy;
DD1_3=(Cex(i+3,j+4)-Cex(i+3,j+3))/dy;

  %i+2 i+1 i
DD2_0=(DD1_1-DD1_0)/(2*dy);
  %i+6 i+5 i+4
DD2_4=(DD1_5-DD1_4)/(2*dy);

  %i+3 i+2 i+1 
DD2_1=(DD1_2-DD1_1)/(2*dy);
  %i+4 i+3 i+2
DD2_2=(DD1_3-DD1_2)/(2*dy);
  %i+5 i+4 i+3 
DD2_3=(DD1_4-DD1_3)/(2*dy);

%i+3 i+2 i+1 i 
DD3_1=(DD2_1-DD2_0)/(3*dy);
%i+4 i+3 i+2 i+1  
DD3_2=(DD2_2-DD2_1)/(3*dy);
%i+5 i+4 i+3 i+2 
DD3_3=(DD2_3-DD2_2)/(3*dy);
%i+6 i+5 i+4 i+3
DD3_4=(DD2_4-DD2_3)/(3*dy);


if Vavg(i,j)>0
      S1=DD1_2+DD2_1*(YEx(j+3)-YEx(j+3)+YEx(j+3)-YEx(j+2))+DD3_1*((YEx(j+3)-YEx(j+2))*((YEx(j+3)-YEx(j+1))));
      S2=DD1_2+DD2_1*(YEx(j+3)-YEx(j+3)+YEx(j+3)-YEx(j+2))+DD3_2*((YEx(j+3)-YEx(j+2))*((YEx(j+3)-YEx(j+1))));
      S3=DD1_2+DD2_2*(YEx(j+3)-YEx(j+3)+YEx(j+3)-YEx(j+2))+DD3_3*((YEx(j+3)-YEx(j+2))*((YEx(j+3)-YEx(j+4))));
 
 
                    gamma1=1/10;
              gamma2=6/10;
              gamma3=3/10;
              
           %i+3 i+2 i+1 i 
              B1=13/12*(Cex(i+3,j+1)-2*Cex(i+3,j+2)+Cex(i+3,j+3))^2+1/4*(Cex(i+3,j+1)-4*Cex(i+3,j+2)+3*Cex(i+3,j+3))^2;
              %i+4 i+3 i+2 i+1  
              B2=13/12*(Cex(i+3,j+2)-2*Cex(i+3,j+3)+Cex(i+3,j+4))^2+1/4*(Cex(i+3,j+2)-Cex(i+3,j+4))^2;
              %i+5 i+3 i+2 i+1
              B3=13/12*(Cex(i+3,j+3)-2*Cex(i+3,j+4)+Cex(i+3,j+5))^2+1/4*(Cex(i+3,j+3)-4*Cex(i+3,j+4)+3*Cex(i+3,j+5))^2;
                
              epsilon=10^-6;
              Wt1=gamma1/(epsilon+B1)^2;
              Wt2=gamma2/(epsilon+B2)^2;
              Wt3=gamma3/(epsilon+B3)^2;
        
        W1=Wt1/(Wt1+Wt2+Wt3);
        W2=Wt2/(Wt1+Wt2+Wt3);
        W3=Wt3/(Wt1+Wt2+Wt3);
        
else
      S1=DD1_3+DD2_2*(YEx(j+3)-YEx(j+4))+DD3_2*((YEx(j+3)-YEx(j+2))*((YEx(j+3)-YEx(j+4))));
      S2=DD1_3+DD2_2*(YEx(j+3)-YEx(j+4))+DD3_3*((YEx(j+3)-YEx(j+2))*((YEx(j+3)-YEx(j+4))));
      S3=DD1_3+DD2_3*(YEx(j+3)-YEx(j+4))+DD3_4*((YEx(j+3)-YEx(j+5))*((YEx(j+3)-YEx(j+4))));
              gamma1=3/10;
              gamma2=6/10;
              gamma3=1/10;
               
                %i +4 i+3 i+2 i+1
              B1=13/12*(Cex(i+3,j+1)-2*Cex(i+3,j+2)+Cex(i+3,j+3))^2+1/4*(Cex(i+3,j+1)-4*Cex(i+3,j+2)+3*Cex(i+3,j+3))^2;
              %i+5 i+4 i+3 i+2  
              B2=13/12*(Cex(i+3,j+2)-2*Cex(i+3,j+3)+Cex(i+3,j+4))^2+1/4*(Cex(i+3,j+2)-Cex(i+3,j+4))^2;
              %i+6 i+5 i+3 i+2 
              B3=13/12*(Cex(i+3,j+3)-2*Cex(i+3,j+4)+Cex(i+3,j+5))^2+1/4*(Cex(i+3,j+3)-4*Cex(i+3,j+4)+3*Cex(i+3,j+5))^2;
                
              epsilon=10^-6;
              Wt1=gamma1/(epsilon+B1)^2;
              Wt2=gamma2/(epsilon+B2)^2;
              Wt3=gamma3/(epsilon+B3)^2;
        
        W1=Wt1/(Wt1+Wt2+Wt3);
        W2=Wt2/(Wt1+Wt2+Wt3);
        W3=Wt3/(Wt1+Wt2+Wt3);

end
Cy(i,j)=S1*W1+S2*W2+S3*W3;
% Adv(i,j)=(Cx*Uavg(i,j)+Cy*Vavg(i,j));


end
end
Cx=Cx';
Cy=Cy';
end
