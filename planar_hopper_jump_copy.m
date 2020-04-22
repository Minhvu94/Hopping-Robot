clear all
clc
load('u0_horizontal.mat','u1');
%%
global control indicator x_td norminal_traj
indicator = 1;
% indicator: indicate flight (indi=1) and contact (indi=-1) phases
% x_td: x_coordinates of the point foot when it first touch the ground (to compute tangential ground force)

lambda_coeff = 0.01;
eps = 0.00001;
dt = 0.01;
T = 1.3;
iter_max = ceil(T/dt);

x0 = [0; 1; deg2rad(-10); deg2rad(10); 1;
      1; 0;     0;          0;      0];
x_target = [1.5; 2; deg2rad(0); deg2rad(10); 1;
              0; 0;      0;            0;      0];
        
% u1 = [u0;zeros(2*25,1)];
option = odeset('RelTol',1e-3);

rejected = false;
continue_iterating = true;

Q1 = [];
for ii=1:iter_max
Q1 = blkdiag(Q1,[1 0;0 1]);
end

% Q2 = [0,0,1,0,0,0,0,0,0,0;
%       0,0,0,1,0,0,0,0,0,0;
%       0,0,0,0,0,0,0,0,0,0];
base = [1,1,0,0,0,0,0,0,0,0];
Q2 = [];
for ii=1:10
Q2 = blkdiag(Q2,base(:,ii));
end
%  rad2deg(x_traj(3:4,end))

M_bound = -[0,0,0,0,1,0,0,0,0,0];
r0 = 0;
V_bound = -r0;
% for repeat=1:3
while continue_iterating      
    if rejected==false
        x_traj = [];
        x_norm = x0;
        M_constraint = [];
        V_constraint = [];
        for iter=1:iter_max
            % This is the norminal trajectory
            norminal_traj = true;
            u_nominal = u1(2*iter-1:2*iter);
            control = u_nominal;
            [t,y] = ode45(@dynamics, [0 dt], x_norm, option);       
            x_norminal = y(end,:)';
            x_traj = [x_traj x_norminal];
            
            % Calculate partial derivative dPhi/du
            norminal_traj = false;
            dPhi_du = zeros(10,2);
            for i=1:length(u_nominal)
                control = u_nominal;
                control(i) = u_nominal(i) + eps;
                [t,y] = ode45(@dynamics, [0 dt], x_norm, option);
                x_perturbed = y(end,:)';
                dPhi_du(:,i) = (x_perturbed - x_norminal)/eps;
            end
            B_store{iter} = dPhi_du;
            
            % Calculate partial derivative dPhi/dx
            norminal_traj = false;
            dPhi_dx = zeros(10,10);
            for i=1:10
                control = u_nominal;
                x_ = x_norm;
                x_(i) = x_norm(i) + eps;
                [t,y] = ode45(@dynamics, [0 dt], x_, option);
                x_perturbed = y(end,:)';
                dPhi_dx(:,i) = (x_perturbed - x_norminal)/eps;
            end
            A_store{iter} = dPhi_dx;
            
            x_norm = x_norminal;
            
            % Calculate constraint matrices
%             q31L=x_norm(1); q32L=x_norm(2); q41L=x_norm(3); q42L=x_norm(4); q1L=x_norm(5); xH=x_norm(6); yH=x_norm(7);
%             foot1 = yH + L_fem*cos(q31L+q1L) + L_tib*cos(q31L+q41L+q1L);
%             Dfoot1_dx = [-L_fem*sin(q31L+q1L)-L_tib*sin(q31L+q41L+q1L), 0, -L_tib*sin(q31L+q41L+q1L), 0, -L_fem*sin(q31L+q1L)-L_tib*sin(q31L+q41L+q1L), 0, 1,zeros(1,9)]; 
%             foot2 =  yH + L_fem*cos(q32L+q1L) + L_tib*cos(q32L+q42L+q1L);
%             Dfoot2_dx = [0, -L_fem*sin(q32L+q1L)-L_tib*sin(q32L+q42L+q1L), 0, -L_tib*sin(q32L+q42L+q1L), -L_fem*sin(q32L+q1L)-L_tib*sin(q32L+q42L+q1L), 0, 1,zeros(1,9)];
%             
%             M = [-Dfoot1_dx;-Dfoot2_dx; M_bound];
%             V = [foot1 - ground_level; foot2 - ground_level; V_bound - M_bound*x_norm];
%             M_constraint = blkdiag(M_constraint,M);
%             V_constraint = [V_constraint; V];
%             M_constraint = blkdiag(M_constraint,M_bound);
%             V_constraint = [V_constraint; V_bound - M_bound*x_norm];

%             M_constraint = blkdiag(M_constraint,M_bound);
%             V_constraint = [V_constraint; V_bound - M_bound*x_norm];

            % keep this
            M = [M_bound;-[0,1,0,0,0,0,0,0,0,0]*x_norm*[0,1,0,0,0,0,0,0,0,0]];
            V = [V_bound - M_bound*x_norm;([0,1,0,0,0,0,0,0,0,0]*x_norm)^2];
            M_constraint = blkdiag(M_constraint,M);
            V_constraint = [V_constraint; V];

        end
        cost_old = (x_norm-x_target)'*Q2*(x_norm-x_target);
        if cost_old <= 0.01 
            cost_old
            break
        end        
        % Calculate H
        H = B_store{1};
        HH = H;
        for iter = 2:iter_max
            H = [A_store{iter}*H, B_store{iter}];
            HH = [HH zeros((iter-1)*10,2); H];
        end
        % Denote switching locations
        indices = find(diff(sign(x_traj(2,:))));
        
%         % Compute Aeq and Beq
%         Aeq = [];
%         Beq = [];
%         for index=indices
%             Aeq = [Aeq; [0,1,0,0,0,0,0,0,0,0]*HH(10*index-9:10*index,:)];
%             Beq = [Beq; -x_traj(2,index)];
%         end

         figure(1)
         hold off 
         scatter(x_traj(1,:),x_traj(2,:),10,'b','filled')
         hold on 
         scatter(x_traj(1,end),x_traj(2,end),13,'r','filled')
         plot([-3,1,2,3],[0,0,2,0],'LineWidth',2)
         plot(x_target(1),x_target(2),'rx','LineWidth',2)
         grid on 
         xlim([-5 8])
         ylim([-1 5])
        
    end
    options = optimset('display','off', 'TolFun',1e-3);
    du = quadprog(H'*Q2*H + (lambda_coeff*cost_old^2)*eye(2*iter_max), (x_norm-x_target)'*Q2*H, M_constraint*HH, V_constraint,[],[],[],[],[],options);
    u_proposal = u1 + du;
    
    figure(2)
    clf
    subplot(2,1,1)
    hold on
    stairs(u1(1:2:end))
    stairs(u_proposal(1:2:end))
    
    subplot(2,1,2)
    hold on
    stairs(u1(2:2:end))
    stairs(u_proposal(2:2:end))
    drawnow; 
     
    % simulate real system using proposed u
    x_sim = x0;
    x_traj_sim = [];
    norminal_traj = true;
    for ite = 1:iter_max
        control = u_proposal(2*ite-1:2*ite);
        [t,y] = ode45(@dynamics, [0 dt], x_sim, option);
        x_sim = y(end,:)';
        x_traj_sim = [x_traj_sim x_sim];
    end
    
    % Lambda adaptive scheme
    cost_proposal_actual = (x_sim-x_target)'*Q2*(x_sim-x_target);
    if cost_proposal_actual < cost_old      
        lambda_coeff = 0.9*lambda_coeff;
        u1 = u_proposal;
        rejected = false;     % accept u_proposal due to cost benefits
    else                                    
        lambda_coeff = 1.05*lambda_coeff;
        rejected = true;      % reject u_proposal
    end
        
    datta = [indices(1) indices(2) x_traj(2,indices(1)) x_traj(2,indices(2)) cost_old cost_proposal_actual lambda_coeff norm(du)];
    formatSpec = 'i1 = %.3f; i2= %.3f; S1 = %.8f; S2 = %.4f; cost_old = %.3f; cost_proposal = %.3f; lambda_coeff = %.4f; norm_u  = %.4f;\n';
    fprintf(formatSpec, datta)  
end
u_step1 = u1;
stopp = ajnvajnv
%    u1=u_step1;
%%      
figure(1)
simulation(x_traj)
%%
figure(2)
hold on 
uuu = u1(2:2:end);
stairs(uuu)
%%
% clear all
clc
% load('u_jump.mat','u1')
global control indicator x_td norminal_traj
norminal_traj = true;
indicator = 1;

llll = 1

dt = 0.01;
T = 1.3;
iter_max = ceil(T/dt);

x0 = [0; 1; deg2rad(-10); deg2rad(10); 1;
      1; 0;     0;          0;      0];
         
u1 = 0*ones(2*iter_max,1);              % 45,100
u1(46*2:2:96*2) = 250*ones((96-46)+1,1);       % dr>0 at 66
u1(82*2:2:95*2) = 300*ones((95-82)+1,1);
u1(98*2:2:109*2) = -30*ones((109-98)+1,1);
u1(2:2:end) = u1(2:2:end)/20;

u1(46*2-1:2:80*2-1) = -15*ones(35,1); 

option = odeset('RelTol',1e-5);
% lllll = 1
save('u0_horizontal.mat','u1');

% u1 = [u0;0*ones(2*iter_max-length(u0),1)];   % 45,89    [46->88]=300;     
% u1(138:2:170) = [300:1:316]';
x_traj_1 = [];
x_norm = x0;

for iter=1:iter_max
    % This is the norminal trajectory
    u_nominal = u1(2*iter-1:2*iter);
    control = u_nominal;
    [t,y] = ode45(@dynamics, [0 dt], x_norm);
    x_norminal = y(end,:)';
    x_traj_1 = [x_traj_1 x_norminal];
    x_norm = x_norminal;
end
% mmm =avkjna
switching = find(diff(sign(x_traj_1(2,:))))
% yy = x_traj_1(2,switching(2))
thrust_phase = find(diff(sign(x_traj_1(10,:))))


figure(1)
simulation(x_traj_1)

%     figure
%     subplot(2,1,1)
%     hold on
%     stairs(u1(1:2:end))
%     stairs(u_step1(1:2:end))
%     grid on
%     
%     subplot(2,1,2)
%     hold on
%     stairs(u1(2:2:end))
%     stairs(u_step1(2:2:end)/20)
%     grid on
   
    
% figure(1)
% grid on 
% hold off
% plot(x_traj_1(2,:),'k','LineWidth',1.5);
% hold on
% plot(x_traj_1(5,:),'b','LineWidth',1.5);
% % plot(x_traj_1(7,:),'r','LineWidth',1.5);
% plot(x_traj_1(10,:),'g','LineWidth',1.5);
% legend('y','r','dr')


% figure(3)
% simulation(x_traj_1(:,[thrust_phase:thrust_phase+1]))

% grid on 
% figure(3)
% plot(store_Fy)
% plot(u1(1:2:end))
% plot(u1(2:2:end))
% grid on 
% x_traj_0 =x_traj_1;

function simulation(x_traj)
g = 9.81;
m = 10;
J = 10;
ml = 1;
Jl = 1;
l1 = 0.5;
l2 = 0.4;

[row,column] = size(x_traj);
for i = 1:column
    x=x_traj(1,i); y=x_traj(2,i); theta=x_traj(3,i); phi=x_traj(4,i); r=x_traj(5,i);
    hold off
    p_foot = [x;y];
    p_hip = p_foot + r*[sin(theta); cos(theta)];
    p_body = p_hip + l1*[sin(phi); cos(phi)];

    plot([p_foot(1),p_hip(1)],[p_foot(2),p_hip(2)],'k','LineWidth',2)
    hold on
    grid on
    plot([p_hip(1),p_body(1)],[p_hip(2),p_body(2)],'r','LineWidth',2)
    plot(p_foot(1),p_foot(2),'ko','MarkerSize',5,'LineWidth',2)
    plot(p_hip(1),p_hip(2),'bo','MarkerSize',5,'LineWidth',2)
    plot(p_body(1),p_body(2),'ro','MarkerSize',5,'LineWidth',2)    
        % draw obstacles
    plot([-3,3],[0,0],'LineWidth',2)
    
    xlim([-5 8])
    ylim([-1 5])
    axis equal 
    grid on
    pause(0.01)
end
end

function dxdt = dynamics(t,X)
global control indicator x_td norminal_traj
scale = 20;
g = 9.81;
m = 10;
J = 10;
ml = 1;
Jl = 1;
l1 = 0.5;
l2 = 0.4;

kg = 1e4;
bg = 50;

x=X(1) ;  y=X(2);  theta=X(3);  phi=X(4);  r=X(5);
dx=X(6); dy=X(7); dtheta=X(8); dphi=X(9); dr=X(10);
R = r-l1;

% if r >= 1 && control(2) > 0
%     control(2) = 0;
% end

if y<=0
    if indicator>0 && norminal_traj == true
%         disp('yessss')
        x_td = x;
    end
    Fx = -kg*(x-x_td) - bg*dx; 
    Fy = -kg*y - bg*dy;

%     torque = 50*(theta/2-phi)-20*dphi;
%     if dr>0
%         control = [torque;control(2)];
%     else
%         control = [torque;control(2)];
%     end



    indicator = -1;
else
    Fx = 0; Fy = 0;
    indicator = 1;
%     control(1) = 0;
end

% if r >= 1 && control(2) > 0
%     control(2) = 0;
% end



W = [-ml*R       0     (Jl-ml*R*l1)*cos(theta)         0               0;
        0      ml*R    (Jl-ml*R*l1)*sin(theta)         0               0;
       m*R       0     (Jl+m*R*r)*cos(theta)     m*R*l2*cos(phi)   m*R*sin(theta);
        0      -m*R    (Jl+m*R*r)*sin(theta)     m*R*l2*sin(phi)  -m*R*cos(theta);
        0        0      Jl*l2*cos(theta-phi)          -J*R              0];
    
V = [cos(theta)*(l1*Fy*sin(theta)-l1*Fx*cos(theta)-control(1)) - R*(Fx-scale*control(2)*sin(theta)-ml*l1*dtheta^2*sin(theta));
     sin(theta)*(l1*Fy*sin(theta)-l1*Fx*cos(theta)-control(1)) + R*(ml*l1*dtheta^2*cos(theta)+Fy-scale*control(2)*cos(theta)-ml*g);
     cos(theta)*(l1*Fy*sin(theta)-l1*Fx*cos(theta)-control(1)) + R*scale*control(2)*sin(theta) + m*R*(r*dtheta^2*sin(theta) + l2*dphi^2*sin(phi) - 2*dr*dtheta*cos(theta));
     sin(theta)*(l1*Fy*sin(theta)-l1*Fx*cos(theta)-control(1)) - R*(scale*control(2)*cos(theta)-m*g) - m*R*(2*dr*dtheta*sin(theta) + r*dtheta^2*cos(theta) + l2*dphi^2*cos(phi));
     l2*cos(theta-phi)*(l1*Fy*sin(theta)-l1*Fx*cos(theta)-control(1)) - R*(l2*scale*control(2)*sin(phi-theta)+control(1))];    

 dxdt = [X(6:10);
        inv(W)*V];

% datta = [Fx Fy X(2) dxdt(2) dxdt(7) X(5) dxdt(5) dxdt(10)];
% formatSpec = 'Fx = %.2f; Fy = %.2f; y = %.4f; dy = %.4f; ddy = %.4f; r = %.1f; dr = %.1f; ddr = %.1f; \n';
% if y<=0
% datta = [y dy Fy -kg*y -bg*dy indicator];
% formatSpec = 'y = %.3f; dy = %.3f; Fy = %.2f; kg = %.2f; bg = %.2f; indicator = %.0f\n';
% fprintf(formatSpec, datta) 
% end
% datta = [Fx Fy control(1) control(2) x dx dxdt(6) r dr dxdt(10)];
% formatSpec = 'Fx = %.4f; Fy = %.4f; u1=%.4f; u2=%.4f; x=%.4f; dx=%.4f; ddx=%.4f; r=%.4f; dr=%.4f; ddr=%.4f;\n';
% fprintf(formatSpec, datta)

end












































