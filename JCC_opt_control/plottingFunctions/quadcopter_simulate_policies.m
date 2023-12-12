figure(14)
fontsize=25;
imagesc(mask)
set(gca,'FontName','Times Roman','Fontsize',fontsize)
colormap(mycolors);
hold on
%quiver(quiverDataMine(:,1),quiverDataMine(:,2),quiverDataMine(:,3),quiverDataMine(:,4),'r')

% Extending our policy to full state space to get indices right (so far it was only defined over safe states).
mu_k_star_our_ext_up = zeros(size(mu_k_star));
mu_k_star_our_ext_low = zeros(size(mu_k_star));
mu_k_star_our_ext_up(indices_of_safe_states,:) = mu_k_star_our_up;
mu_k_star_our_ext_low(indices_of_safe_states,:) = mu_k_star_our_low;

% Record total_attempts many simulations
total_attempts=150;

for attempt=1:total_attempts
    success = 1;
    x_idx = zeros(N+1,1);
    x_idx(1) = x_0_idx;
    x = zeros(2,N+1);
    x(:,1) = transformVectorIdxToMatrixIdx([Mx,My],x_idx(1));

    randomNumber = rand();
    if randomNumber<=p_up
        mu_k_star_our_ext = mu_k_star_our_ext_up;
    else
        mu_k_star_our_ext = mu_k_star_our_ext_low;
    end

    for k=1:N
        if success==1
            u = mu_k_star_our_ext(x_idx(k),k);
        else
            u = mu_k_star(x_idx(k),k);
        end
        
        psi = u/Mu*2*pi; 
        d = [0;0];
        d = mvnrnd([0;0],Sigma_disturbance);
        dx = round(3*Mx/50*cos(psi)) + round(d(1));
        dy = round(3*My/50*sin(psi)) + round(d(2));

        x(:,k+1) = [x(1,k)+dx; x(2,k)+dy];
        if x(1,k+1)<=0
            x(1,k+1)=1;
        end
        if x(1,k+1)>Mx
            x(1,k+1) = Mx;
        end
        if x(2,k+1)<=0
            x(2,k+1)=1;
        end
        if x(2,k+1)>My
            x(2,k+1)=My;
        end    
        x_idx(k+1) = transformMatrixIdxToVectorIdx([Mx,My], x(:,k+1));
        if abs(mask(x(1,k+1),x(2,k+1)))<0.01
            success = 0;
        end
    end

    if success==0
        plot(x(2,:),x(1,:),'LineWidth',1,'Color','r')
    else
        plot(x(2,:),x(1,:),'LineWidth',1,'Color','blue')
    end
    hold on


end
set(gca,'FontName','Times Roman','Fontsize',fontsize)
scatter(25.5,25.5,90,'white','filled','o','LineWidth',3,'MarkerEdgeColor','black','SizeData',400)
scatter(25.5,25.5,90,'black','x','LineWidth',3,'SizeData',400)
scatter(x(1,1),x(2,1),90,'white','filled','o','LineWidth',3,'MarkerEdgeColor','black','SizeData',400)
xlabel('x-coordinate')
ylabel('y-coordinate')
set(gca, 'XTick', [1:10:50,50], 'XTickLabel', [-25:10:25,25]) % 10 ticks 
set(gca, 'YTick', [1:10:50,50], 'YTickLabel', [-25:10:25,25]) % 10 ticks 
set(gca,'FontName','Times Roman','Fontsize',fontsize)
title('Policy Our Approach')
%%
figure(15)
imagesc(mask)
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
set(gca,'FontName','Times Roman','Fontsize',fontsize)
colormap(mycolors);
hold on
%quiver(quiverDataOno(:,1),quiverDataOno(:,2),quiverDataOno(:,3),quiverDataOno(:,4),'r')

total_attempts=200;

for attempt=1:total_attempts
    success = 1;
    x_idx = zeros(N+1,1);
    x_idx(1) = x_0_idx;
    x = zeros(2,N+1);
    x(:,1) = transformVectorIdxToMatrixIdx([Mx,My],x_idx(1));
    for k=1:N
        u = pi_star_ono(x_idx(k),k);
        
        psi = u/Mu*2*pi; 
        d = [0;0];
        d = mvnrnd([0;0],Sigma_disturbance);
        dx = round(3*Mx/50*cos(psi)) + round(d(1));
        dy = round(3*My/50*sin(psi)) + round(d(2));
    
        x(:,k+1) = [x(1,k)+dx; x(2,k)+dy];
        if x(1,k+1)<=0
            x(1,k+1)=1;
        end
        if x(1,k+1)>Mx
            x(1,k+1) = Mx;
        end
        if x(2,k+1)<=0
            x(2,k+1)=1;
        end
        if x(2,k+1)>My
            x(2,k+1)=My;
        end    
        x_idx(k+1) = transformMatrixIdxToVectorIdx([Mx,My], x(:,k+1));
        if abs(mask(x(1,k+1),x(2,k+1)))<0.01
            success = 0;
        end
    end
% 
    if success==0
        plot(x(2,:),x(1,:),'LineWidth',1,'Color','r')
    else
        plot(x(2,:),x(1,:),'LineWidth',1,'Color','blue')
    end
end
scatter(25.5,25.5,90,'white','filled','o','LineWidth',3,'MarkerEdgeColor','black','SizeData',400)
scatter(25.5,25.5,90,'black','x','LineWidth',3,'SizeData',400)
scatter(x(1,1),x(2,1),90,'white','filled','o','LineWidth',3,'MarkerEdgeColor','black','SizeData',400)
xlabel('x-coordinate')
ylabel('y-coordinate')
set(gca, 'XTick', [1:10:50,50], 'XTickLabel', [-25:10:25,25]) % 10 ticks 
set(gca, 'YTick', [1:10:50,50], 'YTickLabel', [-25:10:25,25]) % 10 ticks 
set(gca,'FontName','Times Roman','Fontsize',fontsize)
title('Policy by Ono et al.')