%% Plotting the policy of the approach by Ono
quiverDataOno=zeros(num_vars,4);

for state = 1:num_vars
    u = pi_star_ono(state,1);
    psi = u/Mu*2*pi; % Under disturbances on the input
    dx = round(3*Mx/50*cos(psi));
    dy = round(3*My/50*sin(psi));

    coords = transformVectorIdxToMatrixIdx([Mx,My],state);
    if mod(coords(1),quiver_spacing)==0 && mod(coords(2),quiver_spacing)==0 && coords(1) || mod(coords(2),quiver_spacing)==0 && coords(1)==1 || mod(coords(1),quiver_spacing)==0 &&coords(2)==1 || coords(1)==1&&coords(2)==1
        quiverDataOno(state,:)=[coords(2),coords(1),dy,dx];
    end
end
quiverDataOno = quiverDataOno(any(quiverDataOno,2),:);
figure(11)
colormap(quivercolors);
colormap(mycolors);
imagesc(mask)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
hold on
quiver(quiverDataOno(1:end,1),quiverDataOno(1:end,2),quiverDataOno(1:end,3)*quiver_arrowscaling,quiverDataOno(1:end,4)*quiver_arrowscaling,0.7,'black','LineWidth',3)
 p1 = [quiverDataOno(1:end,1) quiverDataOno(1:end,2)]; % data start point
 u = quiverDataOno(1:end,3); v=quiverDataOno(1:end,4);
 arrow3(p1,p1+LineLength*[u,v],linestyle,headWidth,headLength);
coords = transformVectorIdxToMatrixIdx([Mx,My],x_0_idx);
%scatter(coords(1),coords(2),150,'red','filled','o')
scatter(coords(1),coords(2),90,'white','filled','o','LineWidth',3,'MarkerEdgeColor','black','SizeData',400)
title('Policy k=0 by Ono et al.')

%% Plotting the policy of our approach
quiverDataMine=zeros(num_vars,4);   % We will store the cheapest policy in this variable for some plots
for state = 1:numel(indices_of_safe_states)
    u = mu_k_star_our_up(state,1);
    psi = u/Mu*2*pi; % Under disturbances on the input
    dx = round(3*Mx/50*cos(psi));                   % positive x is movement to the right, which changes the column coordinate
    dy = round(3*My/50*sin(psi));                   % positive y is movement down, which changes the row coordinate
    x = transformVectorIdxToMatrixIdx([Mx,My],indices_of_safe_states(state));   % First state is rows, second state is columns (coordinate system starts in top left, then goes right)

    %quiverDataMine(indices_of_safe_states(state),:)=[x(1),x(2),dx, dy];
    if mod(x(1),quiver_spacing)==0 && mod(x(2),quiver_spacing)==0 && x(1) || mod(x(2),quiver_spacing)==0 && x(1)==1 || mod(x(1),quiver_spacing)==0 &&x(2)==1 || x(1)==1&&x(2)==1
        quiverDataMine(state,:)=[x(2),x(1),dy,dx];
    end
    %quiverDataMine(indices_of_safe_states(state),:)=[x(2),x(1),dy,-dx];          % First state are columns, second state are rows (coordinate system starts in bottom left, then goes right)
end        
quiverDataMine = quiverDataMine(any(quiverDataMine,2),:);                                                                    % third state how to move along horizontaly, fourth state how to move vertically (positive = up)
figure(17)
% Plot mask
imagesc(mask)
hold on
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap(quivercolors);
colormap(mycolors);

quiver(quiverDataMine(1:end,1),quiverDataMine(1:end,2),quiverDataMine(1:end,3)*quiver_arrowscaling,quiverDataMine(1:end,4)*quiver_arrowscaling,0.7,'black','LineWidth',2)
 p1 = [quiverDataMine(1:end,1) quiverDataMine(1:end,2)]; % data start point
 u = quiverDataMine(1:end,3); v=quiverDataMine(1:end,4);
 arrow3(p1,p1+LineLength*[u,v],linestyle,headWidth,headLength);
coords = transformVectorIdxToMatrixIdx([Mx,My],x_0_idx);
%catter(coords(1),coords(2),150,'red','filled','o')
scatter(coords(1),coords(2),90,'white','filled','o','LineWidth',3,'MarkerEdgeColor','black','SizeData',400)
title('Policy k=0 by our approach, b_k=1')
%% Plotting the cheapest policy 

quiverDataCheapest=zeros(num_vars,4);   % We will store the cheapest policy in this variable for some plots

for state = 1:num_vars
    u = mu_k_star(state,1);
    psi = u/Mu*2*pi; % Under disturbances on the input
    dx = round(3*Mx/50*cos(psi));
    dy = round(3*My/50*sin(psi));

    coords = transformVectorIdxToMatrixIdx([Mx,My],state);
    if mod(coords(1),quiver_spacing)==0 && mod(coords(2),quiver_spacing)==0 && coords(1) || mod(coords(2),quiver_spacing)==0 && coords(1)==1 || mod(coords(1),quiver_spacing)==0 &&coords(2)==1 || coords(1)==1&&coords(2)==1
        quiverDataCheapest(state,:)=[coords(2),coords(1),dy,dx];
    end
end
quiverDataCheapest = quiverDataCheapest(any(quiverDataCheapest,2),:);
figure(13)
% Plot mask
imagesc(mask)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap(quivercolors);
colormap(mycolors);
hold on
quiver(quiverDataCheapest(1:end,1),quiverDataCheapest(1:end,2),quiverDataCheapest(1:end,3)*quiver_arrowscaling,quiverDataCheapest(1:end,4)*quiver_arrowscaling,0.7,'black','LineWidth',3)
 p1 = [quiverDataCheapest(1:end,1) quiverDataCheapest(1:end,2)]; % data start point
 u = quiverDataCheapest(1:end,3); v=quiverDataCheapest(1:end,4);
 arrow3(p1,p1+LineLength*[u,v],linestyle,headWidth,headLength);
coords = transformVectorIdxToMatrixIdx([Mx,My],x_0_idx);
%scatter(coords(1),coords(2),150,'red','filled','o')
scatter(coords(1),coords(2),90,'white','filled','o','LineWidth',3,'MarkerEdgeColor','black','SizeData',400)
title('Policy k=0 by our approach, b_k=0')