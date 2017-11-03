close all; clear all; clc;

% for detailed insturctions on bvp4c, check the help document

% initial mesh for bvp
meshinit = 101;

% non-relativistic
% y'' + 2y'/x + y^(3/2) = 0
% y'(0) = 0, y(1) = 0
sol_nonrel = bvp4c(...
    @(x,y)([y(2);-y(1).^(1.5)]),...% equation
    @(ya,yb)([ya(2);yb(1)]),...% boundary value
    bvpinit(linspace(0,1,meshinit),[80;-160]),...% initialize
    bvpset('SingularTerm',[0,0;0,-2])...% singularity
);
% density plot rho/rho0
nicefigure1;
setpaper2([16 9]);
axes('position',[0.1 0.1 0.8 0.8]);
hold on; box on;
plot(sol_nonrel.x,(sol_nonrel.y(1,:)./sol_nonrel.y(1,1)).^(1.5));
xlabel('r / R');
ylabel('\rho / \rho_0');

% ultra-relativistic
% y'' + 2y'/x + y^3 = 0
% y'(0) = 0, y(1) = 0
sol_ultrel = bvp4c(...
    @(x,y)([y(2);-y(1).^(3.0)]),...% equation
    @(ya,yb)([ya(2);yb(1)]),...% boundary value
    bvpinit(linspace(0,1,meshinit),[3;-6]),...% initialize
    bvpset('SingularTerm',[0,0;0,-2])...% singularity
);
% density plot rho/rho0
nicefigure1;
setpaper2([16 9]);
axes('position',[0.1 0.1 0.8 0.8]);
hold on; box on;
plot(sol_ultrel.x,(sol_ultrel.y(1,:)./sol_ultrel.y(1,1)).^3);
xlabel('r / R');
ylabel('\rho / \rho_0');


% relativistic - another path
% y'' + 2y'/x + PARA^2 (y^2-1/yy0^2)^(3/2) = 0
% y(0) = 1, y'(0) = 0, y(1) = 1/yy0
yy0 = [0,50,30,20,15,10,6,4,2.5,1.7,1.4,1.22,1.14,1.1,1.07,1.04];
% yy0(0) is supposed to be +\infty
% in substitution we change our first eqaution to be ultra-relativistic,
% and set yy0(0) as 0 (not used)
sol_rel = cell(1,length(yy0));
MASS = zeros(1,length(yy0));
RADIUS = zeros(1,length(yy0));
PARA = zeros(1,length(yy0));
% solve the ultra-relativistic equation first
% use this solution as the initial guess for the following ones
fprintf('%d\t',1);
sol_rel{1} = bvp4c(...
    @(x,y,p2)([y(2);-p2.*y(1).^3]),...
    @(ya,yb,p2)([ya(2);yb(1);ya(1)-1]),...
    bvpinit(linspace(0,1,meshinit),[0.5;-1],(6.9).^2),...
    bvpset('SingularTerm',[0,0;0,-2])...
);
PARA(1) = real(sqrt(sol_rel{1}.parameters));
RADIUS(1) = 0;
MASS(1) = -PARA(1).*real(sol_rel{1}.y(2,end));
fprintf('%0.2f\t%0.2f\t%0.2f\n',RADIUS(1),MASS(1),PARA(1));
for ii = 2:length(yy0)
    fprintf('%d\t',ii);
    % use former solution as the initial guess for the coming one
    sol_rel{ii} = bvp4c(...
        @(x,y,p2)([y(2);-p2.*(y(1).^2-(1./yy0(ii)).^2).^(1.5)]),...
        @(ya,yb,p2)([ya(2);yb(1)-1./yy0(ii);ya(1)-1]),...
        bvpinit(sol_rel{ii-1}.x,@(x)(deval(sol_rel{ii-1},x)),PARA(ii-1).^2),...
        bvpset('SingularTerm',[0,0;0,-2])...
    );
    PARA(ii) = real(sqrt(sol_rel{ii}.parameters));
    RADIUS(ii) = 2.*PARA(ii)./yy0(ii);
    MASS(ii) = -PARA(ii).*real(sol_rel{ii}.y(2,end));
    fprintf('%0.2f\t%0.2f\t%0.2f\n',RADIUS(ii),MASS(ii),PARA(ii));
end

nicefigure1;
setpaper2([10 10]);
axes('position',[0.1 0.1 0.8 0.8]);
hold on; box on;
phdl = []; chdl = {};
phdl = [phdl,plot(MASS,RADIUS,'k-o','LineWidth',1.5)]; hold on;
chdl = [chdl,{'Relativistic'}];
RADIUSn = 3.9:0.01:RADIUS(end);
MASSn = -real(sol_nonrel.y(2,end)) .* RADIUSn.^(-3);
phdl = [phdl,plot(MASSn,RADIUSn,'--','LineWidth',1.5,'Color',[0.5 0.5 0.5])];
chdl = [chdl,{'Non-Relativistic'}];
lhdl = legend(phdl,chdl);
set(lhdl,'Box','off');
xlabel('\it{M} / \it{M}_{0}','FontName','Times New Roman');
ylabel('\it{R} / \it{R}_{0}','FontName','Times New Roman');
% set(gca, 'XLim', [0 2.2]); 
set(gca, 'XTickLabel', sprintf('%0.1f|',get(gca,'XTick')));
% set(gca, 'YTick', []);

save('result.mat','sol_rel','MASS','RADIUS','PARA','yy0','MASSn','RADIUSn');

save('mr.mat','MASS','RADIUS','PARA','yy0');
save('mrn.mat','MASSn','RADIUSn');




XX = cell(1,length(yy0)+1);
YY = cell(1,length(yy0)+1);
DENSE = cell(1,length(yy0)+1);
XX{1} = sol_ultrel.x; 
YY{1} = real(sol_ultrel.y(1,:));
DENSE{1} = (YY{1}./YY{1}(1)).^3;
for ii = 2:length(yy0)
    XX{ii} = real(sol_rel{ii}.x);
    YY{ii} = real(sol_rel{ii}.y(1,:));
    DENSE{ii} = (((yy0(ii).*YY{ii}).^2-1)./((yy0(ii).*YY{ii}(1)).^2-1)).^(1.5);
end
XX{end} = sol_nonrel.x; 
YY{end} = real(sol_nonrel.y(1,:));
DENSE{end} = (YY{end}./YY{end}(1)).^(1.5);

nicefigure1;
setpaper2([10 10]); axes('position',[0.1 0.1 0.8 0.8]);
hold on; box on; cmap = colormap(hot);
phdl = []; chdl = {}; sel = [1 2 4 6 7 8 9 10 11 13 17];
for ii = sel
    phdl = [phdl, plot(XX{ii},DENSE{ii},'LineWidth',1.5,...
        'Color',cmap(floor((ii-1)./(max(sel).*1.5).*length(cmap))+1,:))]; 
    if (ii==1)
        chdl = [chdl, { '1/\it{\gamma}_c = \rm{0.00}' } ];
    elseif (ii==sel(end))
        chdl = [chdl, { '1/\it{\gamma}_c = \rm{1.00}' } ];
    else
        chdl = [chdl, { strcat('1/\it{\gamma}_c','\rm',...
        sprintf(' = %0.2f',1./yy0(ii))) } ];
    end
    hold on;
end
lhdl = legend(phdl,chdl);
set(lhdl,'Box','off');
xlabel('\it{r} / \it{R}','FontName','Times New Roman');
ylabel('\it{\rho} / \it{\rho}_{c}','FontName','Times New Roman');
set(gca, 'XTickLabel', sprintf('%0.1f|',get(gca,'XTick')));
set(gca, 'YTickLabel', sprintf('%0.1f|',get(gca,'YTick')));




nicefigure1;
setpaper2([10 10]); axes('position',[0.1 0.1 0.8 0.8]);
hold on; box on; cmap = colormap(hot);
phdl = []; chdl = {}; sel = [1 2 4 6 7 8 9 10 11 13 17];
% phdl = [phdl, plot(real(sol_rel{1}.x),real(sol_rel{1}.y(1,:)),'LineWidth',1.5,...
%         'Color',cmap(1,:))]; 
% chdl = [chdl, { '\it{R}/\it{R}_0\rm = 0.00'} ];
for ii = sel(2:(end-1))
    phdl = [phdl, plot(XX{ii},YY{ii},'LineWidth',1.5,...
        'Color',cmap(floor((ii-1)./(max(sel).*1.5).*length(cmap))+1,:))]; 
    chdl = [chdl, { strcat('\it{R}/\it{R}_0','\rm',...
        sprintf(' = %0.2f',RADIUS(ii))) } ];
    hold on;
end
% phdl = [phdl, plot(0:0.01:1,ones(1,101),'LineWidth',1.5,...
%         'Color',cmap(floor((sel(end)-1)./(max(sel).*1.5).*length(cmap))+1,:))]; 
% chdl = [chdl, { '\it{R}/\it{R}_0\rm = \infty'} ];
lhdl = legend(phdl,chdl,'Location','SouthWest');
set(lhdl,'Box','off');
xlabel('\it{r} / \it{R}','FontName','Times New Roman');
ylabel('\it{\gamma} / \it{\gamma}_{c}','FontName','Times New Roman');
set(gca, 'XTickLabel', sprintf('%0.1f|',get(gca,'XTick')));
set(gca, 'YTickLabel', sprintf('%0.1f|',get(gca,'YTick')));



nicefigure1;
setpaper2([10 10]);
axes('position',[0.1 0.1 0.8 0.8]);
hold on; box on;
plot([0,1./yy0(2:end)],RADIUS,'k-o','LineWidth',1.5); hold on;
xlabel('1 / \it{\gamma}_{c}','FontName','Times New Roman');
ylabel('\it{R} / \it{R}_{0}','FontName','Times New Roman');
set(gca, 'XTickLabel', sprintf('%0.1f|',get(gca,'XTick')));
% set(gca, 'YTickLabel', sprintf('%0.1f|',get(gca,'YTick')));
% set(gca, 'YTick', []);