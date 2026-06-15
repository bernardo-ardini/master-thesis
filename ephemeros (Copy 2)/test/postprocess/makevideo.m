%load 110006.mat; % large
%load 110037.mat; % von mises
load 110100.mat; % narrow

figure(3);
clf;
theme(gcf,"light");

times=arrayfun(@(s) md.hist{s}.t,1:length(md.hist));

t=33.6;

[~,s]=min(abs(times-t));
h=md.hist{s};

p=@(T) -1/2*sum(reshape(eye(2,2),2,2,1).*T,[1,2]);

geo.plotfield(h.al{2,1},2,@(al) al,displ=h.displ,scale=10);
shading interp;
pbaspect([1,1,1]);
daspect([1,1,1]);
axis off;
xlim([-1.1*2,1.1*2]);
ylim([-1.1*1,1.1*1]);
colormap parula(30);
colorbar;
title("\alpha");

% geo.plotfield(h.Te,2,p,displ=h.displ,scale=10);
% shading interp;
% pbaspect([1,1,1]);
% daspect([1,1,1]);
% axis off;
% xlim([-1.1*2,1.1*2]);
% ylim([-1.1*1,1.1*1]);
% colormap parula(30);
% colorbar;
% title("p");

%saveas(gcf,"pressure4.pdf");