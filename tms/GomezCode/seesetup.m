function seesetup(loop,jvec,sigma,mu)
%%% plots cerry simulation setup
%%% loop is a cell array loop{number of loops}(number of elements,
%%% direction)
%%% if sigma is <10^-5 it is assumed zero
%%% mu=1 is consided outside of boundary
global dx dy dz nx ny nz
for i=1:length(loop)
    %   l=sqrt(jvec{i}(:,1).^2+jvec{i}(:,2).^2+jvec{i}(:,3).^2);
    l=1;
    hold on;
%quiver3(loop{i}(:,1),loop{i}(:,2),loop{i}(:,3),jvec{i}(:,1)./l,jvec{i}(:,2)./l,jvec{i}(:,3)./l);
    plot3(loop{i}(:,1),loop{i}(:,2),loop{i}(:,3),'b');
end
 axis([-0 0.35 -0 0.35 -0 0.35])
sigma(find(abs(sigma)>10^-5))=1;
mu(find(mu~=1))=2;
mu=mu-1;
for i=1:3:nz
    [B,L] = bwboundaries(sigma(:,:,i));
    for l=1:length(B)
        plot3(dx*B{l}(:,1),dy*B{l}(:,2),i*dz*ones(size(B{l}(:,1))),'r');
    end
end
%  axis([-0 0.35 -0 0.35 -0 0.35])
% for i=1:3:nz
%     [B,L] = bwboundaries(mu(:,:,i));
%     for l=1:length(B)
%         plot3(dx*B{l}(:,1),dy*B{l}(:,2),i*dz*ones(size(B{l}(:,1))),'g');
%     end
% end
 axis([-0 0.35 -0 0.35 -0 0.35])
view([0.6820    0.7314    0.0000 -0.7067; ...
     -0.3206    0.2990    0.8988   -0.4386;...
     -0.6573    0.6130   -0.4384    8.9016;...
           0         0         0    1.0000])
% for i=1:nx
% [B,L] = bwboundaries(reshape(sigma(i,:,:),[ny nz]));
% for l=1:length(B)
%     plot3(i*dx*ones(size(B{l}(:,1))),dy*B{l}(:,1),dz*B{l}(:,2),'r');
% end
% end
%
% for i=1:nx
% [B,L] = bwboundaries(reshape(mu(i,:,:),[ny nz]));
% for l=1:length(B)
%     plot3(i*dx*ones(size(B{l}(:,1))),dy*B{l}(:,1),dz*B{l}(:,2),'g');
% end
% endclear all;
