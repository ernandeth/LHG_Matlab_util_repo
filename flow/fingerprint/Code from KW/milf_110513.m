% im=load_nii('vol_e2264_11_13_113_0191.nii');
% %%Select slice
% slice=10;
% COFF_im=im.img;
% sub_pairs=squeeze(COFF_im(:,:,slice,1:2:end)-COFF_im(:,:,slice,2:2:end));
% figure;for n=1:size(sub_pairs,3),imagesc(abs(sub_pairs(:,:,n))),pause,end
% 
% nxy=size(COFF_im,1);
% nt=size(COFF_im,4);
% % obssignals=double(abs(squeeze(im.img(20,40,4,:)))).';
% slice1=squeeze(double(COFF_im(:,:,slice,:)));
% figure;imagesc(abs(squeeze(slice1(:,:,1))));
% mask=roipoly();
% mask_big=repmat(mask,[1,1,nt]);
% slice1_mask=slice1.*mask_big;
% slice1_mask=reshape(slice1_mask,[nxy*nxy,nt]);
% slice1_mask=mean(slice1_mask(find(slice1_mask(:,1)),:),1);

slice1_mask=[];%%%%% Set this parameter to the ROI data

Dictionary Generation
useFileDictionary=0;

if ~useFileDictionary
    fprintf('\nGenerating Dictionary ... \n')
    [dict, parms] = gen_dictionary;
    fprintf('\n ... done ')
else
    load dictionary.mat
end

nim=size(COFF_im,4)-2;
fprintf('\n Comparing the observed signal at each voxel with each dictionary entry ... \n');
P = size(slice1_mask,1);
index_big=zeros(1,P);
M0mat=zeros(1,P);
IP_big=zeros(1,P);
recon_save =zeros(nim,P);
for m=1:P
    if slice1_mask(m,1)==0
    else
        obssignals=slice1_mask(m,3:end);
        fprintf('\r percent pixels .... %03f', 100*m/P);
        dict_temp=abs(dict).';
        dict_temp_norm = sqrt(sum(dict_temp.*conj(dict_temp)));
        
        inner_product = conj(obssignals)*dict_temp./dict_temp_norm./norm(obssignals);
        [best index] = max(abs(inner_product));
        dictcol=squeeze(dict_temp(:,index));
        coef_vector =  pinv(dictcol)*obssignals.';
        recon_curve = dictcol*coef_vector;
        M0mat(m)=coef_vector;
        index_big(m)=index;
        recon_save(:,m) = recon_curve;
        IP_big(m)=best;
    end
end

index=index_big;
figure;plot(1:nim,squeeze(slice1_mask(index,3:end)),1:nim,squeeze(recon_save(:,index)))
parms(index)
% index2=find(IP_big==min(nonzeros(IP_big)));
% figure;plot(1:nim,squeeze(slice1_mask(index2,3:end)),1:nim,squeeze(recon_save(:,index2)))

% figure;
% subplot(2,1,1)
% % plot(1:270,squeeze(slice1_mask(find(slice1_mask(:,4)>0),3:end)))
% plot(1:nim,squeeze(slice1_mask(index,3:end)))
% subplot(2,1,2)
% plot(1:nim,squeeze(dict(:,:)*M0mat(index)).')

% indexmap=zeros(1,P);
% for n=1:P
%     if index_big(n)>0
%         indexmap(n)=index_big(n);
%     end
% end
% figure;imagesc(reshape(indexmap,[nxy,nxy]));colormap hot;
% 
% IPmap=zeros(1,P);
% for n=1:P
%     if index_big(n)>0
%         IPmap(n)=IP_big(n);
%     end
% end
% figure;imagesc(reshape(IPmap,[nxy,nxy]));colormap hot;

% figure;imagesc((rot90(squeeze(slice1(:,:,1)))));colormap gray; axis image; axis off;
% 
% figure;imagesc(rot90((reshape(M0mat./max(M0mat),[nxy,nxy]))));colormap hot;axis image;axis off;title('M_o Map')
% 
% T1map=zeros(1,P);
% for n=1:P
%     if index_big(n)>0
%         T1map(n)=1./parms(index_big(n)).r1tis;
%     end
% end
% figure;imagesc(rot90((reshape(T1map,[nxy,nxy]))));colormap hot;axis image; axis off;title('T_1 Map (s)')
% 
% fmap=zeros(1,P);
% for n=1:P
%     if index_big(n)>0
%         fmap(n)=parms(index_big(n)).f;
%     end
% end
% figure;imagesc(rot90((reshape(fmap,[nxy,nxy]))));colormap hot;title('CBF Map (ml/s/g)');axis image; axis off;
% 
% ttmap=zeros(1,P);
% for n=1:P
%     if index_big(n)>0
%         ttmap(n)=parms(index_big(n)).transit;
%     end
% end
% figure;imagesc(rot90(reshape(ttmap,[nxy,nxy])));colormap hot;title('Transit Time Map (s)');axis image; axis off;
% 
% 
% 
% 
% 
