function IFD_calc_save(sp,sname,SubID,R1_NGR_path,R2_NGR_path,R1_GSR_path,R2_GSR_path)

%% check ID
% Ensure the correspondings for each subject

% si1=cellfun(@(x) str2double(x(end-38:end-33)),R1_NGR_path);
% si2=cellfun(@(x) str2double(x(end-38:end-33)),R2_NGR_path);
% si3=cellfun(@(x) str2double(x(end-38:end-33)),R1_GSR_path);
% si4=cellfun(@(x) str2double(x(end-38:end-33)),R2_GSR_path);
% if sum(SubID-si1)+sum(SubID-si2)+sum(SubID-si3)+sum(SubID-si4)
%     exit
% end

%% load gradients for asymmetric distance
SubSize=length(SubID);
gmats=zeros(SubSize,9354);
IFdist=repmat({gmats,gmats,gmats,gmats},10,1);
for i=1:SubSize
    % load gradients
    tmpg1=load(R1_NGR_path{i});
    tmpg2=load(R2_NGR_path{i});
    tmpg3=load(R1_GSR_path{i});
    tmpg4=load(R2_GSR_path{i});   
    % calculate the between-hemisphere functional distance in multiple
    % dimensions (from 1 to 10)
    for j =1:10

        gradL1=tmpg1.grads_L(:,[1:j]);
        gradR1=tmpg1.realigned_RtoL(:,[1:j]);
        IFdist{j,1}(i,:)=sqrt(sum((gradL1 - gradR1) .^ 2,2))'; % R1 NGR
        
        gradL2=tmpg2.grads_L(:,[1:j]);
        gradR2=tmpg2.realigned_RtoL(:,[1:j]);
        IFdist{j,2}(i,:)=sqrt(sum((gradL2 - gradR2) .^ 2,2))'; % R2 NGR
        
        gradL3=tmpg3.grads_L(:,[1:j]);
        gradR3=tmpg3.realigned_RtoL(:,[1:j]);
        IFdist{j,3}(i,:)=sqrt(sum((gradL3 - gradR3) .^ 2,2))'; % R1 GSR
        
        gradL4=tmpg4.grads_L(:,[1:j]);
        gradR4=tmpg4.realigned_RtoL(:,[1:j]);
        IFdist{j,4}(i,:)=sqrt(sum((gradL4 - gradR4) .^ 2,2))'; % R2 GSR
    end 
end
%% calculate the ICC between 2 session
IFdist_whole_AVG = cellfun(@(x) mean(x,2), IFdist, 'UniformOutput',false);
ICC_IFdist_whole=zeros(10,2);
ICC_IFdist=cell(10,2);

for i=1:10
	ICC_IFdist_whole(i,1)=IPN_icc([IFdist_whole_AVG{i,1},IFdist_whole_AVG{i,2}],2,'single');
    ICC_IFdist_whole(i,2)=IPN_icc([IFdist_whole_AVG{i,3},IFdist_whole_AVG{i,4}],2,'single');
    ICC_vertices=zeros(2,9354);
    for j = 1:9354
        ICC_vertices(1,j) = IPN_icc([IFdist{i,1}(:,j),IFdist{i,2}(:,j)],2,'single');
        ICC_vertices(2,j) = IPN_icc([IFdist{i,3}(:,j),IFdist{i,4}(:,j)],2,'single');
    end
    ICC_IFdist{i,1} = ICC_vertices(1,:);
    ICC_IFdist{i,2} = ICC_vertices(2,:);
end

%% average across sessions 
IFdist_AVG=cell(10,2);
IFdistAll_AVG={zeros(SubSize,10),zeros(SubSize,10)};
for i=1:10
    IFdist_AVG{i,1}=(IFdist{i,1}+IFdist{i,2})/2; %NSR
    IFdist_AVG{i,2}=(IFdist{i,3}+IFdist{i,4})/2; %GSR
    IFdistAll_AVG{1}(:,i)=(IFdist_whole_AVG{i,1}+IFdist_whole_AVG{i,2})/2; %NSR
    IFdistAll_AVG{2}(:,i)=(IFdist_whole_AVG{i,3}+IFdist_whole_AVG{i,4})/2; %GSR
end

%% save results
save([sp,filesep,sname,'.mat'], 'SubID','IFdist','IFdist_whole_AVG','ICC_IFdist','ICC_IFdist_whole','IFdist_AVG','IFdistAll_AVG','-nocompression', '-v7.3');

end