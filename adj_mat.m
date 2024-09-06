function mat=adj_mat(data)
uniq_node=(unique(data));
idx=grp2idx(uniq_node);
for i=1:length(uniq_node)
    [q1]=find(strcmp(data(:,1),uniq_node(i)));
    data_idx(q1,1)=idx(i);
    [q2]=find(strcmp(data(:,2),uniq_node(i)));
    data_idx(q2,2)=idx(i);
end
for i=1:length(uniq_node)
    mat(idx(i),data_idx(find(data_idx(:,1)==idx(i)),2))=1;
end