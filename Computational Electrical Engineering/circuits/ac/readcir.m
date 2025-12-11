function [nodes,types,labels,vals,err]=readcir (fnam)

disp('Reading netlist...');disp(' ');

nodes=[];
types=[];
labels=[];
vals=[];

err=0;
fid=fopen(fnam,'r');
if fid==-1
    err=fid;
   return
end

F = textscan(fid, '%d %d %s %n %n');
nodes(:,1)=F{1};
nodes(:,2)=F{2};
labels=F{3};
vals(:,1)=F{4};
vals(:,2)=F{5};
types=char(labels);
types=types(:,1);
fclose(fid);

nn=max(max(nodes));
nl=size(nodes,1);

ss=sprintf('Circuit has %d nodes and %d ports',nn,nl);
disp(ss);

nv=numel(find(types=='V'));
ni=numel(find(types=='I'));

ss=sprintf('Circuit has %d voltage source/s and %d current source/s \n',nv,ni);
disp(ss);

