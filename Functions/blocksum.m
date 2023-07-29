function out=blocksum(A,blocksize)
newsize=floor(floor(size(A,1)/blocksize)/1000)*1000;
out=zeros(newsize, size(A,2), size(A,3));

for i=1:size(out,1)
    out(i,:,:)=sum(A(blocksize*(i-1)+1:blocksize*i,:,:),1);
    
    
end
