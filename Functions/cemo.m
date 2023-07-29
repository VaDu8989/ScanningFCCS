function out=cemo(x,a)
  
  %central moment of order a for vector x
  
  out=mean((x-mean(x)).^a);

  
end
