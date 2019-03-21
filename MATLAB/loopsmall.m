for i =1:20
  x =2+1.0/(8^i) 
  y = atan(x)-atan(2)
  z = (8^i)*y %#ok<*NOPTS>
  disp(i);
    
end     