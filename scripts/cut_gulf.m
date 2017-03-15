function field = cut_gulf(field,mask,x_cut,flag_cut,ybc)

[si_x,si_y] = size(field);
if size(field,2)>1
  if flag_cut == 1
    for nx = 1:x_cut
      field(end-x_cut+nx,:) = field(end-x_cut+nx,:).*(1-mask(nx,:)) + field(nx,:).*mask(nx,:);
    end
    field = field(x_cut+1:end,:);
    
  elseif flag_cut == 2
    nx2 = 0;
    nx0 = x_cut+1;
    for nx = nx0: si_x
      nx2 = nx2 + 1;
      field(nx2,:) = field(nx2,:).*(1-mask(nx,:)) + field(nx,:).*mask(nx,:);
    end
    field = field(1:x_cut,:);
  end
else
  if flag_cut == 1
    for nx = 1:x_cut
      field(end-x_cut+nx) = field(end-x_cut+nx).*(1-mask(nx,ybc)) + field(nx).*mask(nx,ybc);
    end
    field = field(x_cut+1:end);
    
  elseif flag_cut == 2
    nx2 = 0;
    nx0 = x_cut+1;
    for nx = nx0: si_x
      nx2 = nx2 + 1;
      field(nx2) = field(nx2).*(1-mask(nx,ybc)) + field(nx).*mask(nx,ybc);
    end
    field = field(1:x_cut);
  end
end  