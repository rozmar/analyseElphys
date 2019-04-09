% createScale calculates values of a range with a given stepsize and with a
% given scale type. Scale type can be linear or logarithmic. If stepsize
% hasn't been given, default values will be used.

function scaleVector = createScale(range, step, scale)

  %% -------------------------
  %  Constant definitions
  %% -------------------------
  LINEAR = 1;
  LOGARITHMIC = 2;
  %% -------------------------

  if nargin<3
    scale = LINEAR;
  end
  
  if nargin<2
    step = 0;
  end
    
%   if mod(diff(range),max(1, step))~=0
%     warning('Stepsize is not whole fraction of the range, some values will be cut from the end.');
%   end
    
  if scale==LINEAR
    % Between 20 and 50 the default value 1
    if step<=0
      scaleVector = range(1):1:min(range(2), 50);  
    else
      scaleVector = range(1):step:min(range(2), 50);
    end
    if range(2)<=50 ; return; end;
    % Above 50 the default value 5
    if step==0
      scaleVector = [ scaleVector(1:end-1) , max(50, range(1)):5:range(2) ];
    else
      scaleVector = [ scaleVector(1:end-1) , max(50, range(1)):step:range(2) ];
    end
  elseif scale==LOGARITHMIC  
%     if range(2)<=50
%       nValue = min(range(2), 50)-min(range(1), 50)+1;
%       scaleVector = logspace(log10(min(range(1), 50)), log10(min(range(2), 50)), nValue);
%     else
        if step<=0
            nValue = ceil(range(2)-range(1))+1;
        else
            nValue = round(ceil(range(2)-range(1))/step)+1;
        end 
      scaleVector = logspace(log10(range(1)), log10(range(2)), nValue);
%     end
  end
    
end