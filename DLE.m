function [error1,error2,error3] = DLE(e_o, e_m, vert, n_dip)
% Comoutes the Dipole localization error
%compares the coordenates of the original dipole localization
%against the method estimated dipole localization
% e_o = original dipole location (head model)
% e_m = estimated dipole location (head model)
% vert = vertex from the head model
% error = euclidean distance between both coordenates.+

[~,a] = max(e_o);
[~,b] = max(e_m);
l_o = vert(a,:);
l_e = vert(b,:);
if n_dip<=1
    error1 = sqrt((1/length(l_e))*sum((l_o-l_e).^2));
    error2 = norm(l_o-l_e);
    l1 = a-round(length(e_o)*0.01);
    l2 = a+round(length(e_o)*0.01);
    l3 = b-round(length(e_m)*0.01);
    l4 = b+round(length(e_m)*0.01);
    
    if l1<=0
        l1=1;
    end
    if l2>length(e_o)
        l2=length(e_o);
    end
    if l3<=0
        l3=1;
    end
    if l4>length(e_m)
        l4=length(e_m);
    end
        
        
%         if b-round(length(e_m)*0.01)<0
%             error3 = abs((var(e_o(1 : a+round(length(e_o)*0.01))) ...
%             - var(e_m(1 : b+round(length(e_m)*0.01))))/ ...
%             var(e_o(1 : a+round(length(e_o)*0.01))))*100;
%         elseif b+round(length(e_m)*0.01)>length(e_m)
%             error3 = abs((var(e_o(1 : a+round(length(e_o)*0.01))) ...
%             - var(e_m(b-round(length(e_m)*0.01) : length(e_m))))/ ...
%             var(e_o(1 : a+round(length(e_o)*0.01))))*100;
%         end
%     elseif a+round(length(e_o)*0.01)>length(e_o)
%         
%         if b-round(length(e_m)*0.01)<0
%             error3 = abs((var(e_o(a-round(length(e_o)*0.01) : length(e_o))) ...
%             - var(e_m(1 : b+round(length(e_m)*0.01))))/ ...
%             var(e_o(a-round(length(e_o)*0.01) : length(e_o))))*100;
%         elseif b+round(length(e_m)*0.01)>length(e_m)
%             error3 = abs((var(e_o(a-round(length(e_o)*0.01) : length(e_o))) ...
%             - var(e_m(b-round(length(e_m)*0.01) : length(e_m))))/ ...
%             var(e_o(a-round(length(e_o)*0.01) : length(e_o))))*100;
%         end
    %else
        error3 = abs((var(e_o(l1 : l2))- var(e_m(l3 : l4)))/var(e_o(l1 : l2)))*100;
    end

    %error1 = 1;
end
