function [firstgoodindex,badframes,height,width]=timelapsesetup_PRIMO(mem_temp,badframes)
%%% determine median cell size for blur detection %%%%%%%%%%%%%%%%%%%%%%%%%
nuc_area=zeros(3,1); numcells=zeros(3,1);
dims=size(mem_temp);
height=dims(1); width=dims(2);
%%% determine first good frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstgoodindex=1;
 %was 0.5
end

