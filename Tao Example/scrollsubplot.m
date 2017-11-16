function theAxis = scrollsubplot(nrows, ncols, thisPlot)
%SCROLLSUBPLOT Create axes in tiled positions.
%   SCROLLSUBPLOT(m,n,p), breaks the Figure window into
%   an m-by-n matrix of small axes, selects the p-th axes for 
%   for the current plot, and returns the axis handle.  The axes 
%   are counted along the top row of the Figure window, then the
%   second row, etc. 
%   Example:
% 
%       scrollsubplot(3,1,-1), plot(income)
%       scrollsubplot(3,1,1), plot(declared_income)
%       scrollsubplot(3,1,2), plot(tax)
%       scrollsubplot(3,1,3), plot(net_income)
%       scrollsubplot(3,1,4), plot(the_little_extra)
% 
%   plots declared_income on the top third of the window, tax in
%   the middle, and the net_income in the bottom third. Above the
%   top of the figure income is ploted and below the lower edge
%   the_little_extra is to be found. To navigate there is a slider
%   along the right figure edge.
% 
%   SCROLLSUBPLOT now also work with less regular subplot-layouts:
%
%   axs3(1) = scrollsubplot(3,3,1);
%   axs3(2) = scrollsubplot(3,3,3);
%   axs3(3) = scrollsubplot(3,3,[4,7]);
%   axs3(4) = scrollsubplot(3,3,5);
%   axs3(5) = scrollsubplot(3,3,[3,6]);
%   axs3(6) = scrollsubplot(3,3,10);
%   axs3(7) = scrollsubplot(3,3,[8,9,11,12]);
%   axs3(8) = scrollsubplot(3,3,[13,14,15]);
%   for i1 = 1:8,
%     if ishandle(axs3(i1))
%       axes(axs3(i1))
%       imagesc(randn(2+3*i1))
%     end
%   end
%
%   The function works well for regular grids where m,n is constant
%   for all p. When m,n varies there is no guarantee that the steps
%   of the slider is nicely adjusted to the sizes of the
%   subplots, further the slider also responds to mouse-wheel scrolling.
%
%   Differences with SUBPLOT: SCROLLSUBPLOT requires 3 input
%   arguments, no compatibility with subplot(323), no handle as
%   input. Further  PERC_OFFSET_L is decreased from 2*0.09 to 0.07
%   and PERC_OFFSET_R is decreased from 2*0.045 to 0.085. This
%   leaves less space for titles and labels, but give a plaid grid
%   of subplots even outside the visible figure area.
%   
%   Bug/feature when the slider is shifted from its initial
%   position and then extra subplots is added, they get
%   mis-positioned.
%   
%   See also SUBPLOT,

%   Copyright © Bjorn Gustavsson 20050526-2014, Modification/extension of
%   Mathworks subplot.
%   Version 2, modified from version 1 in that scroll now is a subfunction,
%   together with 

persistent maxrownr minrownr

%%% This is a matlab bug(?) that blanks axes that have been rescaled to
%%% accomodate a colorbar. I've not tried to fix it. BG
% we will kill all overlapping siblings if we encounter the mnp
% specifier, else we won't bother to check:
narg = nargin;
% kill_siblings = 0;
create_axis = 1;
delay_destroy = 0;
tol = sqrt(eps);
if narg ~= 3 % not compatible with 3.5, i.e. subplot ==
             % subplot(111) errors out
  error('Wrong number of arguments')
end

%check for encoded format
handle = '';
position = '';

kill_siblings = 1;

% if we recovered an identifier earlier, use it:
if(~isempty(handle))
  
  set(get(0,'CurrentFigure'),'CurrentAxes',handle);
  
elseif(isempty(position))
  % if we haven't recovered position yet, generate it from mnp info:
  if (min(thisPlot) < 1)&&0 % negative Thisplot corresponds to
                            % panels above top row, that is OK
    error('Illegal plot number.')
  else
    % This is the percent offset from the subplot grid of the plotbox.
    PERC_OFFSET_L = 0.07;
    PERC_OFFSET_R = 0.085;
    PERC_OFFSET_B = PERC_OFFSET_L;
    PERC_OFFSET_T = PERC_OFFSET_R;
    if nrows > 2
      PERC_OFFSET_T = 0.9*PERC_OFFSET_T;
      PERC_OFFSET_B = 0.9*PERC_OFFSET_B;
    end
    if ncols > 2
      PERC_OFFSET_L = 0.9*PERC_OFFSET_L;
      PERC_OFFSET_R = 0.9*PERC_OFFSET_R;
    end

    % Subplots version:
    % row = (nrows-1) -fix((thisPlot-1)/ncols)
    % col = rem (thisPlot-1, ncols);
    % Slightly modified to allow for having negative thisPlot
    row = (nrows-1) -floor((thisPlot-1)/ncols);
    col = mod (thisPlot-1, ncols);
    
    % From here on to line 190 essentially identical to SUBPLOT (==untouched)
    
    % For this to work the default axes position must be in normalized coordinates
    if ~strcmp(get(gcf,'defaultaxesunits'),'normalized')
      warning('DefaultAxesUnits not normalized.')
      tmp = axes;
      set(axes,'units','normalized')
      def_pos = get(tmp,'position');
      delete(tmp)
    else
      def_pos = get(gcf,'DefaultAxesPosition')+[-.05 -.05 +.1 +.05];
    end
    col_offset = def_pos(3)*(PERC_OFFSET_L+PERC_OFFSET_R)/ ...
	(ncols-PERC_OFFSET_L-PERC_OFFSET_R);
    row_offset = def_pos(4)*(PERC_OFFSET_B+PERC_OFFSET_T)/ ...
	(nrows-PERC_OFFSET_B-PERC_OFFSET_T);
    totalwidth = def_pos(3) + col_offset;
    totalheight = def_pos(4) + row_offset;
    width = totalwidth/ncols*(max(col)-min(col)+1)-col_offset;
    height = totalheight/nrows*(max(row)-min(row)+1)-row_offset;
    position = [def_pos(1)+min(col)*totalwidth/ncols ...
		def_pos(2)+min(row)*totalheight/nrows ...
		width height];
    if width <= 0.5*totalwidth/ncols
      position(1) = def_pos(1)+min(col)*(def_pos(3)/ncols);
      position(3) = 0.7*(def_pos(3)/ncols)*(max(col)-min(col)+1);
    end
    if height <= 0.5*totalheight/nrows
      position(2) = def_pos(2)+min(row)*(def_pos(4)/nrows);
      position(4) = 0.7*(def_pos(4)/nrows)*(max(row)-min(row)+1);
    end
  end
end

% kill overlapping siblings if mnp specifier was used:
nextstate = get(gcf,'nextplot');
if strncmp(nextstate,'replace',7), nextstate = 'add'; end
if(kill_siblings)
  if delay_destroy
    if nargout %#ok<UNRCH>
      error('Function called with too many output arguments')
    else
      set(gcf,'NextPlot','replace'); return,
    end
  end
  sibs = get(gcf, 'Children');
  got_one = 0;
  for i1 = 1:length(sibs)
    if(strcmp(get(sibs(i1),'Type'),'axes'))
      units = get(sibs(i1),'Units');
      set(sibs(i1),'Units','normalized')
      sibpos = get(sibs(i1),'Position');
      set(sibs(i1),'Units',units);
      intersect = 1;
      if(     (position(1) >= sibpos(1) + sibpos(3)-tol) || ...
	      (sibpos(1) >= position(1) + position(3)-tol) || ...
	      (position(2) >= sibpos(2) + sibpos(4)-tol) || ...
	      (sibpos(2) >= position(2) + position(4)-tol))
	intersect = 0;
      end
      if intersect
	if got_one || any(abs(sibpos - position) > tol)
	  delete(sibs(i1));
	else
	  got_one = 1;
	  set(gcf,'CurrentAxes',sibs(i1));
	  if strcmp(nextstate,'new')
	    create_axis = 1;
	  else
	    create_axis = 0;
	  end
	end
      end
    end
  end
  set(gcf,'NextPlot',nextstate);
end

% create the axis:
if create_axis
  if strcmp(nextstate,'new'), figure, end
  ax = axes('units','normal','Position', position);
  set(ax,'units',get(gcf,'defaultaxesunits'))
else 
  ax = gca; 
end


% return identifier, if requested:
if(nargout > 0)
  theAxis = ax;
end

%% SCROLLSUBPLOT modification part
%
% From here on out set up scrollbar if needed
scroll_hndl = findall(gcf,'Type','uicontrol','Tag','scroll');

ax_indx = findall(gcf,'Type','axes');

if length(ax_indx)==1
  maxrownr = -inf;
  minrownr = inf;
end

%% This is setting up the scroll-slider (invisible)
if isempty(scroll_hndl)
  uicontrol('Units','normalized',...
            'Style','Slider',...
            'Position',[.98,0,.02,1],...
            'Min',0,...
            'Max',1,...
            'Value',1,...
            'visible','off',...
            'Tag','scroll',...
            'Callback',@(scr,event) scroll(1));
  set(gcf,'WindowScrollWheelFcn',@wheelScroll)      
end
% making it visible when needed
if ( nrows*ncols < max(thisPlot) || min(thisPlot) < 1 )
  set(scroll_hndl,'visible','on')
end
scroll(1)


maxrownr = max(maxrownr(:),max(-row(:)));
minrownr = min(minrownr(:),min(-row(:)));

% Adjust the slider step-sizes to account for the number of rows of
% subplots: 
set(scroll_hndl,...
    'sliderstep',[1/nrows 1]/(1/((nrows)/max(1,1+maxrownr(:)-minrownr(:)-nrows))))
set(scroll_hndl,...
    'value',1)

   function scroll(old_val)
   % SCROLL - Scroll subplots vertically
   %   Used by scrollsubplot.
   % Calling:
   %   scroll(old_val)
   % Set as callback function handle:
   %   
   %   See also SCROLLSUBPLOT
   
   % Copyright Bjorn Gustavsson 20050526
   % Version 2, modified to become a subfunction of scrollsubplot.
   
   %% Scroll the subplot axeses
   %  That is change their position some number of steps up or down
   
   % Get the handle of the scroll-slider handle and all the handle
   % to all the subplot axes
   clbk_ui_hndl = findall(gcf,'Type','uicontrol','Tag','scroll');
   ax_hndl = findall(gcf,'Type','axes');
   
   for i2 = length(ax_hndl):-1:1,
     a_pos(i2,:) = get(ax_hndl(i2),'position');
   end
   
   pos_y_range = [min(.07,min(a_pos(:,2))) max(a_pos(:,2) + a_pos(:,4) )+.07-.9];
   
   val = get(clbk_ui_hndl,'value');
   step = ( old_val - val) * diff(pos_y_range);
   
   
   for i2 = 1:length(ax_hndl),
     set(ax_hndl(i2),'position',get(ax_hndl(i2),'position') + [0 step 0 0]);
   end
   
   set(clbk_ui_hndl,'callback',@(scr,event) scroll(val));
   
   end % end scroll

   function wheelScroll(src,evnt) %#ok<INUSL>
   % wheelScroll - mouse-wheel wrapper to scroll-function
   % Calling:
   %   wheelScroll(src,evnt)
   % Set as calback
   %   set(gcf,'WindowScrollWheelFcn',@wheelScroll)
   try
     C = findobj(gcf,'type','uicontrol');
     MaxVal = get(C,'Max');%  = [1]
	 minval = get(C,'Min');% = [0]
	 sliderstep = get(C,'SliderStep'); % = [0.0555556 0.222222]
     oldval = get(C,'value');
     if evnt.VerticalScrollCount > 0 
       set(C,'value',max(minval,oldval-sliderstep(1)/10))
     elseif evnt.VerticalScrollCount < 0 
       set(C,'value',min(MaxVal,oldval+sliderstep(1)/10))
     end
     scroll(oldval)
   catch
   end
   end %wheelScroll

end
