function tighten_margin(fig, ax, h)
    % auto tighten the margin of figure 
    % while keeping the ratio
    arguments
        fig;
        ax;
        % default height is 1000 px
        h = 1000;
    end
    fig.Units = 'pixels';
    ax.Units = 'pixels';
    % ax position left, bottom, width, height
    axp_l = ax.Position(1);
    axp_b = ax.Position(2);
    axp_w = ax.Position(3);
    axp_h = ax.Position(4);
    % ax smallest possible margin for labels, etc.
    axm_l = ax.TightInset(1);
    axm_b = ax.TightInset(2);
    axm_r = ax.TightInset(3);
    axm_t = ax.TightInset(4);
    % ax true tight width and height
    axt_w = axm_l + axm_r + axp_w;
    axt_h = axm_t + axm_b + axp_h;
    % true ratio
    ratio = axt_w / axt_h;
    
    % determine fig size in px
    fig_h = h;
    fig_w = fig_h * ratio;
    fig.Position = [0, 0, fig_w, fig_h];
    % resize the ax inside the figure tightly
    new_axp_l = axm_l;
    new_axp_b = axm_b;
    new_axp_w = fig_w - axm_l - axm_r;
    new_axp_h = fig_h - axm_t - axm_b;
    ax.Position = [new_axp_l, new_axp_b, new_axp_w, new_axp_h];  
end

