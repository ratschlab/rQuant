function plot_profiles(profile_weights, CFG)
% plot_profiles(fname)
%
% -- input --
% profile_weights: weights of profile functions
% CFG: configuration struct
%
% -- output --
% visualises profile weights

figure(); hold on;
col = colormap;
%col(find(col(:,3)==0 & col(:,1)>0.8 & col(:,1)>0.8),:) = []; % no yellow
col_bins = floor(size(col,1)/size(profile_weights,2));
col = col(1:col_bins:end,:);

lmt = get_limits(CFG.max_side_len, CFG.num_plifs/2);
max_len = [CFG.transcript_len_ranges(1:end-1,2)' CFG.max_side_len];

% positions that are ignored
min_gap = [0, CFG.num_plifs];
for t = 1:size(profile_weights,2),
  dist = max_len(t)*0.5;
  max_bin1 = find(lmt(1:end-1)<= dist & lmt(2:end)>dist);
  max_bin2 = CFG.num_plifs+1-find(lmt(1:end-1)<= dist & lmt(2:end)>dist);
  min_gap(1) = max(min_gap(1), max_bin1);
  min_gap(2) = min(min_gap(2), max_bin2);
end
offset = min_gap(2)-min_gap(1)-1 - 5;

% plot weights
for t = 1:size(profile_weights,2),
  dist = max_len(t)*0.5;
  max_bin1 = find(lmt(1:end-1)<= dist & lmt(2:end)>dist);
  ph(t) = plot([1:max_bin1], profile_weights(1:max_bin1,t), 'LineWidth', 1, 'Color', col(t,:));
  max_bin2 = CFG.num_plifs+1-find(lmt(1:end-1)<= dist & lmt(2:end)>dist);
  plot([max_bin2:CFG.num_plifs]-offset, profile_weights(max_bin2:end,t), 'LineWidth', 1, 'Color', col(t,:));
  plot([max_bin1:min_gap(1)], profile_weights(max_bin1:min_gap(1),t), ':', 'LineWidth', 1, 'Color', col(t,:));
  plot([min_gap(2):max_bin2]-offset, profile_weights(min_gap(2):max_bin2,t), ':', 'LineWidth', 1, 'Color', col(t,:));
  lg{t} = sprintf('%i<=len<=%i', CFG.transcript_len_ranges(t,1), CFG.transcript_len_ranges(t,2));
end

% axis label and limits
xlb = {};
for p = 1:min_gap(1),
  xlb{p} = sprintf('%i', round(lmt(p)));
end
for p = 1:CFG.num_plifs-min_gap(2)+1,
  xlb{CFG.num_plifs-p+1-min_gap(2)+min_gap(1)+1} = sprintf('-%i', round(lmt(p)));
end
plot_idx1 = 1:5:min_gap(1);
plot_idx2 = length(xlb):-5:min_gap(1)+1;
plot_idx2 = plot_idx2(end:-1:1);
xlim([0 length(xlb)+5]); set(gca, 'XTick', [plot_idx1, plot_idx2+5]);
xlb = {xlb{[plot_idx1, plot_idx2]}};
set(gca,'XTickLabel',xlb);
axis tight;

% labels
xlabel('Distance to transcript start/end');
ylabel('Profiles');
legend(ph, lg, 'Location', 'EastOutside');

% figure size
set(gcf, 'Position', [100 800 1000 700]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.25 2.5 13 5]);

print('-depsc2', '-r300', sprintf('%s/profiles.eps', CFG.out_dir));
