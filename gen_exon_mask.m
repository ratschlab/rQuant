function [exon_mask, exon_mask_part1, exon_mask_part2] = gen_exon_mask(profile_weights, tscp_len_bin, exon_feat, exon_feat_val, exon_feat_row, exon_feat_col)
% [exon_mask, exon_mask_part1, exon_mask_part2] = gen_exon_mask(profile_weights, tscp_len_bin, exon_feat, exon_feat_val, exon_feat_row, exon_feat_col)
%
% -- input --
% profile_weights: weights of profile functions
% tscp_len_bin: vector of length bins for each transcript
% exon_feat: vector of features for P exonic positions
% exon_feat_val: vector of indices to supporting points
% exon_feat_row: vector of row indices for sparse P x T matrix (position index)
% exon_feat_col: vector column indices for sparse P x T matrix (transcript index)
%
% -- output --
% exon_mask: P x T matrix of PLiF values (evaluated for each position and transcript)
% exon_mask_part1: P x T matrix of PLiF part corresponding to left supporting point
% exon_mask_part2: P x T matrix of PLiF part corresponding to right supporting point


P = size(exon_feat, 1);
T = size(exon_feat, 2);
profiles = sparse(profile_weights(:,tscp_len_bin));
if nargout==1
  % theta_f + (theta_f+1-theta_f)*feat(p)
  thf = sparse(exon_feat_row, exon_feat_col, profiles(exon_feat_val), P, T);
  thdiff = sparse(exon_feat_row, exon_feat_col, profiles(exon_feat_val+1)-profiles(exon_feat_val), P, T);
  exon_mask = thf + exon_feat.*thdiff;
  exon_mask_part1 = [];
  exon_mask_part2 = [];
else
  % theta_f*(1-feat(p)) + theta_f+1*feat(p)
  exon_mask_part1 = sparse(exon_feat_row, exon_feat_col, profiles(exon_feat_val), P, T) - exon_feat.*sparse(exon_feat_row, exon_feat_col, profiles(exon_feat_val), P, T);
  exon_mask_part2 = exon_feat.*sparse(exon_feat_row, exon_feat_col, profiles(exon_feat_val+1), P, T);
  exon_mask = exon_mask_part1 + exon_mask_part2;
end