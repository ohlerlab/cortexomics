finalfolder="~/Dropbox/projects/cortexomics/manuscript/October"
figfolder="/fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/"
mkdir -p $finalfolder

# rsync -azhr --update $finalfolder Archive
rsync -r --update $finalfolder Archive
rm -r $finalfolder

mkdir -p "$finalfolder"/Figure_1/1_schematic

# Figure 1
mkdir -p "$finalfolder"/Figure_1
## Figure 1 1
rsync -u -r "Archive"/${finalfolder}/Figure_1/1_schematic/ ${finalfolder}/Figure_1/1_schematic/
## Figure 1 2

## Figure 1 3
mkdir -p "$finalfolder"/Figure_1/3_metagene_all
rsync -u "$figfolder"figure1/fig1c_myribowaltz_allsec_stageov.pdf "$finalfolder"/Figure_1/3_metagene_all
rsync -u "$figfolder"figure1/fig1c_myribowaltz_allsec_stagesep.pdf "$finalfolder"/Figure_1/3_metagene_all
rsync -u "$figfolder"figure1/fig1c_myribowaltz_aug_stageov.pdf "$finalfolder"/Figure_1/3_metagene_all
rsync -u "$figfolder"figure1/fig1c_myribowaltz_stop_stageov.pdf "$finalfolder"/Figure_1/3_metagene_all

## Figure 1 4
mkdir -p "$finalfolder"/Figure_1/4_differential_TE_scatter
rsync -u "$figfolder"figure1/foldchangecomp_limma.pdf "$finalfolder"/Figure_1/4_differential_TE_scatter
rsync -u "$figfolder"figure1/foldchangecomp_xtail.pdf "$finalfolder"/Figure_1/4_differential_TE_scatter
rsync -u "$figfolder"figure1/figures/figure1/foldchange_numbers.pdf "$finalfolder"/Figure_1/4_differential_TE_scatter
rsync -u "$figfolder"figures/figure1/foldchangecomp.pdf "$finalfolder"/Figure_1/4_differential_TE_scatter
rsync -u "$figfolder"figures/figure1/foldchange_numbers.pdf "$finalfolder"/Figure_1/4_differential_TE_scatter


##Figure 1 5
mkdir -p "$finalfolder"/Figure_1/5_GO_for_TE
rsync -u "$figfolder"figure1/ggoplot_MF_Fisher.elim_All_Translational_up_vs_all.pdf "$finalfolder"/Figure_1/5_GO_for_TE
rsync -u "$figfolder"figure1/ggoplot_MF_Fisher.elim_All_Translational_down_vs_all.pdf "$finalfolder"/Figure_1/5_GO_for_TE

## Figure 1 6
mkdir -p "$finalfolder"/Figure_1/6_Traj_for_RP_translation
echo "$figfolder"figure1/RP_traj.pdf "$finalfolder"/Figure_1/6_Traj_for_RP_translation
rsync -u "$figfolder"figure1/RP_traj.pdf "$finalfolder"/Figure_1/6_Traj_for_RP_translation

## Figure 1 - additional analysis
mkdir -p  "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/startlog2foldchange.pdf "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/end_up_vs_TEchange.pdf "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure1 "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure1/fig1c_myribowaltz_allsec_stageov.pdf "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure1/metaplottesribodipa_notup.pdf "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure1/metaplottesribodipaup.pdf "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/start_up_vs_TEchange.pdf "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/tables/startup_poseffect_go.tsv "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/tables/startdown_poseffect_go.tsv "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/tables/position_effect.tsv "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure1/ggoplot_BP_Fisher.elim_start_poseff_up_.pdf "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure1/ggoplot_BP_Fisher.elim_start_poseff_down_.pdf "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure1/ggoplot_CC_Fisher.elim_start_poseff_up_.pdf "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure1/ggoplot_CC_Fisher.elim_start_poseff_down_.pdf "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure1/ggoplot_MF_Fisher.elim_start_poseff_dup_.pdf "$finalfolder"/Figure_1/start_effect_analysis
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure1/ggoplot_MF_Fisher.elim_start_poseff_down_.pdf "$finalfolder"/Figure_1/start_effect_analysis

cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure1/fig1c_myribowaltz_frac_metaplots.pdf  ~/Dropbox/projects/cortexomics/manuscript/October/Figure_1/start_effect_analysis/

# Figure 2
mkdir -p "$finalfolder"/Figure_2
mkdir -p "$finalfolder"/Figure_2/1_UTR_reg
#1_te_utr_distplots.R
rsync -u "$figfolder"figure2/cumdist_TEchange_tplen.pdf "$finalfolder"/Figure_2/1_UTR_reg
rsync -u "$figfolder"figure2/cumdist_TEchange_fplen.pdf "$finalfolder"/Figure_2/1_UTR_reg
mkdir -p "$finalfolder"/Figure_S/horizontal_te_dist
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure2/TEchange_horizontaldist.pdf  "$finalfolder"/Figure_S/horizontal_te_dist/
mkdir -p "$finalfolder"/Figure_S/te_vs_txn_foldchange_density
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/Figure_S/TE_vs_txn_dens.pdf "$finalfolder"/Figure_S/te_vs_txn_foldchange_density
mkdir -p "$finalfolder"/Figure_S/te_vs_tpm_density
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/Figure_S/TPM_vs_TE_dens.pdf "$finalfolder"/Figure_S/te_vs_tpm_density
mkdir -p "$finalfolder"/Figure_S/neurite_foldchange_vs_start_eff_groups
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/zapp_violin_startgrps.pdf "$finalfolder"/Figure_S/neurite_foldchange_vs_start_eff_groups

#2_codon_profiles.R
mkdir -p "$finalfolder"/FIgure_2/2_A-site_var
rsync -u /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure2/trna_codons/fppos_vs_codon_variance.pdf "$finalfolder"/Figure_2/2_A-site_var

#/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure2/8_tRNA_array_analysis.R
mkdir -p "$finalfolder"/Figure_2/3_AA-site_var
rsync -u /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure2/trna_codons/stripplot_aa_codon.pdf "$finalfolder"/Figure_2/3_AA-site_var

mkdir -p "$finalfolder"/codon_gc_plots
cp /fast/AG_Ohler/dharnet/cortexomics/plots/codon_gc_dt_boxplots.pdf "$finalfolder"/codon_gc_plots

rsync -u /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure2/trna_codons/trna_sig_allfrac_mergedecod.pdf "$finalfolder"/Figure_2/4_tRNA
rsync -u /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure2/trna_codons/codon_usage_vs_summed_tRNAab.pdf "$finalfolder"/Figure_2/4_tRNA

mkdir -p "$finalfolder"/Figure_2/5_reg_cors
rsync -u  /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure2/trna_codons/codon_stat_grid_Total_E13-E145-E16-E175-P0.pdf "$finalfolder"/Figure_2/5_reg_cors

mkdir -p "$finalfolder"/Figure_3
mkdir -p "$finalfolder"/Figure_3/3_traj_meodeling
mkdir -p "$finalfolder"/Figure_3/1_isoform_dp_comps
mkdir -p "$finalfolder"/Figure_3/2_Corr_with_Mass_spec
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/cortilesribo.pdf  "$finalfolder"/Figure_3/2_Corr_with_Mass_spec


rm -rf "$finalfolder"/Figure_3/lfc_cors
mkdir -p "$finalfolder"/Figure_3/lfc_cors
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/lfc_cors_*.pdf "$finalfolder"/Figure_3/lfc_cors/
 

mkdir "$finalfolder"/Figure_3/cor_tiles
rsync -u  /fast/work/groups/ag_ohler/dharnet_m/cortexomics/'plots/cortilesribo.pdf' "$finalfolder"/Figure_3/cor_tiles


mkdir "$finalfolder"/Figure_3/Var_explained
rsync -u  /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure3/tp_var_explained.pdf "$finalfolder"/Figure_3/Var_explained
rsync -u /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure3/within_tp_var_explained.pdf $finalfolder/Figure_3/Var_explained


mkdir -p "$finalfolder"/Figure_4
mkdir -p "$finalfolder"/Figure_4/1_Clusters
cp /Users/dharnet/projects/cortexomics/plots/contrasts_heatmapContrasts_pca_t0_jribo_changegenes.pdf "$finalfolder"/Figure_4/1_Clusters
cp /Users/dharnet/projects/cortexomics/plots/clusters/Contrasts_pca_t0_jribo_changegenes_effect.pdf "$finalfolder"/Figure_4/1_Clusters

cp /Users/dharnet/projects/cortexomics/plots/contrasts_heatmapContrasts_allPca_t0_noribosplit..pdf "$finalfolder"/Figure_4/1_Clusters
cp /Users/dharnet/projects/cortexomics/plots/clusters/Contrasts_allPca_t0_noribo_effect.pdf "$finalfolder"/Figure_4/1_Clusters
mkdir -p "$finalfolder"/Figure_4/2_Cluster_GO
cp /Users/dharnet/projects/cortexomics/plots/cluster_go_bpContrasts_pca_t0_jribo_changegenes_effect.pdf "$finalfolder"/Figure_4/2_Cluster_GO
cp /Users/dharnet/projects/cortexomics/plots/cluster_go_bpContrasts_allPca_t0_noribo_effect.pdf "$finalfolder"/Figure_4/2_Cluster_GO
   mkdir -p "$finalfolder"/Figure_4/3_MuSiC
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/musicplots.pdf  "$finalfolder"/Figure_4/3_MuSiC
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/telley_MuSiC_bulk_both_telleycoregenes.pdf "$finalfolder"/Figure_4/3_MuSiC

mkdir -p "$finalfolder"/Figure_4/TE_telleycore_boxes/
cp /fast/AG_Ohler/dharnet/cortexomics/plots/telleyweight_vs_dTE.pdf "$finalfolder"/Figure_4/TE_telleycore_boxes/
cp /fast/AG_Ohler/dharnet/cortexomics/plots/telleyweight_Time_vs_dTE.pdf "$finalfolder"/Figure_4/TE_telleycore_boxes/

mkdir -p "$finalfolder"/Figure_4/4_Single_gene_examples
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure4/trajectory.pdf "$finalfolder"/Figure_4/4_Single_gene_examples
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure4/trajectory_te.pdf "$finalfolder"/Figure_4/4_Single_gene_examples
mkdir -p "$finalfolder"/Figure_S/riverplot
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figureS1/detection_riverplot_isdetected_.pdf "$finalfolder"/Figure_S/riverplot

mkdir -p "$finalfolder"/Figure_S/samplecor_grid
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/cortilesribo.pdf "$finalfolder"/Figure_S/samplecor_grid

mkdir -p "$finalfolder"/Figure_S/PCA_plots
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/hclust_pca_loadings.pdf "$finalfolder"/Figure_S/PCA_plots
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/hclust_pca_loadings_2.pdf "$finalfolder"/Figure_S/PCA_plots

mkdir -p "$finalfolder"/Figure_S/startcheck_plots
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/te_quart_metaplot.pdf "$finalfolder"/Figure_S/startcheck_plots
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/dtxn_quart_metaplot.pdf "$finalfolder"/Figure_S/startcheck_plots
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/mono_quart_metaplot.pdf "$finalfolder"/Figure_S/startcheck_plots
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/figures/figure1/metaplottesribodipa_notup_polymono.pdf "$finalfolder"/Figure_S/startcheck_plots
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/sigpepmetaplot.pdf "$finalfolder"/Figure_S/startcheck_plots
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/te_startup_density.pdf "$finalfolder"/Figure_S/startcheck_plots

mkdir -p "$finalfolder"/Figure_S/peptidecheckplots
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/PPE_motif_Psite_Alignment.pdf "$finalfolder"/Figure_S/peptidecheckplots
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/PPG_motif_Psite_Alignment.pdf "$finalfolder"/Figure_S/peptidecheckplots

find  $finalfolder -type f ! -iname '*DS_Store*' -exec realpath {} \;

find  $finalfolder -type f -mtime -1

mkdir -p "$finalfolder"/Figure_S/opt_vs_dt/
cp /fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/rfreq_vs_dwell_time.pdf "$finalfolder"/Figure_S/opt_vs_dt/


## Kinetics Redo
destdir=~/Dropbox/projects/cortexomics/manuscript/Draft_August/kinetics_redo
mkdir -p $destdir
cp /fast/AG_Ohler/dharnet/cortexomics/plots/kinetics_redo/becker_trajectory_classes.pdf $destdir
cp /fast/AG_Ohler/dharnet/cortexomics/plots/kinetics_redo/demoplot.pdf $destdir
cp /fast/AG_Ohler/dharnet/cortexomics/plots/kinetics_redo/NED_vs_notmsdev_barplots.pdf $destdir
cp /fast/AG_Ohler/dharnet/cortexomics/plots/kinetics_redo/trajectory_example_plots.pdf $destdir
cp /fast/AG_Ohler/dharnet/cortexomics/plots/kinetics_redo/hierarch_tefix_indiv_v_mcshane_pihalf.pdf $destdir
cp /fast/AG_Ohler/dharnet/cortexomics/plots/kinetics_redo/beck_traj_class_goplots.pdf $destdir
## Psites-redo
destdir=/home/zslastman/Dropbox/projects/cortexomics/manuscript/Draft_August/p_a_redux
mkdir -p  $destdir
#
cp   /fast/AG_Ohler/dharnet/cortexomics/plots/p_a_redux/a_occ_ab_vs_dt.pdf $destdir
cp   /fast/AG_Ohler/dharnet/cortexomics/plots/p_a_redux/a_occ_av_vs_dt.pdf $destdir
cp  /fast/AG_Ohler/dharnet/cortexomics/plots/p_a_redux/p_occ_ab_vs_dt.pdf $destdir
cp /fast/AG_Ohler/dharnet/cortexomics/plots/p_a_redux/p_occ_av_vs_dt.pdf $destdir
cp /fast/AG_Ohler/dharnet/cortexomics/plots/a_occ_stage_dtdist.pdf $destdir
cp /fast/AG_Ohler/dharnet/cortexomics/plots/p_occ_stage_dtdist.pdf $destdir
cp /fast/AG_Ohler/dharnet/cortexomics/plots/atimes_adjecency_vs_psite_occ.pdf $destdir
cp  /fast/AG_Ohler/dharnet/cortexomics/plots/p_a_redux/abundance_freq.pdf $destdir 
cp   /fast/AG_Ohler/dharnet/cortexomics/plots/p_a_redux/ab_v_wusage.pdf $destdir
cp  /fast/AG_Ohler/dharnet/cortexomics/plots/p_a_redux/p_site_occ_vs_poly_av_vs_dt.pdf $destdir
cp /fast/AG_Ohler/dharnet/cortexomics/plots/kinetics_redo/hierarch_tefix_indiv_v_mcshane_pihalf.pdf $destdir

cp /fast/AG_Ohler/dharnet/cortexomics/plots/rust_fppos_vs_codon_variance.pdf $destdir

cp /fast/AG_Ohler/dharnet/cortexomics/plots/Figures/Figure4/lfc_cors_TE.pdf /home/zslastman/Dropbox/projects/cortexomics/manuscript/Draft_August/
cp /fast/AG_Ohler/dharnet/cortexomics/plots/Figures/Figure4/lfc_cors_TEvmrna.pdf /home/zslastman/Dropbox/projects/cortexomics/manuscript/Draft_August/
cp /fast/AG_Ohler/dharnet/cortexomics/plots/Figures/Figure4/lfc_cors_all.pdf /home/zslastman/Dropbox/projects/cortexomics/manuscript/Draft_August/
cp /fast/AG_Ohler/dharnet/cortexomics/plots/Figures/Figure4/lfc_cors_ribo.pdf /home/zslastman/Dropbox/projects/cortexomics/manuscript/Draft_August/
cp /fast/AG_Ohler/dharnet/cortexomics/plots/Figures/Figure4/lfc_cors_TEv_msR.pdf /home/zslastman/Dropbox/projects/cortexomics/manuscript/Draft_August/

cp /fast/AG_Ohler/dharnet/cortexomics/plots/broad_dt_aa_stripplot.pdf /home/zslastman/Dropbox/projects/cortexomics/manuscript/Draft_August/
cp /fast/AG_Ohler/dharnet/cortexomics/plots/kinetics_redo/NED_vs_notmsdev_barplots.pdf /home/zslastman/Dropbox/projects/cortexomics/manuscript/Draft_August/
cp //fast/AG_Ohler/dharnet/cortexomics/plots/atimes_adjecency_vs_psite_occ.pdf ~/Dropbox/projects/cortexomics/manuscript/Draft_August/ 
