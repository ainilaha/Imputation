#############################################################

### ENWAS:  Part 3  Quick top hits plot###

#############################################################

### Pull in results, analyze

#library(rlib23)
library(dplyr)
library(ggplot2)

## Color Blind Palette ##
# blacks, colors, grey



cbPal = c("#000000", "#E69F00", "#56B4E9", 
          "#009E73", "#F0E4420", "#0072B2",
          "#D55E00", "#CC79A7","#999999",
          "#009E73", "#F0E4420", "#0072B2")

## load enwas run file first
# plot.res.dir = paste(enwas_run$save.location,enwas_run$pop,"/all_results/Top_Hits/", sep = "")
# 
# ifelse(!dir.exists(file.path(plot.res.dir)), dir.create(file.path(plot.res.dir), recursive = TRUE), FALSE)

res.df = #enwas results file 
  #read.csv(paste(enwas_run$save.location,enwas_run$pop,'/all_results/enwas_results/enwas.results_' ,enwas_run$outcome.var,'.csv', sep = ''), stringsAsFactors = FALSE) %>% distinct(phenotype, .keep_all = TRUE)
#pheno.tags = enwas_run$tags   %>% distinct(pheno, .keep_all = TRUE)

ntags = ncol(pheno.tags)-1


res.df_filt = res.df # %>% left_join(pheno.tags, by = c('phenotype' = 'pheno'))
colnames(res.df_filt)[(ncol(res.df)+1):ncol(res.df_filt)] = paste('tag_',1:ntags, sep = "")

res.df_gen = res.df_filt %>% filter(!is.na(post_stratified_beta)) %>% 
  mutate(post_stratified_standard = post_stratified_beta*SD,
           y_low = post_stratified_standard - 1.96*san_se*SD,
         y_high = post_stratified_standard + 1.96*san_se*SD,
         N_ch = format(N, big.mark = ',', scientific = FALSE,trim = TRUE),
         PSB = paste(as.character(phenotype), " (", N_ch,")", sep = '')
  ) %>% arrange(abs_post_beta_SD) %>% mutate(PSB=factor(PSB, levels=PSB)) %>% 
  rename('Post-Stratified Beta (stnd)' = post_stratified_standard, SE = san_se) %>% 
 arrange(desc(abs_post_beta_SD)) 

## Dot plots
nvars = min(50, nrow(res.df_gen))

plot.df = res.df_gen %>% filter(!is.na(post_stratified_beta)) %>% slice(1:nvars) %>%
  mutate(col = case_when(
    (as.numeric(rownames(.))%%2)==0~1,
    TRUE~0),
    xmin = seq(0.5, nvars - 0.5, by = 1),
    xmax = seq(1.5,nvars + 0.5, by = 1)
  )

ss_leg = ifelse(length(unique(plot.df$sex_specific))==1, "none",'bottom')
#plot.df = res.df_gen  %>% slice(1:nvars)

ylims= c(round(min(plot.df$y_low, na.rm = TRUE),1)-0.1, 
         round(max(plot.df$y_high, na.rm = TRUE),1)+.1)

# col.tags = data.frame(coltag = viridis(n = length(unique(plot.df$tag_1))), tag_1 = c(sort(unique(plot.df$tag_1)))) 
# 
 plot.df$tag_1 = as.factor(plot.df$tag_1)
 plot.df$tag_1 = factor(plot.df$tag_1, levels = c('Exposure', 'Outcome', 'Other'), labels = c('Exposure', 'Outcome', 'Other'))
# # 
# # #col.tags = data.frame(coltag = viridis(n = length(unique(plot.df$tag_1))), tag_1 = c(sort(unique(plot.df$tag_1)))) 
# # 
# # 
# # color.tags = plot.df %>% left_join(col.tags, by = c('tag_1'))
# # ct = color.tags$coltag
library(viridis)

p1 = ggplot(plot.df, 
            aes(x = PSB, y = `Post-Stratified Beta (stnd)`, 
                ymin = y_low, ymax = y_high)) + theme_bw() +
  geom_point(size = 2, position = position_dodge(width = 1)) + 
  geom_errorbar(aes(ymin = y_low, ymax = y_high), 
                position = position_dodge(width = 1), width=0.5,cex=1) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, face = 'bold'),
        legend.position = 'bottom', legend.title = element_blank()) + xlab('') +
  scale_color_manual(values = rev(cbPal[1:3])) + 
  geom_rect(aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = +Inf,
    fill = factor(col)
  ),color = 'black',
  alpha = 0.2) +
  scale_fill_manual(values = c('white','grey20'), guide = 'none') +
   coord_flip(ylim = ylims) + 
  ylab("Standardized Effect Estimate") + 
  ggtitle(paste("Top ",nvars, " Hits: ",enwas_run$pop,"-",enwas_run$outcome.var, sep ="")) 


ggsave(filename = paste(plot.res.dir,"Top_",nvars,"_Hits_",
                        enwas_run$outcome.var,".pdf", sep = ''), width = 12, height = 10, units = 'in')


