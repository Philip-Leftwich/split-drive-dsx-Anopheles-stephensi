source("scripts/functions.R")
source("scripts/custom_theme.R")

### homing and mosaicis, data

# read and input all files ending .csv====
data <- list.files(path="./data/", pattern = "*.csv", full.names = T) %>% 
  read_csv(., na = c("", "-"))  


# Analysis for sgRNA ====

## Clean and tidy====
data_clean <- data %>% 
  fill(`F1 cross`) %>% 
  drop_na(`Female no.`) %>% 
  mutate(win = `A` + `AB`,
         loss = `B` + `WT`) %>% 
  separate(`F1 cross`, into =c("female_parent", "male_parent"), 
         sep="x", extra = "merge", fill = "left", remove = FALSE) %>% 
  mutate(Cas9_parent = if_else(female_parent == "SDA-500", "Male", "Female"),
         .after = `F1 cross`) %>% 
  mutate(Cas9_grandparent = if_else(str_detect(`F1 cross`, ":2288"), "Male", "Female"), 
         .before = Cas9_parent)






model <- glm(cbind(win,loss) ~ `F1 cross`, family = binomial, data = data_clean)

DHARMa_check(model)

# simulated residuals indicate overdispersion, correction with mixed model

model_mixed <- glmmTMB(cbind(win,loss) ~ `F1 cross` + (1|`F1 cross`/`Female no.`), family = binomial, data = data_clean)

DHARMa_check(model_mixed)

# produce model predictions table
total_CI <- emmeans::emmeans(model_mixed, specs = ~ `F1 cross`, type = "response") %>% 
  as_tibble() %>% 
  separate(`F1 cross`, into =c("female_parent", "male_parent"), 
           sep="x", extra = "merge", fill = "left", remove = FALSE) %>% 
  mutate(Cas9_parent = if_else(female_parent == "SDA-500", "Male", "Female"),
         .after = `F1 cross`) %>% 
  mutate(Cas9_grandparent = if_else(str_detect(`F1 cross`, ":2288"), "Male", "Female"), 
         .before = Cas9_parent) %>% 
  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("♀", "♂")))


summary_data <- data_clean %>% 
  group_by(`F1 cross`) %>% 
  filter(`No. of embryos` > 0) %>% 
  summarise(n = n())

total_CI <- full_join(total_CI, summary_data) %>% drop_na()

## Figure 4. ====

custom_labeller <- labeller(
  Line = as_labeller(label_parsed),
  Cas9_grandparent = label_value, # keep the default labeller
  Cas9_parent = label_value,
  genotype = as_labeller(label_parsed)# keep the default labeller
)

data_clean %>% 
  mutate(homing=((win/(win+loss)*100))) %>% 
  filter(`F1 cross` != "(2288:2073)xSDA-500") %>% 
  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  ggplot(aes(x=`F1 cross`, y=homing, colour=Cas9_parent, fill = Cas9_parent))+
  geom_point(aes(size = win+loss), fill = "white", alpha = 0.4, shape = 21,
             position=position_nudge_any(x=0.12, y=0, ggbeeswarm::position_beeswarm(cex = 1.8)))+
 # geom_beeswarm(aes(size=win+loss), fill="white", alpha=0.6, shape=21)+
  geom_errorbar(data=total_CI, aes(min=(asymp.LCL*100), max=(asymp.UCL*100), y=prob, group=NA),width=0, linewidth=2,
                position=position_nudge(x=-0.12))+
  geom_point(data=total_CI, shape=21, aes(y=prob*100, group=NA), size=4, stroke=0.6,
             position=position_nudge(x=-0.12))+
  geom_label(data = total_CI,
             fill = "white",
             aes(label = paste0("n=",n)),
             show.legend = FALSE,
             y = 105) +
  scale_size(range=c(0,4),
             breaks=c(50,100,150))+
  scale_color_manual(values=c("#64310F", "#558BB8"))+
  scale_fill_manual(values = c("#EBB391", "#BFD3E4")) +
  guides(colour=FALSE)+
  labs(x="", 
       y="Percentage of individuals scored",
       size="Number of offspring",
       shape="",
       fill = "Cas9 parent",
       tag = "Cas9 F1\n\n\nCas9 F0")+
  scale_y_continuous(limits=c(50,105),
                     labels=scales::percent_format(scale=1) # automatic percentages
  )+
  theme_custom()+
  theme(legend.position = "right",
        legend.box="vertical",
        axis.text.x = element_blank(),
        strip.placement = "outside",
        ggh4x.facet.nestline = element_line(colour = "grey60"),
        plot.tag.position = c(.65,.15),
        plot.tag = element_text(size = rel(1)))+
  facet_nested_wrap(~Cas9_grandparent+Cas9_parent, 
                    nest_line = element_line(),
                    nrow=1,
                    strip.position="bottom",
                    scales="free_x",
                    labeller = custom_labeller) # reduce line

## Table S11. ====

tbl_regression(model_mixed,
               intercept = TRUE, exponentiate = FALSE,
               label = list(`F1 cross`~ "Cross"),
               pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% 
  modify_header(label = "Coefficient", std.error = "SE",
                estimate = "Estimate", statistic = "z") %>% 
  modify_column_hide(conf.low) %>%  
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  as_gt () %>% 
  gt::tab_source_note(gt::md("*This is a binomial mixed model*")) 

ggsave("homing.png", dpi=900, width=5000, height=5000, units = "px")


# Supplementary analysis Cas9 inheritance ====

## Analysis====


data_clean2 <- data %>% 
  fill(`F1 cross`) %>% 
  drop_na(`Female no.`) %>% 
  mutate(cas9_win = `B` + `AB`,
         cas9_loss = `A` + `WT`) %>% 
  separate(`F1 cross`, into =c("female_parent", "male_parent"), 
           sep="x", extra = "merge", fill = "left", remove = FALSE) %>% 
  mutate(Cas9_parent = if_else(female_parent == "SDA-500", "Male", "Female"),
         .after = `F1 cross`) %>% 
  mutate(Cas9_grandparent = if_else(str_detect(`F1 cross`, ":2288"), "Male", "Female"), 
         .before = Cas9_parent)


model_mixed <- glmmTMB(cbind(cas9_win,cas9_loss) ~ `F1 cross` + (1|`Female no.`), family = binomial, data = data_clean2)

DHARMa_check(model_mixed)

total_CI <- emmeans::emmeans(model_mixed, specs = ~ `F1 cross`, type = "response") %>% 
  as_tibble() %>% 
  separate(`F1 cross`, into =c("female_parent", "male_parent"), 
           sep="x", extra = "merge", fill = "left", remove = FALSE) %>% 
  mutate(Cas9_parent = if_else(female_parent == "SDA-500", "Male", "Female"),
         .after = `F1 cross`) %>% 
  mutate(Cas9_grandparent = if_else(str_detect(`F1 cross`, ":2288"), "Male", "Female"), 
         .before = Cas9_parent) %>% 
  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("♀", "♂")))


summary_data <- data_clean %>% 
  group_by(`F1 cross`) %>% 
  filter(`No. of embryos` > 0) %>% 
  summarise(n = n())

total_CI <- full_join(total_CI, summary_data) %>% drop_na()

## Figure S5. ====

custom_labeller <- labeller(
  Line = as_labeller(label_parsed),
  Cas9_grandparent = label_value, # keep the default labeller
  Cas9_parent = label_value,
  genotype = as_labeller(label_parsed)# keep the default labeller
)

data_clean2 %>% 
  mutate(homing=((cas9_win/(cas9_win+cas9_loss)*100))) %>% 
  filter(`F1 cross` != "(2288:2073)xSDA-500") %>% 
  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  ggplot(aes(x=`F1 cross`, y=homing, colour=Cas9_parent, fill = Cas9_parent))+
  geom_point(aes(size = cas9_win+cas9_loss), fill = "white", alpha = 0.4, shape = 21, stroke =.6,
             position=position_nudge_any(x=0.12, y=0, ggbeeswarm::position_beeswarm(cex = 1.8)))+
  # geom_beeswarm(aes(size=win+loss), fill="white", alpha=0.6, shape=21)+
  geom_errorbar(data=total_CI, aes(min=(asymp.LCL*100), max=(asymp.UCL*100), y=prob, group=NA
                                   ),width=0, linewidth=2,
                position=position_nudge(x=-0.12))+
  geom_point(data=total_CI, shape=21, aes(y=prob*100, group=NA), size=4, stroke=0.6,
             position=position_nudge(x=-0.12))+
  geom_label(data = total_CI,
             fill = "white",
             aes(label = paste0("n=",n)),
             show.legend = FALSE,
             y = 105) +
  scale_size(range=c(0,4),
             breaks=c(50,100,150))+
  scale_color_manual(values=c("#64310F", "#558BB8"))+
  scale_fill_manual(values = c("#EBB391", "#BFD3E4")) +
  guides(colour=FALSE)+
  labs(x="", 
       y="Percentage of individuals scored",
       size="Number of offspring",
       shape="",
       fill = "Cas9 parent",
       tag = "Cas9 F1\n\n\nCas9 F0")+
  scale_y_continuous(limits=c(0,105),
                     labels=scales::percent_format(scale=1) # automatic percentages
  )+
  theme_custom()+
  theme(legend.position = "right",
        legend.box="vertical",
        axis.text.x = element_blank(),
        strip.placement = "outside",
        ggh4x.facet.nestline = element_line(colour = "grey60"),
        plot.tag.position = c(.65,.15),
        plot.tag = element_text(size = rel(1)))+
  facet_nested_wrap(~Cas9_grandparent+Cas9_parent, 
                    nest_line = element_line(),
                    nrow=1,
                    strip.position="bottom",
                    scales="free_x",
                    labeller = custom_labeller) # reduce line

ggsave("cas9.png", dpi=900, width=5000, height=5000, units = "px")




# Fecundity of crosses====

fit_zipoisson <- glmmTMB(`No. of embryos`~ Cas9_parent * Cas9_grandparent +
                           +(1|`F1 cross`/`Female no.`),
                         data=data_clean,
                         ziformula=~Cas9_parent * Cas9_grandparent,
                         family=poisson)

total_CI <- emmeans::emmeans(fit_zipoisson, specs =  ~ Cas9_parent * Cas9_grandparent, type = "response") %>% 
  as_tibble()

summary_data <- data_clean %>% 
  group_by(`Cas9_grandparent`, `Cas9_parent`) %>% 
  summarise(n = n())

total_CI <- full_join(total_CI, summary_data) %>% drop_na()%>% 
  filter(!(Cas9_grandparent == "Female" & Cas9_parent == "Female")) %>% 
  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>% 
  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("♀", "♂")))
  

custom_labeller <- labeller(
  Line = as_labeller(label_parsed),
  Cas9_grandparent = label_value, # keep the default labeller
  Cas9_parent = label_value,
  genotype = as_labeller(label_parsed)# keep the default labeller
)

plot_appendix_1 <- data_clean %>% 
  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  ggplot(aes(x=interaction(Cas9_parent, Cas9_grandparent), y=`No. of embryos`, colour=Cas9_parent, fill = Cas9_parent))+
  geom_point(aes(size = win+loss), fill = "white", alpha = 0.4, shape = 21,
             position=position_nudge_any(x=0.12, y=0, ggbeeswarm::position_beeswarm(cex = 1.8)))+
  # geom_beeswarm(aes(size=win+loss), fill="white", alpha=0.6, shape=21)+
  geom_errorbar(data=total_CI, aes(min=(asymp.LCL), max=(asymp.UCL), y=rate, group=NA),width=0, linewidth=2,
                position=position_nudge(x=-0.12))+
  geom_point(data=total_CI, shape=21, aes(y=rate, group=NA), size=4, stroke=0.6,
             position=position_nudge(x=-0.12))+
  geom_label(data = total_CI,
             fill = "white",
             aes(label = paste0("n=",n)),
             show.legend = FALSE,
             y = 184) +
  scale_size(range=c(0,4),
             breaks=c(50,100,150))+
  scale_color_manual(values=c("#64310F", "#558BB8"))+
  scale_fill_manual(values = c("#EBB391", "#BFD3E4")) +
  guides(colour=FALSE, size = FALSE)+
  labs(x="", 
       y="Percentage of individuals scored",
       size="Number of offspring",
       shape="",
       fill = "Cas9 parent",
       tag = "Cas9 F1\n\n\nCas9 F0")+
  theme_custom()+
  theme(legend.position = "right",
        legend.box="vertical",
        axis.text.x = element_blank(),
        strip.placement = "outside",
        ggh4x.facet.nestline = element_line(colour = "grey60"),
        plot.tag.position = c(1.2,.1),
        plot.tag = element_text(size = rel(1)))+
  facet_nested_wrap(~Cas9_grandparent+Cas9_parent, 
                    nest_line = element_line(),
                    nrow=1,
                    strip.position="bottom",
                    scales="free_x",
                    labeller = custom_labeller) # reduce line

# Fertility of crosses====

## Figure S4. ====

data2 <- data_clean %>% 
  mutate(`Unhatched` = `No. of embryos` - `Total`) %>% 
  mutate(id = row_number())


bin_mod <- glm(cbind(`Total`, `Unhatched`) ~ Cas9_parent * Cas9_grandparent, family = binomial, data = data2)

bin_mix_mod <- glmmTMB(cbind(`Total`, `Unhatched`) ~ Cas9_parent * Cas9_grandparent + (1|id), family = binomial, data = data2)



total_CI <- emmeans::emmeans(bin_mix_mod, specs = ~ Cas9_parent * Cas9_grandparent, type = "response") %>% 
  as_tibble()

summary_data <- data2 %>% 
  group_by(`Cas9_parent`, `Cas9_grandparent`) %>% 
  summarise(n = n())

total_CI <- full_join(total_CI, summary_data) %>% drop_na()%>% 
  filter(!(Cas9_grandparent == "Female" & Cas9_parent == "Female")) %>% 
  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>% 
  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("♀", "♂")))


plot_appendix_2 <- data2 %>% 
  mutate(hatching = ((`Total`/(`No. of embryos`)*100))) %>% 
  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("♀", "♂"))) %>%
  ggplot(aes(x=interaction(Cas9_parent, Cas9_grandparent), y=hatching, colour=Cas9_parent, fill = Cas9_parent)) +
  geom_point(aes(size = `No. of embryos`), fill = "white", alpha = 0.4, shape = 21,
             position=position_nudge_any(x=0.12, y=0, ggbeeswarm::position_beeswarm(cex = 1.8)))+
  # geom_beeswarm(aes(size=win+loss), fill="white", alpha=0.6, shape=21)+
  geom_errorbar(data=total_CI, aes(min=(asymp.LCL*100), max=(asymp.UCL*100), y=prob, group=NA),width=0, linewidth=2,
                position=position_nudge(x=-0.12))+
  geom_point(data=total_CI, shape=21, aes(y=prob*100, group=NA), size=4, stroke=0.6,
             position=position_nudge(x=-0.12))+
  geom_label(data = total_CI,
             fill = "white",
             aes(label = paste0("n=",n)),
             show.legend = FALSE,
             y = 105) +
  scale_size(range=c(0,4),
             breaks=c(50,100,150))+
  scale_color_manual(values=c("#64310F", "#558BB8"))+
  scale_fill_manual(values = c("#EBB391", "#BFD3E4")) +
  guides(colour=FALSE)+
  labs(x="", 
       y="Percentage of individuals scored",
       size="Number of offspring",
       shape="",
       fill = "Cas9 parent",
       tag = "Cas9 F1\n\n\nCas9 F0")+
  scale_y_continuous(limits=c(0,110),
                     labels=scales::percent_format(scale=1) # automatic percentages
  )+
  theme_custom()+
  theme(legend.position = "right",
        legend.box="vertical",
        axis.text.x = element_blank(),
        strip.placement = "outside",
        ggh4x.facet.nestline = element_line(colour = "grey60"),
        plot.tag.position = c(.75,.1),
        plot.tag = element_text(size = rel(1)))+
  facet_nested_wrap(~Cas9_grandparent+Cas9_parent, 
                    nest_line = element_line(),
                    nrow=1,
                    strip.position="bottom",
                    scales="free_x",
                    labeller = custom_labeller)

design <- "AA#BB
           AA#BB
           AA#BB"

plot_appendix_1 + plot_appendix_2 + plot_layout(guides = "collect", design = design)

ggsave("cas_fecundity.png", dpi = 600, width = 12, height = 9)



## Tables S6.& S9.====

egg_count <- tbl_regression(fit_zipoisson,
               intercept = TRUE, exponentiate = FALSE,
         #      label = list(`F1 cross`~ "Cross"),
               pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% 
  modify_header(label = "Coefficient", std.error = "SE",
                estimate = "Estimate", statistic = "z") %>% 
  modify_column_hide(conf.low) %>%  
  bold_labels() %>% 
  bold_p(t = 0.05) 


egg_hatch <- tbl_regression(bin_mix_mod,
               intercept = TRUE, exponentiate = FALSE,
                            pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% 
  modify_header(label = "Coefficient", std.error = "SE",
                estimate = "Estimate", statistic = "z") %>% 
  modify_column_hide(conf.low) %>%  
  bold_labels() %>% 
  bold_p(t = 0.05) 



tbl_merge(
  tbls = list(egg_count, egg_hatch),
  tab_spanner = c("**Egg count**", "**Egg hatching**")
) %>% 
   as_gt () %>% 
    gt::tab_source_note(gt::md("*Egg count is a zero-inflated poisson mixed model. Egg hatching is a binomial mixed model*")) 
