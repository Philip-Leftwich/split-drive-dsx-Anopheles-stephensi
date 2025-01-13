# Load required scripts and functions ====
source("scripts/functions.R")       # Custom functions for analysis
source("scripts/custom_theme.R")    # Custom ggplot theme settings


# Data import====

data <- list.files(path="./data/", pattern = "*.xlsx", full.names = T) %>% 
  readxl::read_excel(., na = c("", "-"))  # read and input all files ending .xlsx

# Select relevant columns and process the data
data <- data %>% 
  dplyr::select(`Cross`, `Replicate`, `Egg number`, `Larvae number`) %>% 
  separate(`Cross`, into =c("male_parent", "female_parent"), # crosses set in reverse
           sep="x", extra = "merge", fill = "left", remove = FALSE) %>% 
  mutate(line = case_when(female_parent == "Het" | male_parent == "Het" ~ "Het",
                          female_parent == "Hom" | male_parent == "Hom" ~ "Hom",
                          female_parent == "WT" | male_parent == "WT" ~ "WT")) %>% 
  mutate(parent = if_else(female_parent == "SDA-500", "male", "female")) %>% 
  mutate(line = factor(line, levels = c("WT", "Het", "Hom"))) %>% 
  mutate(parent = factor(parent, levels = c("male", "female"))) %>% 
  mutate(id = row_number()) %>% 
  mutate(`Egg number` = if_else(id == 15, 188, `Egg number`)) # correction from Mireia's lab book


# Zero-inflated fecundity model====

library(pscl)

# Fit a zero-inflated Poisson model to the data
fit_zipoisson <- glmmTMB(`Egg number`~ line + parent+
                           +(1|Cross/Replicate),
                         data=data,
                         ziformula=~line + parent,
                         family=poisson)

# Estimated marginal means (EMMs) for the model
emmeans::emmeans(fit_zipoisson, ~ line, 
        type = "response", 
        component = "zi")

emmeans::emmeans(fit_zipoisson, ~ line + parent, 
                 type = "response", 
                 component = "zi")

emmeans::emmeans(fit_zipoisson, ~ line + parent, 
                 type = "response")

data %>% group_by(parent) %>% do(model = glmmTMB(`Egg number`~ line
                           +(1|Cross/Replicate),
                         data=.,
                         ziformula=~line,
                         family=poisson) %>%summary())

total_CI <- emmeans::emmeans(fit_zipoisson, specs =  ~ line + parent, type = "response") %>% 
  as_tibble()

summary_data <- data%>% 
  group_by(`line`, `parent`) %>% 
  summarise(n = n())

# Combine EMMs and summary data

total_CI <- full_join(total_CI, summary_data) %>% drop_na()%>% 
  mutate(parent = factor(parent, levels=c("female", "male"), labels=c("♀", "♂")))

## Figure====

custom_labeller <- labeller(
  Line = as_labeller(label_parsed),
  Cas9_grandparent = label_value, # keep the default labeller
  Cas9_parent = label_value,
  genotype = as_labeller(label_parsed)# keep the default labeller
)

plot_1 <- data %>% 
  mutate(parent = factor(parent, levels=c("female", "male"), labels=c("♀", "♂"))) %>%
  ggplot(aes(x=interaction(parent,line), y=`Egg number`, colour=parent, fill = parent))+
  geom_point(fill = "white", alpha = 0.4, shape = 21,
             position=position_nudge_any(x=0.12, y=0, ggbeeswarm::position_beeswarm(cex = .5)))+
  # geom_beeswarm(aes(size=win+loss), fill="white", alpha=0.6, shape=21)+
  geom_errorbar(data=total_CI, aes(min=(asymp.LCL), max=(asymp.UCL), y=rate, group=parent),width=0, linewidth=2,
                position=position_nudge(x=-0.12))+
  geom_point(data=total_CI, shape=21, aes(y=rate, group=NA), size=4, stroke=0.6,
             position=position_nudge(x=-0.12))+
  geom_label(data = total_CI,
             fill = "white",
             aes(label = paste0("n=",n)
                 ),
             y = 200,
             show.legend = FALSE) +
  scale_color_manual(values=c("#64310F", "#558BB8"))+
  scale_fill_manual(values = c("#EBB391", "#BFD3E4")) +
  scale_y_continuous(limits = c(0,210))+
  guides(colour="none")+
  labs(x="", 
       y="Number of eggs laid",
       shape="",
       fill = "gRNAs parent")+
  theme_custom()+
  theme(legend.position = "right",
        legend.box="vertical",
        axis.text.x = element_blank(),
        strip.placement = "outside",
        ggh4x.facet.nestline = element_line(colour = "grey60"),
#        plot.tag.position = c(.65,.15),
        plot.tag = element_text(size = rel(1)))+
  facet_nested_wrap(~ line + parent, 
                    nest_line = element_line(),
                    nrow=1,
                    strip.position="bottom",
                    scales="free_x",
                    labeller = custom_labeller)

## Table S4.====

tbl_regression(fit_zipoisson,
               intercept = TRUE, exponentiate = FALSE,
               #      label = list(`F1 cross`~ "Cross"),
               pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% 
  modify_header(label = "Coefficient", std.error = "SE",
                estimate = "Estimate", statistic = "z") %>% 
  modify_column_hide(conf.low) %>%  
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
  as_gt () %>% 
  gt::tab_source_note(gt::md("*This is a zero-inflated poisson mixed model*")) 

# Binomial model====

data2 <- data %>% 
  mutate(`Larvae number` = replace_na(`Larvae number`, 0)) %>%
  mutate(`Unhatched` = `Egg number` - `Larvae number`) %>% 
  filter(`Unhatched` >= 0) 


bin_mod <- glm(cbind(`Larvae number`, `Unhatched`) ~ line * parent, family = binomial, data = data2)

bin_mix_mod <- glmmTMB(cbind(`Larvae number`, `Unhatched`) ~ line * parent + (1|Cross/Replicate), family = binomial, data = data2)


total_CI <- emmeans::emmeans(bin_mix_mod, specs = ~ line + parent, type = "response") %>% 
  as_tibble()

summary_data <- data2 %>% 
  group_by(`line`, `parent`) %>% 
  summarise(n = n())

total_CI <- full_join(total_CI, summary_data) %>% drop_na() %>% 
  mutate(parent = factor(parent, levels=c("female", "male"), labels=c("♀", "♂")))

##Figure====
plot_2 <- data2 %>% 
  mutate(hatching = ((`Larvae number`/(`Larvae number`+`Unhatched`)*100))) %>% 
  mutate(parent = factor(parent, levels=c("female", "male"), labels=c("♀", "♂"))) %>% 
  ggplot(aes(x=interaction(parent,line), y= hatching, colour=parent, fill = parent))+
  geom_point(aes(size = `Larvae number`+`Unhatched`), fill = "white", alpha = 0.4, shape = 21,
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
       y="Percentage of eggs hatching",
       size="Number of offspring",
       shape="",
       fill = "gRNAs parent")+
  scale_y_continuous(limits=c(0,110),
                     labels=scales::percent_format(scale=1) # automatic percentages
  )+
  theme_custom()+
  theme(legend.position = "right",
        legend.box="vertical",
        axis.text.x = element_blank(),
        strip.placement = "outside",
        ggh4x.facet.nestline = element_line(colour = "grey60"),
#        plot.tag.position = c(.65,.15),
        plot.tag = element_text(size = rel(1)))+
  facet_nested_wrap(~ line + parent, 
                    nest_line = element_line(),
                    nrow=1,
                    strip.position="bottom",
                    scales="free_x",
                    labeller = custom_labeller)

#Count table====

Cross <- factor(c("SDA-500xWT", "SDA-500xHet", "SDA-500xHom",
                  "WTxSDA-500", "HetxSDA-500", "HomxSDA-500"))
Type <- rep(c("WT", "Het", "Hom"), times = 2)
Sex <- rep(c("Female", "Male"), each = 3)

Trials <- rep(50,6)
Engorge <- c(41,35,35,39,43,0)
Eggs <- c(33,24,11,30,20,NA)
Hatch <- c(29,22,8,29,6,NA)

table <- tibble(Cross, Type, Sex, Trials, Engorge, Eggs, Hatch)

table <- table %>% 
  dplyr::mutate(
    'Engorgement' = paste0(Engorge, "/", Trials),
    'Egg laying' = if_else(is.na(Eggs), "-", paste0(Eggs, "/", Engorge)),
    'Hatching' = if_else(is.na(Eggs), "-", paste0(Hatch, "/", Eggs))
  ) %>% 
  dplyr::select(Type, Sex, 'Engorgement':'Hatching') %>% 
  gt(groupname_col = 'Sex', rowname_col = 'Type') %>% 
  tab_style(
    style = cell_text(
      size = "smaller",
      weight = "bold",
      transform = "uppercase"
    ),
    locations =   cells_row_groups()
  ) %>% 
  tab_style(
    style = cell_fill(color = "#EBB391"),
    locations =   cells_row_groups(groups = "Female")
  )  %>% 
  tab_style(
    style = cell_fill(color = "#BFD3E4"),
    locations =   cells_row_groups(groups = "Male")
  ) %>% 
  opt_horizontal_padding(scale = 3) %>% 
  as_gtable()

# Figure 3.====

design <- "#######
           CCCC###
           #######
           AAA#BBB
           AAA#BBB
           AAA#BBB
           AAA#BBB"

plot_1 + plot_2 + table+ plot_layout(guides = "collect", design = design)  



## Table S5.==== 

tbl_regression(bin_mix_mod,
               intercept = TRUE, exponentiate = FALSE,
               pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% 
  modify_header(label = "Coefficient", std.error = "SE",
                estimate = "Estimate", statistic = "z") %>% 
  modify_column_hide(conf.low) %>%  
  bold_labels() %>% 
  bold_p(t = 0.05) %>% 
as_gt () %>% 
  gt::tab_source_note(gt::md("*This is a binomial mixed model*")) 

