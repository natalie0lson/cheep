today <- Sys.Date()
#remotes::install_github("Russel88/MicEco")
#remotes::install_github("vmikk/metagMisc")

pacman::p_load(phyloseq, vegan, ggplot2, ggpubr, microViz, corncob, ggraph, DT, readr, dplyr, tidyverse, gridExtra, microbiome, metagMisc, scales, GUniFrac, kableExtra, optimx, edgeR, xaringan, xaringanExtra, magick, MicEco)
physeq <- readRDS("Data/phyloseq_obj_2024-04-07.rds")
pathseq <- readRDS("Data/pathseq_obj_2024-02-27.rds")
ARGseq <- readRDS("Data/ARGseq_obj_2024-04-22.rds")
clinseq <- readRDS("Data/clinseq_obj_2024-04-18.rds")
phy_meta <- read.csv("Data/phy_meta_2024-01-30.csv")
pathogens <- read.csv("Data/pathogens_2024-01-30.csv")
ARG <- read.csv("Data/ARG_2024-01-30.csv")

## Rename "Indigenous" to "Local" 
physeq@sam_data$breed[physeq@sam_data$breed=="Indigenous"] <- "Local"
pathseq@sam_data$breed[pathseq@sam_data$breed=="Indigenous"] <- "Local"
ARGseq@sam_data$breed[ARGseq@sam_data$breed=="Indigenous"] <- "Local"
clinseq@sam_data$breed[clinseq@sam_data$breed=="Indigenous"] <- "Local"

levels(physeq@sam_data$loc.type.n) <- c(levels(physeq@sam_data$loc.type.n), "Local/Household/Fecal, N=7", "Local/Market/Fecal, N=7", "Local/Market/Carcass, N=7")

physeq@sam_data$loc.type.n[(physeq@sam_data$breed=="Local") & (physeq@sam_data$Group=="Production")] <- "Local/Household/Fecal, N=7"
physeq@sam_data$loc.type.n[(physeq@sam_data$breed=="Local")& (physeq@sam_data$Group=="Market")] <- "Local/Market/Fecal, N=7"
physeq@sam_data$loc.type.n[(physeq@sam_data$breed=="Local") & (physeq@sam_data$Group=="Processing")] <- "Local/Market/Carcass, N=7"
table(physeq@sam_data$loc.type.n)

pathseq@sam_data$loc.type.n[pathseq@sam_data$loc.type.n=="Indigenous/Household/Fecal, N=7"] <- "Local/Household/Fecal, N=7"
pathseq@sam_data$loc.type.n[pathseq@sam_data$loc.type.n=="Indigenous/Market/Fecal, N=7"] <- "Local/Market/Fecal, N=7"
pathseq@sam_data$loc.type.n[pathseq@sam_data$loc.type.n=="Indigenous/Market/Carcass, N=7"] <- "Local/Market/Carcass, N=7"

ARGseq@sam_data$loc.type.n[ARGseq@sam_data$loc.type.n=="Indigenous/Household/Fecal, N=7"] <- "Local/Household/Fecal, N=7"
ARGseq@sam_data$loc.type.n[ARGseq@sam_data$loc.type.n=="Indigenous/Market/Fecal, N=7"] <- "Local/Market/Fecal, N=7"
ARGseq@sam_data$loc.type.n[ARGseq@sam_data$loc.type.n=="Indigenous/Market/Carcass, N=7"] <- "Local/Market/Carcass, N=7"

clinseq@sam_data$loc.type.n[clinseq@sam_data$loc.type.n=="Indigenous/Household/Fecal, N=7"] <- "Local/Household/Fecal, N=7"
clinseq@sam_data$loc.type.n[clinseq@sam_data$loc.type.n=="Indigenous/Market/Fecal, N=7"] <- "Local/Market/Fecal, N=7"
clinseq@sam_data$loc.type.n[clinseq@sam_data$loc.type.n=="Indigenous/Market/Carcass, N=7"] <- "Local/Market/Carcass, N=7"

physeq_fixed <- physeq %>%
  tax_fix(
    min_length = 2,
    unknowns = c("sp. ANNP30", "Ponticoccus"),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

# Custom colors for plotting
custom_colors <- c(
  "Commercial/Farm/Fecal, N=4" = "darkolivegreen1",
  "Commercial/Market/Fecal, N=8" = "darkolivegreen3",
  "Commercial/Market/Carcass, N=42" = "darkolivegreen4",
  "Local/Household/Fecal, N=7" = "lightblue1",
  "Local/Market/Fecal, N=7" = "lightblue3",
  "Local/Market/Carcass, N=7" = "lightblue4",
  "Market Rinse Water, N=26" = "goldenrod"
)


## Update formatting on loc.type.n variable
physeq_fixed@sam_data$loc.type.n <- gsub(", ", "\n", physeq_fixed@sam_data$loc.type.n)
physeq_fixed@sam_data$loc.type.n <- gsub("/Farm/Fecal", " ", physeq_fixed@sam_data$loc.type.n)
physeq_fixed@sam_data$loc.type.n <- gsub("/Household/Fecal", " ", physeq_fixed@sam_data$loc.type.n)

pathseq@sam_data$loc.type.n <- gsub(", ", "\n", pathseq@sam_data$loc.type.n)
pathseq@sam_data$loc.type.n <- gsub("/Farm/Fecal", " ", pathseq@sam_data$loc.type.n)
pathseq@sam_data$loc.type.n <- gsub("/Household/Fecal", " ", pathseq@sam_data$loc.type.n)

ARGseq@sam_data$loc.type.n <- gsub(", ", "\n", ARGseq@sam_data$loc.type.n)
ARGseq@sam_data$loc.type.n <- gsub("/Farm/Fecal", " ", ARGseq@sam_data$loc.type.n)
ARGseq@sam_data$loc.type.n <- gsub("/Household/Fecal", " ", ARGseq@sam_data$loc.type.n)

clinseq@sam_data$loc.type.n <- gsub(", ", "\n", clinseq@sam_data$loc.type.n)
clinseq@sam_data$loc.type.n <- gsub("/Farm/Fecal", " ", clinseq@sam_data$loc.type.n)
clinseq@sam_data$loc.type.n <- gsub("/Household/Fecal", " ", clinseq@sam_data$loc.type.n)

#update names of custom colors 
names(custom_colors) <- gsub(", ", "\n", names(custom_colors))
names(custom_colors) <- gsub("/Farm/Fecal", " ", names(custom_colors))
names(custom_colors) <- gsub("/Household/Fecal", " ", names(custom_colors))

pathseq_fec <- subset_samples(pathseq, Feces == 1)
pathseq_fec_prod <- subset_samples(pathseq_fec, Group=="Production")
#install.packages("ggside")
library(ggside)
physeq_fec <- subset_samples(physeq_fixed, Feces == 1)
physeq_fec_prod <- subset_samples(physeq_fec, Group=="Production")

table(physeq_fec_prod@sam_data$breed)

# Transform phyloseq object
physeq.rel <- microbiome::transform(physeq_fec_prod, "compositional")
otu <- abundances(physeq.rel)
meta <- meta(physeq.rel)
# Significance Test (PERMANOVA)
formula_str <- paste("t(otu) ~", "breed")  # Create the formula string
permanova <- adonis2(as.formula(formula_str), data = meta, permutations = 999999, method = "bray")
p_value <- sprintf("%.3f", permanova$`Pr(>F)`[1])  # Format p-value to four decimal places
p_value
bac_nmds <- physeq_fec_prod %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc("NMDS") %>%
  ord_plot(color = "loc.type.n", shape = "loc.type.n", size = 16) +
  stat_ellipse(aes(group = loc.type.n)) + ## add ellipses
  scale_color_manual(values=custom_colors)+
  theme_bw() + ggtitle(paste0("Bacterial Composition \n P-value =", p_value)) + 
  ggside::geom_xsidedensity(aes(fill = loc.type.n), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = loc.type.n), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void() + scale_fill_manual(values=custom_colors) + theme(
    text = element_text(size = 50),           # General text size
    legend.title = element_blank(),           # No legend title
    plot.title = element_text(size = 50),     # Plot title size
    axis.title = element_text(size = 50),     # Axis titles size
    axis.text = element_text(size = 50),      # Axis text size
    legend.text = element_text(size = 50)     # Legend text size
  )
bac_nmds

# Transform phyloseq object
physeq.rel <- microbiome::transform(pathseq_fec_prod, "compositional")
otu <- abundances(physeq.rel)
meta <- meta(physeq.rel)

# Significance Test (PERMANOVA)
formula_str <- paste("t(otu) ~", "breed")  # Create the formula string
permanova <- adonis2(as.formula(formula_str), data = meta, permutations = 999999, method = "bray")
p_value <- sprintf("%.3f", permanova$`Pr(>F)`[1])  # Format p-value to four decimal places
p_value
path_nmds <- pathseq_fec_prod %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc("NMDS") %>%
  ord_plot(color = "loc.type.n", shape = "loc.type.n", size = 16) +
  stat_ellipse(aes(group = loc.type.n)) + ## add ellipses
  scale_color_manual(values=custom_colors)+
  theme_bw() + ggtitle(paste0("Pathogen Composition \n P-value =", p_value)) + 
  ggside::geom_xsidedensity(aes(fill = loc.type.n), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = loc.type.n), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void() + scale_fill_manual(values=custom_colors) + theme(
      text = element_text(size = 50),           # General text size
      legend.title = element_blank(),           # No legend title
      plot.title = element_text(size = 50),     # Plot title size
      axis.title = element_text(size = 50),     # Axis titles size
      axis.text = element_text(size = 50),      # Axis text size
      legend.text = element_text(size = 50)     # Legend text size
    )
path_nmds

ARGseq_fec <- subset_samples(ARGseq, Feces == 1)
ARGseq_fec_prod <- subset_samples(ARGseq_fec, Group=="Production" )
# Transform phyloseq object
physeq.rel <- microbiome::transform(ARGseq_fec, "compositional")
otu <- abundances(physeq.rel)
meta <- meta(physeq.rel)

# Significance Test (PERMANOVA)
formula_str <- paste("t(otu) ~", "breed")  # Create the formula string
permanova <- adonis2(as.formula(formula_str), data = meta, permutations = 999999, method = "bray")
p_value <- sprintf("%.5f", permanova$`Pr(>F)`[1])  # Format p-value to four decimal places
p_value
arg_nmds <- ARGseq_fec_prod %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc("NMDS") %>%
  ord_plot(color = "loc.type.n", shape = "loc.type.n", size = 16) +
  stat_ellipse(aes(group = loc.type.n)) + ## add ellipses
  scale_color_manual(values=custom_colors)+
  theme_bw() + ggtitle(paste0("ARG Composition \n P-value =", p_value)) + 
  ggside::geom_xsidedensity(aes(fill = loc.type.n), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = loc.type.n), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void() + scale_fill_manual(values=custom_colors) + theme(
    text = element_text(size = 50),           # General text size
    legend.title = element_blank(),           # No legend title
    plot.title = element_text(size = 50),     # Plot title size
    axis.title = element_text(size = 50),     # Axis titles size
    axis.text = element_text(size = 50),      # Axis text size
    legend.text = element_text(size = 50)     # Legend text size
  )
arg_nmds

## C-args
clinseq_fec <- subset_samples(clinseq, Feces == 1)
clinseq_fec_prod <- subset_samples(clinseq_fec, Group=="Production" )# Transform phyloseq object
physeq.rel <- microbiome::transform(clinseq_fec_prod, "compositional")
otu <- abundances(physeq.rel)
meta <- meta(physeq.rel)

# Significance Test (PERMANOVA)
formula_str <- paste("t(otu) ~", "breed")  # Create the formula string
permanova <- adonis2(as.formula(formula_str), data = meta, permutations = 999999, method = "bray")
p_value <- sprintf("%.3f", permanova$`Pr(>F)`[1])  # Format p-value to four decimal places
p_value

hr_arg_nmds <- clinseq_fec_prod %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc("NMDS") %>%
  ord_plot(color = "loc.type.n", shape = "loc.type.n", size = 16) +
  stat_ellipse(aes(group = loc.type.n)) + ## add ellipses
  scale_color_manual(values=custom_colors)+
  theme_bw() + ggtitle(paste0("HR-ARG Composition \n P-value =", p_value)) + 
  ggside::geom_xsidedensity(aes(fill = loc.type.n), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = loc.type.n), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void() + scale_fill_manual(values=custom_colors) + theme(
    text = element_text(size = 50),           # General text size
    legend.title = element_blank(),           # No legend title
    plot.title = element_text(size = 50),     # Plot title size
    axis.title = element_text(size = 50),     # Axis titles size
    axis.text = element_text(size = 50),      # Axis text size
    legend.text = element_text(size = 50)     # Legend text size
  )
hr_arg_nmds

# Save the multi-panel figure to a PDF file
#ggexport(multi_panel_plot, filename = "figures/path_nmds_plots.jpg", height = 15, width = 25)

#######################################
#### Differential abundance 
##########################################
source("R/daa_function.R")

bac_prod <- subset_samples(physeq, Group=="Production") %>% tax_glom("Phylum")

#DIFFERENTIAL ABUNDANCE BETWEEN BREEDS (AMONG PRODUCTION SITE)
daa_prod <- differentialTest(formula = ~ breed,
                             phi.formula = ~ breed,
                             formula_null = ~ 1,
                             phi.formula_null = ~ breed,
                             test = "Wald", boot = FALSE,
                             full_output=FALSE,
                             data = bac_prod,
                             fdr_cutoff = 0.05)

bac_daa <- plot.differentialTest(daa_prod, level= "Phylum")
bac_daa

path_prod <- subset_samples(pathseq, Group=="Production") %>% tax_glom("Genus")

#DIFFERENTIAL ABUNDANCE BETWEEN BREEDS (AMONG PRODUCTION SITE)
daa_prod <- differentialTest(formula = ~ breed,
                             phi.formula = ~ breed,
                             formula_null = ~ 1,
                             phi.formula_null = ~ breed,
                             test = "Wald", boot = FALSE,
                             full_output=FALSE,
                             data = path_prod,
                             fdr_cutoff = 0.05)

path_daa <- plot.differentialTest(daa_prod, level= "Genus")
path_daa
### DAA ARGS

argprod <- subset_samples(ARGseq, Group=="Production")

daa_argseq <- differentialTest(formula = ~ breed,
                               phi.formula = ~ breed,
                               formula_null = ~ 1,
                               phi.formula_null = ~ breed,
                               test = "Wald", boot = FALSE,
                               full_output=FALSE,
                               data = argprod,
                               fdr_cutoff = 0.05)

arg_daa<- plot.differentialTest(daa_argseq, level= "Genus")
arg_daa

### DAA c-args 

clinprod <- subset_samples(clinseq, Group=="Production")

daa_clinseq <- differentialTest(formula = ~ breed,
                               phi.formula = ~ breed,
                               formula_null = ~ 1,
                               phi.formula_null = ~ breed,
                               test = "Wald", boot = FALSE,
                               full_output=FALSE,
                               data = clinprod,
                               fdr_cutoff = 0.05)

clin_daa <- plot.differentialTest(daa_clinseq, level= "Genus")
clin_daa

### Venn Diagram

argprod <- subset_samples(ARGseq, Group=="Production")
argvenn <- ps_venn(
  argprod,
  group = "breed",
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE, 
  labels=list(cex=5),
 #main = "Production ARGs",
  
  quantities = list(cex=5, type = c('counts',"percent")),
  
  fill=c("Commercial" = "darkolivegreen1", 
         "Local" = "lightblue1")
  
)
argvenn

### Venn Diagram
clinprod <- subset_samples(clinseq, Group=="Production")
clinvenn <- ps_venn(
  clinprod,
  group = "breed",
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE, 
  labels=list(cex=5),
 #main = "Production HR-ARGs",
  
  quantities = list(cex=5, type = c('counts',"percent")),
  
  fill=c("Commercial" = "darkolivegreen1", 
         "Local" = "lightblue1")
  
)
clinvenn
60+51+18
51/129
60/129
18/129

bac_prod <- subset_samples(physeq, Group=="Production")

## Venn diagram 
bacvenn <- ps_venn(
  bac_prod,
  group = "breed",
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE,
  labels=list(cex=5),
 # main = "Production Bacteria",
  
  quantities = list(cex=5, type = c('counts',"percent")),
  
  fill=c("Commercial" = "darkolivegreen1", 
         "Local" = "lightblue1")
)
bacvenn

path_prod <- subset_samples(pathseq, Group=="Production")

## Venn diagram 
pathvenn <- ps_venn(
  path_prod,
  group = "breed",
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE,
  labels=list(cex=5),
  #main = "Production Pathogens",
  
  quantities = list(cex=5, type = c('counts',"percent")),
  
  fill=c("Commercial" = "darkolivegreen1", 
         "Local" = "lightblue1")
)
pathvenn

### Most Abundant 

phy_df <- phyloseq_to_df(physeq_fec_prod, addtax=T)

head(phy_df)

df_long <- phy_df%>% 
  gather(key="ID", value="count", -c("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
head(df_long)
# Genus & species name
df_long$g.species <- paste0(substr(df_long$Genus, 1, 1), ".", " ", df_long$Species)

df_stack <- merge(df_long, phy_meta, by="ID")
phyla = aggregate(df_stack$count,
                  by = list(df_stack$Phylum),
                  FUN = sum)
phyla$num <- phyla$x
phyla$Phylum <- phyla$Group.1
phyla <- phyla[with(phyla,order(-num)),]

phyla <-  na.omit(phyla[-7, ][1:5, ])

phylist <- (phyla$Phylum)
phylist

df_stack$phystack <- NA

phydat1 <- df_stack[df_stack$Phylum %in% phylist,]
phydat1$phystack <- phydat1$Phylum

phydat2 <- df_stack[!df_stack$Phylum %in% phylist,]
phydat2$phystack <- "Other"

phydat <- rbind(phydat1, phydat2)
# Reorder the 'Category' variable based on the 'Value' variable
phydat$phystack <- factor(phydat$phystack, levels=c( 'Other', 'Fusobacteria',  'Bacteroidetes'  ,    'Firmicutes' , "Actinobacteria",  'Proteobacteria'  ))

# Order the levels consistently across all data frames
#phydat$Group <- factor(phydat$Group, levels = desired_order)
# Define the desired order of x-axis labels
unique_categories <- unique(phy_df$loc.type.n)

# Define the desired order of the levels for the "Group" factor
desired_order <- c("Production", "Market", "Processing")

# Order the levels consistently across all data frames
phydat$Group <- factor(phydat$Group, levels = desired_order)
#stacked bar chart
bac_stack <- phydat %>% 
  ggplot() + 
  geom_bar(aes(x=labels, y=count, fill=phystack, color=phystack), stat="identity", position="fill") + 
  ggtitle("Proportional Abundance of Bacterial Phyla") + 
  ylab("Proportional Abundance") + theme_classic() +  theme(#text = element_text(size = 20),
    legend.title=element_blank(), text=element_text(size=50)) + 
  xlab("") + scale_y_continuous(labels = scales::percent_format(scale = 100),  # Format y-axis as percentages
                                limits = c(0, 1)) 
bac_stack

## Path Stack 

phy_df <- phyloseq_to_df(pathseq_fec_prod, addtax=T)

df_long <- phy_df%>% 
  gather(key="ID", value="count", -c("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
head(df_long)
# Genus & species name
df_long$g.species <- paste0(substr(df_long$Genus, 1, 1), ".", " ", df_long$Species)

df_stack <- merge(df_long, phy_meta, by="ID")
phyla = aggregate(df_stack$count,
                  by = list(df_stack$g.species),
                  FUN = sum)
phyla$num <- phyla$x
phyla$g.species <- phyla$Group.1
phyla <- phyla[with(phyla,order(-num)),]

phyla <-  na.omit(phyla[1:5, ])

phylist <- (phyla$g.species)
phylist
df_stack$g.species

df_stack$phystack <- NA

phydat1 <- df_stack[df_stack$g.species %in% phylist,]
phydat1$phystack <- phydat1$g.species

phydat2 <- df_stack[!df_stack$g.species %in% phylist,]
phydat2$phystack <- "Other"

phydat <- rbind(phydat1, phydat2)
# Reorder the 'Category' variable based on the 'Value' variable
phydat$phystack <- factor(phydat$phystack, levels=c( "Other", "E. coli"     ,  "E. faecalis"  , "S. flexneri",   "K. pneumoniae" ,"S. enterica"     ))

# Define the desired order of x-axis labels
unique_categories <- unique(phy_df$loc.type.n)

# Check for missing values
sum(is.na(phydat$labels))
sum(is.na(phydat$count))
sum(is.na(phydat$phystack))

# Check data range
summary(phydat$count)

#stacked bar chart

phydat$phystack
path_stack <- phydat %>% 
  ggplot() + 
  geom_bar(aes(x=labels, y=count, fill=phystack, color=phystack), stat="identity", position="fill") + 
  ggtitle("Proportional Abundance of Pathogens") + 
  ylab("Proportional Abundance") + theme_classic() +  theme(#text = element_text(size = 20),
    legend.title=element_blank(), text=element_text(size=50)) + 
  xlab("") + scale_y_continuous(labels = scales::percent_format(scale = 100)) 
path_stack

### ARGs

## Most abundant ARGs
prod_args <- subset(ARG, ARG$Group=="Production")
top.args = aggregate(prod_args$TotalIDs,
                     by = list(prod_args$name2),
                     FUN = sum)
top.args$num <- top.args$x
top.args$name2 <- top.args$Group.1
top.args <- top.args[with(top.args,order(-num)),]
top.args <-  na.omit(top.args[1:5, ])

arglist <- (top.args$name2)
arglist
table(prod_args$clinical)
length(prod_args$name2)

ARG$argstack <- NA 
argdat1 <- prod_args[prod_args$name2 %in% arglist,]

argdat1 <- prod_args[prod_args$name2 %in% arglist,]
argdat1$argstack <- argdat1$name2

argdat2 <- prod_args[!prod_args$name2 %in% arglist,]
argdat2$argstack <- "Other"

argdat <- rbind(argdat1, argdat2)

# Reorder the 'Category' variable based on the 'Value' variable

argdat$argstack <- factor(argdat$argstack, levels=c( "Other","tet(Q)" ,   "tet(W)"  ,  "ant(6)-Ia" ,"erm(F)"  ,  "mef(A)"   ))

# Order the levels consistently across all data frames

arg_stack <-  argdat %>% ggplot(aes(x=labels,y=TotalIDs)) + 
  geom_bar(aes(fill=argstack, color=argstack), stat="identity", position="fill") + 
  ggtitle("Proportional Abundance of ARGs") + 
  ylab("Proportional Abundance") +xlab("")+ theme_classic() + 
  theme( legend.title=element_blank(), text=element_text(size=50))

arg_stack
#ggsave(stacked_bar, file=paste0(fig_dir, "AMR_stacked_bar_", today,".png"), 
#   width = 44.45, height = 27.78, units = "cm", dpi=300)
### hr-ARGs

## Most abundant ARGs
HR_ARG <- subset(prod_args, prod_args$clinical==1)
top.args = aggregate(HR_ARG$TotalIDs,
                     by = list(HR_ARG$name2),
                     FUN = sum)
top.args$num <- top.args$x
top.args$name2 <- top.args$Group.1
top.args <- top.args[with(top.args,order(-num)),]
top.args <-  na.omit(top.args[1:5, ])

arglist <- (top.args$name2)
arglist

ARG$argstack <- NA 
argdat1 <- HR_ARG[HR_ARG$name2 %in% arglist,]

argdat1 <- HR_ARG[HR_ARG$name2 %in% arglist,]
argdat1$argstack <- argdat1$name2

argdat2 <- HR_ARG[!HR_ARG$name2 %in% arglist,]
argdat2$argstack <- "Other"

argdat <- rbind(argdat1, argdat2)

# Reorder the 'Category' variable based on the 'Value' variable

argdat$argstack <- factor(argdat$argstack, levels=c( "Other","tet(Q)"  ,  "tet(W)" ,   "ant(6)-Ia" ,"erm(F)"    ,"tet(X)"     ))

# Order the levels consistently across all data frames

hr_arg_stack <-  argdat %>% ggplot(aes(x=labels,y=TotalIDs)) + 
  geom_bar(aes(fill=argstack, color=argstack), stat="identity", position="fill") + 
  ggtitle("Proportional Abundance of HR-ARGs") + 
  ylab("Proportional Abundance") +xlab("")+ theme_classic() +
  theme( legend.title=element_blank(), text=element_text(size=50))

hr_arg_stack

### Bacteria Heatmap
# Glom the taxa at the Phylum level
physeq_g <- tax_glom(physeq, "Phylum")

# Get unique taxa names
taxaorder <- unique(taxa_names(physeq_g))

physeq_prod <- subset_samples(physeq_g, Group=="Production")
physeq_prod_com <- subset_samples(physeq_prod, breed=="Commercial")
physeq_prod_ind <- subset_samples(physeq_prod, breed=="Local")

hm_bac_prod_com <-plot_heatmap(physeq_prod_com, taxa.label= "Phylum", taxa.order=taxaorder, title = "Commercial", trans = log_trans(10),
                               low="#e3fa9e", high="#59730d", na.value="white")  + theme_classic() +   theme(
                                 text = element_text(size = 30),        # General text size
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_blank(),
                                 axis.text = element_text(size = 30),   # Axis text size
                                 axis.title = element_text(size = 30),  # Axis title size
                                 plot.title = element_text(size = 50),  # Plot title size
                                 legend.text = element_text(size = 30), # Legend text size
                                 legend.title = element_text(size = 30) # Legend title size
                               ) + guides(fill=guide_legend(
                                 keywidth=0.1,
                                 keyheight=0.1,
                                 default.unit="inch")
                               ) + scale_x_discrete(guide = guide_axis(angle = 30)) 

hm_bac_prod_com

hm_bac_prod_ind <-plot_heatmap(physeq_prod_ind, taxa.label= "Phylum", taxa.order=taxaorder,  title = "Local",   trans = log_trans(10), low="#d6f5f4", high="#074e54", na.value="white") +
  theme_classic() +   theme(
    text = element_text(size = 30),        # General text size
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 30),   # Axis text size
    axis.title = element_text(size = 30),  # Axis title size
    plot.title = element_text(size = 50),  # Plot title size
    legend.text = element_text(size = 30), # Legend text size
    legend.title = element_text(size = 30) # Legend title size
  ) + guides(fill=guide_legend(
    keywidth=0.1,
    keyheight=0.1,
    default.unit="inch")
  ) + scale_x_discrete(guide = guide_axis(angle = 30)) 
hm_bac_prod_ind

ggarrange(hm_bac_prod_com, hm_bac_prod_ind)


### Pathogen HEatmaps

pathseq_g <- tax_glom(pathseq, "Genus")

taxaorder <- unique(taxa_names(pathseq_g))

pathseq_prod <- subset_samples(pathseq_g, Group=="Production")
pathseq_prod_com <- subset_samples(pathseq_prod, breed=="Commercial")
pathseq_prod_ind <- subset_samples(pathseq_prod, breed=="Local")

hm_path_prod_com <-plot_heatmap(pathseq_prod_com, taxa.label= "Genus", taxa.order=taxaorder, title = "Commercial", trans = log_trans(10), low="#e3fa9e", high="#59750d", na.value="white")  +
  theme_classic() +   theme(
    text = element_text(size = 50),        # General text size
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 50),   # Axis text size
    axis.title = element_text(size = 50),  # Axis title size
    plot.title = element_text(size = 50),  # Plot title size
    legend.text = element_text(size = 50), # Legend text size
    legend.title = element_text(size = 50) # Legend title size
  ) + guides(fill=guide_legend(
    keywidth=0.1,
    keyheight=0.1,
    default.unit="inch")
  ) + scale_x_discrete(guide = guide_axis(angle = 30)) 
hm_path_prod_com

hm_path_prod_ind <-plot_heatmap(pathseq_prod_ind, taxa.label= "Genus", taxa.order=taxaorder,  title = "Local",   trans = log_trans(10), low="#d6f5f4", high="#074e54", na.value="white") +
  theme_classic() +   theme(
    text = element_text(size = 50),        # General text size
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 50),   # Axis text size
    axis.title = element_text(size = 50),  # Axis title size
    plot.title = element_text(size = 50),  # Plot title size
    legend.text = element_text(size = 50), # Legend text size
    legend.title = element_text(size = 50) # Legend title size
  ) + guides(fill=guide_legend(
    keywidth=0.1,
    keyheight=0.1,
    default.unit="inch")
  ) + scale_x_discrete(guide = guide_axis(angle = 50)) 

hm_path_prod_ind

ggarrange(hm_path_prod_com, hm_path_prod_ind)

### ARG heatmap 
ARGseq@tax_table
ARGseq_g <- tax_glom(argprod, "Domain")

taxaorder <- unique(taxa_names(ARGseq_g))

ARGseq_prod <- subset_samples(ARGseq_g, Group=="Production")
ARGseq_prod_com <- subset_samples(ARGseq_prod, breed=="Commercial")
ARGseq_prod_ind <- subset_samples(ARGseq_prod, breed=="Local")

ARGseq_com <- subset_samples(ARGseq_g, breed=="Commercial")
ARGseq_com_fec <- subset_samples(ARGseq_com, Feces==1)

hm_arg_prod_com <-plot_heatmap(ARGseq_prod_com, taxa.label= "Domain", taxa.order=taxaorder, title = "Commercial", trans = log_trans(10), low="#e3fa9e", high="#59720d", na.value="white")  + theme_classic() + 
  theme(
    text = element_text(size = 50),        # General text size
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 50),   # Axis text size
    axis.title = element_text(size = 50),  # Axis title size
    plot.title = element_text(size = 50),  # Plot title size
    legend.text = element_text(size = 50), # Legend text size
    legend.title = element_text(size = 50) # Legend title size
  ) + guides(fill=guide_legend(
    keywidth=0.1,
    keyheight=0.1,
    default.unit="inch")
  ) + scale_x_discrete(guide = guide_axis(angle = 50)) 
hm_arg_prod_com

hm_arg_prod_ind <-plot_heatmap(ARGseq_prod_ind, taxa.label= "Domain", taxa.order=taxaorder,  title = "Local",   trans = log_trans(10), low="#d6f5f4", high="#074e54", na.value="white") + theme_classic() +
  theme(
    text = element_text(size = 50),        # General text size
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 50),   # Axis text size
    axis.title = element_text(size = 50),  # Axis title size
    plot.title = element_text(size = 50),  # Plot title size
    legend.text = element_text(size = 50), # Legend text size
    legend.title = element_text(size = 50) # Legend title size
  ) + guides(fill=guide_legend(
    keywidth=0.1,
    keyheight=0.1,
    default.unit="inch")
  ) + scale_x_discrete(guide = guide_axis(angle = 30)) 
hm_arg_prod_ind


### HR-ARG heatmap 

clinseq_g <- tax_glom(clinprod, "Genus")

taxaorder <- unique(taxa_names(clinseq_g))

clinseq_prod <- subset_samples(clinseq_g, Group=="Production")
clinseq_prod_com <- subset_samples(clinseq_prod, breed=="Commercial")
clinseq_prod_ind <- subset_samples(clinseq_prod, breed=="Local")

clinseq_com <- subset_samples(clinseq_g, breed=="Commercial")
clinseq_com_fec <- subset_samples(clinseq_com, Feces==1)

hm_clin_prod_com <-plot_heatmap(clinseq_prod_com, taxa.label= "Genus", taxa.order=taxaorder, title = "Commercial", trans = log_trans(10), low="#e3fa9e", high="#59730d", na.value="white")  + theme_classic() + 
  theme(
    text = element_text(size = 40),        # General text size
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 40),   # Axis text size
    axis.title = element_text(size = 40),  # Axis title size
    plot.title = element_text(size = 40),  # Plot title size
    legend.text = element_text(size = 40), # Legend text size
    legend.title = element_text(size = 40) # Legend title size
  ) + guides(fill=guide_legend(
    keywidth=0.1,
    keyheight=0.1,
    default.unit="inch")
  ) + scale_x_discrete(guide = guide_axis(angle = 30)) 
hm_clin_prod_com

hm_clin_prod_ind <-plot_heatmap(clinseq_prod_ind, taxa.label= "Genus", taxa.order=taxaorder,  title = "Local",   trans = log_trans(10), low="#d6f5f4", high="#074e54", na.value="white") + theme_classic() +
  theme(
    text = element_text(size = 40),        # General text size
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 40),   # Axis text size
    axis.title = element_text(size = 40),  # Axis title size
    plot.title = element_text(size = 40),  # Plot title size
    legend.text = element_text(size = 40), # Legend text size
    legend.title = element_text(size = 40) # Legend title size
  ) + guides(fill=guide_legend(
    keywidth=0.1,
    keyheight=0.1,
    default.unit="inch")
  ) + scale_x_discrete(guide = guide_axis(angle = 30)) 
hm_clin_prod_ind

##### Arrange Figures 

# Arrange and save the multi-panel figure
sup_bac <- ggarrange( bacvenn, bac_stack, bac_nmds, bac_daa, hm_bac_prod_com, hm_bac_prod_ind, ncol=3, nrow=2, labels=c("A", "B", "C", "D", "E", " "), font.label = list(size = 70))
ggsave(filename = "figures/sup1_bac.png", plot = sup_bac, width =150, height = 150, units = "cm", limitsize = FALSE)

# Arrange and save the multi-panel figure
path_1 <- ggarrange(pathvenn,path_stack, path_nmds, path_daa, hm_path_prod_com, hm_path_prod_ind, ncol=3, nrow=2, labels=c("A", "B", "C", "D", "E", " "), font.label = list(size = 70))
# Arrange and save the multi-panel figure as a PNG file
ggsave(filename = "figures/sup1_path.png", plot = path_1, width =150, height = 150, units = "cm", limitsize = FALSE)

# Arrange and save the multi-panel figure
arg_1<-ggarrange(argvenn,arg_stack, arg_nmds, arg_daa, hm_arg_prod_com, hm_arg_prod_ind, ncol=3, nrow=2, labels=c("A", "B", "C", "D", "E", " "), font.label = list(size = 70))
# Arrange and save the multi-panel figure as a PNG file
ggsave(filename = "figures/sup1_arg.png", plot = arg_1, width =150, height = 150, units = "cm", limitsize=FALSE)

# Arrange and save the multi-panel figure
clin_1 <- ggarrange(clinvenn, hr_arg_stack, hr_arg_nmds, clin_daa, hm_clin_prod_com, hm_clin_prod_ind, ncol=3, nrow=2, labels=c("A", "B", "C", "D", "E", " "), font.label = list(size = 70))
ggsave(filename = "figures/sup1_clin.png", plot = clin_1, width =150, height = 150, units = "cm", limitsize=FALSE)

### Arrange and save fig 1 path & ARGs 
FIG_1 <- ggarrange( path_nmds, path_daa, hm_clin_prod_com, hm_clin_prod_ind, ncol=2, nrow=2,  labels = c("A", "B", "C", " "), 
                    font.label = list(size = 70)) 
FIG_1
ggsave(filename = "figures/fig_1.png", plot = FIG_1, width =100, height = 100, units = "cm")
