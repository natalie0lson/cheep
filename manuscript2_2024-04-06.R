today <- Sys.Date()
pacman::p_load(phyloseq, vegan, ggplot2, ggpubr, microViz, corncob, ggraph, DT, readr, dplyr, tidyverse, gridExtra, breakaway, microbiome, metagMisc, scales, GUniFrac, kableExtra, optimx, edgeR, xaringan, xaringanExtra, magick)

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


## Pathogens 
phyfec <- subset_samples(pathseq, Feces==1)
path_mkt <- subset_samples(phyfec, site=="Market")

bac_fec <- subset_samples(physeq, Feces == 1)
bac_fec_mkt <- subset_samples(bac_fec, Group=="Market")

# Custom colors for plotting
custom_colors <- c(
  "Commercial/Farm/Fecal, N=4" = "darkolivegreen1",
  "Commercial/Market/Fecal, N=8" = "darkolivegreen3",
  "Commercial/Market/Carcass, N=42" = "darkolivegreen4",
  "Indigenous/Household/Fecal, N=7" = "lightblue1",
  "Indigenous/Market/Fecal, N=7" = "lightblue3",
  "Indigenous/Market/Carcass, N=7" = "lightblue4",
  "Market Rinse Water, N=26" = "goldenrod"
)

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


physeq@sam_data$loc.type.n2 <- gsub(", ", "\n", physeq@sam_data$loc.type.n)
physeq@sam_data$loc.type.n2 <- gsub("/Market/Fecal", " ", physeq@sam_data$loc.type.n2)

pathseq@sam_data$loc.type.n2 <- gsub(", ", "\n", pathseq@sam_data$loc.type.n)
pathseq@sam_data$loc.type.n2 <- gsub("/Market/Fecal", " ", pathseq@sam_data$loc.type.n2)

ARGseq@sam_data$loc.type.n2 <- gsub(", ", "\n", ARGseq@sam_data$loc.type.n)
ARGseq@sam_data$loc.type.n2 <- gsub("/Market/Fecal", " ", ARGseq@sam_data$loc.type.n2)

clinseq@sam_data$loc.type.n2 <- gsub(", ", "\n", clinseq@sam_data$loc.type.n)
clinseq@sam_data$loc.type.n2 <- gsub("/Market/Fecal", " ", clinseq@sam_data$loc.type.n2)

#update names of custom colors 
names(custom_colors) <- gsub(", ", "\n", names(custom_colors))
names(custom_colors) <- gsub("Indigenous", "Local", names(custom_colors))
names(custom_colors) <- gsub("/Market/Fecal", " ", names(custom_colors))


custom_colors
table(physeq@sam_data$loc.type.n2)
### Stacked Bars 
### Most Abundant 
phy_df <- phyloseq_to_df(bac_fec_mkt, addtax=T)

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

phydat$labels[phydat$labels=="IND\n(N=7)"] <- "LOC\n(N=7)"
table(phydat$labels)

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
pathseq_fec <- subset_samples(pathseq, Feces==1)
pathseq_fec_mkt <- subset_samples(pathseq_fec, Group=="Market")
phy_df <- phyloseq_to_df(pathseq_fec_mkt, addtax=T)

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
phydat$phystack <- factor(phydat$phystack, levels=c( "Other", "E. coli"   ,  "S. flexneri" ,"C. coli" ,    "S. sonnei",   "S. enterica"    ))

# Define the desired order of x-axis labels
unique_categories <- unique(phy_df$loc.type.n)


phydat$labels[phydat$labels=="IND\n(N=7)"] <- "LOC\n(N=7)"
table(phydat$labels)


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
mkt_args <- subset(ARG, ARG$Group=="Market")
top.args = aggregate(mkt_args$TotalIDs,
                     by = list(mkt_args$name2),
                     FUN = sum)
top.args$num <- top.args$x
top.args$name2 <- top.args$Group.1
top.args <- top.args[with(top.args,order(-num)),]
top.args <-  na.omit(top.args[1:5, ])

arglist <- (top.args$name2)
arglist
table(mkt_args$clinical)
length(mkt_args$name2)

ARG$argstack <- NA 
argdat1 <- mkt_args[mkt_args$name2 %in% arglist,]

argdat1 <- mkt_args[mkt_args$name2 %in% arglist,]
argdat1$argstack <- argdat1$name2

argdat2 <- mkt_args[!mkt_args$name2 %in% arglist,]
argdat2$argstack <- "Other"

argdat <- rbind(argdat1, argdat2)

# Reorder the 'Category' variable based on the 'Value' variable

argdat$argstack <- factor(argdat$argstack, levels=c( "Other","erm(B)" , "tet(Q)" , "erm(F)","tet(W)", "ant(6)-Ia"  ))



argdat$labels[argdat$labels=="IND\n(N=7)"] <- "LOC\n(N=7)"

# Order the levels consistently across all data frames

arg_stack <-  argdat %>% ggplot(aes(x=labels,y=TotalIDs)) + 
  geom_bar(aes(fill=argstack, color=argstack), stat="identity", position="fill") + 
  ggtitle("Proportional Abundance of AMR Genes") + 
  ylab("Proportional Abundance") +xlab("")+ theme_classic() +  theme( legend.title=element_blank(), text=element_text(size=50))

arg_stack
#ggsave(stacked_bar, file=paste0(fig_dir, "AMR_stacked_bar_", today,".png"), 
#   width = 44.45, height = 27.78, units = "cm", dpi=300)




### hr-ARGs

## Most abundant ARGs
HR_ARG <- subset(mkt_args, mkt_args$clinical==1)
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

argdat$argstack <- factor(argdat$argstack, levels=c( "Other", "erm(B)"  ,  "tet(Q)"  ,  "erm(F)"   , "tet(W)" ,   "ant(6)-Ia"   ))

argdat$labels[argdat$labels=="IND\n(N=7)"] <- "LOC\n(N=7)"
# Order the levels consistently across all data frames

hr_arg_stack <-  argdat %>% ggplot(aes(x=labels,y=TotalIDs)) + 
  geom_bar(aes(fill=argstack, color=argstack), stat="identity", position="fill") + 
  ggtitle("Proportional Abundance of HR-ARGs") + 
  ylab("Proportional Abundance") +xlab("")+ theme_classic() +  theme( legend.title=element_blank(), text=element_text(size=50))

hr_arg_stack


##### NMDS 

## Bacteria NMDS 

physeq_fixed <- physeq %>%
  tax_fix(
    min_length = 2,
    unknowns = c("sp. ANNP30", "Ponticoccus"),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

bac_fec <- subset_samples(physeq_fixed, Feces == 1)
bac_fec_mkt <- subset_samples(bac_fec, Group=="Market")

# Transform phyloseq object

physeq.rel <- microbiome::transform(bac_fec_mkt, "compositional")
otu <- abundances(physeq.rel)
meta <- meta(physeq.rel)
# Significance Test (PERMANOVA)
formula_str <- paste("t(otu) ~", "breed")  # Create the formula string
permanova <- adonis2(as.formula(formula_str), data = meta, permutations = 999999, method = "bray")
p_value <- sprintf("%.4f", permanova$`Pr(>F)`[1])  # Format p-value to four decimal places
p_value

### Side Density Plot 
bac_nmds <- bac_fec_mkt %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc("NMDS") %>%
  ord_plot(color = "loc.type.n2", shape = "loc.type.n2", size = 16) +
  stat_ellipse(aes(group = loc.type.n2)) + ## add ellipses
  scale_color_manual(values=custom_colors)+
  theme_bw() + ggtitle(paste0("Bacterial Composition \n P-Value=", p_value)) + 
  ggside::geom_xsidedensity(aes(fill = loc.type.n2), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = loc.type.n2), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void() + scale_fill_manual(values=custom_colors) +
  theme(
    text = element_text(size = 50),           # General text size
    legend.title = element_blank(),           # No legend title
    plot.title = element_text(size = 50),     # Plot title size
    axis.title = element_text(size = 50),     # Axis titles size
    axis.text = element_text(size = 50),      # Axis text size
    legend.text = element_text(size = 50)     # Legend text size
  )
bac_nmds

## Pathogen NMDS 

pathseq_fec <- subset_samples(pathseq, Feces == 1)
pathseq_fec_mkt <- subset_samples(pathseq_fec, Group=="Market" )

# Transform phyloseq object
physeq.rel <- microbiome::transform(pathseq_fec_mkt, "compositional")
otu <- abundances(physeq.rel)
meta <- meta(physeq.rel)
# Significance Test (PERMANOVA)
formula_str <- paste("t(otu) ~", "breed")  # Create the formula string
permanova <- adonis2(as.formula(formula_str), data = meta, permutations = 999999, method = "bray")
p_value <- sprintf("%.3f", permanova$`Pr(>F)`[1])  # Format p-value to four decimal places
p_value

### Side Density Plot 
path_nmds <- pathseq_fec_mkt %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc("NMDS") %>%
  ord_plot(color = "loc.type.n2", shape = "loc.type.n2", size = 16) +
  stat_ellipse(aes(group = loc.type.n2)) + ## add ellipses
  scale_color_manual(values=custom_colors)+
  theme_bw() + ggtitle(paste0("Pathogen Composition \n P-Value=", p_value)) + 
  ggside::geom_xsidedensity(aes(fill = loc.type.n2), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = loc.type.n2), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void() + scale_fill_manual(values=custom_colors) +
  theme(
    text = element_text(size = 50),           # General text size
    legend.title = element_blank(),           # No legend title
    plot.title = element_text(size = 50),     # Plot title size
    axis.title = element_text(size = 50),     # Axis titles size
    axis.text = element_text(size = 50),      # Axis text size
    legend.text = element_text(size = 50)     # Legend text size
  )
path_nmds

## ARG NMDS 

ARGseq_fec <- subset_samples(ARGseq, Feces == 1)
ARGseq_fec_mkt <- subset_samples(ARGseq_fec, Group=="Market" )
# Transform phyloseq object
physeq.rel <- microbiome::transform(ARGseq_fec_mkt, "compositional")
otu <- abundances(physeq.rel)
meta <- meta(physeq.rel)

# Significance Test (PERMANOVA)
formula_str <- paste("t(otu) ~", "breed")  # Create the formula string
permanova <- adonis2(as.formula(formula_str), data = meta, permutations = 999999, method = "bray")
p_value <- sprintf("%.4f", permanova$`Pr(>F)`[1])  # Format p-value to four decimal places
p_value
arg_nmds <- ARGseq_fec_mkt %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc("NMDS") %>%
  ord_plot(color = "loc.type.n2", shape = "loc.type.n2", size = 16) +
  stat_ellipse(aes(group = loc.type.n2)) + ## add ellipses
  scale_color_manual(values=custom_colors)+
  theme_bw() + ggtitle(paste0("ARG Composition \n P-Value=", p_value)) + 
  ggside::geom_xsidedensity(aes(fill = loc.type.n2), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = loc.type.n2), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void() + scale_fill_manual(values=custom_colors) +
  theme(
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
clinseq_fec_mkt <- subset_samples(clinseq_fec, Group=="Market" )# Transform phyloseq object
physeq.rel <- microbiome::transform(clinseq_fec_mkt, "compositional")
otu <- abundances(physeq.rel)
meta <- meta(physeq.rel)

# Significance Test (PERMANOVA)
formula_str <- paste("t(otu) ~", "breed")  # Create the formula string
permanova <- adonis2(as.formula(formula_str), data = meta, permutations = 999999, method = "bray")
p_value <- sprintf("%.3f", permanova$`Pr(>F)`[1])  # Format p-value to four decimal places
p_value
clin_nmds <- clinseq_fec_mkt %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc("NMDS") %>%
  ord_plot(color = "loc.type.n2", shape = "loc.type.n2", size = 16) +
  stat_ellipse(aes(group = loc.type.n2)) + ## add ellipses
  scale_color_manual(values=custom_colors)+
  theme_bw() + ggtitle(paste0("HR-ARG Composition \n P-Value=", p_value)) + 
  ggside::geom_xsidedensity(aes(fill = loc.type.n2), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = loc.type.n2), alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void() + scale_fill_manual(values=custom_colors) +
  theme(
    text = element_text(size = 50),           # General text size
    legend.title = element_blank(),           # No legend title
    plot.title = element_text(size = 50),     # Plot title size
    axis.title = element_text(size = 50),     # Axis titles size
    axis.text = element_text(size = 50),      # Axis text size
    legend.text = element_text(size = 50)     # Legend text size
  )
clin_nmds
#######################################
#### Differential abundance 
##########################################

source("R/daa_function.R")

### DAA BACTERIA 
bacmkt <- subset_samples(physeq, Group=="Market") %>% tax_glom("Phylum")

#DIFFERENTIAL ABUNDANCE BETWEEN BREEDS (AMONG PRODUCTION SITE)
daa_mkt <- differentialTest(formula = ~ breed,
                            phi.formula = ~ breed,
                            formula_null = ~ 1,
                            phi.formula_null = ~ breed,
                            test = "Wald", boot = FALSE,
                            full_output=FALSE,
                            data = bacmkt,
                            fdr_cutoff = 0.05)

bac_daa <- plot.differentialTest(daa_mkt, level= "Phylum")
bac_daa

### DAA PATH

mkt_prod <- subset_samples(pathseq, Group=="Market") %>% tax_glom("Genus")

#DAA ARGS
daa_mkt <- differentialTest(formula = ~ breed,
                             phi.formula = ~ breed,
                             formula_null = ~ 1,
                             phi.formula_null = ~ breed,
                             test = "Wald", boot = FALSE,
                             full_output=FALSE,
                             data = mkt_prod,
                             fdr_cutoff = 0.05)

path_daa <- plot.differentialTest(daa_mkt, level= "Genus")
path_daa

### DAA ARGS

argmkt <- subset_samples(ARGseq, Group=="Market")  %>% tax_glom("Genus")

daa_argseq <- differentialTest(formula = ~ breed,
                               phi.formula = ~ breed,
                               formula_null = ~ 1,
                               phi.formula_null = ~ breed,
                               test = "Wald", boot = FALSE,
                               full_output=FALSE,
                               data = argmkt,
                               fdr_cutoff = 0.05)

arg_daa <- plot.differentialTest(daa_argseq, level= "Genus")
arg_daa


### DAA HR-ARGS 

clinmkt <- subset_samples(clinseq, Group=="Market")  %>% tax_glom("Genus")

daa_clinseq <- differentialTest(formula = ~ breed,
                               phi.formula = ~ breed,
                               formula_null = ~ 1,
                               phi.formula_null = ~ breed,
                               test = "Wald", boot = FALSE,
                               full_output=FALSE,
                               data = clinmkt,
                               fdr_cutoff = 0.05)

clin_daa <- plot.differentialTest(daa_clinseq, level= "Genus")
clin_daa


library(MicEco)
### Venn Diagrams
### Bacteria 
bacmkt <- subset_samples(physeq, Group=="Market")

### Venn Diagram
bacvenn <- ps_venn(
  bacmkt,
  group = "breed",
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE, 
  labels=list(cex=5),
  #main = "Bacteria",
  
  quantities = list(cex=5, type = c('counts',"percent")),
  
  fill=c("Commercial" = "darkolivegreen3", 
         "Indigenous" = "lightblue3")
  
)


## Pathogen venn diagram 
pathmkt <- subset_samples(physeq, Group=="Market")

### Venn Diagram
pathvenn <- ps_venn(
  pathmkt,
  group = "breed",
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE, 
  labels=list(cex=5),
  #main = "Pathogens",
  
  quantities = list(cex=5, type = c('counts',"percent")),
  
  fill=c("Commercial" = "darkolivegreen3", 
         "Indigenous" = "lightblue3")
  
)



###ARGs 
argvenn <- ps_venn(
  argmkt,
  group = "breed",
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE, 
  labels=list(cex=5),
#main = "Market Site ARGs",
  
  quantities = list(cex=5, type = c('counts',"percent")),
  
  fill=c("Commercial" = "darkolivegreen3", 
         "Indigenous" = "lightblue3")
  
)


## HR-ARG venn 


### Venn Diagram
clinvenn <- ps_venn(
  clinmkt,
  group = "breed",
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE, 
  labels=list(cex=5),
  #main = "Market Site HR-ARGs",
  
  quantities = list(cex=5, type = c('counts',"percent")),
  
  fill=c("Commercial" = "darkolivegreen3", 
         "Indigenous" = "lightblue3")
  
)
#################
### Heat Maps 
#################
### Bacteria Heatmaps 


physeq_g <- tax_glom(physeq, "Phylum")

taxaorder <- unique(taxa_names(physeq_g))

physeq_mkt <- subset_samples(physeq_g, Group=="Market")
physeq_mkt_com <- subset_samples(physeq_mkt, breed=="Commercial")
physeq_mkt_ind <- subset_samples(physeq_mkt, breed=="Local")

hm_phy_mkt_com <-plot_heatmap(physeq_mkt_com, taxa.label= "Phylum", taxa.order=taxaorder, title = "Commercial", trans = log_trans(10), low="#e3fa9e", high="#59730d", na.value="white")  +
  theme_classic() + theme(text=element_text(size=50)) + scale_x_discrete(guide = guide_axis(angle = 40)) #+  scale_fill_gradient(limits=range(67,11661349),  low="#deebe4", high="#2c3831", na.value="white")
hm_phy_mkt_com

hm_phy_mkt_ind <-plot_heatmap(physeq_mkt_ind, taxa.label= "Phylum", taxa.order=taxaorder,  title = "Local",   trans = log_trans(10), low="#d6f5f4", high="#074e54", na.value="white") +
  theme_classic() + theme(text=element_text(size=50))+ scale_x_discrete(guide = guide_axis(angle = 40)) #+  scale_fill_gradient(limits=range(67,11661349),  low="#deebe4", high="#2c3831", na.value="white")
hm_phy_mkt_ind

ggarrange(hm_phy_mkt_com, hm_phy_mkt_ind)
### Pathogen HEatmaps

pathseq_g <- tax_glom(pathseq, "Genus")

taxaorder <- unique(taxa_names(pathseq_g))

pathseq_mkt <- subset_samples(pathseq_g, Group=="Market")
pathseq_mkt_com <- subset_samples(pathseq_mkt, breed=="Commercial")
pathseq_mkt_ind <- subset_samples(pathseq_mkt, breed=="Local")

hm_path_mkt_com <-plot_heatmap(pathseq_mkt_com, taxa.label= "Genus", taxa.order=taxaorder, title = "Commercial", trans = log_trans(10), low="#e3fa9e", high="#59730d", na.value="white")  + theme_classic() +
  theme(text=element_text(size=50), axis.title.x=element_blank(), axis.title.y=element_blank())+ scale_x_discrete(guide = guide_axis(angle = 40))  #+  scale_fill_gradient(limits=range(67,11661349),  low="#deebe4", high="#2c3831", na.value="white")
hm_path_mkt_com

hm_path_mkt_ind <-plot_heatmap(pathseq_mkt_ind, taxa.label= "Genus", taxa.order=taxaorder,  title = "Local",   trans = log_trans(10), low="#d6f5f4", high="#074e54", na.value="white") +
  theme_classic() + theme(text=element_text(size=50), axis.title.x=element_blank(), axis.title.y=element_blank())+ scale_x_discrete(guide = guide_axis(angle = 40)) #+  scale_fill_gradient(limits=range(67,11661349),  low="#deebe4", high="#2c3831", na.value="white")
hm_path_mkt_ind

ggarrange(hm_path_mkt_com, hm_path_mkt_ind)


#### ARG Heat map 


argseq_g <- tax_glom(argmkt, "Domain")

taxaorder <- unique(taxa_names(argseq_g))

argseq_mkt <- subset_samples(argseq_g, Group=="Market")
argseq_mkt_com <- subset_samples(argseq_mkt, breed=="Commercial")
argseq_mkt_ind <- subset_samples(argseq_mkt, breed=="Local")

argseq_com <- subset_samples(argseq_g, breed=="Commercial")
argseq_com_fec <- subset_samples(argseq_com, Feces==1)

hm_arg_mkt_com <-plot_heatmap(argseq_mkt_com, taxa.label= "Domain", taxa.order=taxaorder, title = "Commercial", trans = log_trans(10), low="#e3fa9e", high="#59730d", na.value="white")  + 
  theme_classic() +   theme(
    text = element_text(size = 50), 
    axis.title.x = element_blank() ,  
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
  ) + scale_x_discrete(guide = guide_axis(angle = 40)) 
hm_arg_mkt_com

hm_arg_mkt_ind <-plot_heatmap(argseq_mkt_ind, taxa.label= "Domain", taxa.order=taxaorder,  title = "Local",   trans = log_trans(10), low="#d6f5f4", high="#074e54", na.value="white") +
  theme_classic() +  theme(
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
  ) + scale_x_discrete(guide = guide_axis(angle = 40)) 
hm_arg_mkt_ind



### HR-ARG heatmap 

clinseq_g <- tax_glom(clinmkt, "Genus")

taxaorder <- unique(taxa_names(clinseq_g))

clinseq_mkt <- subset_samples(clinseq_g, Group=="Market")
clinseq_mkt_com <- subset_samples(clinseq_mkt, breed=="Commercial")
clinseq_mkt_ind <- subset_samples(clinseq_mkt, breed=="Local")

clinseq_com <- subset_samples(clinseq_g, breed=="Commercial")
clinseq_com_fec <- subset_samples(clinseq_com, Feces==1)

hm_clin_mkt_com <-plot_heatmap(clinseq_mkt_com, taxa.label= "Genus", taxa.order=taxaorder, title = "Commercial", trans = log_trans(10), low="#e3fa9e", high="#59730d", na.value="white")  + 
  theme_classic() +   theme(
    text = element_text(size = 50), 
    axis.title.x = element_blank() ,  
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
  ) + scale_x_discrete(guide = guide_axis(angle = 40)) 
hm_clin_mkt_com

hm_clin_mkt_ind <-plot_heatmap(clinseq_mkt_ind, taxa.label= "Genus", taxa.order=taxaorder,  title = "Local",   trans = log_trans(10), low="#d6f5f4", high="#074e54", na.value="white") +
  theme_classic() +  theme(
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
  ) + scale_x_discrete(guide = guide_axis(angle = 40)) 
hm_clin_mkt_ind


##### Arrange Figures 

# Arrange and save the multi-panel figure
bac_2 <- ggarrange(bacvenn,bac_stack, bac_nmds, bac_daa, hm_phy_mkt_com, hm_phy_mkt_ind, ncol=3, nrow=2, labels=c("A", "B", "C", "D", "E", " "), font.label = list(size = 70)) 
ggsave(filename = "figures/sup2_bac.png", plot = bac_2, width =150, height = 150, units = "cm", limitsize = FALSE)


# Arrange and save the multi-panel figure
path_2 <- ggarrange(pathvenn,path_stack, path_nmds, path_daa, hm_path_mkt_com, hm_path_mkt_ind, ncol=3, nrow=2, labels=c("A", "B", "C", "D", "E", " "), font.label = list(size = 70)) 
# Arrange and save the multi-panel figure as a PNG file
ggsave(filename = "figures/sup2_path.png", plot = path_2, width =150, height = 150, units = "cm", limitsize=FALSE)

# Arrange and save the multi-panel figure
arg_2<-ggarrange(argvenn,arg_stack, arg_nmds, arg_daa, hm_arg_mkt_com, hm_arg_mkt_ind, ncol=3, nrow=2, labels=c("A", "B", "C", "D", "E", " "), font.label = list(size = 70)) 
# Arrange and save the multi-panel figure as a PNG file
ggsave(filename = "figures/sup2_arg.png", plot = arg_2, width =150, height = 150, units = "cm", limitsize=FALSE)

# Arrange and save the multi-panel figure

clin_2 <- ggarrange(clinvenn, hr_arg_stack, clin_nmds, clin_daa, hm_clin_mkt_com, hm_clin_mkt_ind,ncol=3, nrow=2, labels=c("A", "B", "C", "D", "E", " "), font.label = list(size = 70)) 

ggsave(filename = "figures/sup2_clin.png", plot = clin_2, width =150, height = 150, units = "cm", limitsize=FALSE)



### Arrange and save fig 1 path & ARGs 
FIG_2 <- ggarrange( path_nmds, path_daa, hm_clin_mkt_com, hm_clin_mkt_ind, ncol=2, nrow=2,  labels = c("A", "B", "C", " "), 
                    font.label = list(size = 70)) 
FIG_2
ggsave(filename = "figures/fig_2.png", plot = FIG_2, width =100, height = 100, units = "cm")

