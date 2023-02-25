#!/usr/bin/Rscript
library(optparse)
library(data.table)
library(plyr)
library(ggplot2)
library(showtext)

# apt-get install libfreetype6-dev 
# install.packages("showtext")

parser <- OptionParser()
parser <- add_option(parser, c("-e", "--eigenvector"), type="character", 
                help="eigenvector from plink")
parser <- add_option(parser, c("-v", "--eigenvalue"), type="character", 
                help="eigenvalue from plink", default = NULL)
parser <- add_option(parser, c("-f", "--fam"), type="character",  
                help="fam file from plink")                
parser <- add_option(parser, c("-i", "--info"), type="character", default = "/volume/prsdata/GWAS_REF/relationships_w_pops_121708.hapmap3",  
                help="population info file")      
parser <- add_option(parser, c("-p", "--pop"), type="character", default = "CHD-CHB-JPT",  
                help="population to use")    
parser <- add_option(parser, c("-m", "--method"), type="character", default = "DIST",  
                help="SD or DIST")
parser <- add_option(parser, c("-s", "--sd_multi"), type="integer", default = 3,  
                help="SD multiplier; used when method")              
parser <- add_option(parser, c("-n", "--num_pc"), type="integer", default = 10,  
                help="number of PC")              
parser <- add_option(parser, c("-o", "--output"), type="character",
                default="PCA", help="output path")
opt = parse_args(parser)

cat("PCA.R: Analyze population stratification for merge of hapmap and target\n")

eigenvector <- opt$eigenvector
eigenvalue <- opt$eigenvalue
Input_fam <- opt$fam
pop_info <- opt$info
selected_pop <- opt$pop
METHOD <- opt$method
OUT_PATH <- opt$output
SD_mutiplier = opt$sd_multi
num_pc = opt$num_pc

selected_pop = strsplit(selected_pop,"-")[[1]]
col_to_cal = sapply(1:num_pc,function(x) paste0("PC",x))
#col_to_cal = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

##############################################
# Step 1: process eigen value
##############################################
if(is.null(eigenvalue)){
  eigenvalue = sub("vec$","val",eigenvector)
  cat("PCA.R: read ",eigenvalue, " according to path of ", eigenvector,"\n")
}
if(file.exists(eigenvalue)){
  EIGENVAL = fread(eigenvalue, header = FALSE)
  EIGENVAL = EIGENVAL$V1/sum(EIGENVAL$V1)
  EIGENVAL = EIGENVAL[1:length(col_to_cal)]
}else{
  EIGENVAL = rep(1,length(col_to_cal))
}

##############################################
# Step 2: Read example datasets
##############################################
mds <- fread(eigenvector, header=TRUE)
target <- fread(Input_fam, header = FALSE)
names(target) = c("#FID","IID","P","M","SEX","PHENO1")
hapmap <- fread(pop_info, header = TRUE)

target = target[,c("#FID","IID")]
hapmap = hapmap[,c("FID","IID", "population")]
names(hapmap) = c("#FID","IID", "population")
DA = merge(target, hapmap, by = c("#FID","IID"), all.x = TRUE, sort = FALSE)
DA[ is.na(DA$population), population := "Target"]
ll = length(DA$IID)

# need to modify this to pop_list, but since the color is fixed, it may take some time
population= c("Target", "CHB", "CHD", "JPT", "CEU", "TSI",
                "GIH", "MEX",  "ASW", "MKK", "LWK", "YRI")
list_colormap_race <- c("#f24405", "#f99600", "#a0060d", "#ffc700", "#d5d7e3", "#7e84a3",
                        "#21d59b", "#106a4d", "#0058ff", "#acc8ff", "#7a0dff", "#e5eeff")

DA$population = factor(DA$population, levels = c("Target", "CHB", "CHD", "JPT", "CEU", "TSI",
                "GIH", "MEX",  "ASW", "MKK", "LWK", "YRI"))

all_population = factor(as.character(DA$population))


#############################################
# Calculate the distance between samples and races
#############################################
# yilun 20210715: Add METHOD == "SD"; Fix bug that population is not selected

if ( ! all(selected_pop %in% levels(DA$population))){
  cat("selected_pop is invalid or None, not doing population selection according mds\n")
}else if(METHOD=="SD"){
  # Standardized EIGENVAL since not all PCs were used but we need the sum of weight as 1
  EIGENVAL_SUM1 = EIGENVAL/sum(EIGENVAL)

  #Calculate pop center for "CHD", "CHB", "JPT"
  pos <- which(all_population %in% selected_pop)
  M_hapmap = as.matrix(mds[pos,..col_to_cal])
  pop_center <- colMeans(M_hapmap)

  # Calcuate stdev for each componenet
  pc_sd = apply(M_hapmap, 2, sd)

  # Get target sample
  pos <- which(all_population == "Target")
  M_target = as.matrix(mds[pos,..col_to_cal])

  # Calculate how many sd is for each sample and component 
  deviation = abs(sweep(M_target,2,pop_center,FUN="-"))
  dist_over_sd = sweep(deviation,2,pc_sd,FUN="/")
  # sum and weighted for each componenet 
  w_dist_over_sd = sweep(dist_over_sd,2,EIGENVAL_SUM1,FUN="*")
  overall_SD = rowSums(w_dist_over_sd)
  # Get valid sample
  II = which(overall_SD <= SD_mutiplier)

  if (length(II) == 0){
    cat("PCA.R: Error! No sample left\n")
    stop()
  }

  cat(paste0("PCA.R: ",length(pos) - length(II)," out of ",length(pos)," were filtered out based on population stratification SD method\n"))
  mds_pass = mds[pos[II], c("#FID","IID")] 
  write.table(mds_pass, paste0(OUT_PATH, ".pca_qc.ind_list"), col.names = F, row.names = F, quote = F, sep = "\t")
  cat("PCA.R: Save as ",paste0(OUT_PATH, ".pca_qc.ind_list"),"\n")

  DA[pos[-II], population := "failed_target" ]
  write.table(DA, paste0(OUT_PATH, ".pca_qc.csv"), row.names = F, quote = F, sep = ",")
  cat("PCA.R: Save as ",paste0(OUT_PATH, ".pca_qc.csv"),"\n")

}else if(METHOD=="DIST"){

  # population center
  col_to_cal = c("PC1","PC2","PC3")
  all_population = levels(DA$population)
  all_population = all_population[which(all_population!="Target")]
  pop_center <- c()
  for(race in all_population){
    if(race != "target"){
      M = as.matrix(mds[all_population == race,..col_to_cal])
      pop_center <- rbind(pop_center, colMeans(M))
    }
  }

  # closest population from samples
  II <- c()
  target_pos = which(DA$population == "target")
  for(target_row in target_pos){
    dist <- pop_center
    dist[,1] <- dist[,1] - mds[target_row][[col_to_cal[1]]]
    dist[,2] <- dist[,2] - mds[target_row][[col_to_cal[2]]]
    dist[,3] <- dist[,3] - mds[target_row][[col_to_cal[3]]]
    euclidean_dist <- rowSums(dist^2)
    closest_pop <- all_population[which(euclidean_dist == min(euclidean_dist))]
    if(!closest_pop %in% selected_pop){
      DA[target_row, "population"] = paste0("target_",closest_pop)
    }else{
      II = c(II, target_row)
    }
  }
  
  if (length(II) == 0){
    cat("PCA.R: Error! No sample left\n")
    stop()
  }

  cat(paste0("PCA.R: ",length(target_pos) - length(II)," out of ",length(target_pos)," were filtered out based on population stratification DIST method\n"))
  mds_pass = mds[pos[II], c("#FID","IID")] 
  write.table(mds_pass, paste0(OUT_PATH, ".pca_qc.ind_list"), col.names = F, row.names = F, quote = F, sep = "\t")

  write.table(DA, paste0(OUT_PATH, ".pca_qc.csv"), row.names = F, quote = F, sep = ",")

}else{
  print("METHOD should be SD or DIST")
}

# Get File Path
GetFilePath <- function(){
  comman.args = commandArgs() 
  my.file.path = grep("--file=", comman.args, value = T)
  file.path = strsplit(my.file.path,"=")[[1]][2]
  file.dir = dirname(file.path)
  return(file.dir)
}

# Deal with font
GetFont <- function(file.path){
  default_font = names(pdfFonts())[1]
  # We assume that font located in file.path + /Report/fonts/
  # <<- just assign to to upper env, please dont use this in nested function
  font_path = paste0(file.path,"/Report/fonts/SFProText-Regular.ttf")
  SFPro <<- default_font
  if(file.exists(font_path)){
    SFPro <<- "SFPro"
    font_add(family = "SFPro", regular = font_path)
  }
  font_path = paste0(file.path,"/Report/fonts/SFProText-Light.ttf")
  SFProLight <<- default_font
  if(file.exists(font_path)){
    SFProLight <<- "SFProLight"
    font_add(family = "SFProLight", regular = font_path)
  }
  font_path = paste0(file.path,"/Report/fonts/Quicksand-Medium.ttf")
  Quicksand <<- default_font
  if(file.exists(font_path)){
    Quicksand <<- "Quicksand"
    font_add(family = "Quicksand", regular = font_path)
  }
  showtext_auto()
}
GetFont(GetFilePath())

# Plot 
size_list <- rep(0.6, ll)
size_list[which(all_population == "target")] <- 1.2

EIGENVAL2 = EIGENVAL*100
my_theme = theme(
      axis.title = element_text(size=32, family = Quicksand, lineheight = 0),
      axis.title.x = element_text(margin = margin(t =  12)),
      axis.title.y = element_text(margin = margin(r =  12)),
      axis.text.y = element_text(margin = margin(r =  6)),
      axis.text.x = element_text(margin = margin(t =  6)),
      axis.text = element_text(size=40, color = "#7e84a3",family = Quicksand),
      axis.ticks = element_blank(),
      axis.line.y = element_blank(),
      axis.line.x.top = element_line(colour = "#f1f1f5", size = 0.8),
      axis.line.x.bottom = element_line(colour = "#f1f1f5", size = 0.8),
      legend.text = element_text(size = 32, family = Quicksand, margin = margin(r = 8, unit = "pt")),
      legend.text.align = 0,
      legend.justification = 'bottom',
      legend.box.margin = margin(),
      legend.key = element_blank(),
      plot.title = element_text(size = 54, margin = margin(b =  30, l = 0), family = SFPro, hjust=-0.03),
      plot.title.position = 'plot',
      plot.margin = margin(b = 10, t =  10, l = 20, r = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#f1f1f5", size = 0.8),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      )

mds$pop = factor(all_population, levels = c("Target", "CHB", "CHD", "JPT", "CEU", "TSI",
                "GIH", "MEX",  "ASW", "MKK", "LWK", "YRI"))

p <- ggplot(data = mds) +
  geom_point(aes(x = PC1, y = PC2, color = pop), size = size_list) + 
  scale_color_manual(breaks=population, values=list_colormap_race) + my_theme +
  labs(x = paste0("PC",1,": ",round(EIGENVAL2[1],2),"%"),
    y = paste0("PC",2,": ",round(EIGENVAL2[2],2),"%"), color = "", title = "All Population") +
  guides(color = guide_legend(label.position = "left")) 

ggsave(paste0(OUT_PATH, ".C1_C2.png"),p, height = 7, width = 8.5, dpi = 300)
cat("PCA.R: Save as ",paste0(OUT_PATH, ".C1_C2.png"),"\n")

p <- ggplot(data = mds) +
  geom_point(aes(x = PC1, y = PC3, color = pop), size = size_list) + 
  scale_color_manual(breaks=population, values=list_colormap_race) + my_theme +
  labs(x = paste0("PC",1,": ",round(EIGENVAL2[1],2),"%"),
    y = paste0("PC",3,": ",round(EIGENVAL2[3],2),"%"), color = "", title = "All Population") +
  guides(color = guide_legend(label.position = "left")) 

ggsave(paste0(OUT_PATH, ".C1_C3.png"),p, height = 7, width = 8.5, dpi = 300) 
cat("PCA.R: Save as ",paste0(OUT_PATH, ".C1_C3.png"),"\n")

p <- ggplot(data = mds) +
  geom_point(aes(x = PC2, y = PC3, color = pop), size = size_list) + 
  scale_color_manual(breaks=population, values=list_colormap_race) + my_theme +
  labs(x = paste0("PC",2,": ",round(EIGENVAL2[2],2),"%"),
    y = paste0("PC",3,": ",round(EIGENVAL2[3],2),"%"), color = "", title = "All Population") +
  guides(color = guide_legend(label.position = "left")) 

ggsave(paste0(OUT_PATH, ".C2_C3.png"),p, height = 7, width = 8.5, dpi = 300)
cat("PCA.R: Save as ",paste0(OUT_PATH, ".C2_C3.png"),"\n")
