
## this function plots rna and atac profile of genes in a table of interest
plot_rna_atac <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw() +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, size = 20, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, size = 20, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 20, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(. ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=20, face="bold", hjust = 1),
      axis.title.y = element_text(size=20, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  return(p)
  
}

## this function extracts the expression and accessibility profile of genes of interest
## need to give rna and atac splines and genes of interest as input
## set the scale T/F
get_rna_atac_profile <- function(rna.splines, atac.splines, genes.tab, scale = T) {
  
  
  sc.rna.dtw.wide <- rna.splines %>% 
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
    mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = scale)) %>%
    as.data.frame()
  
  sc.atac.dtw.wide <- atac.splines %>% 
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
    mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = scale)) %>%
    as.data.frame()
  
  
  sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
    pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')
  
  sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
    pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')
  
  
  ## table of genes to plot their expression 
  tab.genes <- data.frame(TGME49 = gsub("_", "-", genes.tab$gene_name), 
                          Name = genes.tab$ProductDescription)
  
  
  sc.rna.long <- inner_join(sc.rna.mu.scale, tab.genes, by = c('GeneID' = 'TGME49')) 
  sc.rna.long <- sc.rna.long %>% 
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Name)  %>% distinct()
  
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab.genes, by = c('GeneID' = 'TGME49')) 
  sc.atac.long <- sc.atac.long %>% 
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Name) %>% distinct()
  
  
  sc.rna.sc.atac.joint <- inner_join(sc.rna.long, sc.atac.long, 
                                     by = c("time", "GeneID", "Name"))
  colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "scATAC")
  
  sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>% 
    pivot_longer(-c('time', "GeneID", "Name"), 
                 names_to = 'data', values_to = 'normExpr') 
  
  return(sc.rna.sc.atac.joint.long)
  
}

## This function performs dynamic time warping clustering for both rna and atac profiles
## it gets a table of genes, 
## the table should include the product description (gene ID and gene Name) 
## you can specify the number of clusters you are looking for (cannot be less than 2)

clust.df <- function(tab , num.clust) {
  
  k <- num.clust
  sc.rna <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in% tab$gene_name ]
  sc.atac <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tab$gene_name]
  sc.atac<- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% colnames(sc.rna)]
  
  sc.rna.markers.hc_dtw <- dtwClustCurves(sc.rna, nclust = k)
  sc.atac.markers.hc_dtw <- dtwClustCurves(sc.atac, nclust = k)
  
  tab <- tab[tab$gene_name %in% colnames(sc.rna),]
  
  sc.rna.long <- inner_join(sc.rna.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.rna.long <- sc.rna.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription)  %>% distinct()
  
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.atac.long <- sc.atac.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription) %>% distinct()
  
  
  sc.rna.clust.info <- data.frame(GeneID = colnames(sc.rna), cluster = cutree(sc.rna.markers.hc_dtw, k = k))
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac), cluster = cutree(sc.atac.markers.hc_dtw, k = k))
  
  sc.rna.long.clust <- inner_join(sc.rna.long, sc.rna.clust.info, by = 'GeneID')
  sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')
  
  sc.rna.sc.atac.joint <- inner_join(sc.rna.long.clust, sc.atac.long.clust,
                                     by = c("time", "GeneID", "Name"))
  colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "cluster.RNA", "scATAC", "cluster.ATAC")
  
  sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>%
    pivot_longer(-c('time', "GeneID", "Name", "cluster.RNA", "cluster.ATAC"),
                 names_to = 'data', values_to = 'normExpr')
  
  sc.rna.sc.atac.joint.long$cluster.RNA <- paste('C', sc.rna.sc.atac.joint.long$cluster.RNA)
  sc.rna.sc.atac.joint.long$cluster.ATAC <- paste('C', sc.rna.sc.atac.joint.long$cluster.ATAC)
  
  return(sc.rna.sc.atac.joint.long)
  
}


clust.atac.df <- function(tab, num.clust = num.clust){
  
  tab <- tab
  k <- num.clust
  
  sc.atac <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tab$gene_name]
  sc.atac.markers.hc_dtw <- dtwClustCurves(sc.atac, nclust = k)
  tab <- tab[tab$gene_name %in% colnames(sc.atac),]
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.atac.long <- sc.atac.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, group = group, trans.cluster.rna) %>% 
    distinct()
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac), cluster = cutree(sc.atac.markers.hc_dtw, k = k))
  sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')
  colnames(sc.atac.long.clust) <- c("time", "GeneID", "normExpr", "group", "trans.cluster.rna", "cluster.ATAC")
  sc.atac.long.clust$cluster.ATAC <- paste('C', sc.atac.long.clust$cluster.ATAC)
  
  return(sc.atac.long.clust)
}


plot_atac_trand <- function(sc.atac.long.clust){
  
  p  <- ggplot(sc.atac.long.clust, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 15, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 15, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.ATAC ~ ., scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=18, face="bold", hjust = 1),
      axis.title.y = element_text(size=18, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold")) +
    theme(axis.ticks =  element_blank())
  
  return(p)
}


clust.atac.df <- function(tab, num.clust = num.clust){
  
  tab <- tab
  k <- num.clust
  
  sc.atac <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tab$gene_name]
  sc.atac.markers.hc_dtw <- dtwClustCurves(sc.atac, nclust = k)
  tab <- tab[tab$gene_name %in% colnames(sc.atac),]
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.atac.long <- sc.atac.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, group = group, trans.cluster.rna) %>% 
    distinct()
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac), cluster = cutree(sc.atac.markers.hc_dtw, k = k))
  sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')
  colnames(sc.atac.long.clust) <- c("time", "GeneID", "normExpr", "group", "trans.cluster.rna", "cluster.ATAC")
  sc.atac.long.clust$cluster.ATAC <- paste('C', sc.atac.long.clust$cluster.ATAC)
  
  return(sc.atac.long.clust)
}


plot_atac_trand <- function(sc.atac.long.clust){
  
  p  <- ggplot(sc.atac.long.clust, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 15, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 15, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.ATAC ~ ., scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=18, face="bold", hjust = 1),
      axis.title.y = element_text(size=18, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold")) +
    theme(axis.ticks =  element_blank())
  
  return(p)
}


plot_rna_atac_trends.ord <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw() +
    #theme_bw(base_size = 16) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 22, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 22, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.RNA.ordered ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=20, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=22, face="bold", hjust = 1),
      axis.title.y = element_text(size=22, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  return(p)
  
}

