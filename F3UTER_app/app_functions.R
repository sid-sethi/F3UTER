
tissue_format <- c(     "amygdala" =	"brain_amygdala",
                        "anteriorcingulatecortexba24" =	"brain_anterior_cingulate_cortex_ba24",
                        "caudatebasalganglia" =	"brain_caudate_basal_ganglia",
                        "frontalcortexba9" =	"brain_frontal_cortex_ba9",
                        "hippocampus" =	"brain_hippocampus",
                        "hypothalamus"	= "brain_hypothalamus",
                        "nucleusaccumbensbasalganglia" =	"brain_nucleus_accumbens_basal_ganglia",
                        "putamenbasalganglia" =	"brain_putamen_basal_ganglia",
                        "spinalcordcervicalc-1" =	"brain_spinal_cord_cervical_c_1",
                        "substantianigra" =	"brain_substantia_nigra"
)



extract_gene_er_data <- function(gene_id, tissue_input, datadb){
  
  if(tissue_input %in% "All") {
    
      x<- glue_sql("SELECT seqnames, start, end, width, tissue, [predicted.response], [predicted.prob], ER_tissue_specificity_category, distance_from_geneEnd, annotationType_split_read_annot, uniq_genes_split_read_annot, associated_gene, id FROM er_raw_all WHERE associated_gene = {gene_id}", .con = datadb) %>% 
        dbGetQuery(datadb, .) 
  }else{
    
    if(sum(tissue_format %in% tissue_input) > 0){
      tissue_formatted <- tissue_format[tissue_format %in% tissue_input] %>% names()
    }
    else{
      tissue_formatted <- tissue_input
    }
    
      x <- x<- glue_sql("SELECT seqnames, start, end, width, tissue, [predicted.response], [predicted.prob], ER_tissue_specificity_category, distance_from_geneEnd, annotationType_split_read_annot, uniq_genes_split_read_annot, associated_gene, id FROM er_raw_all WHERE associated_gene = {gene_id} AND tissue = {tissue_formatted}", .con = datadb) %>% 
        dbGetQuery(datadb, .)
      
  }
  
  return(x)
  
}


extract_er_data <- function(er_id, datadb){
  
  y <- glue_sql("SELECT * FROM erData WHERE id = {er_id}", .con = datadb) %>% 
    dbGetQuery(datadb, .) %>% 
    mutate(class = "ER") %>%
    dplyr::select(-id, -polyA_signal) %>% 
    pivot_longer(cols = meanPd:RepeatsOverlap, names_to = "feature", values_to = "mean") %>%
    dplyr::mutate(sd = NA) %>% 
    bind_rows(dbGetQuery(datadb, 'SELECT * FROM gsData')) %>% 
    dplyr::mutate(dummy = "dummy")
  
  return(y)
  
}




plot_features <- function(df){
  
    # df$feature <- factor(df$feature, levels = c(
    #   "mean_phastCons7way", "meanPd", "meanEE", "RepeatsOverlap", "GA", "CG", "bdnatwist", "AG", "TC", "TG",
    #   "bendability", "proteindnatwist", "AC", "G", "cpgislands", "CA", "CT", "AA", "A", "CC", "cpnpgislands",
    #   "TA", "TT", "GT", "AT", "basestacking", "GG", "C", "aphilicity", "GC", "proteindeformation",
    #   "cpnpgcpgislands", "nucleosomeposition", "duplexstabilitydisruptenergy", "T", "propellortwist",
    #   "bendingstiffness", "zdna", "dnadenaturation", "duplexstabilityfreeenergy"
    #   )
    # )
  
    df <- df %>% mutate(feature = fct_relevel(feature, 
                                             "mean_phastCons7way", "meanPd", "meanEE", "RepeatsOverlap", "GA", "CG", "bdnatwist", "AG", "TC", "TG",
                                             "bendability", "proteindnatwist", "AC", "G", "cpgislands", "CA", "CT", "AA", "A", "CC", "cpnpgislands",
                                             "TA", "TT", "GT", "AT", "basestacking", "GG", "C", "aphilicity", "GC", "proteindeformation",
                                             "cpnpgcpgislands", "nucleosomeposition", "duplexstabilitydisruptenergy", "T", "propellortwist",
                                             "bendingstiffness", "zdna", "dnadenaturation", "duplexstabilityfreeenergy"))  
    
    z <- ggplot(df, aes(x=factor(dummy), y=mean, color = class, fill= class)) +
    geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd),lwd = 1, width = 0.2, position = position_dodge(0.7)) +
    geom_point(size = 1, shape = 15, fill="white", stroke = 2, position = position_dodge(0.7)) +
    scale_fill_manual(values= c("#D40636", "#00AFBB", "#E7B800")) +
    scale_color_manual(values= c("#D40636", "#00AFBB", "#E7B800")) +
    theme_bw() +
    facet_wrap(~ feature, scales = "free", ncol=6) +
    theme(plot.title = element_text(size=13, hjust=0.5),
          legend.key.size = unit(0.50,"cm"),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position= c(0.80,0.07),
          legend.background = element_rect(fill = "transparent",colour = NA),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 9),
          #panel.border = element_rect(colour="BLACK",size=0.4),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill="transparent"),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill="#FFF6EE"),
          strip.text = element_text(size=10, face="bold")
    )
    
    return(z)
}




.get_coords_to_plot <- function(gene_tx_to_plot, er_id, geneName, zoom_factor){
  
  coords_to_plot <-
    list(
      chr = gene_tx_to_plot %>% seqnames() %>% as.character() %>% unique(),
      strand = gene_tx_to_plot %>% strand() %>% as.character() %>% unique(),
      gene_start = gene_tx_to_plot %>% start(),
      gene_end = gene_tx_to_plot %>% end(),
      gene_name = geneName,
      er_id = er_id
    )
  
  coords_to_plot[["er_start"]] <- as.integer(str_split(er_id, ":")[[1]][3])
  coords_to_plot[["er_end"]] <- as.integer(str_split(er_id, ":")[[1]][4])
  coords_to_plot[["er_tissue"]] <- str_split(er_id, ":")[[1]][1]
  coords_to_plot[["ensg_id"]] <- gene_tx_to_plot$gene_id
  
  #coords_to_plot[["range_exon_start_end"]] <- ifelse(IRanges::width(gene_tx_to_plot)/zoom_factor > 3000, IRanges::width(gene_tx_to_plot)/zoom_factor, 3000)
  
  coords_to_plot[["range_exon_start_end"]] <- IRanges::width(gene_tx_to_plot)/zoom_factor
  
  if(coords_to_plot[["strand"]] %in% "+"){
    coords_to_plot[["x_min"]] <- coords_to_plot[["gene_end"]] - coords_to_plot[["range_exon_start_end"]]
    coords_to_plot[["x_max"]] <- coords_to_plot[["er_end"]] + (coords_to_plot[["range_exon_start_end"]]/5)
    coords_to_plot[["segment_start"]] <- coords_to_plot[["x_min"]]
    coords_to_plot[["segment_end"]] <- coords_to_plot[["gene_end"]]
  }else{
    coords_to_plot[["x_min"]] <- coords_to_plot[["er_start"]] - (coords_to_plot[["range_exon_start_end"]]/5)
    coords_to_plot[["x_max"]] <- coords_to_plot[["gene_start"]] + coords_to_plot[["range_exon_start_end"]]
    coords_to_plot[["segment_start"]] <- coords_to_plot[["gene_start"]]
    coords_to_plot[["segment_end"]] <- coords_to_plot[["x_max"]]
  }
  
  return(coords_to_plot)
}




.get_gene_exons_to_plot <- function(gene_tx_id_list, coords_to_plot_gr, ref){
  
  # filter for exons of gene/tx of interest
  gene_exons_to_plot <- GenomicFeatures::exons(ref, filter = gene_tx_id_list)
  # check
  if (length(gene_exons_to_plot) == 0) {stop("No exons were retrieved for the gene from TxDb object")}
  # if exons overlap (e.g. for gene inputs), disjoin for plotting
  gene_exons_to_plot <- gene_exons_to_plot %>% GenomicRanges::disjoin()
  # filter for exons which overlap the plot coordinates
  gene_exons_hits <- GenomicRanges::findOverlaps(gene_exons_to_plot, coords_to_plot_gr)
  gene_exons_to_plot <- gene_exons_to_plot[S4Vectors::queryHits(gene_exons_hits)] %>% 
    as.data.frame()
  #mutate(predicted.response = "exons")
  
  return(gene_exons_to_plot)
}



.plot_gene_track <- function(coords_to_plot, gene_exons_to_plot) {
  
  # plot gene line
  gene_track <- ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(
      x = coords_to_plot[["segment_start"]],
      xend = coords_to_plot[["segment_end"]],
      y = 0, yend = 0
    ),
    size = 2
    )
  
  # plot exons
  gene_track <- gene_track +
    ggplot2::geom_rect(
      data = gene_exons_to_plot,
      ggplot2::aes(
        xmin = start, xmax = end,
        ymin = -0.25, ymax = 0.25
      ),
      colour = "black",
      fill = ggpubr::get_palette("Greys", 10)[4]
    )
  
  # plot strand arrow
  gene_track <- gene_track +
    ggplot2::geom_segment(ggplot2::aes(
      x = ifelse(coords_to_plot[["strand"]] == "+",
                 coords_to_plot[["segment_start"]],
                 coords_to_plot[["segment_end"]]
      ),
      xend = ifelse(coords_to_plot[["strand"]] == "+",
                    coords_to_plot[["segment_start"]] + coords_to_plot[["range_exon_start_end"]] / 8,
                    coords_to_plot[["segment_end"]] - coords_to_plot[["range_exon_start_end"]] / 8
      ),
      y = 0.50, yend = 0.50
    ),
    size = 0.75,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.3, units = "cm"))
    ) +
    # add gene strand label
    ggplot2::annotate(geom="text",
      x=ifelse(coords_to_plot[["strand"]] == "+",
               coords_to_plot[["segment_start"]],
               coords_to_plot[["segment_end"]] - coords_to_plot[["range_exon_start_end"]] / 8
      ), 
      y=0.58, 
      label="gene strand", hjust=0
    )+
    # add gene name
    ggplot2::annotate(geom="text",
                      x = coords_to_plot[["segment_start"]] ,
                      y=0.30, 
                      label= coords_to_plot[["gene_name"]], hjust = 0
    ) +
    # add track label
    ggplot2::annotate(geom="text",
                      x = (coords_to_plot[["x_max"]] + coords_to_plot[["x_min"]])/2 ,
                      y=0.40, 
                      label="Gene and ER track", col = "blue", fontface = "bold", hjust = 0
    )
  
  # add scale/theme aesthetic tweaks
  gene_track <- gene_track +
    ggplot2::scale_y_continuous(limits = c(-1, 1.3)) +
    ggplot2::scale_x_continuous(
      name = stringr::str_c("Chromosome: ", coords_to_plot[["chr"]] )
    ) +
    coord_cartesian(xlim = c(coords_to_plot[["x_min"]], coords_to_plot[["x_max"]])) +
    ggpubr::theme_pubclean(flip = TRUE) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )
  
  return(gene_track)
  
}



.plot_add_ers <- function(coords_to_plot, gene_track_plot, ers_to_plot){
  
  class <- ers_to_plot$predicted.response %>% unique()
  
  colors = list(
    three_utrs = c("#00AFBB"),
    non_utrs = c("#E7B800"),
    both = c("#00AFBB", "#E7B800")
  )
  
  if(length(class) > 1){
    selected_col <- colors[["both"]]
  }else if(class %in% "3_UTR"){
    selected_col <- colors[["three_utrs"]]
  }else{
    selected_col <- colors[["non_utrs"]]
  }
  
  
  # plot ER exons
  gene_track <- gene_track_plot +
    ggplot2::geom_rect(
      data = ers_to_plot,
      ggplot2::aes(
        xmin = start, xmax = end,
        ymin = -0.25, ymax = 0.25,
        fill = predicted.response
      ),
      colour = "black"
    ) +
    scale_fill_manual(values= selected_col) +
    # add ER highlight
    ggplot2::annotate(geom="text",
                      x = (coords_to_plot[["er_start"]] + coords_to_plot[["er_end"]])/2 ,
                      y=0.28, 
                      label="*", col = "black", fontface = "bold", hjust = 0.5, size=8
    )
}




.plot_add_pas <- function(coords_to_plot, gene_er_track_plot, pas_to_plot){
  # plot PAS exons
  gene_track <- gene_er_track_plot +
    ggplot2::geom_rect(
      data = pas_to_plot,
      ggplot2::aes(
        xmin = start, xmax = end,
        ymin = -0.75, ymax = -0.50
      ),
      colour = "black", fill = "cornflowerblue"
    ) +
    # add track label
    ggplot2::annotate(geom="text",
                      x = (coords_to_plot[["x_max"]] + coords_to_plot[["x_min"]])/2 ,
                      y=-0.40, 
                      label="Poly(A) atlas track", col = "blue", fontface = "bold", hjust = 0
    )
}



# mark the midpoints to label with the count
.mark_mid <- function(x) {
  # the point with the lowest absolute diff from the mid
  mid_index <- which.min(abs(x - mean(x)))
  mid_marked <- logical(length = length(x))
  mid_marked[mid_index] <- TRUE
  return(mid_marked)
}




get_junction_points <- function(er_id, datadb, coords_to_plot, ncp=25){
  
  junc_type <- dbGetQuery(datadb, 'SELECT * FROM er_raw_all WHERE id = ?', params = c(er_id))
  
  if(junc_type$annotationType_split_read_annot %in% "none"){
    return()
  }else{
    
    p_annot_junc_ids = junc_type$p_annot_junc_ids_split_read_annot %>% 
      str_split(";") %>%
      unlist() %>%
      unique() %>% 
      na.omit()
    
    unannot_junc_ids = junc_type$unannot_junc_ids_split_read_annot %>%
      str_split(";") %>%
      unlist() %>%
      unique() %>% 
      na.omit()
    
    # make junID vector
    junc_ids <- c(p_annot_junc_ids, unannot_junc_ids)
    
    if(coords_to_plot[["er_tissue"]] %in% names(tissue_format)){
      tissue <- tissue_format[coords_to_plot[["er_tissue"]]]
    }
    else{
      tissue <- coords_to_plot[["er_tissue"]]
    }
    
    
    # load junction data
    
    sr <- glue_sql("SELECT * FROM splitReadData WHERE tissue = {tissue} AND junID IN ({junc_ids*})", .con = datadb) %>% 
      dbGetQuery(datadb, .)
    
    sr$type = c(rep("annotated", times = length(p_annot_junc_ids)), rep("unannotated", times = length(unannot_junc_ids)))
    sr$index <- 1:nrow(sr)
    #gc()
    
    # calculate the points to plot the curve for each junction
    junctions_points <- grid:::calcControlPoints(
      x1 = sr$start, x2 = sr$stop,
      y1 = 0, y2 = 0,
      angle = 90,
      curvature = -0.5,
      ncp = ncp
    ) 
    
    junctions_points <-
      dplyr::tibble(
        x = junctions_points[["x"]],
        y = junctions_points[["y"]],
        index = sr$index %>%
          rep(times = ncp) %>% # repeat these indexes the same number of ctrl points
          sort()
      )
    
    # add the intial start/end positions, since these are not return from the ctrl points
    junctions_points <-
      dplyr::tibble(
        x = c(sr$start, sr$stop),
        y = 0, # start from middle of exons
        index = rep(sr$index, 2)
      ) %>% 
      dplyr::bind_rows(junctions_points)
    
    # normalise y so values should always sit between -1 and 1
    # change, even junctions y values to negative - to be plotted on bottom
    junctions_points <- junctions_points %>%
      dplyr::mutate(
        y = y / max(y)
        #y = ifelse(index %% 2 == 0, -y, y)
      )
    
    # mark midpoints of junction curves to add label
    junctions_points <- junctions_points %>%
      dplyr::group_by(index) %>%
      dplyr::mutate(mid_point = .mark_mid(x)) %>%
      dplyr::ungroup()
    
    
    # add junction categories to colour by
    # add counts and tidy from wide into a long format for plotting
    junctions_points <- sr %>%
      dplyr::select(junID, countsSamples, type, index) %>%
      dplyr::left_join(junctions_points, ., by = "index") 
    
    junctions_points <- junctions_points %>%
      #dplyr::mutate(linetype = ifelse(count == 0, 2, 1)) %>%
      dplyr::arrange(index, x)
    
    
    return(junctions_points)
    
  }
  
}




.plot_add_junctions <- function(coords_to_plot, gene_er_pas_track_plot, junctions_points){
  
  if(is.null(junctions_points)){
    
    return(gene_er_pas_track_plot)
    
  }else{
    
    gene_er_pas_jun_plot <- gene_er_pas_track_plot +
      ggplot2::geom_path(
        data = junctions_points,
        ggplot2::aes(
          x = x, y = y,
          group = as.factor(index),
          colour = type
        ),
        lineend = "round", size=1
      ) +
      ggplot2::scale_size_continuous(
        range = c(0.2, 1.5),
        limits = c(0, 1),
        guide = "none"
      ) +
      ggplot2::scale_colour_manual(
        name = "Junction type",
        values = c("gray40", "#FF6666"),
        breaks = c("annotated", "unannotated"),
        labels = c("Partially annotated", "Unannotated"),
        drop = FALSE
      ) +
      ggplot2::theme(
        legend.key = ggplot2::element_rect(colour = "black", fill = "white")
      ) +
      # add track label
      ggplot2::annotate(geom="text",
                        x = (coords_to_plot[["x_max"]] + coords_to_plot[["x_min"]])/2 ,
                        y=1.2, 
                        label="Junction-read track", col = "blue", fontface = "bold", hjust = 0
      ) +
      # add label to the junction
      ggrepel::geom_label_repel(
        data = junctions_points %>% dplyr::filter(mid_point),
        ggplot2::aes(
          x = x, y = y,
          label = countsSamples,
          colour = type
        ),
        min.segment.length = 0,
        seed = 32,
        show.legend = FALSE,
        size = 3.5
      )
    
    
    return(gene_er_pas_jun_plot)
    
  }
  
}



.plot_annotation <- function(gene_er_pas_jun_plot, coords_to_plot) {
  gene_er_pas_jun_plot_annotated <- gene_er_pas_jun_plot %>%
    ggpubr::annotate_figure(
      top = ggpubr::text_grob(stringr::str_c(
        "Tissue: ", coords_to_plot[["er_tissue"]], ", Gene: ",
        coords_to_plot[["ensg_id"]], ", strand: ", coords_to_plot[["strand"]] ,
        ", ER-id: ", coords_to_plot[["er_id"]]
      ),
      x = 0.98, face = "italic", size = 10, just = "right"
      )
    )
  
  return(gene_er_pas_jun_plot_annotated)
  
}




