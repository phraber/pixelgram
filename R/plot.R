#' Renders the pixgram object
#'
#' This layout-compatible method validates the pixgram object, establishes plot coordinates, and draws the pixel and phylogram plots.
#'
#' The x-axis of the plot coordinate system is positive for the tree-plotting region, and negative for the pixel-plotting region.  Setting \code{raster_width} rescales the minimum x-axis value by this factor.  The default is 1, so the pixel and raster plots are equal width.  The \code{raster_margin} parameter rescales the right (least negative) edge of the pixel plot, to leave a gap between pixel and tree, e.g. to include a heatmap of phyenotypic assay results.
#'
#' The tree can be decorated by passing \code{leaf_colors}, \code{tip_labels}, \code{edge_widths}, and \code{edge_colors}.  See the last vignette example.
#'
#' If defined, data in the \code{pheno_matrix} are plotted in the margin between tree and pixel plots.  Entry values should be positive (non-zero) integers, which are interpreted as indices to \code{palette()}, color names, or color values.  Matrix rows must have names that match tree tip-label (i.e. \code{pg$tre$tip.label}  strings, or they will not be plotted.  This is to facilitate their proper ordering, which follows that of the leaves on the tree.  All columns are plotted, and in the order given.  However, column widths are scaled to fit together into the marginal area provided within the pixel-plot area (x-values < 0) by \code{raster_margin}.
#' 
#' The x- and y-axis lim arguments are provided to facilitate plotting multiple trees on the same scale in different panels.
#' 
#' The option \code{color_lut_type} is used to specify alternative columns from the defined color scheme (see \code{help(pixmap_colors)}).  For example, 'taylor' in the vignette uses the colors described in Wittgensteinean fashion by W. Taylor, Protein Engineering, Vol 10 , 743-746 (1997).
#' 
#' Entries in the \code{scale_bar} list argument should be named \code{'pos'} for plot-position coordinates, \code{'lwd'}, \code{'len'} for bar length, in terms of the plot coordinates, \code{'cex'}, and \code{'text'}.
#'
#' Entries in the \code{legend_attributes} list argument should be named \code{'location'} for legend positioning as defined in \code{?legend}, \code{'pchs'}, \code{'cex'}, \code{'pt.cex'}, \code{'lwds'}, \code{'bg'}, \code{'bty'}, \code{'box.col'}, \code{'cols'} for colors, \code{'title'}, and \code{'text'}.
#' 
#' @param pg pixgram object to be plotted
#' @param xform_type see \code{help(pixgramr::pixgram())}.
#' @param xform_master see \code{help(pixgramr::pixgram())}.
#' @param leaf_colors vector of tree leaf colors; can be colors or integers
#' @param tip_labels vector of tree tip labels
#' @param legend_attributes Named list of entries that go into the legend
#' @param edge_widths tree edge widths
#' @param edge_colors tree edge colors
#' @param scale_bar Named list of entries used to draw the scale bar
#' @param pheno_matrix heatmap to draw in raster_margin
#' @param show_top_axis Draw axis along top of plot?
#' @param show_tip_label Show taxon labels on tree tips?
#' @param selected_rows Vector of selected taxa
#' @param no_margin No margin?  This is passed to \code{ape::plot.phylo()}.
#' @param plot_margin_points Plot margin points?
#' @param vbars vertical bars?
#' @param raster_width
#' @param raster_margin
#' @param color_lut_type
#' @param annotate_env If true, draw Env landmarks.
#' @export
plot.pixgram <- function(P,
    xform_type=NULL,
    xform_master=NULL,
    leaf_colors=NULL,
    tip_labels=NULL,
    legend_attributes=NULL,
    edge_widths=NULL,
    edge_colors=NULL,
    scale_bar=NULL,
    pheno_matrix=NULL,
    show_top_axis=F,
    show_tip_label=F,
    selected_rows=NULL,
    no_margin=F,
    plot_margin_points=F,
    vbars=NULL,
    raster_width=NULL,
    raster_margin=NULL,
    color_lut_type=NULL,
#    x_lim=NULL,
#    y_lim=NULL,
    annotate_env=F,
    ...) { #, hbars=NULL) {

    if (class(P) != "pixgram")
	stop("plot.pixgram ERROR: Please specify pixgram object")

    P <- pixgram.validate(P)
    P <- pixgram.reorder(P)

    if (is.null(P$nt_master) & is.null(P$aa_master))
        P <- set.master(P)

    if (!is.null(raster_margin))
	P$raster_margin = raster_margin

    if (!is.null(raster_width))
	P$raster_width = raster_width

    if (!is.null(xform_type))
	P$xform_type = xform_type

    if (!is.null(xform_master))
	P$xform_master = xform_master

    if (!is.null(P$nt_master) | !is.null(P$aa_master))
        P <- pixgram.xform(P)

#    if (!is.null(x_lim)) P$x_lim = x_lim # should sanitize
#    if (!is.null(y_lim)) P$y_lim = y_lim

    R <- P
    class(R) <- "pixgram"

#    old.par = par()
#    old.pal = palette()

#    if (!is.null(cex))
#        par(cex=cex)

#    par(lend=1, ljoin=1)

#    if (no_margin) {
#        par(oma=c(0,0,0,0))#1,1,1,1))
#        par(oma=c(0,0,0,0))
#        par(mar=c(1,0,0/2,0))
#        par(xaxs='i', yaxs='i')
#        par(xaxs='r', yaxs='r')
#    } else {
#        par(oma=c(1/4,0,0,0))
#        par(mar=c(0,0,1,0))
#        par(xaxs='r', yaxs='r')
#    }

#    if (!is.null(R$palette)) palette(R$palette)

    R <- pixgram.tree.plot(R,
	leaf_colors=leaf_colors,
	tip_labels=tip_labels,
	legend_attributes=legend_attributes,
	edge_widths=edge_widths,
	edge_colors=edge_colors,
	no_margin=no_margin,
	pheno_matrix=pheno_matrix,
	plot_margin_points=plot_margin_points,
	show_tip_label=show_tip_label,
	selected_rows=selected_rows,
	scale_bar=scale_bar,
	...)#, hbars=hbars)

    if (is.null(color_lut_type)) {
	if (!is.null(R$aas)) 
	    R$color_lut_type = "aa" 
    } else {
	if (color_lut_type %in% colnames(pixmap_colors))
	    R$color_lut_type = color_lut_type
    }

    if (!is.null(R$aas))
        R <- pixgram.raster.aa(R, vbars=vbars, show_top_axis=show_top_axis)

    if (!is.null(R$nts))
	R$color_lut_type = "nt" 

    if (!is.null(R$nts) & !R$is_orf)
        R <- pixgram.raster.nt(R, vbars=vbars, show_top_axis=show_top_axis)
# to do: add axis options or ...

     if (annotate_env & !is.null(R$refseq_lut))
 	annotate.env(refseq_lut=R$refseq_lut,
 	    show_marginal_annotation=T,
 	    y_lim=c(R$y_lim[1] + ifelse(R$invert_y==T, -1/2, 1/2),
 		R$y_lim[2] - ifelse(R$invert_y==T, -1/2, 1/2)))

#    par(old.par)
#    palette(old.pal)

    R
}

#' @keywords internal
color.node <- function(edge_colors, this_tree, node) {

    ancestor = this_tree$edge[which(this_tree$edge[,2]==node), 1]

    if (length(this_tree$edge[which(this_tree$edge[,1] == ancestor), 2]) > 0) {

        mycolor = edge_colors[which(this_tree$edge[,2] == node)]

#	if (is.na(mycolor) | length(mycolor) == 0)
#	    stop(paste("pixgram ERROR in pixgram::color.node()!  No color specified for", ))

        sibs = this_tree$edge[which(this_tree$edge[,1] == ancestor & this_tree$edge[,2] != node), 2]

        # compare colors
        is_same_color = T

        if (length(sibs) > 0) {

            for (sib in c(1:length(sibs))) {

                if (edge_colors[which(this_tree$edge[,2] == sibs[sib])] != mycolor) {
                    if (edge_colors[which(this_tree$edge[,2] == sibs[sib])] != "#444444") { # default.edge.color
                        is_same_color = F
                    }
                }
            }
        }

        if (is_same_color) {
            edge_colors[which(this_tree$edge[,2] == ancestor)] = mycolor
        } else {
            edge_colors[which(this_tree$edge[,2] == ancestor)] = "#444444" # default.edge.color
        }
    }

    if (length(this_tree$edge[which(this_tree$edge[,2] == ancestor), 1]) > 0)
        edge_colors = color.node(edge_colors, this_tree, ancestor)

    edge_colors
}

#' @keywords internal
color.edges <- function(T, tip_colors) {

    edge_colors <- rep("#444444", length(T$edge)) # default.edge.color

    for (i in 1:length(T$tip.label))
        edge_colors[which(T$edge[,2] == i)] <- tip_colors[i]

    for (i in c(1:length(T$tip.label)))
        edge_colors <- color.node(edge_colors, T, i)

    edge_colors
}

#' @keywords internal
add.legend <- function(L) {

    if (!is.null(L$text)) {

        if (is.null(L$location)) L$location='topright'
	if (is.null(L$bty)) L$bty = 'n'
#	if (is.null(L$bg)) L$bg = 'white'
	if (is.null(L$pt.bg)) L$pt.bg = 'white'
	if (is.null(L$box_col)) L$box_col = NA
	if (is.null(L$cex)) L$cex = 1

        if (length(L$location) == 2) {

	    if (is.null(L$lwds)) {

		legend(L$location[1], L$location[2],
		    L$text, col=L$cols, cex=L$cex,
		    pch=L$pchs, bty=L$bty, bg=L$bg,
		    pt.bg=L$pt.bg,
		    box.col=L$box_col,
		    title=L$title, pt.cex=L$pt.cex)

	    } else if (is.null(L$pchs)) {

		legend(L$location[1], L$location[2],
		    L$text, col=L$cols, cex=L$cex,
                    lwd=L$lwds, bty=L$bty, bg=L$bg,
		    box.col=L$box_col,
		    pt.bg=L$pt.bg,
		    title=L$title, pt.cex=L$pt.cex)
	    } else {

		legend(L$location[1], L$location[2],
		   L$text, col=L$cols, cex=L$cex,
                   pch=L$pchs, lwd=L$lwds, bty=L$bty, bg=L$bg,
		   pt.bg=L$pt.bg,
		   box.col=L$box_col, title=L$title,
		   pt.cex=L$pt.cex,
		   merge=T)
	    }

	} else { # assume location is a character string?

	    if (is.null(L$lwds)) {

		legend(as.character(L$location), L$text, col=L$cols, cex=L$cex,
                    pch=L$pchs, bty=L$bty, bg=L$bg,
		    box.col=L$box_col,
		    title=L$title,
		    pt.bg=L$pt.bg,
		    pt.cex=L$pt.cex)
	    } else if (is.null(L$pchs)) {

		legend(as.character(L$location), L$text, col=L$cols, cex=L$cex,
                   lwd=L$lwds, bty=L$bty, bg=L$bg, box.col=L$box_col,
		   title=L$title, pt.cex=L$pt.cex)
	    } else {

		legend(as.character(L$location), L$text, col=L$cols, cex=L$cex,
                   pch=L$pchs, lwd=L$lwds, bty=L$bty, bg=L$bg,
		   pt.bg=L$pt.bg,
		   box.col=L$box_col, title=L$title,
		   pt.cex=L$pt.cex,
		   merge=T)
	    }
	}
    }
}

#' @keywords internal
pixgram.tree.plot <- function(P,
    leaf_colors=leaf_colors,
    tip_labels=tip_labels,
    edge_widths=edge_widths,
    edge_colors=edge_colors,
    legend_attributes=legend_attributes,
    no_margin=no_margin,
    pheno_matrix=pheno_matrix,
    plot_margin_points=plot_margin_points,
    show_tip_label=show_tip_label,
    selected_rows=selected_rows,
    scale_bar=scale_bar,
    ...) {

    if (class(P) != "pixgram")
	stop("pixgram.tree.plot ERROR: Please specify pixgram object")

    R <- P
    message("*** pixgram.tree.plot ***\n")
    class(R) <- "pixgram"

    if (is.null(R$x_lim) & is.null(R$y_lim)) {
	tree_out <- ape::plot.phylo(R$tre, font=1, 
	    no.margin=no_margin, plot=F, ...)
    } else if (!is.null(R$x_lim) &  is.null(R$y_lim)) {
	tree_out <- ape::plot.phylo(R$tre, font=1, 
	    no.margin=no_margin, x.lim=R$x_lim, plot=F, ...)
    } else if ( is.null(R$x_lim) & !is.null(R$y_lim)) {
	tree_out <- ape::plot.phylo(R$tre, font=1, 
	    no.margin=no_margin, y.lim=R$y_lim, plot=F, ...)
    } else {
	tree_out <- ape::plot.phylo(R$tre, font=1, 
	    no.margin=no_margin, x.lim=R$x_lim, y.lim=R$y_lim, plot=F, ...)
    }

    ### why not [1]?
#    R$x_lim <- c(-1.025* R$raster_width * tree_out$x.lim[2], tree_out$x.lim[2])
    R$x_lim <- c(-R$raster_width * tree_out$x.lim[2], tree_out$x.lim[2])
    R$y_lim <- c(max(tree_out$y.lim), min(tree_out$y.lim))
#    R$y_lim <- c(length(R$tre$tip.label), 1)

    if (no_margin)
	R$y_lim <- c(max(tree_out$y.lim) + 1, min(tree_out$y.lim) - 1)
#	R$y_lim <- c(1+length(R$tre$tip.label), 0) # or 1/2?
#	R$y_lim <- c(1+length(R$tre$tip.label), 0) # or 1/2?

    if (R$invert_y == T)
	R$y_lim <- rev(R$y_lim)

    message(paste("tree x_lim: ", R$x_lim[1], R$x_lim[2]))
    message(paste("tree y_lim: ", R$y_lim[1], R$y_lim[2]))

    plot_size <- dev.size("px")
    R$my_cex <- 5/4*(plot_size[2]/length(R$tre$tip.label))/12

    if (!is.null(leaf_colors)) {

        if (length(which(is.na(leaf_colors))) > 0)
            stop(paste("pixgram.tree.plot() ERROR! Undefined leaf_colors:",
		paste(R$tre$tip.label[which(is.na(leaf_colors))], collapse=" ")))

	if (is.null(edge_colors))
            edge_colors <- color.edges(R$tre, leaf_colors)

    } else {
	if (is.null(edge_colors))
            edge_colors = "black" # default.edge.color
    }

    if (is.null(edge_widths))
        edge_widths = 1

  if (!is.null(selected_rows)) {
    if (R$invert_y == T) {

	for (i in 1:length(selected_rows))
	  rect(R$x_lim[1],
	    selected_rows[i]+1/2,
            R$x_lim[2],
	    selected_rows[i]-1/2,
	    border=NA, col=6)

    } else {
# the default is not to invert
	if (no_margin) {
#	R$y_lim <- c(1+length(R$tre$tip.label), 0)

#	    stop("Sorry, this is not working correctly")
	    x_factor <- (1 + R$y_lim[1] - R$y_lim[2]) / length(R$tre$tip.label)

	  for (i in 1:length(selected_rows))
	    rect(R$x_lim[1],
		length(R$tre$tip.label)-selected_rows[i]/x_factor + 1/2,
		R$x_lim[2],
		length(R$tre$tip.label)-selected_rows[i]/x_factor - 1/2,
		border=NA, col=6)
	} else {
#       R$y_lim <- c(length(R$tre$tip.label), 1)
	  for (i in 1:length(selected_rows))
	    rect(R$x_lim[1],
		length(R$tre$tip.label)-selected_rows[i]+1/2,
		R$x_lim[2],
		length(R$tre$tip.label)-selected_rows[i]+3/2,
		border=NA, col=6)
	}
    }
  }
#    show_tip_label = T #ifelse(is.null(tip_labels), T, F)
    safe.tiplabels = R$tre$tip.label
    R$tre$tip.label = sapply(1:length(R$tre$tip.label), function(i)
        paste("  ", R$tre$tip.label[i]))
    par(new=T)

    tree_out <- ape::plot.phylo(R$tre,
                                font=1,
                                no.margin=no_margin,
                                cex=R$my_cex,
                                x.lim=R$x_lim,
                                y.lim=R$y_lim,
                                plot=R$show_tree,
				edge.width=edge_widths,
                                edge.color=edge_colors,
                                tip.color=tip_labels$col, # NOTE THIS IS NOT ROBUST! assumes tree order is unchanged
                                show.tip.label=show_tip_label, ...) # 

    R$tre$tip.label = safe.tiplabels

        # design issue here is how to map leaf names to { pch, bg, col }
    if (!is.null(tip_labels) & !is.null(tip_labels$pchs)) {

        if (!is.null(tip_labels) & length(which(is.na(tip_labels$pchs))) > 1)
	    stop(paste("pixgram ERROR! Undefined tip_labels$pchs:",
		paste(R$tre$tip.label[which(is.na(tip_labels$pchs))],
		    collapse=" ")))

	if (!is.null(tip_labels) & length(which(is.na(tip_labels$col))) > 1)
	    stop(paste("pixgram ERROR! Undefined tip_labels$col:",
		paste(R$tre$tip.label[which(is.na(tip_labels$col))],
		    collapse=" ")))

        ape::tiplabels(pch=tip_labels$pchs,
                  col=tip_labels$col,
                  bg=tip_labels$bgs,
                  cex=R$my_cex, lwd=1/2)

# to do: position properly with raster_margin considered:
	if (plot_margin_points) {

	    if ((!is.null(leaf_colors) & is.null(tip_labels$col)) |
		(!is.null(tip_labels$col) & length(which(tip_labels$col != leaf_colors)) > 0))
		    for (i in 1:length(R$tre$tip.label))
                    points(-1.01 * R$raster_width * R$x_lim[2], i, col=leaf_colors[i], cex=R$my_cex, lwd=1/2)

	    if (!is.null(tip_labels$col)) {

	        if (is.null(tip_labels$pchs)) {
		    for (i in 1:length(R$tre$tip.label))
		        points(-0.02 * R$raster_width * R$x_lim[2], i, col=tip_labels$col[i], cex=R$my_cex, lwd=1/2)
		} else {
		    if (length(tip_labels$bgs)==1)
			tip_labels$bgs = rep(tip_labels$bgs, length(tip_labels$pchs))

		    for (i in 1:length(R$tre$tip.label))
		        points(-0.01 * R$raster_width * R$x_lim[2], i,
			    col=tip_labels$col[i], pch=tip_labels$pchs[i], bg=tip_labels$bgs[i], cex=R$my_cex, lwd=1/2)
			    # bgs could be a scalar or vector
		}
	    }
        }

        ### HERE we plot the phenotypic/binding/neutralization data
        # should this really be nested within if (!is.null(tip_labels))?

	if (!is.null(pheno_matrix)) {

	    x_a <- R$raster_width * R$x_lim[2] * R$raster_margin
	    x_b <- 0

	    for (col_num in c(1:ncol(pheno_matrix))) {

                # transform col_num onto raster_margin interval
		L_pos = -R$raster_width * R$x_lim[2] * R$raster_margin * (1 - (col_num-1)/(1+ncol(pheno_matrix)))
		R_pos = -R$raster_width * R$x_lim[2] * R$raster_margin * (1 - col_num/(1+ncol(pheno_matrix)))

#		if (R$invert_y == T) {

		myordering <- R$tre$edge[which(R$tre$edge[, 2] <= length(R$tre$tip.label)), 2]

		rect(rep(L_pos, length(R$tre$tip.label)),
			c(length(R$tre$tip.label):1)-1/2,
			rep(R_pos, length(R$tre$tip.label)),
			c(length(R$tre$tip.label):1)+1/2,
#			col=pheno_matrix[, col_num],
		    col=pheno_matrix[myordering[c(length(R$tre$tip.label):1)], col_num],
			border=NA)

#		} else {

#		    rect(rep(L_pos, length(R$tre$tip.label)),
#			c(length(R$tre$tip.label):1)-1/2,
#			rep(R_pos, length(R$tre$tip.label)),
#			c(length(R$tre$tip.label):1)+1/2,
#			col=pheno_matrix[, col_num],
##		    col=pheno_matrix[myordering[c(length(R$tre$tip.label):1)], col_num],
#			border=NA)
#		}
#    mtext(letters[col_num], 1, at=myxpos, line=-1/3, cex=cexS, adj=0.5, padj=1/2)
	    }
	} else {

	    if (plot_margin_points) {

	        if ((!is.null(leaf_colors) & is.null(tip_labels$col)) |
		    (!is.null(tip_labels$col) & length(which(tip_labels$col != leaf_colors)) > 0))
	                for (i in 1:length(R$tre$tip.label))
                            points(-0.01 * R$raster_width * R$x_lim[2], i, col=leaf_colors[i], cex=R$my_cex, lwd=1/2)

	        if (!is.null(tip_labels$col)) {

	            if (is.null(tip_labels$pchs)) {
			for (i in 1:length(R$tre$tip.label))
		            points(-0.02 * R$raster_width * R$x_lim[2], i, col=tip_labels$col[i], cex=R$my_cex, lwd=1/2)
		    } else {
			if (length(tip_labels$bgs)==1)
			    tip_labels$bgs = rep(tip_labels$bgs, length(tip_labels$pchs))

			for (i in 1:length(R$tre$tip.label))
		            points(-0.01 * R$raster_width * R$x_lim[2], i,
				col=tip_labels$col[i], pch=tip_labels$pchs[i], bg=tip_labels$bgs[i], cex=R$my_cex, lwd=1/2)
			    # bgs could be a scalar or vector
		    }
		}
	    }
        }
    }

    if (!is.null(legend_attributes)) {

	if (is.null(legend_attributes$cex))
	    legend_attributes$cex = 1

	if (is.null(legend_attributes$pt.cex))
	    legend_attributes$pt.cex = 4/5*R$my_cex

        add.legend(legend_attributes)
    }

    if (!is.null(scale_bar)) {

	sb.xpos = ifelse(is.null(scale_bar$pos[1]), NULL, scale_bar$pos[1])
	sb.ypos = ifelse(is.null(scale_bar$pos[2]), NULL, scale_bar$pos[2])

	ape::add.scale.bar(sb.xpos, sb.ypos, length=scale_bar$len,
	    lwd=scale_bar$lwd, cex=scale_bar$cex/100)

	if (R$invert_y == T) {
	    text(sb.xpos+scale_bar$len/2, sb.ypos-1, scale_bar$text,
		cex=scale_bar$cex, adj=c(1/2, 1))
	} else {
	    text(sb.xpos+scale_bar$len/2, sb.ypos+1, scale_bar$text,
		cex=scale_bar$cex, adj=c(1/2, 1))
	}
    }

    R$x_lim <- c(tree_out$x.lim[1], tree_out$x.lim[2])

    R
}

#' @keywords internal
pixgram.raster.nt <- function(P, vbars, show_top_axis) {

    if (class(P) != "pixgram")
	stop("pixgram.raster.nt ERROR: Please specify pixgram object")

    R <- P
    class(R) <- "pixgram"

    if (is.null(R$color_lut_type))
	R$color_lut_type = "aa"

    if (!is.null(R$nt_rast)) {

        message("*** Making nt raster map ***")

        # is it really necessary to use nested 'for' loops here?
        for (j in 1:ncol(R$nt_rast)) {

            col_alphabet <- unique(c(R$nt_rast[, j]))

            for (i in 1:length(col_alphabet))
                R$nt_rast[, j] = replace(R$nt_rast[, j],
                    which(R$nt_rast[, j] == col_alphabet[i]),
		        pixmap_colors[col_alphabet[i], R$color_lut_type])
        }

        rasterImage(as.raster(R$nt_rast),
            -R$raster_width * R$x_lim[2],
              1/2+nrow(R$nts),
            -R$raster_width * R$x_lim[2] * R$raster_margin,
              1/2, interpolate=F)
    }

#par(yaxt='s', xaxt='s')
#	axis(2, line=-5, yaxt='s')
#	axis(2, line=1, yaxt='s')
#	axis(2, line=0, yaxt='s')
#	axis(2, line=-1, yaxt='s')

    if (is.null(R$aas))
        rect(-R$raster_width * R$x_lim[2], 1/2+nrow(R$nt_rast),
             -R$raster_width * R$x_lim[2] * R$raster_margin, 1/2,
             border='black', col=NA, lwd=1/2)

    R
}


#' @keywords internal
annotate.env <- function(refseq_lut=NULL,
    show_marginal_annotation=FALSE,
    y_lim=NULL) {

    message("*** annotate.env ***\n")

#    if (!is.null(plot_annotation)) {
#	Lannotation <- plot_annotation$Lannotation
#	Rannotation <- plot_annotation$Rannotation
#	 annotation <- plot_annotation$annotation
#    } else {

#	if (!is.null(region) & grepl("^gp", region)) {
	    Lannotation <- c(132, 185, 276, 363, 396, 459)
	    Rannotation <- c(152, 190, 283, 373, 410, 466)
	    Nannotation <- c('V1', 'V2', 'Loop D', 'CD4 Loop', 'V4', 'V5')
#	} else if (!is.null(region) & grepl("^V", region)) {
#	    Lannotation <- c(26, 50, 98)
#	    Rannotation <- c(35, 65, 110)
#	     annotation <- c('H1', 'H2', 'H3')
#	} else {
#	    Lannotation <- NULL
#	    Rannotation <- NULL
#	     annotation <- NULL
#	}
#    }

    if (is.null(y_lim))
	y_lim = (usr <- par('usr'))[c(3,4)]

    if (!is.null(refseq_lut)) {

 	for (x in 1:length(Nannotation)) {

 	    x_1 = NULL
	    x_1 <- refseq_lut$aln[min(which(refseq_lut$l==Lannotation[x]))]

 	    x_2 = NULL
	    x_2 <- refseq_lut$aln[max(which(refseq_lut$r==Rannotation[x]))]

 	    if (!is.null(x_1) & !is.null(x_2) & !is.na(x_1) & !is.na(x_2))
 		rect(x_1, y_lim[1], x_2, y_lim[2], border=NA, col="grey")
 	}

        # gp120/gp41 boundary
	x_pos = NULL
	x_pos <- mean(refseq_lut$aln[min(which(refseq_lut$l==511))],
	    refseq_lut$aln[max(which(refseq_lut$r==512))])

	if (!is.null(x_pos) & !is.na(x_pos))
	    segments(x_pos, y_lim[1], x_pos, y_lim[2], col="grey", lwd=1, lty=2)

        # signal peptide
	x_pos = NULL
	x_pos <- mean(refseq_lut$aln[which(refseq_lut$l==30&refseq_lut$r==30)],
	              refseq_lut$aln[which(refseq_lut$l==31&refseq_lut$r==31)])

	if (!is.null(x_pos) & !is.na(x_pos))
	    segments(x_pos, y_lim[1], x_pos, y_lim[2], col="grey", lwd=1, lty=2)

#	if (!is.null(x_pos) & !is.na(x_pos)) 
#	    abline(v=x_pos, col=5, lwd=1, lty=3)

        ### add margin text and legend if top plot
 	if (show_marginal_annotation) {

 	    for (X in 1:length(Nannotation)) {

 		l.xpos=NULL
		l.xpos <- refseq_lut$aln[min(which(refseq_lut$l==Lannotation[X]))]

 		r.xpos=NULL
		r.xpos <- refseq_lut$aln[max(which(refseq_lut$r==Rannotation[X]))]

 		my_padj <- ifelse(X < 4, 1/2, ifelse(X==4, 1, 0))

 		if (!is.null(r.xpos) & !is.null(r.xpos))
 		    mtext(Nannotation[X], 3, at=(l.xpos + r.xpos)/2,
 			line=1, adj=1/2, padj=my_padj, cex=8/12) # not really marginal, is it?

 	    }
        }
    }

    mtext("Site", side=1, outer=T, line=6/4)

}

#' @keywords internal
pixgram.raster.aa <- function(P, vbars, show_top_axis) {

    if (class(P) != "pixgram")
	stop("pixgram.raster.aa ERROR: Please specify pixgram object")

    R <- P
    class(R) <- "pixgram"

    if (!is.null(R$refseq_lut)) {

        axis_xlim = c(-R$raster_width * R$x_lim[2], -R$raster_width * R$x_lim[2] * R$raster_margin)

        my_slope = (ncol(R$aas)-1)/(axis_xlim[2]-axis_xlim[1])

 	my_refseq_lut = R$refseq_lut
 	my_refseq_lut$aln = (my_refseq_lut$aln - 1)/my_slope - R$raster_width * R$x_lim[2]

# #	if (!is.null(vbars)) {
# #	    axis_xlim = c(-R$raster_width * R$x_lim[2], -R$raster_width * R$x_lim[2] * 0.025)
# #            my_slope = (ncol(R$aas)-1)/(axis_xlim[2]-axis_xlim[1])
# #	    x_locs = (vbars-1)/my_slope - R$raster_width * R$x_lim[2]
# #	    abline(v=vbars, col='#66666666', lwd=1/2)
# #	}
    }

    if (!is.null(R$aa_rast)) {

        message("*** Making aa raster map ***")

	if (is.null(R$color_lut_type))
	    R$color_lut_type = "aa"

	if (R$xform_type==3)
	    R$color_lut_type = "charge"

#        # is it really necessary to use nested 'for' loops here?
        for (j in 1:ncol(R$aa_rast)) {

            col_alphabet <- unique(c(R$aa_rast[, j]))

#	    R$aa_rast[i, ] <- sapply(1:length(row_aa_alphabet), function(j)
            for (i in 1:length(col_alphabet))
		R$aa_rast[, j] = replace(R$aa_rast[, j],
                    which(R$aa_rast[, j] == col_alphabet[i]),
                        pixmap_colors[col_alphabet[i], R$color_lut_type])
        }

        rasterImage(as.raster(R$aa_rast),
            -R$raster_width * R$x_lim[2],
             1/2+nrow(R$aas),
            -R$raster_width * R$x_lim[2] * R$raster_margin, #0.025,
             1/2, interpolate=F)
    }

    rect(-R$raster_width * R$x_lim[2], 1/2 + nrow(R$aas),
         -R$raster_width * R$x_lim[2] * R$raster_margin, # 0.025,
	 1/2,
         border='black', col=NA, lwd=1/2)

  # new code: FOR NOW ONLY WORKS WITH AA ALIGNMENTS
#-R$raster_width * R$x_lim[2] * raster_margin
# ncol(R$aas) is number of aligned sites -> min(R$refseq_lut$aln),

    if (!is.null(R$refseq_lut)) {

         axis_xlim = c(-R$raster_width * R$x_lim[2], -R$raster_width * R$x_lim[2] * R$raster_margin)#0.025)

         my_slope = (ncol(R$aas)-1)/(axis_xlim[2]-axis_xlim[1])


     tick_interval = ifelse(ncol(R$aas) >= 3000, 1000,
 	ifelse(ncol(R$aas) >= 300, 100,
 	    ifelse(ncol(R$aas) >= 200, 50,
 		ifelse(ncol(R$aas) >= 150, 30,
 		    ifelse(ncol(R$aas) >= 80, 20,
 			ifelse(ncol(R$aas) >= 40, 10,
 			    ifelse(ncol(R$aas) >= 20, 5, 1)))))))

 	x_locs <- R$refseq_lut$aln[which(R$refseq_lut$aln < ncol(R$aas) &
 	    R$refseq_lut$l == R$refseq_lut$r & R$refseq_lut$l %% tick_interval == 0)]

 	x_locs = (x_locs-1)/my_slope - R$raster_width * R$x_lim[2]
# #	axis_xlim = c(1, ncol(R$aas)) / my_slope# - R$raster_width * R$x_lim[2]

         x_labs <- R$refseq_lut$l[which(R$refseq_lut$l == R$refseq_lut$r &
 		R$refseq_lut$l %% tick_interval == 0)]

### to do: think about what to do when my_lim %in% x_locs
 	my_cexl=4/5
 	my_cexa=2/3

#If either line or pos is set, they (rather than par("mgp")[3])
#determine the position of the axis line and tick marks, and the tick
#labels are placed par("mgp")[2] further lines into (or towards for
#pos) the margin.

 	my_mgp=c(2/3, 0/5, 0)
 	if (length(x_locs) == length(x_labs)) {

# 	message(paste("* xlabs=", paste(x_labs, collapse=" ")))
# 	message(paste("* xlocs=", paste(x_locs, collapse=" ")))
# 	message(paste("* ylims=", paste(R$y_lim, collapse=" ")))

#        axis(1, at=c(axis_xlim[1], x_locs, axis_xlim[2]), labels=c("", x_labs, ""),
         axis(1, at=x_locs, labels=x_labs,
 		xaxt='s', cex.lab=my_cexl, cex.axis=my_cexa, mgp=my_mgp,
#		pos=1/2,
 		pos=ifelse(R$invert_y==T, R$y_lim[1]-1/2, R$y_lim[1]+1/2),
 		padj=0/2,
 #		hadj=0/2,
 		lwd=1/2,
 		tck=-1/nrow(R$aas))

 	my_mgp=c(2/3, 0/5, 0)
#        axis(3, at=c(axis_xlim[1], x_locs, axis_xlim[2]), labels=c("", x_labs, ""),
 	if (R$invert_y == F & show_top_axis)
             axis(3, at=x_locs, labels=x_labs,
 		xaxt='s', cex.lab=my_cexl, cex.axis=my_cexa, mgp=my_mgp,
 		pos=1/2, padj=0, lwd=1/2, tck=-1/nrow(R$aas))
 	} else {
 	    message(paste("* NO xlabs=", paste(x_labs, collapse=" ")))
 	    message(paste("* NO xlocs=", paste(x_locs, collapse=" ")))
 	}
     }

    #if (!is.null(R$ptid))
	#mtext(paste0(R$ptid, " "), 3, at=0, line=0, cex=1, outer=F, adj=1)

#    if (!is.null(R$region))
#	mtext(R$region, 1, at=mean(axis_xlim), line=1/2, cex=my_cexl)
#	mtext(R$region, 3, at=0, line=0, cex=1, outer=F, adj=1/2, font=3)

    R
}