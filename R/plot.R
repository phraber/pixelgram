#' Renders the pixelgram object
#'
#' This layout-compatible method validates the pixelgram object, establishes plot coordinates, and draws the pixel and phylogram plots.
#'
#' The x-axis of the plot coordinate system is positive for the tree-plotting region, and negative for the pixel-plotting region.  Setting \code{raster_width} rescales the minimum x-axis value by this factor.  The default is 1, so the pixel and raster plots are equal width.  The \code{raster_margin} parameter rescales the right (least negative) edge of the pixel plot, to leave a gap between pixel and tree, e.g. to include a heatmap of phyenotypic assay results.
#'
#' Tree tip labels can be added by passing a list called \code{tip_labels}, with elements named "pch", "bgs", and "col".  See vignette.  The tree can be further decorated by passing \code{leaf_colors}, \code{edge_widths}, and \code{edge_colors}.  See the last vignette example.
#'
#' If defined, data in the \code{pheno_matrix} are plotted in the margin between tree and pixel plots.  Entry values should be positive (non-zero) integers, which are interpreted as indices to \code{palette()}, color names, or color values.  Matrix rows must have names that match tree tip-label (i.e. \code{pg$tre$tip.label}  strings, or they will not be plotted.  This is to facilitate their proper ordering, which follows that of the leaves on the tree.  All columns are plotted, and in the order given.  However, column widths are scaled to fit together into the marginal area provided within the pixel-plot area (x-values < 0) by \code{raster_margin}.
#'
#' The x- and y-axis lim arguments are provided to facilitate plotting multiple trees on the same scale in different panels.
#'
#' The option \code{color_lut_type} is used to specify alternative columns from the defined color scheme (see \code{help(pixmap_colors)}).  For example, 'taylor' in the vignette uses the colors described in Wittgensteinean fashion by W. Taylor, Protein Engineering, Vol 10 , 743-746 (1997).
#'
#' Entries in the \code{scale_bar} list argument should be named \code{'pos'} for plot-position coordinates, \code{'lwd'}, \code{'len'} for bar length, in terms of the plot coordinates, \code{'cex'}, and \code{'text'}.
#'
#' Entries in the \code{legend_attributes} list argument should be named \code{'location'} for legend positioning as defined in \code{?legend}, \code{'pch'}, \code{'cex'}, \code{'pt.cex'}, \code{'lwd'}, \code{'bg'}, \code{'bty'}, \code{'box.col'}, \code{'col'} for colors, \code{'title'}, and \code{'text'}.
#'
#' @param x pixelgram object to be plotted
#' @param xform_type see \code{help(pixelgram::pixelgram())}.
#' @param xform_master see \code{help(pixelgram::pixelgram())}.
#' @param leaf_colors vector of tree leaf colors; can be colors or integers
#' @param legend_attributes Named list of entries that go into the legend
#' @param edge_widths tree edge widths
#' @param edge_colors tree edge colors
#' @param scale_bar Named list of entries used to draw the scale bar
#' @param pheno_matrix Heatmap of RGB values to draw in raster_margin.
#' @param pheno_colnames A vector of names to label pheno_matrix.
#' @param pheno_letters A vector of letters to label pheno_matrix.
#' @param show_top_axis Draw axis along top of plot?
#' @param show_tip_label Show taxon labels on tree tips?
#' @param selected_rows Vector of selected taxa.
#' @param selected_row_colors Colors to depict vector of selected taxa.
#' @param no_margin No margin?  This is passed to \code{ape::plot.phylo()}.
#' @param tip_labels List for adding tree tip labels; see details.
#' @param plot_margin_points Plot margin points?
#' @param vbars vertical bars?
#' @param raster_width Proportionate width of the pixel plot, e.g. 1 indicates equal width.
#' @param raster_margin Offset in raster portion of the layout, e.g. for heatmap.
#' @param color_lut_type Optional string for non-default color lookup table, currently only implemented for amino acids as 'charge' and 'taylor'.
#' @param x_labs Vector of values to plot as labels along x axis of raster plot, e.g. to identify specific sites.
#' @param notes Is a named list used to draw landmarks for the region sequenced, with elements "Lhs", "Rhs", "col" and "txt".  Using numbering from the reference sequence lookup-table, rectangles are drawn bounded by Lhs and Rhs for each entry and filled with the color specified by col (using alpha transparency gives better results).  If present, labels with values in "txt" are plotted along the margin.  If annotate_env is true, these values are populated internally for the HIV-1 env.
#' @param annotate_env If true, draw Env landmarks.
#' @param show_tree Show the tree? If not, only show the pixel plot.
#' @param main Plot title, e.g. coded subject id.  If given, it is plotted at the top of the pixelgram output.
#' @param sub Plot subtitle, i.e. region sequenced.  If specified, it appears with "site" below the pixel plot.
#' @param ... dots, passed to ape::plot.phylo()
#'
#' @examples
#' \dontrun{
#'   plot( pg <- pixelgram::pixelgram(nts=pixelgram::hiv.ref$nts) )
#' }
#'
#' @export
plot.pixelgram <- function(x,
    xform_type=NULL,
    xform_master=NULL,
    leaf_colors=NULL,
    tip_labels=NULL,
    legend_attributes=NULL,
    edge_widths=NULL,
    edge_colors=NULL,
    scale_bar=NULL,
    pheno_matrix=NULL,
    pheno_colnames=NULL,
    pheno_letters=NULL,
    show_top_axis=F,
    show_tip_label=F,
    selected_rows=NULL,
    selected_row_colors=NULL,
    no_margin=F,
    plot_margin_points=F,
    vbars=NULL,
    raster_width=NULL,
    raster_margin=NULL,
    color_lut_type=NULL,
x_labs=NULL,
#    x_lim=NULL,
#    y_lim=NULL,
    notes=NULL,
    annotate_env=F,
    show_tree=NULL,
    main=NULL,
    sub=NULL,
    ...) { #, hbars=NULL) {

    dots = list(...)

    if (class(x) != "pixelgram")
	stop("plot.pixelgram ERROR: Please specify pixelgram object")

    x <- try(pixelgram.validate(x))
    x <- try(pixelgram.reorder(x))

    if (is.null(x$nt_master) & is.null(x$aa_master))
        x <- set.master(x)

    if (!is.null(raster_margin))
	x$raster_margin = raster_margin

    if (!is.null(raster_width))
	x$raster_width = raster_width

    if (!is.null(xform_type))
	x$xform_type = xform_type

    if (!is.null(xform_master))
	x$xform_master = xform_master

    if (!is.null(x$nt_master) | !is.null(x$aa_master))
        x <- pixelgram.xform(x)

#    if (!is.null(x_lim)) x$x_lim = x_lim # should sanitize
#    if (!is.null(y_lim)) x$y_lim = y_lim

    R <- x
    class(R) <- "pixelgram"

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
    if (!is.null(show_tree)) R$show_tree = show_tree
    if (!is.null(main)) R$main = main
    if (!is.null(sub)) R$sub = sub
    if (!is.null(selected_rows) & is.null(selected_row_colors))
        selected_row_colors = rep("#66666666", length(selected_rows))

      R <- pixelgram.tree.plot(R,
                             leaf_colors=leaf_colors,
                             tip_labels=tip_labels,
                             legend_attributes=legend_attributes,
                             edge_widths=edge_widths,
                             edge_colors=edge_colors,
                             no_margin=no_margin,
                             pheno_matrix=pheno_matrix,
			     pheno_colnames=pheno_colnames,
			     pheno_letters=pheno_letters,
                             plot_margin_points=plot_margin_points,
                             show_tip_label=show_tip_label,
                             selected_rows=selected_rows,
                             selected_row_colors=selected_row_colors,
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
        R <- pixelgram.raster.aa(R, vbars=vbars, show_top_axis=show_top_axis, x_labs=x_labs)

    if (!is.null(R$nts))
	R$color_lut_type = "nt"

    if (!is.null(R$nts) & !R$is_orf)
        R <- pixelgram.raster.nt(R, vbars=vbars, show_top_axis=show_top_axis)
# to do: add axis options or ...

# the second condition/s implicitly assume env with hxb2 present
    if ((!is.null(notes) & is.list(notes)) |
	(annotate_env & !is.null(R$refseq_lut)))
 	    annotate.region(R, notes=notes,
 	        y_lim=c(R$y_lim[1] + ifelse(R$invert_y, -1/2, 1/2),
 		    R$y_lim[2] - ifelse(R$invert_y, -1/2, 1/2)))

#    par(old.par)
#    palette(old.pal)

    R
}

#' @keywords internal
color.node <- function(edge_colors, this_tree, node) {

    default.color = ifelse(is.numeric(edge_colors), 1, "#444444")

    ancestor = this_tree$edge[which(this_tree$edge[,2]==node), 1]

    if (length(this_tree$edge[which(this_tree$edge[,1] == ancestor), 2]) > 0) {

        mycolor = edge_colors[which(this_tree$edge[,2] == node)]

#	if (is.na(mycolor) | length(mycolor) == 0)
#	    stop(paste("pixelgram ERROR in pixelgram::color.node()!  No color specified for", ))

        sibs = this_tree$edge[which(this_tree$edge[,1] == ancestor & this_tree$edge[,2] != node), 2]

        # compare colors
        is_same_color = T

        if (length(sibs) > 0) {

            for (sib in c(1:length(sibs))) {

                if (edge_colors[which(this_tree$edge[,2] == sibs[sib])] != mycolor) {
                    if (edge_colors[which(this_tree$edge[,2] == sibs[sib])] != default.color) { # default.edge.color
                        is_same_color = F
                    }
                }
            }
        }

        if (is_same_color) {
            edge_colors[which(this_tree$edge[,2] == ancestor)] = mycolor
        } else {
            edge_colors[which(this_tree$edge[,2] == ancestor)] = default.color # default.edge.color
        }
    }

    if (length(this_tree$edge[which(this_tree$edge[,2] == ancestor), 1]) > 0)
        edge_colors = color.node(edge_colors, this_tree, ancestor)

    edge_colors
}

#' Recursively propagate tip colors internally to color interior branches
#'
#' Branches with descendents that have different colors will be represented with the default color.
#'
#' NB: the default color is currently hard-coded, which seems unwise.
#' Some trees may not be amenable to the recursion scheme used here.
#' Alternatively, starting to label on the root node may be a bad idea.
#' One tree was reported to have caused a problem, so it is worth keeping in mind.
#'
#' @param T ape::phylo() tree object.
#' @param tip_colors vector of colors, whether rgb hex or integers, whose order matches the names in the tree tip.label slot.
#'
#' @return A vector of colors of length equal to the number of tree branches.
#'
#' @export
color.edges <- function(T, tip_colors) {

    default.color = ifelse(is.numeric(tip_colors), 1, "#444444")
    edge_colors <- rep(default.color, length(T$edge)) # default.edge.color

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
        if (is.null(L$box.col)) L$box.col = NA
        if (is.null(L$cex)) L$cex = 1

        if (length(L$location) == 2) {

            if (is.null(L$lwd)) {

                legend(L$location[1], L$location[2],
                       L$text, col=L$col, cex=L$cex,
                       pch=L$pch, bty=L$bty, bg=L$bg,
                       pt.bg=L$pt.bg,
                       box.col=L$box_col,
                       title=L$title, pt.cex=L$pt.cex)

            } else if (is.null(L$pch)) {

                legend(L$location[1], L$location[2],
                       L$text, col=L$col, cex=L$cex,
                       lwd=L$lwd, bty=L$bty, bg=L$bg,
                       box.col=L$box.col,
                       pt.bg=L$pt.bg,
                       title=L$title, pt.cex=L$pt.cex)
            } else {

                legend(L$location[1], L$location[2],
                       L$text, col=L$col, cex=L$cex,
                       pch=L$pch, lwd=L$lwd, bty=L$bty, bg=L$bg,
                       pt.bg=L$pt.bg,
                       box.col=L$box.col, title=L$title,
                       pt.cex=L$pt.cex,
                       merge=T)
            }

        } else { # assume location is a character string?

            if (is.null(L$lwd)) {

                legend(as.character(L$location), L$text, col=L$col, cex=L$cex,
                       pch=L$pch, bty=L$bty, bg=L$bg,
                       box.col=L$box.col,
                       title=L$title,
                       pt.bg=L$pt.bg,
                       pt.cex=L$pt.cex)
            } else if (is.null(L$pch)) {

                legend(as.character(L$location), L$text, col=L$col, cex=L$cex,
                       lwd=L$lwd, bty=L$bty, bg=L$bg, box.col=L$box.col,
                       title=L$title, pt.cex=L$pt.cex)
            } else {

                legend(as.character(L$location), L$text, col=L$col, cex=L$cex,
                       pch=L$pch, lwd=L$lwd, bty=L$bty, bg=L$bg,
                       pt.bg=L$pt.bg,
                       box.col=L$box.col, title=L$title,
                       pt.cex=L$pt.cex,
                       merge=T)
            }
        }
    }
}

#' @keywords internal
pixelgram.tree.plot <- function(x,
    leaf_colors=leaf_colors,
    tip_labels=tip_labels,
    edge_widths=edge_widths,
    edge_colors=edge_colors,
    legend_attributes=legend_attributes,
    no_margin=no_margin,
    pheno_matrix=pheno_matrix,
    pheno_colnames=pheno_colnames,
    pheno_letters=pheno_letters,
    plot_margin_points=plot_margin_points,
    show_tip_label=show_tip_label,
    selected_rows=selected_rows,
    selected_row_colors=selected_row_colors,
    scale_bar=scale_bar,
    ...) {

    if (class(x) != "pixelgram")
        stop("pixelgram.tree.plot ERROR: Please specify pixelgram object")

     R <- x
     message("*** pixelgram.tree.plot ***\n")
     class(R) <- "pixelgram"

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

    if (R$show_tree & !is.null(leaf_colors)) { # pth 02102016 not sure these should be combined

      if (any(is.na(leaf_colors)))
        stop(paste("pixelgram.tree.plot() ERROR! Undefined leaf_colors:",
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
	    selected_rows[i]+1/4,
            R$x_lim[2],
	    selected_rows[i]-1/4,
	    border=selected_row_colors[i], col=NA, lwd=1/2)

    } else {
# the default is not to invert
	if (no_margin) {
#	R$y_lim <- c(1+length(R$tre$tip.label), 0)

#	    stop("Sorry, this is not working correctly")
	    x_factor <- (1 + R$y_lim[1] - R$y_lim[2]) / length(R$tre$tip.label)

	  for (i in 1:length(selected_rows))
	    rect(R$x_lim[1],
		length(R$tre$tip.label)-selected_rows[i]/x_factor + 1/4,
		R$x_lim[2],
		length(R$tre$tip.label)-selected_rows[i]/x_factor - 1/4,
		border=selected_row_colors[i], col=NA, lwd=1/2)
	} else {
#       R$y_lim <- c(length(R$tre$tip.label), 1)
	  for (i in 1:length(selected_rows))
	    rect(R$x_lim[1],
		length(R$tre$tip.label)-selected_rows[i]+3/4,
		R$x_lim[2],
		length(R$tre$tip.label)-selected_rows[i]+5/4,
		border=selected_row_colors[i], col=NA, lwd=1/4)
	}
    }
  }

#    show_tip_label = T #ifelse(is.null(tip_labels), T, F)
    safe.tiplabels = R$tre$tip.label
    if (is.character(tip_labels$pch))
        R$tre$tip.label = sapply(1:length(R$tre$tip.label), function(i)
            paste0(" ", tip_labels$pch[i]))
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

    if (!R$show_tree & show_tip_label)
      for (i in 1:length(R$tre$tip.label))
        text(0, i, #-1.01 * R$raster_width * R$x_lim[2],
             adj=0, R$tre$tip.label[i], cex=R$my_cex, col=tip_labels$col[i])

    # design issue here is how to map leaf names to { pch, bg, col }
    if (!is.null(tip_labels) & !is.null(tip_labels$pch)) {

        if (!is.null(tip_labels) & length(which(is.na(tip_labels$pch))) > 1)
	    stop(paste("pixelgram ERROR! Undefined tip_labels$pch:",
		paste(R$tre$tip.label[which(is.na(tip_labels$pch))],
		    collapse=" ")))

	if (!is.null(tip_labels) & length(which(is.na(tip_labels$col))) > 1)
	    stop(paste("pixelgram ERROR! Undefined tip_labels$col:",
		paste(R$tre$tip.label[which(is.na(tip_labels$col))],
		    collapse=" ")))

	if (R$show_tree & is.numeric(tip_labels$pch))
            ape::tiplabels(pch=tip_labels$pch,
                  col=tip_labels$col,
                  bg=tip_labels$bgs, #adj=0,
                  cex=R$my_cex, lwd=1/2)


# to do: position properly with raster_margin considered:
	if (plot_margin_points) {

	    if ((!is.null(leaf_colors) & is.null(tip_labels$col)) |
		(!is.null(tip_labels$col) & any(tip_labels$col != leaf_colors)))
		    for (i in 1:length(R$tre$tip.label))
		        points(-1.01 * R$raster_width * R$x_lim[2], i, col=leaf_colors[i], cex=R$my_cex, lwd=1/2, pch=15)

##		if (!is.null(tip_labels$col)) {

##	        if (is.null(tip_labels$pch) | any(is.na(tip_labels$pch))) {

#		      for (i in 1:length(R$tre$tip.label))
#		        points(-0.02 * R$raster_width * R$x_lim[2], i, col=tip_labels$col[i], cex=R$my_cex, lwd=1/2, pch=15)

# unsure why a second column of black tiplabels appears immediately next to pheno_matrix
#		} else {
#		    if (length(tip_labels$bgs)==1)
#			tip_labels$bgs = rep(tip_labels$bgs, length(tip_labels$pch))
#		    for (i in 1:length(R$tre$tip.label))
#		        points(-0.01 * R$raster_width * R$x_lim[2], i,
#			      col=tip_labels$col[i], pch=tip_labels$pch[i], bg=tip_labels$bgs[i], cex=R$my_cex, lwd=1/2)
			    # bgs could be a scalar or vector
##		}
##	    }
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

		myordering <- R$tre$edge[which(R$tre$edge[, 2] <= length(R$tre$tip.label)), 2]

		rect(rep(L_pos, length(R$tre$tip.label)),
			c(length(R$tre$tip.label):1)-1/2,
			rep(R_pos, length(R$tre$tip.label)),
			c(length(R$tre$tip.label):1)+1/2,
			col=pheno_matrix[myordering[c(length(R$tre$tip.label):1)], col_num],
			border=NA)

		if (!is.null(pheno_colnames)) {

		    myxpos = mean(c(L_pos, R_pos))

		    usr=par('usr')
		    ypos.1 <- usr[3] + (usr[4] - usr[3])*0.03
		    ypos.2 <- usr[3] + (usr[4] - usr[3])*0.97

		    text(myxpos, ypos.1, pheno_colnames[col_num], cex=1/2, adj=c(1, 1/2), srt=90)
		    text(myxpos, ypos.2, pheno_colnames[col_num], cex=1/2, adj=c(0, 1/2), srt=90)
		}

		if (!is.null(pheno_letters)) {
		    myxpos = mean(c(L_pos, R_pos))

		    usr=par('usr')
		    ypos.1 <- usr[3] + (usr[4] - usr[3])*0.03
		    ypos.2 <- usr[3] + (usr[4] - usr[3])*0.97

		    text(myxpos, ypos.1, pheno_letters[col_num],
		    		 cex=2/3, adj=c(1/2, 2/2))
		    text(myxpos, ypos.2, pheno_letters[col_num],
		    		 cex=2/3, adj=c(1/2, 0/2))
		}
	    }
	}

    if (!is.null(tip_labels) & !is.null(tip_labels$pch)) {

	if (plot_margin_points) {

	    if ((!is.null(leaf_colors) & is.null(tip_labels$col)) |
		(!is.null(tip_labels$col) & any(tip_labels$col != leaf_colors)))
	        for (i in 1:length(R$tre$tip.label))
                    points(-0.01 * R$raster_width * R$x_lim[2], i, col=leaf_colors[i], cex=R$my_cex, lwd=1/2, pch=15)

	    if (!is.null(tip_labels$col)) {

		if (is.null(tip_labels$pch) | any(is.na(tip_labels$pch))) {

#		    for (i in 1:length(R$tre$tip.label))
#		        points(-0.02 * R$raster_width * R$x_lim[2], i, col=tip_labels$col[i], cex=R$my_cex, lwd=1/2, pch=15)

# unsure why a second column of black tiplabels appears immediately next to pheno_matrix
#		} else {
#		    if (length(tip_labels$bgs)==1)
#			tip_labels$bgs = rep(tip_labels$bgs, length(tip_labels$pch))
#		    for (i in 1:length(R$tre$tip.label))
#		        points(-0.01 * R$raster_width * R$x_lim[2], i,
#			    col=tip_labels$col[i], pch=tip_labels$pch[i], bg=tip_labels$bgs[i],
#			    cex=R$my_cex, lwd=1/2)
			    # bgs could be a scalar or vector
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

    if (R$show_tree & !is.null(scale_bar)) {

	sb.xpos = ifelse(is.null(scale_bar$pos[1]), NULL, scale_bar$pos[1])
	sb.ypos = ifelse(is.null(scale_bar$pos[2]), NULL, scale_bar$pos[2])


	  ape::add.scale.bar(sb.xpos, sb.ypos, length=scale_bar$len,
	                     lwd=scale_bar$lwd, cex=1/1000)

	if (R$invert_y == T) {
	    text(sb.xpos+scale_bar$len/2, sb.ypos-1, scale_bar$text,
		cex=scale_bar$cex, adj=c(1/2, 1))
	} else {
	    text(sb.xpos+scale_bar$len/2, sb.ypos+1, scale_bar$text,
		cex=scale_bar$cex, adj=c(1/2, 1))
	}
    }

    R$x_lim <- c(tree_out$x.lim[1], tree_out$x.lim[2])

    return ( R )
}

#' @keywords internal
pixelgram.raster.nt <- function(P, vbars, show_top_axis) {

    if (class(P) != "pixelgram")
	stop("pixelgram.raster.nt ERROR: Please specify pixelgram object")

    R <- P
    class(R) <- "pixelgram"

    if (is.null(R$color_lut_type))
	R$color_lut_type = "aa"

    if (!is.null(R$nt_rast)) {

        message("*** Making nt raster map ***")

	if (0) {

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
	} else {

	  col_alphabet <- unique(c(R$nt_rast))

	  for (j in 1:ncol(R$nt_rast))
	      for (i in 1:length(col_alphabet))
	          if (any(R$nt_rast[, j] == col_alphabet[i]))
	              R$nt_rast[, j] = replace(R$nt_rast[, j],
	                  which(R$nt_rast[, j] == col_alphabet[i]), i)

	      raster.colors <- pixmap_colors[col_alphabet, R$color_lut_type]

	      my.xlim <- sort(c(-R$raster_width * R$x_lim[2], -R$raster_width * R$x_lim[2] * R$raster_margin))
	      my.ylim <- c(1, nrow(R$nts))

	      nt_num <- matrix(as.numeric(R$nt_rast), ncol=ncol(R$nt_rast), nrow=nrow(R$nt_rast))

              # adjust grid spacing to enable drawing a rect() surrounds *outer* part of image grid
	      x.grid.spacing <- 0.5 * (my.xlim[2] - my.xlim[1]) / (ncol(nt_num))
              x.grid <- seq(my.xlim[1]+x.grid.spacing, my.xlim[2]-x.grid.spacing, length.out=ncol(nt_num))
              y.grid <- seq(my.ylim[1], my.ylim[2], length.out=nrow(nt_num))

	      if (R$invert_y)
	          nt_num <- nt_num[nrow(nt_num):1, ]

	      image(sort(x.grid), sort(y.grid),
	            t(nt_num),  #[ncol(nt_num):1, ],
	            col=raster.colors,
	            useRaster=F,
		    breaks=c(0:length(col_alphabet)) + 1/2,
	            add=T,
		    xlab='', ylab='',
	            ylim=my.ylim,
	            xlim=my.xlim)
        }
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

    return ( R )
}


#' @keywords internal
annotate.region <- function(R, notes=NULL, y_lim=NULL) {

    message("*** annotate.region ***\n")

    clr <- "#88888844"

    if (is.null(notes)) {
	notes <- list()
	notes$Lhs <- c(132, 185, 298, 396, 459, 511, 30)
	notes$Rhs <- c(152, 190, 331, 410, 466, 512, 31)
	notes$txt <- c('V1', 'V2', 'V3', 'V4', 'V5', '', '')
#	notes$Lhs <- c(132, 185, 276, 396, 459, 511, 30)
#	notes$Rhs <- c(152, 190, 283, 410, 466, 512, 31)
#	notes$txt <- c('V1', 'V2', 'Loop D', 'V4', 'V5', '', '')
    } else {
	if (is.null(notes$Lhs) | is.null(notes$Rhs))
	    return (R)
	if (!is.numeric(notes$Lhs) | !is.numeric(notes$Rhs))
	    return (R)

	if (!is.null(notes$txt))
	    if (!is.character(notes$col))
		return (R)

	if (!is.null(notes$col))
	    clr <- notes$col
    }

#	} else if (!is.null(region) & grepl("^V", region)) {
#	    Lnotes <- c(26, 50, 98)
#	    Rnotes <- c(35, 65, 110)
#	     notes <- c('H1', 'H2', 'H3')
#	} else {
#	    Lnotes <- NULL
#	    Rnotes <- NULL
#	     notes <- NULL
#	}
#    }

    if (is.null(y_lim)) y_lim = (usr <- par('usr'))[c(3,4)]

### need an efficient way to transform alignment columns to x-axis coordinates

    axis_xlim = c(-R$raster_width*R$x_lim[2], -R$raster_width*R$x_lim[2]*R$raster_margin)

    my_slope = ifelse(!is.null(R$aas), ncol(R$aas)-1, ncol(R$nts)-1)/
            (axis_xlim[2]-axis_xlim[1])

    refseq_lut = R$refseq_lut
    refseq_lut$aln = (refseq_lut$aln - 1)/my_slope - R$raster_width * R$x_lim[2]

    if (!is.null(refseq_lut)) {

 	for (x in 1:length(notes$txt)) {

 	    x_1 = NULL
# this causes an error message when using gp120s:
	    x_1 <- try(refseq_lut$aln[min(which(refseq_lut$l==notes$Lhs[x]))])

 	    x_2 = NULL
# this causes an error message when using gp120s:
	    x_2 <- try(refseq_lut$aln[max(which(refseq_lut$r==notes$Rhs[x]))])

 	    if (!is.null(x_1) & !is.null(x_2) & !is.na(x_1) & !is.na(x_2))
 		rect(x_1, y_lim[1], x_2, y_lim[2], border=NA, col=clr)
 	}

        ### annotate text at top of plot here so we don't draw boxes over text
 	if (!is.null(notes$txt)) {

 	    for (i in 1:length(notes$txt)) {

 		l.xpos <- NULL
		l.xpos = try(refseq_lut$aln[min(which(refseq_lut$l==notes$Lhs[i]))])

 		r.xpos <- NULL
		r.xpos = try(refseq_lut$aln[max(which(refseq_lut$r==notes$Rhs[i]))])

		usr <- par('usr')

                #  consider inverted y-axis values
		y.pos <- ifelse(usr[3] < usr[4],
		    usr[4],# - (usr[4] - usr[3])*0.995,
		    usr[3])# - (usr[4] - usr[3])*0.995)

		n_s <- ifelse(usr[3] < usr[4], 3, 1) # north or south?
		my.adj=ifelse(usr[3] < usr[4], 1, 0)

 		if (!is.null(l.xpos) & !is.null(r.xpos) &
		      !is.na(l.xpos) & !is.na(r.xpos)) {

		    if (R$xform_master) # use blank space at top of pixel plot
 			text(x=mean(c(l.xpos, r.xpos)),
			    y=1, # could figure out which row is master sequence row but whatever
			    notes$txt[i], #pos=n_s,
			    cex=2/3)

		    if (!R$xform_master) # use margin
 			mtext(notes$txt[i], 3, #pos=n_s,
			    at=mean(c(l.xpos, r.xpos)), line=-1/2, #y.pos,
			    cex=2/3)
		}
 	    }
        }
    }

    R
}

#' @keywords internal
pixelgram.raster.aa <- function(P, vbars, show_top_axis, x_labs=NULL) {

    if (class(P) != "pixelgram")
	stop("pixelgram.raster.aa ERROR: Please specify pixelgram object")

    R <- P
    class(R) <- "pixelgram"

    if (!is.null(R$refseq_lut)) {

        axis_xlim = c(-R$raster_width * R$x_lim[2], -R$raster_width * R$x_lim[2] * R$raster_margin)

        my_slope = (ncol(R$aas)-1)/(axis_xlim[2]-axis_xlim[1])

 	my_refseq_lut = R$refseq_lut
 	my_refseq_lut$aln = (my_refseq_lut$aln - 1)/my_slope - R$raster_width * R$x_lim[2]
    }

    if (!is.null(R$aa_rast)) {

        message("*** Making aa raster map ***")

	if (is.null(R$color_lut_type))
	    R$color_lut_type = "aa"

	if (R$xform_type==3)
	    R$color_lut_type = "charge"

	if (0) { # 10-02-2015 replace rasterImage() with image(..., useRaster=F) to fix fugly antialiasing

            # is it really necessary to use nested 'for' loops here?
	    for (j in 1:ncol(R$aa_rast)) {

		col_alphabet <- unique(c(R$aa_rast[, j]))

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

	} else {

	  col_alphabet <- unique(c(R$aa_rast))

	  for (j in 1:ncol(R$aa_rast))
	    for (i in 1:length(col_alphabet))
	      if (any(R$aa_rast[, j] == col_alphabet[i]))
	        R$aa_rast[, j] = replace(R$aa_rast[, j],
	                                 which(R$aa_rast[, j] == col_alphabet[i]), i)

	      raster.colors <- pixmap_colors[col_alphabet, R$color_lut_type]

	      my.xlim <- sort(c(-R$raster_width * R$x_lim[2], -R$raster_width * R$x_lim[2] * R$raster_margin))
	      my.ylim <- c(1, nrow(R$aas))

	      aa_num <- matrix(as.numeric(R$aa_rast), ncol=ncol(R$aa_rast), nrow=nrow(R$aa_rast))

              # adjust grid spacing to enable drawing a rect() surrounds *outer* part of image grid
	      x.grid.spacing <- 0.5 * (my.xlim[2] - my.xlim[1]) / (ncol(aa_num))
              x.grid <- seq(my.xlim[1]+x.grid.spacing, my.xlim[2]-x.grid.spacing, length.out=ncol(aa_num))
              y.grid <- seq(my.ylim[1], my.ylim[2], length.out=nrow(aa_num))

	      if (R$invert_y)
	          aa_num <- aa_num[nrow(aa_num):1, ]

	      image(sort(x.grid), sort(y.grid),
	            t(aa_num),#[ncol(aa_num):1, ],
	            col=raster.colors,
	            useRaster=F,
		    breaks=c(0:length(raster.colors)) + 1/2,
	            add=T,
		    xlab='', ylab='',
	            ylim=my.ylim,
	            xlim=my.xlim)
	}
    }

    rect(-R$raster_width * R$x_lim[2], 1/2 + nrow(R$aas),
         -R$raster_width * R$x_lim[2] * R$raster_margin, # 0.025,
	 1/2,
         border='black', col=NA, lwd=1/2)

  # new code: FOR NOW ONLY WORKS WITH AA ALIGNMENTS
#-R$raster_width * R$x_lim[2] * raster_margin
# ncol(R$aas) is number of aligned sites -> min(R$refseq_lut$aln),

    if (is.null(x_labs)) {
	 tick_interval = ifelse(ncol(R$aas) >= 3000, 1000,
 	     ifelse(ncol(R$aas) >= 300, 100,
 		 ifelse(ncol(R$aas) >= 200, 50,
 		     ifelse(ncol(R$aas) >= 150, 30,
 			 ifelse(ncol(R$aas) >= 80, 20,
 			     ifelse(ncol(R$aas) >= 40, 10,
 				 ifelse(ncol(R$aas) >= 20, 5, 1)))))))

         axis_xlim = c(-R$raster_width * R$x_lim[2], -R$raster_width * R$x_lim[2] * R$raster_margin)#0.025)

         my_slope = (ncol(R$aas)-1)/(axis_xlim[2]-axis_xlim[1])

        if (!is.null(R$refseq_lut))
 	    x_locs <- R$refseq_lut$aln[which(R$refseq_lut$aln < ncol(R$aas) &
 	        R$refseq_lut$l == R$refseq_lut$r & R$refseq_lut$l %% tick_interval == 0)]
        else
 	    x_locs <- 1:ncol(R$aas)

 	x_locs = (x_locs-1)/my_slope - R$raster_width * R$x_lim[2]
# #	axis_xlim = c(1, ncol(R$aas)) / my_slope# - R$raster_width * R$x_lim[2]

         x_labs <- R$refseq_lut$l[which(R$refseq_lut$l == R$refseq_lut$r &
 		R$refseq_lut$l %% tick_interval == 0)]
### to do: think about what to do when my_lim %in% x_locs
 	my_cexl=4/5
 	my_cexa=2/3
    } else {

         axis_xlim = c(-R$raster_width * R$x_lim[2], -R$raster_width * R$x_lim[2] * R$raster_margin)#0.025)

         my_slope = (ncol(R$aas))/(axis_xlim[2]-axis_xlim[1])
#         my_slope = (ncol(R$aas)-1)/(axis_xlim[2]-axis_xlim[1])

	x_locs = 1:length(x_labs)
 	x_locs = (x_locs - 1/2)/my_slope - R$raster_width * R$x_lim[2]
### to do: think about what to do when my_lim %in% x_locs
 	my_cexl=1/5
 	my_cexa=1/2
    }

    if (length(x_labs) != length(x_locs))
	message("WARNING: Number of raster x-axis labels differs from the number of locations.")


#If either line or pos is set, they (rather than par("mgp")[3])
#determine the position of the axis line and tick marks, and the tick
#labels are placed par("mgp")[2] further lines into (or towards for
#pos) the margin.

 	my_mgp=c(2/3, 0/5, 0)

	usr=par('usr')
#	my.tcl = ifelse(usr[3]<usr[4], (usr[4]-usr[3])/500, (usr[3]-usr[4])/500)

 	if (length(x_locs) <= length(x_labs)) {

            axis(1, at=x_locs[1:length(x_labs)], labels=x_labs,
 		xaxt='s', cex.lab=my_cexl, cex.axis=my_cexa, mgp=my_mgp,
#		pos=1/2,
 		pos=ifelse(R$invert_y==T, R$y_lim[1]-1/2, R$y_lim[1]+1/2),
 		padj=0/2,
 #		hadj=0/2,
 		lwd=1/2)#,
# 		tcl=-my.tcl)#k=-0.5/nrow(R$aas))

 	my_mgp=c(2/3, 0/5, 0)
#        axis(3, at=c(axis_xlim[1], x_locs, axis_xlim[2]), labels=c("", x_labs, ""),
 	if (show_top_axis) {

            axis(3, at=x_locs[1:length(x_labs)],
		labels=x_labs,
 		xaxt='s',
		cex.lab=my_cexl,
		cex.axis=my_cexa,
		mgp=my_mgp,
 		pos=ifelse(R$invert_y, R$y_lim[2]+1/2, R$y_lim[2]-1/2),
 		padj=ifelse(R$invert_y, 2/2, 0/2),
		lwd=1/2)#,
#		tcl=-my.tcl)#k=-0.5/nrow(R$aas))

# 	} else {
# 	    message(paste("* NO xlabs=", paste(x_labs, collapse=" ")))
# 	    message(paste("* NO xlocs=", paste(x_locs, collapse=" ")))
# 	}
      }
    }

    if (!is.null(R$main))
      mtext(paste0(R$main, " "), 3, at=0, line=0, cex=1, outer=F, adj=1)

    if (!is.null(R$sub))
    	mtext(paste(R$sub, "site"), 1, at=mean(axis_xlim), line=1/2, cex=my_cexl)
	    # mtext(R$sub, 3, at=0, line=0, cex=1, outer=F, adj=1/2, font=3)

    R
}
