## ---- fig.show='asis', tidy=F, warning=F, message=F, fig.width=6, fig.height=6----
par(cex=1/2)
library(pixelgram, warn.conflicts=F, quietly=T)

n <- pixelgram(tre=hiv.ref$tre, nts=hiv.ref$nts,
                       master_name="__consensus__", excise_refseq=F,
                       main="Subtype reference", sub='gp41 ectodomain')

par(mar=c(0,0,1,0), oma=c(0,0,0,0))
plot(n, xform_type=2, show_tip_label=T, raster_width=1/4)
plot(n, xform_type=0, show_tree=F, raster_width=1000)

## ---- fig.show='asis', tidy=F, warning=F, message=F, fig.width=6, fig.height=6----

clade.colors <- c("#FFC20E", "#008FD4", "#ED1C24", "#00A99D", "#F7941D", 
                  "#235192", "#8DC73F", "#00A651", "#EC008C")
names(clade.colors) = sort(unique(hiv.ref$clade))

leaf.colors = sapply(seq_along(hiv.ref$tre$tip.label), function(i)
  clade.colors[hiv.ref$clade[i]])

legend.attributes <- list(title='Subtype', text=names(clade.colors),
                          cols=clade.colors, 
                          bg="white", bty='o', lwds=4)

p <- pixelgram(tre=hiv.ref$tre, nts=hiv.ref$nts, is_orf=T,
                       refseq_name="Ref-B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
                       master_name="__consensus__", excise_refseq=F)

hiv.ref.tip.labels <- list()
hiv.ref.tip.labels$pch = hiv.ref$clade
hiv.ref.tip.labels$bgs = NA
hiv.ref.tip.labels$col = "black"

par(mar=c(1,0,1,0), oma=c(0,0,0,0))
plot(p, xform_type=2, xform_master=T, edge_widths=4,
     leaf_colors=leaf.colors, legend_attributes=legend.attributes,
     tip_labels=hiv.ref.tip.labels, show_tip_label = T)

## ---- fig.show='asis', tidy=F, warning=F, message=F, fig.width=6, fig.height=8----
library(pixelgram, warn.conflicts=F, quietly=T)
p <- pixelgram(tre=CH505$tre, aas=CH505$aas, excise_refseq=T)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(p, xform_type=1, xform_master=F, annotate_env=T)
plot(p, xform_type=0, xform_master=F)
plot(p, xform_type=2, xform_master=F, color_lut_type="taylor")

## ---- fig.show='asis', tidy=F, warning=F, message=F, fig.width=8, fig.height=10----
library(pixelgram)

wpi.colors=c( "#EC008C", "#ED1C24", "#F7941D", "#FFC20E", "#CADB2A",
              "#8DC73F", "#00A651", "#00A99D", "#31B6E9", "#008FD4",
              "#235192", "#662D91",  "black",  "#666666", "#888888")

names(wpi.colors) = sort(unique(CH505$wpi))

CH505.colors = sapply(seq_along(CH505$tre$tip.label), function(i)
  wpi.colors[CH505$wpi[i]])

legend.attributes <- list(text=names(wpi.colors), cols=wpi.colors, 
                          bg="white", bty='o', lwds=3/2, title='WPI')

scale.bar <- list(pos=c(1/20,7), lwd=1, len=5/953, cex=4/5, text="5 aas")
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(CH505.pixelgram <- pixelgram(aas=CH505$aas,
                                        tre=CH505$tre, xform_type=3, 
                                        raster_width=5/4,
                                        raster_margin=1/100, main="CH505", 
                                        sub="Env gp160"),
     leaf_colors=CH505.colors, 
     legend_attributes=legend.attributes,
     scale_bar=scale.bar, no_margin=T, xform_master=T, annotate_env=T)

