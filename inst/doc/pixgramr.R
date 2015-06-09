## ---- fig.show='asis', tidy=F, warning=F, message=F, fig.width=11, fig.height=8.5----
par(cex=1/2)
library(pixgramr, warn.conflicts=F, quietly=T)

n <- pixgramr::pixgram(tre=hiv.ref$tre, nts=hiv.ref$nts,
                       raster_width=1/2,
                       master_name="__consensus__", excise_refseq=F)
plot(n, show_tip_label=T)
plot(n, xform_type=2, show_tip_label=T)

## ---- fig.show='asis', tidy=F, warning=F, message=F, fig.width=11, fig.height=8.5----

clade.colors <- c("#FFC20E", "#008FD4", "#ED1C24", "#00A99D", "#F7941D", 
                  "#235192", "#8DC73F", "#00A651", "#EC008C")
names(clade.colors) = sort(unique(hiv.ref$clade))

leaf.colors = sapply(1:length(hiv.ref$tre$tip.label), function(i)
  clade.colors[hiv.ref$clade[i]])

legend.attributes <- list(text=names(clade.colors), cols=clade.colors, 
                          bg="white", bty='o', lwds=3/2)

p <- pixgramr::pixgram(tre=hiv.ref$tre, nts=hiv.ref$nts, is_orf=T,
                       refseq_name="Ref-B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
                       master_name="__consensus__", excise_refseq=F)

plot(p, xform_type=2, xform_master=T, 
     leaf_colors=leaf.colors, legend_attributes=legend.attributes)

## ---- fig.show='asis', tidy=F, warning=F, message=F, fig.width=11, fig.height=8.5----
library(pixgramr, warn.conflicts=F, quietly=T)
p <- pixgramr::pixgram(tre=CH505$tre, aas=CH505$aas, excise_refseq=T)
plot(p, xform_type=1, xform_master=F, annotate_env=T)
plot(p, xform_type=2, xform_master=F)
plot(p, xform_type=0, xform_master=F, color_lut_type="taylor")

## ---- fig.show='asis', tidy=F, warning=F, message=F, fig.width=11, fig.height=8.5----
library(pixgramr)

wpi.colors=c( "#EC008C", "#ED1C24", "#F7941D", "#FFC20E", "#CADB2A",
              "#8DC73F", "#00A651", "#00A99D", "#31B6E9", "#008FD4",
              "#235192",
              "#662D91", "black", "#666666", "#888888")

names(wpi.colors) = sort(unique(CH505$wpi))

CH505.colors = sapply(1:length(CH505$tre$tip.label), function(i)
  wpi.colors[CH505$wpi[i]])

legend.attributes <- list(text=names(wpi.colors), cols=wpi.colors, 
                          bg="white", bty='o', lwds=3/2)

scale.bar <- list(pos=c(1/20,7), lwd=1, len=5/953, cex=4/5, text="5 aas")

plot(CH505.pixgram <- pixgramr::pixgram(aas=CH505$aas,
                                        tre=CH505$tre, xform_type=3, 
                                        raster_width=5/4,
                                        raster_margin=1/100),
     leaf_colors=CH505.colors, 
     legend_attributes=legend.attributes,
     scale_bar=scale.bar, no_margin=T, xform_master=T, annotate_env=T)

