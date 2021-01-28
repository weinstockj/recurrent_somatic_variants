manhattan = function(summary_stats, output_prefix) {

    flog.info("now plotting")

    chrom_lookup = tibble::tibble(
        CHROM = c(glue::glue("chr{1:22}"), "chrX"),
        CHROM_index = 1:23
    )

    scaling = 1e8

    bonferroni_pvalue = .05 / nrow(summary_stats)
    flog.info(glue::glue("bonferonni at {nrow(summary_stats)} tests is {bonferroni_pvalue}"))

    prepared = summary_stats %>% 
        dplyr::inner_join(chrom_lookup, by = "CHROM") %>%
        # Compute chromosome size
        dplyr::group_by(CHROM_index, CHROM) %>% 
        dplyr::summarise(chr_len = (max(POS) - min(POS)) / scaling) %>% 
        dplyr::ungroup(.) %>%
        # Calculate cumulative start position of each chromosome
        dplyr::mutate(tot = cumsum(chr_len) - chr_len) %>%
        # Add this info to the initial dataset
        dplyr::left_join(summary_stats, ., by = "CHROM") %>%
    #     # Add a cumulative position of each SNP
        dplyr::arrange(CHROM_index, POS) %>%
        dplyr::mutate(POScum = POS / scaling + tot) %>%
        # Filter SNP to make the plot lighter
        dplyr::filter(-log10(p.value) > 0.5) # pvalue < .316

    axis = prepared %>% 
           dplyr::group_by(CHROM_index, CHROM) %>%
           dplyr::summarize(center = (max(POScum) + min(POScum)) / 2) %>% # integer overflow concern
           dplyr::ungroup(.)

    p = ggplot2::ggplot(prepared, aes(x=POScum, y=-log10(p.value))) +
# Show all points
            # ggrastr::geom_point_rast(aes(color=as.factor(CHROM_index)), alpha=0.8, size=.3, raster.dpi = 300) +
            ggplot2::geom_point(aes(color=as.factor(CHROM_index)), alpha=0.8, size=.3) +
            ggplot2::scale_color_manual(values = c(rep(c("#7173C9", "#01035F"), 11), "#7173C9")) +
# custom X axis:
            ggplot2::scale_x_continuous(
                    label = stringr::str_replace(axis$CHROM, "chr", ""),
                    breaks= axis$center, expand = c(0, 0)
                    ) +
            ggplot2::scale_y_continuous(
                    breaks = seq(0, max(-log10(prepared$p.value)), by = 2),
                    expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 1))
                    ) +     # remove space between plot area and x axis
# Custom the theme:
            # ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey") +
            ggplot2::geom_hline(yintercept = -log10(bonferroni_pvalue), linetype = "dashed", color = "grey") +
            ggplot2::labs(y = expression(-log[10](pvalue))) + 
            ggplot2::theme_bw(base_size = 7, base_family = "Helvetica") +
            ggplot2::theme( 
                    legend.position="none",
                    axis.title.x = element_blank(),
                    axis.text.x = element_text(size = 3),
                    axis.text.y = element_text(size = 5),
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    panel.grid.major.y = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    panel.border = element_blank(),
                    axis.line.y = element_line(size = .8)
                 )


    fname_pdf = file.path(output_prefix, "manhattan.pdf")
    fname_tiff = file.path(output_prefix, "manhattan.tiff")
    ggplot2::ggsave(fname_pdf, p, units = "mm", width = 89, height = 70, dpi = 300)
    ggplot2::ggsave(fname_tiff, p, units = "mm", width = 89, height = 70, dpi = 300, type = "cairo")

    flog.info("done plotting.")
}

qqunif_plot = function(pvalues, 
                should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                xlab=expression(paste("Expected (",-log[10], " p-value)")),
                ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                par.settings=list(superpose.symbol=list(pch=pch)), ...) {


#error checking
            if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
                if(!(class(pvalues)=="numeric" || 
                            (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric")))) stop("pvalue vector is not numeric, can't draw plot")
            if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
                if (already.transformed==FALSE) {
                    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
                } else {
                    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
                }


            grp<-NULL
            n<-1
            exp.x<-c()
            if(is.list(pvalues)) {
                nn<-sapply(pvalues, length)
                rs<-cumsum(nn)
                re<-rs-nn+1
                n<-min(nn)
                if (!is.null(names(pvalues))) {
                    grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
                        names(pvalues)<-NULL
                } else {
                    grp=factor(rep(1:length(pvalues), nn))
                }
                pvo<-pvalues
                    pvalues<-numeric(sum(nn))
                    exp.x<-numeric(sum(nn))
                    for(i in 1:length(pvo)) {
                        if (!already.transformed) {
                            pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
                                exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
                        } else {
                            pvalues[rs[i]:re[i]] <- pvo[[i]]
                                exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
                        }
                    }
            } else {
                n <- length(pvalues)+1
                    if (!already.transformed) {
                        exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
                            pvalues <- -log10(pvalues)
                    } else {
                        exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
                    }
            }


#this is a helper function to draw the confidence interval
            panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {

                conf.points = min(conf.points, n-1);
                mpts<-matrix(nrow=conf.points*2, ncol=2)
                    for(i in seq(from=1, to=conf.points)) {
                        mpts[i,1]<- -log10((i-.5)/n)
                            mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
                            mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
                            mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
                    }
                grid::grid.polygon(x=mpts[,1],y=mpts[,2], gp=grid::gpar(fill=conf.col, lty=0), default.units="native")
            }

#reduce number of points to plot
            if (should.thin==T) {
                if (!is.null(grp)) {
                    thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
                        grp = thin$grp
                } else {
                    thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
                }
                pvalues <- thin$pvalues
                    exp.x <- thin$exp.x
            }
            gc()

            prepanel.qqunif= function(x,y,...) {
                A = list()
                    A$xlim = range(x, y)*1.02
                    A$xlim[1]=0
                    A$ylim = A$xlim
                    return(A)
            }

#draw the plot
        lattice::xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
                prepanel=prepanel, scales=list(axs="i"), pch=pch,
                panel = function(x, y, ...) {
                if (draw.conf) {
                panel.qqconf(n, conf.points=conf.points, 
                        conf.col=conf.col, conf.alpha=conf.alpha)
                };
                lattice::panel.xyplot(x,y, ...);
                lattice::panel.abline(0,1);
                }, par.settings=par.settings, ...
              )

}

qqunif_plot_save = function(summary_stats, output_prefix) {
    fname = file.path(output_prefix, "qqplot.pdf")
    pdf(file = fname, width = 5, height = 5)
    # png(file = fname, type = "cairo")
    # lattice::trellis.device(file = fname)
    print(qqunif_plot(summary_stats$p.value))
    dev.off()
}



