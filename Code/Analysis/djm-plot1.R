sources <- c("Ground Truth", "Estimated"); types <- c("True (NHCS)", "Flatter")
dfa <- data.frame(t=t0+rep(1:n_ests, 4)-1,
                 HFRs=c(chging_hfrs[est_idx], flatter_hfrs[est_idx], est_chging, est_flat),
                 HFR=factor(rep(sources, each=n_ests*2), levels=sources),
                 Curve=factor(rep(rep(types, each=n_ests), 2), levels=types))
sources <- c("Ground Truth", "Estimated"); types <- c("12", "20", "28")

dfb <- data.frame(
  t = t0+rep(1:n_ests, 4)-1,
  HFRs = c(chging_hfrs[est_idx], est_short, est_med, est_long),
  HFR = factor(c(rep("Ground Truth", n_ests), rep("Estimated", n_ests*3)), levels = sources),
  Type = factor(c(rep("Ground Truth", n_ests), rep(types, each=n_ests)), levels = c("Ground Truth", types))
)
sources <- c("Ground Truth", "Estimated"); types <- c("True", "Rising", "Falling")

dfc <- data.frame(
  t = t0+rep(1:n_ests, 4)-1,
  HFRs = c(chging_hfrs[est_idx], hfrs_true_hosps, hfrs_alt_hosps1, hfrs_alt_hosps2),
  HFR = factor(c(rep("Ground Truth", n_ests), rep("Estimated", n_ests*3)), levels = sources),
  Type = factor(c(rep("Ground Truth", n_ests), rep(types, each=n_ests)), levels = c("Ground Truth", types))
)

cols <- c(
  `Ground Truth` = "black",
  `True (NHCS)` = "black",
  Flatter = '#DDAA33',
  `12` = '#EE7733',
  `20` = '#0077BB',
  `28` = '#EE3377',
  Rising = '#984ea3',
  Falling = '#CC3311',
  True = '#009988'
)

Plot <- ggplot(dfa, aes(x=t, y=HFRs, color=Curve, linetype=HFR)) +
  geom_line() +
  ggtitle("HFR estimates by true HFR") +
  xlab("Date") + theme_bw() +
  coord_cartesian(ylim = c(.1, .16)) +
  guides(linetype = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  scale_color_manual(values = cols)

delayPlot <- ggplot(dfb, aes(x = t, y = HFRs, linetype = HFR, color = Type)) +
  geom_line(aes(group = interaction(HFR, Type))) +  # Group by HFR and Type
  # Manually set color, exclude Ground Truth from legend
  coord_cartesian(ylim = c(.1, .16)) +
  scale_color_manual(
    breaks = c("12", "20", "28"),
    values = cols) +
  ggtitle("HFR estimates by true delay distribution") +
  xlab("Date") + ylab(NULL) +
  theme_bw() +
  guides(linetype = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme(legend.position = "bottom", axis.text.y = element_blank(),
        legend.box = "vertical", plot.margin = margin(6,0,6,0))

delay_inset <- ggplot(dfDistr, aes(x=t, y=Delay, color=Type)) + geom_line() + theme_bw() +
  scale_color_manual(breaks = c("12", "20", "28"), values = cols) +
  theme(legend.position="none", plot.background = element_blank(),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  xlab("Days after primary event") + ylab("") + ggtitle("Delay distribution")

delayPlotFinal <- delayPlot +
  annotation_custom(ggplotGrob(delay_inset), xmin=t0, xmax=t0+60, ymax=.14)

hospPlot <- ggplot(dfc, aes(x = t, y = HFRs, linetype = HFR, color = Type)) +
  coord_cartesian(ylim = c(.1, .16)) +
  geom_line(aes(group = interaction(HFR, Type))) +  # Group by HFR and Type
  scale_color_manual(
    breaks = c("True", "Rising", "Falling"), values = cols) +
  ggtitle("HFR estimates by primary incidence") +
  labs(color="Primary incidence trajectory") + xlab("Date") + ylab(NULL) +
  theme_bw() +
  guides(linetype = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme(legend.position = "bottom", plot.margin = margin(6,0,6,0),
        axis.text.y = element_blank(),
        legend.box = "vertical")

hosp_inset <- ggplot(dfHosp, aes(x=t, y=Hosp, color=Type, group=Type)) + geom_line() + theme_bw() +
  scale_color_manual(breaks = c("True", "Rising", "Falling"), values = cols) +
  theme(legend.position="none", plot.background = element_blank(),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  xlab("") + ylab("") + ggtitle("Hospitalizations")

hospPlotFinal <- hospPlot +
  annotation_custom(ggplotGrob(hosp_inset), xmin=t0, xmax=t0+60, ymax=.14)


fakelegend_df <- data.frame(
  cols = rep(cols[-(1:2)], 2),
  lty = rep(c(1, 2), each = length(cols) - 2),
  fig = rep(c(1, 2, 2, 2, 3, 3, 3), times = 2),
  x = 1:7,
  y = 1:7
)

ggplot(fakelegend_df, aes(x, y, color = cols, linetype = lty)) +
  geom_line() +
  scale_color_manual(values = cols[4:6], name = "Type") +
  scale_color_manual(values = cols[7:9], name = "Primary incidence trajectory") +
  scale_color_manual(values = cols[2:3], name = "Curve") +
  scale_linetype_manual(values = c(1L, 2L), breaks = c("Ground Truth", "Estimated"))
