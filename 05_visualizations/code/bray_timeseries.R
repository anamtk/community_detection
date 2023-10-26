subsamples <- readRDS(here("05_visualizations",
                           "viz_data",
                           "ABUR1_bray_samples.RDS"))

summary <- readRDS(here("05_visualizations",
                        "viz_data",
                        "ABUR1_bray_summary.RDS"))

raw_bray <- readRDS(here("05_visualizations",
                       "viz_data",
                       "sbc_ABUR1_raw_bray.RDS"))

modeled_col <- "#E88C23"
observed_col <- "#438AA8"

ggplot() +
  geom_line(data = subsamples, aes(x = year, y = bray, group = iter), color = "#E88C23", alpha = 0.1) +
  geom_line(data = raw_bray, aes(x = year, y= raw_bray), color = "#438AA8") +
  geom_line(data = summary, aes(x = year, y = mean_bray), color = "#E88C23") +
  theme_bw()
