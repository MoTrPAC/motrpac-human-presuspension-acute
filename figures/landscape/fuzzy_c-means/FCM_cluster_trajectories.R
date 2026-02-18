## Author: Tyler Sagendorf
## Date: 2025-05-07
##
## Purpose: Create line plots to display the trajectories of individual features
## in each cluster from the fuzzy c-means clustering results.

library(MotrpacHumanPreSuspension) # also loads data package

## Folder to save plots, relative to precovid-analyses/
folder <- file.path("figures", "landscape", "fuzzy_c-means",
                    "plots_cluster_trajectories")

if (!dir.exists(folder))
  dir.create(folder)

## All clusters ----------------------------------------------------------------
plot_cmeans(FCM = FCM_CLUSTERS$adipose,
            ncol = 3L, # number of clusters per column
            # Features with cluster membership probabilities below 0.3 will not
            # be plotted.
            min_membership = 0.3,
            filename = file.path(folder,
                                 "FCM_cluster_trajectories_adipose.pdf"))

plot_cmeans(FCM = FCM_CLUSTERS$blood,
            ncol = 3L,
            min_membership = 0.3,
            filename = file.path(folder,
                                 "FCM_cluster_trajectories_blood.pdf"))

plot_cmeans(FCM = FCM_CLUSTERS$muscle,
            ncol = 3L,
            min_membership = 0.3,
            filename = file.path(folder,
                                 "FCM_cluster_trajectories_muscle.pdf"))

## All clusters, single wide strip ---------------------------------------------

plot_cmeans(FCM = FCM_CLUSTERS$adipose,
            ncol = 20L,
            min_membership = 0.3,
            filename = file.path(folder,
                                 "FCM_cluster_trajectories_wide_adipose.pdf"))

plot_cmeans(FCM = FCM_CLUSTERS$blood,
            ncol = 20L,
            min_membership = 0.3,
            filename = file.path(folder,
                                 "FCM_cluster_trajectories_wide_blood.pdf"))

plot_cmeans(FCM = FCM_CLUSTERS$muscle,
            ncol = 20L,
            min_membership = 0.3,
            filename = file.path(folder,
                                 "FCM_cluster_trajectories_wide_muscle.pdf"))
