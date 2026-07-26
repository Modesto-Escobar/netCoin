# Examples: Using external clustering (LCA, LPA) with surScat
# These show how to compute clusters externally and integrate with surScat

# ============================================================================
# EXAMPLE 1: Latent Class Analysis (LCA) with poLCA
# ============================================================================
# LCA is ideal for categorical data
# Install: install.packages("poLCA")

if (require("poLCA", quietly = TRUE)) {
  # Create example data: categorical responses
  set.seed(123)
  n <- 200
  df_lca <- data.frame(
    item1 = factor(sample(1:3, n, replace=TRUE, prob=c(0.5, 0.3, 0.2))),
    item2 = factor(sample(1:3, n, replace=TRUE, prob=c(0.4, 0.4, 0.2))),
    item3 = factor(sample(1:2, n, replace=TRUE, prob=c(0.6, 0.4))),
    age = sample(c("Young", "Old"), n, replace=TRUE),
    sex = sample(c("M", "F"), n, replace=TRUE)
  )

  # Step 1: Create surScat object
  scatter <- surScat(df_lca,
                     variables = c("item1", "item2", "item3", "age", "sex"),
                     active = c("item1", "item2", "item3"),
                     type = "mca",
                     nclusters = 2,
                     vPatterns = c("age", "sex"))

  # Step 2: Prepare data for poLCA (must be numeric 1,2,3...)
  # poLCA requires responses as integers starting at 1
  df_for_lca <- df_lca[, c("item1", "item2", "item3")]
  for (col in names(df_for_lca)) {
    df_for_lca[[col]] <- as.integer(df_for_lca[[col]])
  }

  # Step 3: Fit LCA model using poLCA
  # f = formula of items to use
  lca_formula <- cbind(item1, item2, item3) ~ 1
  lca_fit <- poLCA::poLCA(lca_formula,
                          data = df_for_lca,
                          nclass = 3,      # number of classes
                          verbose = FALSE)

  # Step 4: Extract cluster assignments
  lca_clusters <- lca_fit$predclass  # posterior class assignments

  # Step 5: Add to surScat object
  # NOTE: If surScat was called with vPatterns, case-level clusters are automatically
  # collapsed to pattern level using the mode
  scatter_lca <- addClusters(scatter,
                             clusters = lca_clusters,
                             name = "LCA(3)")

  # Now visualize: scatter_lca will show both k-means and LCA groupings
  # Use color = "LCA(3)" to color by LCA classes


  # ========== Alternative: View LCA model diagnostics ==========
  # print(lca_fit)                    # Model summary
  # lca_fit$probs                     # Item probabilities by class
  # lca_fit$posterior                 # Posterior class probabilities (softer assignments)
}


# ============================================================================
# EXAMPLE 2: Latent Profile Analysis (LPA) with tidyLPA
# ============================================================================
# LPA is ideal for continuous data (numeric variables)
# Install: install.packages("tidyLPA")

if (require("tidyLPA", quietly = TRUE)) {
  # Create example data: continuous measurements
  set.seed(456)
  n <- 200
  df_lpa <- data.frame(
    score1 = rnorm(n, mean=50, sd=10),
    score2 = rnorm(n, mean=60, sd=12),
    score3 = rnorm(n, mean=55, sd=8),
    group = sample(c("Control", "Treatment"), n, replace=TRUE),
    time = sample(c("T1", "T2"), n, replace=TRUE)
  )

  # Step 1: Create surScat object
  scatter <- surScat(df_lpa,
                     variables = c("score1", "score2", "score3", "group", "time"),
                     active = c("score1", "score2", "score3"),
                     type = "pca",
                     nclusters = 2,
                     vPatterns = c("group", "time"))

  # Step 2: Fit LPA models and select best
  # tidyLPA estimates multiple models and returns model with best fit
  lpa_results <- tidyLPA::estimate_profiles(
    df_lpa[, c("score1", "score2", "score3")],
    n_profiles = 3,          # try 3 profiles
    models = 1,              # model specification (1 = Mclust VEE)
    verbose = FALSE)

  # Step 3: Extract cluster assignments
  lpa_clusters <- tidyLPA::get_data(lpa_results)$Class  # posterior assignments

  # Step 4: Add to surScat object
  scatter_lpa <- addClusters(scatter,
                             clusters = lpa_clusters,
                             name = "LPA(3)")

  # Now visualize both k-means and LPA groupings


  # ========== Alternative: View LPA diagnostics ==========
  # tidyLPA::get_summary(lpa_results)  # Model fit indices
  # tidyLPA::plot_profiles(lpa_results) # Profile plot
}


# ============================================================================
# EXAMPLE 3: Hierarchical Clustering (base R)
# ============================================================================
# No dependencies needed! Works with any distance metric

# Create example data
set.seed(789)
n <- 150
df_hc <- data.frame(
  x = rnorm(n, mean=0, sd=2),
  y = rnorm(n, mean=0, sd=2),
  z = rnorm(n, mean=0, sd=2),
  category = sample(c("A", "B", "C"), n, replace=TRUE)
)

# Step 1: Create surScat object
scatter <- surScat(df_hc,
                   variables = c("x", "y", "z", "category"),
                   active = c("x", "y", "z"),
                   type = "pca",
                   nclusters = 2,
                   vPatterns = "category")

# Step 2: Compute hierarchical clustering on the active variables
hc_data <- scale(df_hc[, c("x", "y", "z")])  # standardize
hc <- hclust(dist(hc_data), method = "ward.D2")
hc_clusters <- cutree(hc, k = 3)  # cut into 3 clusters

# Step 3: Add to surScat object
scatter_hc <- addClusters(scatter,
                          clusters = hc_clusters,
                          name = "Hierarchical(3)")


# ============================================================================
# EXAMPLE 4: Comparing multiple clustering methods
# ============================================================================

# Start with one surScat object
scatter_base <- surScat(df_lpa,
                        variables = c("score1", "score2", "score3", "group", "time"),
                        active = c("score1", "score2", "score3"),
                        type = "pca",
                        nclusters = 2)

# Add different clustering results sequentially
scatter_kmeans2 <- addClusters(scatter_base, kmeans(df_lpa[,1:3], centers=2)$cluster, name="kmeans(2)")
scatter_kmeans3 <- addClusters(scatter_base, kmeans(df_lpa[,1:3], centers=3)$cluster, name="kmeans(3)")

# If you have LPA installed:
if (require("tidyLPA", quietly = TRUE)) {
  lpa_2 <- tidyLPA::estimate_profiles(df_lpa[, c("score1", "score2", "score3")],
                                      n_profiles = 2, verbose = FALSE)
  lpa_clusters_2 <- tidyLPA::get_data(lpa_2)$Class
  scatter_lpa2 <- addClusters(scatter_base, lpa_clusters_2, name="LPA(2)")

  lpa_3 <- tidyLPA::estimate_profiles(df_lpa[, c("score1", "score2", "score3")],
                                      n_profiles = 3, verbose = FALSE)
  lpa_clusters_3 <- tidyLPA::get_data(lpa_3)$Class
  scatter_lpa3 <- addClusters(scatter_base, lpa_clusters_3, name="LPA(3)")

  # Now you can compare all methods visually by switching the color attribute
}


# ============================================================================
# Useful tips:
# ============================================================================
# 1. If surScat was called with vPatterns, case-level clusters are automatically
#    collapsed to pattern level. No manual collapsing needed!
#    scatter <- surScat(df, vPatterns=c("var1", "var2"))
#    scatter_lca <- addClusters(scatter, case_level_lca_clusters)  # auto-collapsed
#
# 2. If your clustering result is soft (probabilities), convert to hard assignments:
#    clusters <- apply(posterior_probs, 1, which.max)
#
# 3. Use replaceClusters() if you want to replace the k-means result entirely:
#    scatter <- replaceClusters(scatter, my_clusters, name="MyMethod")
#
# 4. Use weight parameter if your clustering accounts for case weights:
#    scatter <- addClusters(scatter, clusters, weight = my_weights)
#
# 5. Set sort=FALSE to keep the original cluster numbering from your external method:
#    scatter <- addClusters(scatter, clusters, sort=FALSE)
