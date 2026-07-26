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

  # Step 4: Add to surScat object
  # NOTE: If surScat was called with vPatterns, case-level clusters are automatically
  # collapsed to pattern level using the mode (see collapse report for details)
  # Pass the poLCA object directly, addClusters extracts clusters automatically
  scatter_lca <- addClusters(scatter,
                             clusters = lca_fit,
                             name = "LCA(3)",
                             sourceData = df_lca)  # optional: validates cluster order

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

  # Step 3: Add to surScat object
  # Pass the tidyLPA result directly, addClusters extracts clusters automatically
  scatter_lpa <- addClusters(scatter,
                             clusters = lpa_results,
                             name = "LPA(3)",
                             sourceData = df_lpa)  # optional: validates cluster order

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
# Pass the cluster vector directly (hierarchical clustering needs manual cutree)
scatter_hc <- addClusters(scatter,
                          clusters = hc_clusters,
                          name = "Hierarchical(3)",
                          sourceData = df_hc)  # optional: validates cluster order


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
# Pass kmeans objects directly; addClusters extracts clusters automatically
kmeans_2 <- kmeans(df_lpa[,1:3], centers=2)
scatter_kmeans2 <- addClusters(scatter_base, kmeans_2, name="kmeans(2)", sourceData=df_lpa)

kmeans_3 <- kmeans(df_lpa[,1:3], centers=3)
scatter_kmeans3 <- addClusters(scatter_base, kmeans_3, name="kmeans(3)", sourceData=df_lpa)

# If you have LPA installed:
if (require("tidyLPA", quietly = TRUE)) {
  lpa_2 <- tidyLPA::estimate_profiles(df_lpa[, c("score1", "score2", "score3")],
                                      n_profiles = 2, verbose = FALSE)
  scatter_lpa2 <- addClusters(scatter_base, lpa_2, name="LPA(2)", sourceData=df_lpa)

  lpa_3 <- tidyLPA::estimate_profiles(df_lpa[, c("score1", "score2", "score3")],
                                      n_profiles = 3, verbose = FALSE)
  scatter_lpa3 <- addClusters(scatter_base, lpa_3, name="LPA(3)", sourceData=df_lpa)

  # Now you can compare all methods visually by switching the color attribute
}


# ============================================================================
# Useful tips:
# ============================================================================
# 1. Pass clustering objects directly to addClusters()
#    addClusters() automatically extracts clusters from:
#    - kmeans objects (extracts $cluster)
#    - poLCA objects (extracts $predclass)
#    - tidyLPA results (extracts Class from get_data())
#    - vectors/factors (passed as-is)
#
#    scatter_lca <- addClusters(scatter, lca_fit, name="LCA(3)")  # poLCA object
#    scatter_km <- addClusters(scatter, kmeans_result, name="kmeans(3)")  # kmeans object
#
# 2. Validate cluster order with sourceData parameter
#    Pass the original data frame to validate that clusters are in the correct order:
#    scatter <- addClusters(scatter, clusters, sourceData=df)
#    This checks that the number of rows and order match, preventing silent errors.
#
# 3. If surScat was called with vPatterns, case-level clusters are automatically
#    collapsed to pattern level. Check the collapse report for details:
#    scatter <- surScat(df, vPatterns=c("var1", "var2"))
#    scatter <- addClusters(scatter, lca_fit)  # auto-collapsed
#    report <- attr(scatter, "collapse_LCA(3)")  # access the collapse statistics
#
# 4. Use replaceClusters() to replace all existing cluster columns:
#    scatter <- replaceClusters(scatter, new_clusters, name="NewMethod")
#
# 5. If your clustering result is soft (probabilities), convert to hard assignments first:
#    clusters <- apply(posterior_probs, 1, which.max)
#    scatter <- addClusters(scatter, clusters)
#
# 6. Use weight parameter if your clustering accounts for case weights:
#    scatter <- addClusters(scatter, clusters, weight = my_weights)
#
# 7. Set sort=FALSE to keep the original cluster numbering from your external method:
#    scatter <- addClusters(scatter, clusters, sort=FALSE)
