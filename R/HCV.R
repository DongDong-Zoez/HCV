dissimilarity <- function(optimization_domain, dist_method){

  # calculate the dissimilarity matrix
  if(dist_method == 'correlation'){
    dist_matrix <- 1 - cor(t(optimization_domain))
    return(dist_matrix)
  }

  if(dist_method == 'abscor'){
    dist_matrix <- 1 - abs(cor(t(optimization_domain)))
    return(dist_matrix)
  }

  if(dist_method == 'abscor'){
    dist_matrix <- 1 - abs(cor(t(optimization_domain)))
    return(dist_matrix)
  }
  else{
    dist_matrix <- as.matrix(dist(optimization_domain, method = dist_method))
  }
  return(dist_matrix)
}

LanceWilliams_algorithm <- function(n_i, n_j, n_h){

  # Dynamic programming update the dissimilarity matrix

  n_k <- n_i + n_j # merge i-cluster and j-cluster, denoted by k-cluster

  dict <- vector(mode = 'list', length = 8) # Build the dictionary to store the LanceWilliams coefficients
  names(dict) <- c('single', 'complete', 'average', 'centroid', 'ward.D', 'ward.D2','median', 'weight')
  dict[[1]] <- c(0.5,0.5,0,-0.5)
  dict[[2]] <- c(0.5,0.5,0,0.5)
  dict[[3]] <- c(n_i/n_k, n_j/n_k, 0, 0)
  dict[[4]] <- c(n_i/n_k, n_j/n_k, -n_i*n_j/n_k/n_k, 0)
  dict[[5]] <- c((n_i+n_h)/(n_i+n_j+n_h), (n_j+n_h)/(n_i+n_j+n_h), -n_h/(n_i+n_j+n_h), 0)
  dict[[6]] <- c((n_i+n_h)/(n_i+n_j+n_h), (n_j+n_h)/(n_i+n_j+n_h), -n_h/(n_i+n_j+n_h), 0)
  dict[[7]] <- c(0.5, 0.5, -0.25, 0)
  dict[[8]] <- c(0.5, 0.5, 0, 0)
  return(dict)
}

DistanceBetweenCluster <- function(linkage, alpha_i, alpha_j, beta, gamma, d_hi, d_hj, d_ij){

  # return the distance between cluster by given linkage

  if(linkage == 'ward.D2'){
    return((alpha_i * (d_hi ** 2) + alpha_j * (d_hj ** 2) + beta * (d_ij ** 2) + gamma * (d_hi - d_hj) ** 2) ** 0.5)
  }
  else{
    return(alpha_i * d_hi + alpha_j * d_hj + beta * d_ij + gamma * abs(d_hi - d_hj))
  }
}

Voronoi_adjacency_matrix <- function(constraint_domain, boundless){

  # used for building the Delaunay triangulation and Voronoi diagram
  n <- nrow(constraint_domain)
  adj_mat <- matrix(0, nrow = n, ncol = n)
  vor <- alphahull::delvor(constraint_domain)
  if(boundless==T){
    for(i in 1:nrow(vor$mesh)){

      # let adj_mat{xy} = adj_mat{yx} = 1 if x and y have a common Voronoi edge
      # The attribute "mesh" carry the spatial information of constraint_domain

      if(vor$mesh[i,11]==0 & vor$mesh[i,12]==0){
        x <- vor$mesh[,1:2][i,][1] # The spatial coordinate of point x
        y <- vor$mesh[,1:2][i,][2] # The spatial coordinate of point y
        adj_mat[x, y] <- 1
        adj_mat[y, x] <- 1
      }
    }
  }
  if(boundless==F){
    for(i in 1:nrow(vor$mesh)){

      # let adj_mat{xy} = adj_mat{yx} = 1 if x and y have a common Voronoi edge
      # The attribute "mesh" carry the spatial information of constraint_domain

      x <- vor$mesh[,1:2][i,][1] # The spatial coordinate of point x
      y <- vor$mesh[,1:2][i,][2] # The spatial coordinate of point y
      adj_mat[x, y] <- 1
      adj_mat[y, x] <- 1
    }
  }
  return(adj_mat)
}

FindMinimum <- function(dist_matrix, adj_mat){

  # Find the minimum value among the lower triangle section in dissimilarity matrix

  min_dist = c(-1,-1, Inf) # The first two are the (i, j) terns and third is the minimum value
  weighted_matrix <- dist_matrix / adj_mat
  weighted_matrix[adj_mat==0]=Inf
  n <- nrow(weighted_matrix)
  for(i in 2:n){
    for(j in 1:(i-1)){
      dist <- weighted_matrix[i, j]
      if(dist <= min_dist[3]){
        min_dist[3] = dist
        min_dist[1] = i
        min_dist[2] = j
      }
    }
  }
  return(min_dist)
}

AGNES <- function(dist_matrix, adj_mat, linkage, iterate){
  count <- 0
  n <- nrow(dist_matrix)
  dslist <- list() # data structure list
  height <- rep(0, n - 1)
  merge <- matrix(0, ncol = 2, nrow = n - 1)
  n_stat <- nrow(dist_matrix)
  cuttree <- matrix(-1, nrow = n - iterate + 2, ncol = n_stat)
  cluster <- vector(mode = 'list', length = n)
  key <- which(c('single', 'complete', 'average', 'centroid', 'ward.D', 'ward.D2',
                 'median', 'weight') == linkage)
  for(i in 1:n){
    cluster[[i]] <- i
  }
  while(iterate <= n + 1){
    min_dist <- FindMinimum(dist_matrix, adj_mat)
    n_i <- length(cluster[[min_dist[1]]])
    n_j <- length(cluster[[min_dist[2]]])
    for(h in 1:n){
      n_h <- length(cluster[[h]])
      coef <- LanceWilliams_algorithm(n_i, n_j, n_h)[[key]]

      if(h != min_dist[1]){
        dist_matrix[h, min_dist[1]] <- DistanceBetweenCluster(
          linkage, coef[1], coef[2],coef[3],coef[4],
          dist_matrix[h, min_dist[1]], dist_matrix[h, min_dist[2]],
          dist_matrix[min_dist[1], min_dist[2]]
        )

        dist_matrix[min_dist[1], h] <- dist_matrix[h, min_dist[1]]

        adj_mat[h, min_dist[1]] <- (adj_mat[h, min_dist[1]] | adj_mat[h, min_dist[2]])
        adj_mat[min_dist[1], h] <- adj_mat[h, min_dist[1]]
      }
    }


    if(n != 2){
      dist_matrix <- dist_matrix[,-min_dist[2]]
      dist_matrix <- dist_matrix[-min_dist[2],]
      adj_mat <- adj_mat[,-min_dist[2]]
      adj_mat <- adj_mat[-min_dist[2],]
    }
    else if(n == 2){
      dist_matrix <- as.matrix(dist_matrix[,-min_dist[2]])
      dist_matrix <- as.matrix(dist_matrix[-min_dist[2],])
      adj_mat <- as.matrix(adj_mat[,-min_dist[2]])
      adj_mat <- as.matrix(adj_mat[-min_dist[2],])
    }

    n <- n - 1
    count <- count + 1

    dslist[[count]] <- c(cluster[[min_dist[1]]], cluster[[min_dist[2]]])
    height[count] <- min_dist[[3]]


    if(length(cluster[[min_dist[1]]]) == 1){
      merge[count, 1] <- -cluster[[min_dist[1]]]
    }
    if(length(cluster[[min_dist[2]]]) == 1){
      merge[count, 2] <- -cluster[[min_dist[2]]]
    }

    if(length(cluster[[min_dist[1]]]) > 1){
      for(s in count:1){
        if(all(cluster[[min_dist[1]]] %in% dslist[[s]])){
          merge[count, 1] <- s
        }
      }
    }
    if(length(cluster[[min_dist[2]]]) > 1){
      for(s in count:1){
        if(all(cluster[[min_dist[2]]] %in% dslist[[s]])){
          merge[count, 2] <- s
        }
      }
    }


    cluster[[min_dist[1]]] <- c(cluster[[min_dist[1]]], cluster[[min_dist[2]]])
    cluster <- cluster[-min_dist[2]]
    cuttree[n - iterate + 3,] <- Clusterlabels(cluster, n_stat)

  }
  return(list(cuttree, merge, height))
}

Clusterlabels <- function(cluster, n){

  # return the label for the iteration

  labels <- rep(-1, n)
  for(i in 1:length(cluster)){
    for(item in cluster[[i]]){
      labels[item] <- i
    }
  }
  return(labels)
}

HierarchicalVoronoi <- function(constraint_domain, optimization_domain,
                                linkage, iterate=2, diss = 'none', adjacency=F,
                                dist_method = 'euclidean', weighted=F, boundless=F,
                                clusterGain = F){
  n <- nrow(constraint_domain)
  if(diss == 'precomputed'){
    dist_matrix <- optimization_domain
  }
  else{
    dist_matrix <- dissimilarity(optimization_domain, dist_method)
  }
  dist_matrix <- as.matrix(dist_matrix)

  if(adjacency == T){
    adj_mat <- constraint_domain
  }
  else{
    adj_mat <- Voronoi_adjacency_matrix(constraint_domain, boundless)
  }
  adj_mat <- as.matrix(adj_mat)

  if(weighted){
    adj_bool <- delaunayEdge(constraint_domain)
    adj_mat <- (adj_mat & adj_bool) * 1
  }

  result <- AGNES(dist_matrix, adj_mat, linkage, iterate + 1)
  clust <- list()

  clust$label_matrix <- result[[1]]
  clust$merge <- result[[2]]
  clust$height <- result[[3]]
  clust$method <- linkage
  clust$dist.method <- dist_method
  clust <- data_structure(clust, n)
  class(clust) <- 'hclust'

  if(clusterGain){
    clusterGain <- vector(length = (n-2))
    clusterMember <- TreeStructure(clust)$fatherNode
    for(i in 1:(n-2)){
      clusterLabels <- as.numeric(cutree(clust, i+1))
      Group <- cbind(optimization_domain, clusterLabels)
      Group <- as.data.frame(Group)
      adjList <- adjacencyCluster(clusterLabels, adj_mat)
      localMean <- LocalMean(adjList, clusterLabels, optimization_domain)
      perGroupMean <- aggregate(. ~ clusterLabels, Group, mean)[,-1]
      perGroupCount <- aggregate(. ~ clusterLabels, Group, NROW)[,2]
      clusterGain[i] <- sum((perGroupCount - 1) *
                              ((perGroupMean - localMean)**2))
    }
    clust$clusterGain <- clusterGain
  }

  return(clust)
}

data_structure <- function(ds, n){
  ds$labels <- as.character(1:n)
  class(ds) <- 'hclust'
  ds$order <- order.dendrogram(as.dendrogram(ds))
  return(ds)
}

synthetic_data <- function(k, f, r, n, attribute, spatial, homogeneity = T){
  constraint_domain <- matrix(0, ncol = spatial, nrow = n)
  optimization_domain <- matrix(0, ncol = attribute, nrow = n)
  labels <- rep(-1, n)
  constraint_domain_center = runif(spatial * k, 0, 1)
  constraint_domain_center = matrix(constraint_domain_center, ncol = spatial)
  optimization_domain_center = runif(attribute * k, 0, 1)
  optimization_domain_center = matrix(optimization_domain_center, ncol = attribute)
  if(homogeneity){
    optimization_domain_center <- constraint_domain_center
  }
  for(z in 1:n){
    prob <- rep(0, k)
    p <- runif(spatial, 0, 1)
    constraint_domain[z,] <- p
    for(i in 1:k){
      q <- (1 / (sum((p - constraint_domain_center[i,]) ** 2))**0.5) ** f
      sums <- 0
      for(j in 1:k){
        sums <- sums + (1 / (sum((p - constraint_domain_center[j,]) ** 2))**0.5) ** f
      }
      prob[i] <- q / sums
    }
    labels[z] <- sample(1:k, replace = T, size = 1, prob = prob)
    mean <- optimization_domain_center[labels[z],]
    cov <- diag(attribute) * r
    optimization_domain[z,] <- MASS::mvrnorm(1, mean, cov)
  }
  list <- list(geo = constraint_domain, feat = optimization_domain, labels = labels,
               constraint_domain_center = constraint_domain_center,
               optimization_domain_center = optimization_domain_center)
  return(list)
}

delaunayEdge <- function(constraint_domain){
  vor <- delvor(constraint_domain)
  n_k <- nrow(constraint_domain)
  adj_bool <- matrix(1, ncol=n_k, nrow=n_k)
  idx <- vor$mesh[,1:2]
  edgelength <- 0
  for(j in 1:nrow(vor$mesh)){
    len <- (sum((constraint_domain[idx[j,1],] - constraint_domain[idx[j,2],])**2))**0.5
    edgelength <- edgelength + len
  }

  avelen <- edgelength / nrow(vor$mesh)

  for(j in 1:nrow(vor$mesh)){
    len <- (sum((constraint_domain[idx[j,1],] - constraint_domain[idx[j,2],])**2))**0.5
    if(avelen < len){
      adj_bool[idx[j,1], idx[j,2]] <- adj_bool[idx[j,2], idx[j,1]] <- 0
    }
  }
  edge <- list()
  edge$avelen <- avelen
  edge$count <- nrow(vor$mesh)
  edge$bool <- adj_bool
  return(edge)
}

plotMap <- function(map, feat, n = 10, color = rainbow(303), main = "",
                    bar_title = "rank", bar_step = 2000) {
  layout(t(1:2),widths = c(6,1))
  par(mar = c(4,4,1,0.5))
  plot(map, col = color[feat], main = main)
  par(mar = c(5,1,5,2.5))
  by <- (max(feat) - min(feat)) / (bar_step - 1)
  z <- seq(from = min(feat), to = max(feat), by = by)
  image(y = z,
        z = t(z),
        col = color[z],
        axes = FALSE,
        main = bar_title,
        cex.main = 1)
  axis(4, cex.axis = 0.8, mgp = c(0,.5,0))
}

TreeStructure <- function(hclust)
{
  merge <- hclust$merge
  n <- nrow(merge) + 1
  leftNode <- list()
  rightNode <- list()
  fatherNode <- list()
  height <- hclust$height

  for(i in 1:(n-1)){
    node <- merge[i,1]

    if(node < 0)
      leftNode[[i]] <- -node
    else
      leftNode[[i]] <- fatherNode[[node]]

    node <- merge[i,2]

    if(node < 0)
      rightNode[[i]] <- -node
    else
      rightNode[[i]] <- fatherNode[[node]]

    fatherNode[[i]] <- sort(c(leftNode[[i]], rightNode[[i]]))
  }

  return(list(leftNode=leftNode, rightNode=rightNode, fatherNode=fatherNode,
              height=height))
}

adjacencyCluster <- function(clusterLabels, adj_mat){

  n <- nrow(adj_mat)
  k <- max(clusterLabels)
  adjList <- vector(mode = 'list', length = k)

  for(i in 1:k){
    perGroupadj <- adj_mat[clusterLabels == i,]
    perGroupadj <- matrix(perGroupadj, ncol = n)
    for(j in 1:nrow(perGroupadj)){
      adjList[[i]] <- unique(clusterLabels[perGroupadj[j,] == 1])

    }

  }
  return(adjList)
}

LocalMean <- function(adjList, clusterLabels, optimization_domain){

  n_feat <- ncol(optimization_domain)
  k <- length(adjList)
  meanMatrix <- matrix(0, nrow = k, ncol = n_feat)

  for(i in 1:k){

    localGroup <- unique(c(i, adjList[[i]]))
    meanMatrix[i,] <- colMeans(optimization_domain[clusterLabels %in% localGroup, ])
  }

  return(meanMatrix)
}

newCrit <- function(constraint_domain, optimization_domain, hclust, k, standardize = T){
  require(alphahull)
  n <- nrow(optimization_domain)
  edgelength <- rep(0, k)
  edgecount <- rep(0, k)
  labels <- as.numeric(cutree(hclust, k))
  if(standardize){
    constraint_domain <- scale(constraint_domain)
    optimization_domain <- scale(optimization_domain)
  }
  constraint_domain_kcenter <- matrix(0, nrow = k, ncol = ncol(constraint_domain))
  optimization_domain_kcenter <- matrix(0, nrow = k, ncol = ncol(optimization_domain))
  constraint_domain_center <- colMeans(constraint_domain)
  optimization_domain_center <- colMeans(optimization_domain)
  idx <- edgelength(constraint_domain, hclust, k)
  WSS <- rep(0, k)
  weightedWSS <- rep(0, k)
  for(i in 1:k){
    constraint_domain_kcenter[i,] <- colMeans(constraint_domain[labels == i,])
    optimization_domain_kcenter[i,] <- colMeans(optimization_domain[labels == i,])
    WG <- sweep(optimization_domain[labels == i,], 2, optimization_domain_kcenter[i,])
    WWG <- sweep(constraint_domain[labels == i,], 2, constraint_domain_kcenter[i,])
    WSS[i] <- sum(WG**2)
    weightedWSS[i] <- sum(WWG**2)
    edgecount[i] <- table(idx[,4]==i)['TRUE']
    edgelength[i] <- mean(idx[idx[,4]==i,3])
  }
  edgesum <- sum(edgelength * edgecount) / sum(edgecount)
  alpha <- (edgelength - edgesum) / edgesum
  alpha <- 1 / (1 + exp(-alpha))
  total <- sum(WSS * alpha + weightedWSS * (1 - alpha)) / k

  return(total)
}

SMI <- function(constraint_domain, optimization_domain, hclust, max_k){
  index <- vector(length = max_k-1)
  for(i in 2:max_k){
    index[i-1] <- newCrit(constraint_domain, optimization_domain, hclust, i)
  }
  diff <- -diff(index)
  ratio <- diff / index[-(max_k-1)]
  bestCluster <- which.max(ratio)+2
  res <- list()
  res$bestCluster <- bestCluster
  res$index <- index
  res$ratio <- ratio
  return(res)
}

edgelength <- function(constraint_domain, hclust_obj, k){
  require(alphahull)
  vor <- delvor(constraint_domain)
  n_k <- nrow(constraint_domain)
  idx <- vor$mesh[,c(1,2,11,12)]
  idx <- cbind(idx, matrix(0, ncol=2, nrow=nrow(idx)))
  group <- rep(0, nrow(idx))
  clusterLabel <- as.numeric(cutree(hclust_obj, k))
  edgelength <- 0
  for(j in 1:nrow(vor$mesh)){
    idx[j,5] <- sum((constraint_domain[idx[j,1],] - constraint_domain[idx[j,2],])**2)
    if(clusterLabel[idx[j,1]] == clusterLabel[idx[j,2]]
       & (idx[j,3]==0 & idx[j,4]==0)){
      idx[j,6] <- clusterLabel[idx[j,1]]
    }
  }

  idx <- idx[,c(1,2,5,6)]
  colnames(idx) <- c('idx1', 'idx2', 'len', 'label')

  return(idx)
}
