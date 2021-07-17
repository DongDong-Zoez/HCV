dissimilarity <- function(optimization_domain, dist_method){
  
  # calculate the dissimilarity matrix
  
  dist_matrix <- as.matrix(dist(optimization_domain, method = dist_method))
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

Voronoi_adjacency_matrix <- function(constraint_domain){
  
  # used for building the Delaunay triangulation and Voronoi diagram
  
  n <- nrow(constraint_domain)
  adj_mat <- matrix(0, nrow = n, ncol = n)
  vor <- alphahull::delvor(constraint_domain)
  for(i in 1:nrow(vor$mesh)){
    
    # let adj_mat{xy} = adj_mat{yx} = 1 if x and y have a common Voronoi edge
    # The attribute "mesh" carry the spatial information of constraint_domain
    
    x <- vor$mesh[,1:2][i,][1] # The spatial coordinate of point x
    y <- vor$mesh[,1:2][i,][2] # The spatial coordinate of point y
    adj_mat[x, y] <- 1
    adj_mat[y, x] <- 1
  }
  return(adj_mat)
}

FindMinimum <- function(dist_matrix, adj_mat){
  
  # Find the minimum value among the lower triangle section in dissimilarity matrix
  
  min_dist = c(-1,-1, Inf) # The first two are the (i, j) terns and third is the minimum value
  weighted_matrix <- dist_matrix / adj_mat
  n <- nrow(weighted_matrix)
  for(i in 2:n){
    for(j in 1:(i-1)){
      dist <- weighted_matrix[i, j]
      if(dist <= min_dist[3]){
        dist <- weighted_matrix[i, j]
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
  dslist <- list()
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
      n_h = length(cluster[[h]])
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
                                linkage, iterate, diss = 'none',
                                dist_method = 'euclidean'){
  n <- nrow(constraint_domain)
  if(diss == 'precomputed'){
    dist_matrix <- optimization_domain
  }
  else{
    dist_matrix <- dissimilarity(optimization_domain, dist_method)
  }
  dist_matrix <- as.matrix(dist_matrix)
  adj_mat <- Voronoi_adjacency_matrix(constraint_domain)
  
  result <- AGNES(dist_matrix, adj_mat, linkage, iterate + 1)
  clust <- list()
  clust$label_matrix <- result[[1]]
  clust$merge <- result[[2]]
  clust$height <- result[[3]]
  clust$method <- linkage
  clust$dist.method <- dist_method
  clust <- data_structure(clust, n)
  class(clust) <- 'hclust'
  return(clust)
}

data_structure <- function(ds, n){
  ds$labels <- as.character(1:n)
  class(ds) <- 'hclust'
  ds$order <- order.dendrogram(as.dendrogram(ds))
  return(ds)
}

synthetic_data <- function(k, f, r, n, attribute){
  constraint_domain <- matrix(0, ncol = attribute, nrow = n)
  optimization_domain <- matrix(0, ncol = attribute, nrow = n)
  labels <- rep(-1, n)
  constraint_domain_center = runif(attribute * k, 0, 1)
  constraint_domain_center = matrix(constraint_domain_center, ncol = attribute)
  optimization_domain_center = constraint_domain_center
  for(z in 1:n){
    prob <- rep(0, k)
    p <- runif(attribute, 0, 1)
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

usethis::use_data(t)
usethis::create_package('C:\\Users\\ASUS\\Documents\\R package\\HCV')
usethis::use_vignette("introduction")
