#' Global definitions needed for the R-check to pass
#' 

#' Visualize the parallel traversal of the nodes in the tree from the tips to the root.
#' @param tree a phylo object
#' @param x.lim limits on the time (x) -axis
#' @param draw_prune_ranges logical indicating if prune-ranges should be
#' denoted on each plot.
PlotParallelTraversal <- function(tree, x.lim = c(-2, max(node.depth.edgelength(tree)) + 1), draw_prune_ranges = TRUE) {
  # tree must have node-labels
  tree$tip.label <- as.character(tree$tip.label)
  tree$node.label <- as.character(tree$node.label)
  all_nodes <- c(tree$tip.label, tree$node.label)
  ppt <- SPLITT__OrderedTreeStringNodes$new(
    all_nodes[tree$edge[, 1]], 
    all_nodes[tree$edge[, 2]], 
    tree$edge.length)
  
  tree$tip.label <- as.character(sapply(tree$tip.label, function(a) ppt$FindIdOfNode(a)))
  tree$node.label <- as.character(sapply(tree$node.label, function(a) ppt$FindIdOfNode(a)))
  
  all_nodes <- c(tree$tip.label, tree$node.label)
  names(all_nodes) <- NULL
  edge_node_names <- cbind(all_nodes[tree$edge[, 1]], 
                           all_nodes[tree$edge[, 2]])
  
  ending_at <- match(all_nodes, edge_node_names[, 2])
  names(ending_at) <- all_nodes
  
  num_parallel_ranges_prune <- ppt$num_parallel_ranges_prune
  num_levels <- ppt$num_levels
  num_nodes <- ppt$num_nodes
  num_tips <- ppt$num_tips
  
  ranges_id_prune <- ppt$ranges_id_prune
  ranges_id_visit <- ppt$ranges_id_visit
  
  i_prune <- 1
  
  edge_color <- rep("black", nrow(edge_node_names))
  edge_width <- rep(2, nrow(edge_node_names))
  edge_lty <- rep(1, nrow(edge_node_names))
  names(edge_color) <- names(edge_width) <- names(edge_lty) <- edge_node_names[, 2]
  
  tip_bg <- rep("palegreen1", num_tips)
  tip_color <- rep("springgreen4", num_tips)
  tip_frame <- rep("circle", num_tips)
  tip_frame_color <- rep("springgreen4", num_tips)
  names(tip_bg) <- names(tip_color) <- names(tip_frame) <- names(tip_frame_color) <-  tree$tip.label
  tip_id <- match(sort(as.numeric(tree$tip.label)), 
                  as.numeric(tree$tip.label))
  names(tip_id) <- as.character(sort(as.numeric(tree$tip.label)))
  
  node_bg <- rep("palegreen1", length(tree$node.label))
  node_color <- rep("springgreen4", length(tree$node.label))
  node_frame <- rep("circle", length(tree$node.label));
  node_frame_color <- rep("springgreen4", length(tree$node.label));
  node_id <- match(sort(as.numeric(tree$node.label)),
                   as.numeric(tree$node.label))
  names(node_id) <- as.character(sort(as.numeric(tree$node.label)))
  names(node_bg) <- names(node_color) <- names(node_frame) <- names(node_frame_color) <- tree$node.label
  
  DrawNode <- function(node_name) {
    if(1+as.integer(node_name) < 10) {
      cex = 1.24*.85
    } else {
      cex = 1.04*.85
    }
    if(as.numeric(node_name) < num_tips) {
      color <- tip_frame_color[node_name]
      names(color) <- NULL
      par_col_default <- par(col=color)
      
      tiplabels(text = 1 + as.integer(node_name),
                tip = tip_id[node_name],
                frame = tip_frame[node_name],
                bg = tip_bg[node_name],
                col= tip_color[node_name],
                cex = cex)
      par(par_col_default)
    } else {
      color <- node_frame_color[node_name]
      names(color) <- NULL
      par_col_default <- par(col=color)
      
      nodelabels(text = 1 + as.integer(node_name),
                 node = node_id[node_name] + num_tips,
                 frame = node_frame[node_name],
                 bg = node_bg[node_name],
                 col=node_color[node_name],
                 cex = cex)
      par(col=par_col_default)
    }
  }
  
  plot(tree,
       show.tip.label=FALSE, 
       edge.color = edge_color, 
       edge.width = edge_width, 
       edge.lty = edge_lty, 
       x.lim = x.lim,
       no.margin = TRUE)
  
  edgelabels(text = tree$edge.length, 
             #cex = 0.8,
             frame = "none",
             adj = c(0.7, -0.3),
             col = edge_color)
  
  for(node_name in all_nodes) {
    DrawNode(node_name)
  }
  
  # after initialization, set the background color for all nodes back to white
  tip_bg[] <- "white"
  node_bg[] <- "white"
  
  for(i in 1:num_levels) {
    if(i <= num_levels) {
      es <- ending_at[as.character(ranges_id_visit[i]:(ranges_id_visit[i+1] - 1))]
    } else {
      es <- ppt$num_nodes
    }
    
    if(es[1] != num_nodes) {
      num_branches_done <- 0
      
      nodes_to_visit <- as.character((ranges_id_visit[i]:(ranges_id_visit[i+1] - 1)))
      
      if(draw_prune_ranges) {
        nodes_to_prune <- as.character((ranges_id_prune[i_prune]:(ranges_id_prune[i_prune+1] - 1)))
      } else {
        nodes_to_prune <- nodes_to_visit
      }
      
      if(i == 1) {
        # tips
        if(draw_prune_ranges) {
          tip_bg[nodes_to_prune] <- "lightpink1"
          tip_color[nodes_to_prune] <- "indianred4"
          tip_frame_color[nodes_to_prune] <- "indianred4"
        } else {
          tip_bg[nodes_to_visit] <- "lightpink1"
          tip_color[nodes_to_visit] <- "indianred4"
          tip_frame_color[nodes_to_visit] <- "indianred4"
        }
        
      } else {
        if(draw_prune_ranges) {
          node_bg[nodes_to_prune] <- "indianred1"
          node_color[nodes_to_prune] <- "indianred4"
          node_frame_color[nodes_to_prune] <- "indianred4"
        } else {
          node_bg[nodes_to_visit] <- "indianred1"
          node_color[nodes_to_visit] <- "indianred4"
          node_frame_color[nodes_to_visit] <- "indianred4"
        }
      }
      
      while(num_branches_done != length(es)) {
        if(draw_prune_ranges) {
          nodes_to_prune <- as.character((ranges_id_prune[i_prune]:(ranges_id_prune[i_prune+1] - 1)))
        } else {
          nodes_to_prune <- nodes_to_visit
        }
        
        num_branches_done <- num_branches_done + length(nodes_to_prune)
        i_prune <- i_prune+1
        
        if(draw_prune_ranges) {
          edge_color[nodes_to_prune] <- "royalblue4"
          edge_lty[nodes_to_prune] <- 5
        }
        
        plot(tree, #type="cladogram", 
             show.tip.label=FALSE, 
             edge.color = edge_color, 
             edge.width=edge_width, 
             edge.lty = edge_lty, 
             x.lim = x.lim, 
             no.margin = TRUE)
        
        edgelabels(text = tree$edge.length, 
                   #cex = 0.8,
                   frame = "none", 
                   adj = c(0.7, -0.3),
                   col = edge_color)
        
        edge_color[nodes_to_prune] <- "grey94"
        #edge_lty[nodes_to_prune] <- 1
        
        if(i == 1) {
          # level of tips
          tip_bg[nodes_to_prune] <- "lightpink1"
          tip_color[nodes_to_prune] <- "indianred4"
          tip_frame_color[nodes_to_prune] <- "indianred4"
        } else {
          node_bg[nodes_to_prune] <- "lightpink1"
          node_color[nodes_to_prune] <- "indianred4"
          node_frame_color[nodes_to_prune] <- "indianred4"
        }
        
        if(draw_prune_ranges) {
          for(child_node in nodes_to_prune) {
            parent_node <- as.character(ppt$FindIdOfParent(as.integer(child_node)))
            node_bg[parent_node] <- "skyblue1"
            node_color[parent_node] <- "royalblue4"
            node_frame_color[parent_node] <- "royalblue4"
          }  
        }
        
        for(node_name in all_nodes) {
          DrawNode(node_name)
        }
        
        if(draw_prune_ranges) {
          for(child_node in nodes_to_prune) {
            parent_node <- as.character(ppt$FindIdOfParent(as.integer(child_node)))
            node_bg[parent_node] <- "white"
          }  
        }
        
        if(i == 1) {
          # level of tips
          tip_bg[nodes_to_prune] <- "white"
          tip_color[nodes_to_prune] <- "grey90"
          tip_frame_color[nodes_to_prune] <- "grey90"
        } else {
          node_bg[nodes_to_prune] <- "white"
          node_color[nodes_to_prune] <- "grey90"
          node_frame_color[nodes_to_prune] <- "grey90"
        }
      }  
      
      if(i == 1) {
        # tips
        tip_bg[nodes_to_visit] <- "white"
        tip_color[nodes_to_visit] <- "grey90"
        tip_frame_color[nodes_to_visit] <- "grey90"
      } else {
        node_bg[nodes_to_visit] <- "white"
        node_color[nodes_to_visit] <- "grey90"
        node_frame_color[nodes_to_visit] <- "grey90"
      }
    }
  }
  
  ppt
}