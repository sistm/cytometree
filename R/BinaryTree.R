#' Builds a binary tree.
#'
#' @keywords internal
#'
#' @importFrom stats density
BinaryTree <- function(M, minleaf = 1, t = .1, verbose = TRUE)
{
  n <- nrow(M)
  p <- ncol(M)
  if(verbose)
  {
    flush.console()
    pb <- txtProgressBar(min = 0, max = p, style = 3)
    on.exit(close(pb))
    ipbar <- 1
  }
  combinations <- matrix(NA, n, p)
  if(is.null(colnames(M)))
  {
    colnames(M) <- paste0(rep("M",p), 1:p)
  }
  col_names <- colnames(M)
  colnames(combinations) <- col_names
  root <- tree <- mark_tree <- marks_left <- rootmarks <- list()
  labels <- rep(0, n)
  label_counter <- label_graph <- level <- 1
  CytEMRes <- CytEM(M, 1:n, minleaf, level, t)
  if(is.null(CytEMRes$ind))
  {
    if(verbose){setTxtProgressBar(pb, p)}
    return(list("labels"= rep(1, n)))
  }
  root_ind <- CytEMRes$ind[1]
  mark_left <- c(CytEMRes$ind[-c(1)], CytEMRes$mark_not_dis)
  root[[level]] <- M[,root_ind]
  tree[[level]] <- root
  mark_tree[[level]] <- paste0(col_names[root_ind],".",label_graph)
  rootmarks[[level]] <- mark_left
  marks_left[[level]] <- rootmarks
  len_data_plot <- 512
  pl_list <- list()
  KDE <- stats::density(M[,root_ind], n = len_data_plot)
  GMM <- list()
  GMM$x <- KDE$x
  GMM$y <-   GaussMix(GMM$x, CytEMRes$mu1, CytEMRes$mu2,
                      sqrt(CytEMRes$Var1), sqrt(CytEMRes$Var2),
                      CytEMRes$pi1, CytEMRes$pi2)
  pl_list[[1]] <-  pl_list[[2]]<- pl_list[[3]] <- pl_list[[4]] <- list()
  pl_list[[1]][[label_graph]] <- KDE
  pl_list[[2]][[label_graph]] <- paste0(col_names[root_ind],".",label_graph)
  pl_list[[3]][[label_graph]] <- GMM
  pl_list[[4]][[label_graph]] <- paste(paste0("n",label_graph), "=", n,"D =",
                                       round(CytEMRes$nAIC[1],2), sep =" ")
  while(level < (p + 2))
  {
    c_level <- level + 1
    tree[[c_level]] <- marks_left[[c_level]] <- list()
    n_nodes_at_level <- length(tree[[level]])
    stopping_flag <- 0
    node_counter <- 1
    for (i in 1:n_nodes_at_level)
    {
      flag_child <- 0
      temp_node <- tree[[level]][[i]]
      mark_left <- marks_left[[level]][[i]]
      flag_mark_left <- length(mark_left)
      if(level==1)
      {
        flag_child <- 1
        L_child <- CytEMRes$child$L
        R_child <- CytEMRes$child$R
        combinations[L_child,root_ind] <- 0
        combinations[R_child,root_ind] <- 1
      }
      else
      {
        if(verbose)
        {
          ipbar <- p-length(unique(unlist((marks_left[[level]]))))
          setTxtProgressBar(pb, ipbar)
        }
        if(!length(mark_left))
        {
          stopping_flag <- stopping_flag + 1
          labels[temp_node] <- label_counter
          mark_tree[[level]][[i]] <- as.character(label_counter)
          label_counter <- label_counter + 1
          if(stopping_flag == n_nodes_at_level)
          {
            if(verbose){setTxtProgressBar(pb, p)}
            return(list("combinations"=combinations,"labels"=labels,
                        "mark_tree" = mark_tree, "pl_list"= pl_list)
            )
          }
        }
        else
        {
          CytEMRes <- CytEM(M[temp_node,mark_left], temp_node,
                            minleaf, level, t)
          if(is.null(CytEMRes$ind))
          {
            stopping_flag <- stopping_flag + 1
            labels[temp_node] <- label_counter
            mark_tree[[level]][[i]] <- as.character(label_counter)
            label_counter <- label_counter + 1
            if(stopping_flag==n_nodes_at_level)
            {
              if(verbose){setTxtProgressBar(pb, p)}
              return(list("combinations"=combinations,"labels"=labels,
                          "mark_tree" = mark_tree, "pl_list"= pl_list)
              )
            }
          }
          else
          {
            label_graph <- label_graph + 1
            ind <- mark_left[CytEMRes$ind[1]]
            mark_tree[[level]][[i]] <- paste0(col_names[ind],".",label_graph)
            temp_mar_res <- c(CytEMRes$ind[-c(1)],
                              CytEMRes$mark_not_dis)
            mark_left <- mark_left[temp_mar_res]
            flag_mark_left <- length(mark_left)
            flag_child <- 1
            L_child <- CytEMRes$child$L
            R_child <- CytEMRes$child$R
            combinations[L_child,ind] <- 0
            combinations[R_child,ind] <- 1
            KDE <- stats::density(M[temp_node,ind], n = len_data_plot)
            GMM <- list()
            GMM$x <- KDE$x
            GMM$y <-   GaussMix(GMM$x, CytEMRes$mu1, CytEMRes$mu2,
                                sqrt(CytEMRes$Var1), sqrt(CytEMRes$Var2),
                                CytEMRes$pi1, CytEMRes$pi2)
            pl_list[[1]][[label_graph]] <- KDE
            pl_list[[2]][[label_graph]] <- paste0(col_names[ind],".",
                                                  label_graph)
            pl_list[[3]][[label_graph]] <- GMM
            pl_list[[4]][[label_graph]] <- paste(paste0("n",label_graph), "=",
                                                 length(temp_node),"D =",
                                                 round(CytEMRes$nAIC[1],2),
                                                 sep =" ")
          }
        }
      }
      if(flag_child)
      {
        temp_list_lc <- temp_list_rc <- list()
        temp_list_lc[[1]] <- L_child
        temp_list_rc[[1]] <- R_child
        tree[[c_level]][node_counter] <- temp_list_lc
        tree[[c_level]][node_counter + 1] <- temp_list_rc
        if(flag_mark_left)
        {
          marks_left[[c_level]][[node_counter]] <- mark_left
          marks_left[[c_level]][[node_counter + 1]] <- mark_left
        }
        else
        {
          marks_left[[c_level]][[node_counter]] <- numeric()
          marks_left[[c_level]][[node_counter + 1]] <- numeric()
        }
        node_counter <- node_counter + 2
      }
    }
    level <- c_level
    mark_tree[[c_level]] <- list()
  }
}
