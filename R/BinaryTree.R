#' Build a binary tree.
#'
#' @keywords internal
#'
#' @importFrom stats density
BinaryTree <- function(M, minleaf = 1, t = .1)
{
  n <- nrow(M)
  p <- ncol(M)
  combinations <- matrix(NA, n, p)
  if(is.null(colnames(M)))
  {
    colnames(M) <- paste0(rep("M",p), 1:p)
  }
  col_names <- colnames(M)
  colnames(combinations) <- col_names
  root <- tree <- Signtree <- mark_tree <- marks_left <- rootmarks <- list()
  labels <- rep(0, n)
  label_counter <- label_graph <- level <- 1
  S_Comp1 <- S_Comp2 <- Mu_Comp1 <- Mu_Comp2 <- vecMar <- vecMarNode <- c()
  CytEMRes <- CytEM(M, 1:n, minleaf, level, t)
  if(is.null(CytEMRes$ind))
  {
    return(list("labels"= rep(1, n), "Signtree" = Signtree))
  }
  root_ind <- CytEMRes$ind[1]
  mark_left <- c(CytEMRes$ind[-c(1)], CytEMRes$mark_not_dis)
  root[[level]] <- M[,root_ind]
  tree[[level]] <- root
  Signtree[[level]] <- "root"
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
  cste1maxlevel <- log(n + 1)
  cste2maxlevel <- log(2)
  while(cste1maxlevel >= (level + 1) *  cste2maxlevel)
  {
    c_level <- level + 1
    tree[[c_level]] <- Signtree[[c_level]] <- marks_left[[c_level]] <- list()
    n_nodes_at_level <- length(tree[[level]])
    stopping_flag <- 0
    node_counter <- 1
    for (i in 1:n_nodes_at_level)
    {
      flag_pop <- flag_child <- 0
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
        if(!length(mark_left))
        {
          stopping_flag <- stopping_flag + 1
          if(!is.na(temp_node[1]) & length(temp_node)){
            flag_pop <- 1
            labels[temp_node] <- label_counter
          }
          if(flag_pop)
          {
            mark_tree[[level]][[i]] <- as.character(label_counter)
            label_counter <- label_counter + 1
          }
          else
          {
            mark_tree[[level]][[i]] <- NA
          }
          if(stopping_flag == n_nodes_at_level)
          {
            return(list("combinations"=combinations,"labels"=labels,
                        "mark_tree" = mark_tree,
                        "pl_list"= pl_list, "Signtree" = Signtree)
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
            flag_pop <- 1
            labels[temp_node] <- label_counter
            if(flag_pop)
            {
              mark_tree[[level]][[i]] <- as.character(label_counter)
              label_counter <- label_counter + 1
            }
            else
            {
              mark_tree[[level]][[i]] <- NA
            }
            if(stopping_flag==n_nodes_at_level)
            {
              return(list("combinations"=combinations,"labels"=labels,
                          "mark_tree" = mark_tree,
                          "pl_list"= pl_list, "Signtree" = Signtree)
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
        temp_list_lc_sign <- temp_list_rc_sign <- list()
        temp_list_lc_sign[[1]] <- "-"
        temp_list_rc_sign[[1]] <- "+"
        Signtree[[c_level]][node_counter] <- temp_list_lc_sign
        Signtree[[c_level]][node_counter + 1] <- temp_list_rc_sign
        if(level == 1)
        {
          vecMar <- append(vecMar, col_names[root_ind])
        }
        else
        {
          vecMar <- append(vecMar, col_names[ind])
        }
        vecMarNode <- append(vecMarNode, pl_list[[2]][[label_graph]])
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
      else
      {
        temp_list_lc <- temp_list_rc <- list()
        temp_list_lc[[1]] <- temp_list_rc[[1]] <- NA
        tree[[c_level]][node_counter] <- temp_list_lc
        tree[[c_level]][node_counter + 1] <- temp_list_rc
        temp_list_lc_sign <- temp_list_rc_sign <- list()
        temp_list_lc_sign[[1]] <-  temp_list_rc_sign[[1]] <- NA
        Signtree[[c_level]][node_counter] <- temp_list_lc_sign[[1]]
        Signtree[[c_level]][node_counter + 1] <- temp_list_rc_sign[[1]]
        marks_left[[c_level]][[node_counter]] <- numeric()
        marks_left[[c_level]][[node_counter + 1]] <- numeric()
        node_counter <- node_counter + 2
      }
    }
    level <- c_level
    mark_tree[[c_level]] <- list()
  }
  return(list("combinations"=combinations,"labels"=labels,
              "mark_tree" = mark_tree,
              "pl_list"= pl_list, "Signtree" = Signtree)
  )
}


