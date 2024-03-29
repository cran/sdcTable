#' Create a hierarchy
#'
#' create_node() is defunct, please use sdcHierarchies::hier_create()
#'
#' @keywords internal
#' @rdname defunct-sdcTable
create_node <- function(...) {
  msg <- "Please use `hier_create()` instead."
  .Defunct(new = "hier_create", package = "sdcHierarchies", msg = msg)
}

#' Add nodes to a hierarchy
#'
#' add_nodes() is defunct, please use sdcHierarchies::hier_add()
#'
#' @keywords internal
#' @rdname defunct-sdcTable
add_nodes <- function(...) {
  msg <- "Please use `hier_add()` instead."
  .Defunct(new = "hier_add", package = "sdcHierarchies", msg = msg)
}

#' Delete nodes from a hierarchy
#'
#' delete_nodes() is defunct, please use sdcHierarchies::hier_delete()
#'
#' @keywords internal
#' @rdname defunct-sdcTable
delete_nodes <- function(...) {
  msg <- "Please use `hier_delete()`"
  .Defunct(new = "hier_delete", package = "sdcHierarchies", msg = msg)
}

#' Rename a node in a hierarchy
#'
#' rename_node() is defunct, please use sdcHierarchies::hier_rename()
#'
#' @keywords internal
#' @rdname defunct-sdcTable
rename_node <- function(...) {
  msg <- "Please use `hier_rename()` instead."
  .Defunct(new = "hier_rename", package = "sdcHierarchies", msg = msg)
}


#' Extract information about a cell in a protected table
#'
#' cellInfo() is defunct, please use [cell_info()]
#'
#' @keywords internal
#' @rdname defunct-sdcTable
cellInfo <- function(...) {
  .Defunct(
    new = "cell_info",
    package = "sdcTable",
    msg = "please use sdcTable::cell_info();
      this is a more flexible way to extract the required information."
  )
}

#' Extract information about a cell in a protected table
#'
#' changeCellStatus() is defunct, please use [change_cellstatus()]
#'
#' @keywords internal
#' @rdname defunct-sdcTable
changeCellStatus <- function(...) {
  .Defunct(
    new = "change_cellstatus",
    package = "sdcTable",
    msg = "please use sdcTable::change_cellstatus() instead"
  )
}
