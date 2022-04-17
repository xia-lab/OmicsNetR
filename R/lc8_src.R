## LC8 package source
## source: https://github.com/LiChenPU/Formula_manipulation/
## Reference: https://www.nature.com/articles/s41592-021-01303-3
##



my_calculate_formula = function(Formula1,Formula2,sign = 1, Valid_formula = FALSE){
  formula2_ls = list()
  for(i in 1:length(Formula2)){
    formula2 = Formula2[i]
    {
      formula2 <- gsub("D", "[2]H", formula2)
      ende2 <- nchar(formula2)
      element2 <- c()
      number2 <- c()
      j <- c(1)
      while (j <= ende2) {
        if (substr(formula2, j, j) == c("[")) {
          b <- j
          while (any(substr(formula2, j, j) == c("]")) !=
                 TRUE) {
            j <- c(j + 1)
          }
          k <- j
          while (any(substr(formula2, j, j) == c("-", ".",  "0", "1",
                                                 "2", "3", "4", "5", "6", "7", "8", "9")) !=
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          element2 <- c(element2, substr(formula2, b, m))
        }
        if (any(substr(formula2, j, j) == c("-", ".",  "0", "1", "2", "3",
                                            "4", "5", "6", "7", "8", "9")) != TRUE) {
          k <- j
          while (any(substr(formula2, j, j) == c("-", ".",  "0", "1",
                                                 "2", "3", "4", "5", "6", "7", "8", "9")) !=
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          element2 <- c(element2, substr(formula2, k, m))
        }
        if (any(substr(formula2, j, j) == c("-", ".",  "0", "1", "2", "3",
                                            "4", "5", "6", "7", "8", "9")) == TRUE) {
          k <- j
          while (any(substr(formula2, j, j) == c("-", ".",  "0", "1",
                                                 "2", "3", "4", "5", "6", "7", "8", "9")) ==
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          number2 <- c(number2, as.numeric(substr(formula2,
                                                  k, m)))
        }
        j <- j + 1

      }

    }
    number2 = number2 * sign
    formula2_ls[[length(formula2_ls)+1]] = list(element2, number2)
  }

  formula1_ls = list()
  for(i in 1:length(Formula1)){
    formula1 = Formula1[i]
    {
      formula1 <- gsub("D", "[2]H", formula1)
      ende2 <- nchar(formula1)
      element1 <- c()
      number1 <- c()
      j <- c(1)
      while (j <= ende2) {
        if (substr(formula1, j, j) == c("[")) {
          b <- j
          while (any(substr(formula1, j, j) == c("]")) !=
                 TRUE) {
            j <- c(j + 1)
          }
          k <- j
          while (any(substr(formula1, j, j) == c("-", ".",  "0", "1",
                                                 "2", "3", "4", "5", "6", "7", "8", "9")) !=
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          element1 <- c(element1, substr(formula1, b, m))
        }
        if (any(substr(formula1, j, j) == c("-", ".",  "0", "1", "2", "3",
                                            "4", "5", "6", "7", "8", "9")) != TRUE) {
          k <- j
          while (any(substr(formula1, j, j) == c("-", ".",  "0", "1",
                                                 "2", "3", "4", "5", "6", "7", "8", "9")) !=
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          element1 <- c(element1, substr(formula1, k, m))
        }
        if (any(substr(formula1, j, j) == c("-", ".",  "0", "1", "2", "3",
                                            "4", "5", "6", "7", "8", "9")) == TRUE) {
          k <- j
          while (any(substr(formula1, j, j) == c("-", ".",  "0", "1",
                                                 "2", "3", "4", "5", "6", "7", "8", "9")) ==
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          number1 <- c(number1, as.numeric(substr(formula1,
                                                  k, m)))
        }
        j <- j + 1
      }
    }
    formula1_ls[[length(formula1_ls)+1]] = list(element1, number1)
  }

  formula_mat = matrix("", nrow = length(formula1_ls), ncol = length(formula2_ls))
  valid_mat = matrix(F, nrow = length(formula1_ls), ncol = length(formula2_ls))
  for(row_i in 1:length(Formula1)){
    element1 = formula1_ls[[row_i]][[1]]
    number1 = formula1_ls[[row_i]][[2]]
    for(col_j in 1:length(Formula2)){
      element2 = formula2_ls[[col_j]][[1]]
      number2 = formula2_ls[[col_j]][[2]]
      count_all = number1
      elem_all = element1
      for(i in 1:length(element2)){
        if(any(element2[i]==elem_all) == TRUE){
          count_all[elem_all==element2[i]] = count_all[elem_all==element2[i]]+number2[i]
        }else{
          count_all = c(count_all,number2[i])
          elem_all = c(elem_all,element2[i])
        }
      }

      if(any(count_all<0)){
        valid_mat[row_i,col_j]=F
      }

      elem_order = order(elem_all)
      count_all = count_all[elem_order]
      elem_all = elem_all[elem_order]

      elem_all=elem_all[count_all!=0]
      count_all=count_all[count_all!=0]
      formula_all=c()
      for (i in 1:length(count_all)) {
        formula_all <- c(formula_all, elem_all[i],
                         count_all[i])
      }
      formula_mat[row_i,col_j]=paste0(formula_all, collapse = "")
    }
  }

  if(length(Formula1)==1 & length(Formula2)==1){
    formula_mat = as.vector(formula_mat)
    valid_mat = as.vector(valid_mat)
  }

  if(Valid_formula){
    return(list(formula_mat, valid_mat))
  } else {
    return(formula_mat)
  }
}

formula_mz <- function(formula = "C2H4O1", charge = 0, elem_table = OmicsNetR::elem_table){

  mz = numeric()
  for(i in 1:length(formula)){
    temp_formula = formula[i]
    if(is.na(temp_formula) | is.null(temp_formula)){
      mz[i] = NA
      next
    }
    if(temp_formula == ""){
      mz[i] = 0
      next
    }

    temp_formula <- gsub("D", "[2]H", temp_formula)
    ende2 <- nchar(temp_formula)
    element2 <- c()
    number2 <- c()
    j <- c(1)
    while (j <= ende2) {
      if (substr(temp_formula, j, j) == c("[")) {
        b <- j
        while (any(substr(temp_formula, j, j) == c("]")) !=
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(temp_formula, j, j) == c("-", ".","0", "1",
                                                   "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element2 <- c(element2, substr(temp_formula, b, m))
      }
      if (any(substr(temp_formula, j, j) == c("-", ".","0", "1", "2", "3",
                                              "4", "5", "6", "7", "8", "9")) != TRUE) {
        k <- j
        while (any(substr(temp_formula, j, j) == c("-", ".","0", "1",
                                                   "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element2 <- c(element2, substr(temp_formula, k, m))
      }
      if (any(substr(temp_formula, j, j) == c("-", ".","0", "1", "2", "3",
                                              "4", "5", "6", "7", "8", "9")) == TRUE) {
        k <- j
        while (any(substr(temp_formula, j, j) == c("-", ".","0", "1",
                                                   "2", "3", "4", "5", "6", "7", "8", "9")) ==
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number2 <- c(number2, as.numeric(substr(temp_formula,
                                                k, m)))
      }
      j <- j + 1

    }

    mz[i]=0

    for(j in 1:length(element2)){
      mz[i] = mz[i] + elem_table$mass[element2[j] == elem_table$element] * number2[j]
    }


    e_mass = 0.00054857990943
    mz[i] = mz[i] - e_mass*charge

  }
  return(mz)
}

#' formula_rdbe
#'
#' @param formula e.g. "C2H4O1"
#' @param elem_table a table records unsaturation
#' @noRd
#' @return the ring and double bond number
#'
#' @examples formula_rdbe(formula = "C2H4O1", elem_table = OmicsNetR::elem_table)
formula_rdbe = function(formula = "C2H4O1", elem_table = OmicsNetR::elem_table){

  rdbe = numeric()
  for(i in 1:length(formula)){

    temp_formula = formula[i]
    temp_formula <- gsub("D", "[2]H", temp_formula)
    ende2 <- nchar(temp_formula)
    element2 <- c()
    number2 <- c()
    j <- c(1)
    while (j <= ende2) {
      if (substr(temp_formula, j, j) == c("[")) {
        b <- j
        while (any(substr(temp_formula, j, j) == c("]")) !=
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(temp_formula, j, j) == c("-", ".", "0", "1",
                                                   "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element2 <- c(element2, substr(temp_formula, b, m))
      }
      if (any(substr(temp_formula, j, j) == c("-", ".", "0", "1", "2", "3",
                                              "4", "5", "6", "7", "8", "9")) != TRUE) {
        k <- j
        while (any(substr(temp_formula, j, j) == c("-", ".", "0", "1",
                                                   "2", "3", "4", "5", "6", "7", "8", "9")) !=
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element2 <- c(element2, substr(temp_formula, k, m))
      }
      if (any(substr(temp_formula, j, j) == c("-", ".", "0", "1", "2", "3",
                                              "4", "5", "6", "7", "8", "9")) == TRUE) {
        k <- j
        while (any(substr(temp_formula, j, j) == c("-", ".", "0", "1",
                                                   "2", "3", "4", "5", "6", "7", "8", "9")) ==
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number2 <- c(number2, as.numeric(substr(temp_formula,
                                                k, m)))
      }
      j <- j + 1

    }

    rdbe[i]=1
    for (j in 1:length(element2)) {
      rdbe[i] = rdbe[i] + elem_table$unsaturation[element2[j] ==
                                                    elem_table$element] * number2[j]
    }
  }
  return(rdbe)
}

#' isotopic_abundance
#'
#' @param formula e.g. "C2H4O1"
#' @param elem_table a table records unsaturation
#' @param isotope e.g. [13]C, or [13]C2 for M+2
#' @noRd
#' @return the ratio of given formula and its isotopic peak
#'
#' @examples isotopic_abundance(formula = "C2H4O1", elem_table = OmicsNetR::elem_table)
isotopic_abundance = function(formula = "C2H4O1", isotope = "[13]C1", elem_table = OmicsNetR::elem_table){
  formula <- gsub("D", "[2]H", formula)
  ende2 <- nchar(formula)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {

    if (substr(formula, j, j) == c("[")) {
      b <- j
      while (any(substr(formula, j, j) == c("]")) !=
             TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formula, j, j) == c("-", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formula, b, m))
    }
    if (any(substr(formula, j, j) == c("-", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formula, k, m))
    }
    if (any(substr(formula, j, j) == c("-", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) ==
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formula,
                                              k, m)))
    }
    j <- j + 1

  }


  number1 = number2
  element1 = element2


  formula <- gsub("D", "[2]H", isotope)
  ende2 <- nchar(formula)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {

    if (substr(formula, j, j) == c("[")) {
      b <- j
      while (any(substr(formula, j, j) == c("]")) !=
             TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formula, j, j) == c("-", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
        if(j>ende2){   break      }
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formula, b, m))
    }
    if(j>ende2){
      number2 = c(number2, 1)
      break
    }
    if (any(substr(formula, j, j) == c("-", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
        if(j>ende2){   break      }
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formula, k, m))
    }
    if(j>ende2){
      number2 = c(number2, 1)
      break
    }
    if (any(substr(formula, j, j) == c("-", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) ==
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formula,
                                              k, m)))
    }
    j <- j + 1

  }

  selected_isotopes = grepl("\\[\\d+\\]",element2)
  number_isotopes = number2[selected_isotopes]
  element_isotopes = element2[selected_isotopes]

  ratio = 1
  for(i in 1:length(element_isotopes)){
    element_isotope = element_isotopes[i]
    element_parent = gsub("\\[\\d+\\]", "", element_isotope)
    num_isotope_formula2 = number_isotopes[i]


    num_parent_formula1 = 0
    num_isotpe_formula1 = 0
    if(any(element1 == element_parent)){
      num_parent_formula1 = number1[element1 == element_parent]
    }
    if(any(element1 == element_isotope)){
      num_isotpe_formula1 = number1[element1 == element_isotope]
    }

    natural_abundance = elem_table$abundance[elem_table$element == element_isotope]

    inten_formula1 = choose((num_parent_formula1+num_isotpe_formula1),num_isotpe_formula1)*natural_abundance^(num_isotpe_formula1)
    inten_formula2 = choose((num_parent_formula1+num_isotpe_formula1),(num_isotpe_formula1+num_isotope_formula2))*natural_abundance^(num_isotpe_formula1+num_isotope_formula2)

    ratio = ratio * (inten_formula2/inten_formula1)

  }

  return(ratio)
}

#' my_break_formula
#'
#' @param formula any chemical formula
#' @noRd
#' @return a list containing elem and count
#'
#' @examples my_break_formula("[13]C1C-1H4N2")
my_break_formula = function(formula = "[13]C1C-1H4N2"){
  formula <- gsub("D", "[2]H", formula)
  ende2 <- nchar(formula)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {
    if (substr(formula, j, j) == c("[")) {
      b <- j
      while (any(substr(formula, j, j) == c("]")) !=
             TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formula, j, j) == c("-", ".", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formula, b, m))
    }
    if (any(substr(formula, j, j) == c("-", ".", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", ".", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formula, k, m))
    }
    if (any(substr(formula, j, j) == c("-", ".", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", ".", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) ==
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formula,
                                              k, m)))
    }
    j <- j + 1

  }

  return(list(elem=element2, count=number2))

}

#my function to process formula
table_to_formula <- function(temp_formula){
    r = c()
    elem = temp_formula$elem
    count = temp_formula$count
    len = length(elem)
    for(i in 1:len){
      if(count[i]!=0){
        r=paste(r, elem[i],  count[i], sep="")
      }
    }
    if(sum(count<0)!=0){
      r=paste(r, "Illegal_formula")
    }
    return (r)
}

#' elem_num_query
#'
#' @param formula for example "C2H4O2"
#' @param elem_query it refers to element in query
#'
#' @return It returns the number of element in formula
#' @noRd
#'
#' @examples elem_num_query("C2H4O2","C")
elem_num_query = function(formula,elem_query){


  if(!is.character(formula)|is.na(formula)){return(NA)}

  formula <- gsub("D", "[2]H", formula)
  ende2 <- nchar(formula)
  element2 <- c()
  number2 <- c()
  j <- c(1)
  while (j <= ende2) {
    #browser()
    if (substr(formula, j, j) == c("[")) {
      b <- j
      while (any(substr(formula, j, j) == c("]")) !=
             TRUE) {
        j <- c(j + 1)
      }
      k <- j
      while (any(substr(formula, j, j) == c("-", ".", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      element2 <- c(element2, substr(formula, b, m))
    }
    if (any(substr(formula, j, j) == c("-", ".", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", ".", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) !=
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      element2 <- c(element2, substr(formula, k, m))
    }
    if (any(substr(formula, j, j) == c("-", ".", "0", "1", "2", "3",
                                       "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- j
      while (any(substr(formula, j, j) == c("-", ".", "0", "1",
                                            "2", "3", "4", "5", "6", "7", "8", "9")) ==
             TRUE) {
        j <- c(j + 1)
      }
      m <- c(j - 1)
      j <- c(j - 1)
      number2 <- c(number2, as.numeric(substr(formula,
                                              k, m)))
    }
    j <- j + 1
  }

  if(any(element2 == elem_query )){
    num_query = number2[elem_query == element2]

  }
  else{
    num_query=0

  }
  return(num_query)
}





